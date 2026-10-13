"""Synthetic-only tests for public reference joins and local MVStore lookup."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("prepare_v3_refs", ROOT / "scripts/prepare_track1_v3_references.py")
refs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(refs)


class FakeFasta:
    references = ("chr1",)
    sequence = "CAAAAACTGACGT"

    def fetch(self, chrom, start, end):
        assert chrom == "chr1"
        return self.sequence[start:end]


def tsv(path, columns, rows, preamble=""):
    with gzip.open(path, "wt") as handle:
        handle.write(preamble)
        handle.write("#" + "\t".join(columns) + "\n")
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


class PublicReferenceTests(unittest.TestCase):
    def test_offline_import_after_network_guard_in_fresh_process(self):
        program = (
            'import sys, socket; '
            'sys.path.insert(0, "scripts"); '
            'from track1_v2 import install_network_guard; '
            'install_network_guard(); '
            'from prepare_track1_v3_references import annotate_keys; '
            'assert "httpx" not in sys.modules; '
            '\ntry:\n socket.socket()\nexcept RuntimeError:\n pass\nelse:\n raise AssertionError("network guard disabled")\n'
        )
        result = subprocess.run([sys.executable, '-c', program], cwd=ROOT,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_repeat_left_alignment_and_reference_check(self):
        self.assertEqual(refs.normalize_allele("chr1", 4, "AA", "A", FakeFasta()), ("1", 1, "CA", "C"))
        with self.assertRaisesRegex(ValueError, "REF mismatch"):
            refs.normalize_allele("chr1", 4, "C", "T", FakeFasta())

    def test_padding_trim_shifts_position_and_preserves_haplotype(self):
        fasta = FakeFasta()
        for original in [("chr1", 7, "CTG", "CAG"), ("1", 4, "AA", "A"),
                         ("1", 4, "A", "AAA"), ("1", 7, "CT", "CGT")]:
            chrom, pos, ref, alt = refs.normalize_allele(*original, fasta)
            before = fasta.sequence[:original[1]-1] + original[3] + fasta.sequence[original[1]-1+len(original[2]):]
            after = fasta.sequence[:pos-1] + alt + fasta.sequence[pos-1+len(ref):]
            self.assertEqual(before, after)
            self.assertEqual((chrom, pos, ref, alt), refs.normalize_allele(chrom, pos, ref, alt, fasta))
        self.assertEqual(refs.normalize_allele("chr1", 7, "CTG", "CAG", fasta), ("1", 8, "T", "A"))

    def test_source_receipt_rejects_wrong_url_and_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aliases"; path.write_text("public source")
            receipt = path.with_name(path.name + ".source.json")
            payload = {"url": refs.SOURCES["disease_aliases"]["url"], "bytes": path.stat().st_size,
                       "retrieved_utc": "2026-10-07T00:00:00Z", "sha256": refs.sha256(path)}
            receipt.write_text(json.dumps(payload))
            self.assertEqual(refs.verify_source_receipt("disease_aliases", path), payload)
            for change in [{"url": "https://example.com/different"}, {"bytes": 1}, {"sha256": "bad"}]:
                receipt.write_text(json.dumps(dict(payload, **change)))
                with self.assertRaises(ValueError):
                    refs.verify_source_receipt("disease_aliases", path)

    def test_symbolic_and_unsupported_chromosome_are_rejected(self):
        with self.assertRaises(ValueError):
            refs.normalize_allele("1", 4, "A", "<DEL>")
        with self.assertRaises(ValueError):
            refs.normalize_allele("GL000001", 4, "A", "C")
        self.assertEqual(refs.canonical_chrom("chrMT"), "M")

    def test_header_comments_are_not_data_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "submissions.gz"
            tsv(path, ["VariationID", "SCV"], [{"VariationID": "1", "SCV": "SCV1.1"}],
                "##Comment\n#VariationID: explanation\n")
            self.assertEqual(list(refs.public_tsv(path)), [{"VariationID": "1", "SCV": "SCV1.1"}])

    def test_long_public_fields_are_retained_without_truncation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "submissions.gz"
            description = "synthetic public description " * 6000
            tsv(path, ["VariationID", "ReportedPhenotypeInfo"],
                [{"VariationID": "1", "ReportedPhenotypeInfo": description}])
            self.assertEqual(list(refs.public_tsv(path))[0]["ReportedPhenotypeInfo"], description)

    def test_literal_quotes_do_not_swallow_following_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "submissions.gz"
            with gzip.open(path, "wt") as handle:
                handle.write('#VariationID\tDescription\tSCV\n1\t"A quoted name\tSCV1.1\n2\tOther name"\tSCV2.1\n')
            rows = list(refs.public_tsv(path))
            self.assertEqual([row['SCV'] for row in rows], ['SCV1.1', 'SCV2.1'])
            self.assertEqual(rows[0]['Description'], '"A quoted name')

    def test_malformed_tabular_row_fails_instead_of_losing_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "submissions.gz"
            with gzip.open(path, "wt") as handle:
                handle.write('#VariationID\tSCV\n1\n')
            with self.assertRaisesRegex(ValueError, 'column count'):
                list(refs.public_tsv(path))

    def test_disease_aliases_exclude_gene_entries(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aliases"
            path.write_text("#MIM number\tGeneID\ttype\tSource\tMedGenCUI\tComment\n100001\t1\tphenotype\t-\tC0000001\t-\n100002\t2\tgene\t-\tC0000002\t-\n")
            db = refs.connect(Path(tmp) / "public.sqlite")
            refs.build_condition_aliases(db, path)
            aliases = json.loads(db.execute("SELECT value FROM metadata WHERE key='condition_aliases'").fetchone()[0])
            self.assertEqual(aliases, {"MEDGEN:C0000001": ["OMIM:100001"], "OMIM:100001": ["MEDGEN:C0000001"]})
            db.close()

    def test_scv_is_condition_specific_and_composite_is_not_single_allele(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            columns = ["AlleleID", "VariationID", "Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "PhenotypeIDS"]
            base = dict(AlleleID="1", Assembly="GRCh38", Chromosome="1", ReferenceAlleleVCF="A", AlternateAlleleVCF="C", PhenotypeIDS="OMIM:999999")
            tsv(p / "variants.gz", columns, [dict(base, VariationID="1", PositionVCF="2"),
                dict(base, VariationID="2", PositionVCF="3"), dict(base, VariationID="2", PositionVCF="4"),
                dict(base, VariationID="3", PositionVCF="3"), dict(base, AlleleID="2", VariationID="3", PositionVCF="4", AlternateAlleleVCF="<DEL>"),
                dict(base, VariationID="4", PositionVCF="3"), dict(base, AlleleID="2", VariationID="4", PositionVCF="4", Assembly="GRCh37")])
            fields = ["VariationID", "ClinicalSignificance", "ReportedPhenotypeInfo", "SubmittedPhenotypeInfo", "ReviewStatus", "SCV", "ContributesToAggregateClassification"]
            tsv(p / "submissions.gz", fields, [
                dict(VariationID="1", ClinicalSignificance="Pathogenic", ReportedPhenotypeInfo="C0000001:condition A", SubmittedPhenotypeInfo="OMIM:100001", ReviewStatus="criteria provided, single submitter", SCV="SCV1.1", ContributesToAggregateClassification="yes"),
                dict(VariationID="1", ClinicalSignificance="Benign", ReportedPhenotypeInfo="C0000002:condition B", SubmittedPhenotypeInfo="condition B", ReviewStatus="reviewed by expert panel", SCV="SCV2.1", ContributesToAggregateClassification="no"),
                dict(VariationID="2", ClinicalSignificance="Pathogenic", ReportedPhenotypeInfo="C0000001:condition A", SubmittedPhenotypeInfo="condition A", ReviewStatus="criteria provided, single submitter", SCV="SCV3.1", ContributesToAggregateClassification="yes")])
            db = refs.connect(p / "public.sqlite")
            result = refs.build_clinvar_scv(db, p / "variants.gz", p / "submissions.gz", FakeFasta())
            rows = db.execute("SELECT payload FROM clinvar_rcv WHERE chrom='1' AND pos=2 AND ref='A' AND alt='C'").fetchall()
            self.assertEqual(len(rows), 1)
            assertions = json.loads(rows[0][0])
            self.assertEqual([x["assertion_type"] for x in assertions], ["SCV", "SCV"])
            self.assertEqual(assertions[0]["disease_ids"], ["MEDGEN:C0000001", "OMIM:100001"])
            self.assertNotIn("OMIM:999999", str(assertions))
            self.assertEqual(assertions[1]["disease_ids"], ["MEDGEN:C0000002"])
            self.assertFalse(assertions[1]["contributes_to_aggregate"])
            self.assertEqual(assertions[1]["review_stars"], 3)
            self.assertEqual(result["counts"]["submission_assertions_indexed"], 2)
            self.assertEqual(result["counts"]["composite_variation_ids_excluded"], 2)
            self.assertEqual(db.execute("SELECT DISTINCT variation_id FROM cv_variant").fetchall(), [("1",)])
            self.assertTrue(result["complete"])
            db.close()

    def test_scv_alias_parser_never_uses_variant_conditions(self):
        result = refs.scv_condition({"ReportedPhenotypeInfo": "na:unknown|C0123456:example", "SubmittedPhenotypeInfo": "MONDO:0000001;Orphanet:123", "PhenotypeIDS": "OMIM:999999"})
        self.assertEqual(result, ["MEDGEN:C0123456", "MONDO:0000001", "ORPHA:123"])

    def test_repeated_scv_keeps_distinct_conditions_and_deduplicates_identical_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            fields = ["AlleleID", "VariationID", "Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
            tsv(p / "variants.gz", fields, [dict(AlleleID="1", VariationID="1", Assembly="GRCh38",
                Chromosome="1", PositionVCF="2", ReferenceAlleleVCF="A", AlternateAlleleVCF="C")])
            fields = ["VariationID", "ClinicalSignificance", "ReportedPhenotypeInfo", "SubmittedPhenotypeInfo", "ReviewStatus", "SCV"]
            first = dict(VariationID="1", ClinicalSignificance="Pathogenic", ReportedPhenotypeInfo="C0000001:condition A",
                         SubmittedPhenotypeInfo="OMIM:100001", ReviewStatus="criteria provided, single submitter", SCV="SCV1.1")
            second = dict(first, ClinicalSignificance="Benign", ReportedPhenotypeInfo="C0000002:condition B", SubmittedPhenotypeInfo="OMIM:100002")
            tsv(p / "submissions.gz", fields, [first, first, second])
            db = refs.connect(p / "public.sqlite")
            result = refs.build_clinvar_scv(db, p / "variants.gz", p / "submissions.gz", FakeFasta())
            rows = json.loads(db.execute("SELECT payload FROM clinvar_scv").fetchone()[0])
            self.assertEqual(len(rows), 2)
            self.assertEqual({row['classification'] for row in rows}, {'Pathogenic', 'Benign'})
            self.assertEqual(result['counts']['identical_assertion_rows_deduplicated'], 1)
            self.assertEqual(result['counts']['distinct_assertion_records_retained'], 2)
            db.close()


SYNTHETIC_STORE_JAVA = r'''
import org.h2.mvstore.*;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
public class SyntheticStore {
 public static void main(String[] args) {
  try(MVStore s=new MVStore.Builder().fileName(args[0]).open()) {
   var map=s.openMap("alleles",MvStoreUtil.alleleMapBuilder());
   var key=AlleleProto.AlleleKey.newBuilder().setChr(1).setPosition(101).setRef("A").setAlt("C").build();
   var counts=AlleleProto.Frequency.newBuilder().setFrequencySource(AlleleProto.FrequencySource.GNOMAD_G_AFR).setAc(1).setAn(10000).build();
   var percent=AlleleProto.Frequency.newBuilder().setFrequencySource(AlleleProto.FrequencySource.TOPMED).setFrequency(2.5f).build();
   var missing=AlleleProto.Frequency.newBuilder().setFrequencySource(AlleleProto.FrequencySource.KG).build();
   var score=AlleleProto.PathogenicityScore.newBuilder().setPathogenicitySource(AlleleProto.PathogenicitySource.ALPHA_MISSENSE).setScore(0.8f).build();
   map.put(key,AlleleProto.AlleleProperties.newBuilder().addFrequencies(counts).addFrequencies(percent).addFrequencies(missing).addPathogenicityScores(score).build());
   map.put(AlleleProto.AlleleKey.newBuilder().setChr(25).setPosition(101).setRef("A").setAlt("C").build(),AlleleProto.AlleleProperties.newBuilder().addPathogenicityScores(score).build());
   s.commit();
  }
 }
}
'''


class LocalJavaAdapterTests(unittest.TestCase):
    def test_counts_percentages_missingness_and_mitochondrial_key(self):
        java, lib = refs.java_paths()
        adapter = refs.DEFAULT_OUT / "mvstore_adapter"
        compiler = adapter / "ecj-3.39.0.jar"
        if not java.exists() or not compiler.exists() or not (adapter / "ExomiserLookup.class").exists():
            self.skipTest("Public Java adapter not yet prepared")
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp); source = p / "SyntheticStore.java"; source.write_text(SYNTHETIC_STORE_JAVA)
            classpath = ":".join(str(f) for f in sorted(lib.glob("*.jar")))
            subprocess.run([str(java), "-jar", str(compiler), "-21", "-proc:none", "-cp", classpath, "-d", str(p), str(source)], check=True, capture_output=True)
            db = p / "synthetic.mv.db"
            subprocess.run([str(java), "-cp", str(lib / "*") + ":" + str(p), "SyntheticStore", str(db)], check=True, capture_output=True)
            result = refs.annotate_keys([("1", 101, "A", "C"), ("1", 102, "A", "C"), ("MT", 101, "A", "C")], db, p / "local")
            seen = result["1:101:A:C"]
            self.assertAlmostEqual(seen["population"][0]["af"], .0001)
            self.assertAlmostEqual(seen["population"][1]["af"], .025)
            self.assertIsNone(seen["population"][1]["an"])
            self.assertIsNone(seen["population"][0]["filters"])
            self.assertIsNone(seen["population"][0]["homozygotes"])
            self.assertTrue(seen["population"][0]["homozygotes_zero_or_missing"])
            self.assertIsNone(seen["population"][2]["af"])
            self.assertIsNone(seen["population"][2]["observed"])
            self.assertTrue(seen["other_predictors"]["ALPHA_MISSENSE"]["transcript_unresolved"])
            self.assertEqual(seen["alphamissense_annotations"], [])
            self.assertFalse(result["1:102:A:C"]["database_match"])
            self.assertEqual(result["1:102:A:C"]["population"], [])
            self.assertTrue(result["M:101:A:C"]["database_match"])


if __name__ == "__main__":
    unittest.main()
