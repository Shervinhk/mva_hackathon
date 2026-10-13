"""Download and index complete PUBLIC annotation releases, never patient inputs.

All network requests are to fixed whole-release URLs. The SQLite file exposes
allele-keyed source views; an offline consumer must check metadata completeness.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "references/track1_v3"
SOURCES = {
    "alphamissense": {
        "url": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz",
        "filename": "AlphaMissense_hg38.tsv.gz", "release": "2023 major transcripts hg38",
    },
    "clinvar_rcv": {
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_2026-10.xml.gz",
        "filename": "ClinVarRCVRelease_2026-10.xml.gz", "release": "2026-10",
    },
    "clinvar_variants": {
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
        "filename": "variant_summary.txt.gz", "release": "weekly snapshot recorded in receipt",
    },
    "clinvar_submissions": {
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz",
        "filename": "submission_summary.txt.gz", "release": "weekly snapshot recorded in receipt",
    },
    "disease_aliases": {
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/mim2gene_medgen",
        "filename": "mim2gene_medgen", "release": "public NCBI snapshot recorded in receipt",
    },
}
VALID_CHROMS = {str(n) for n in range(1, 23)} | {"X", "Y", "M"}
SCV_INDEX_VERSION = 4


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def verify_source_receipt(source, path):
    """Verify a complete fixed public source; file existence is insufficient."""
    path = Path(path)
    receipt = json.loads(path.with_name(path.name + ".source.json").read_text())
    if (receipt.get("url") != SOURCES[source]["url"]
            or receipt.get("bytes") != path.stat().st_size
            or not receipt.get("retrieved_utc")
            or receipt.get("sha256") != sha256(path)):
        raise ValueError("Public source receipt does not match the complete cached source")
    return receipt


def free_guard(path, extra=0):
    if shutil.disk_usage(path).free < 6 * 1024**3 + extra:
        raise RuntimeError("Public reference operation would violate the 6 GiB free-space reserve")


def download(source, output):
    # The offline consumer imports this module after disabling sockets.
    # Load network dependencies only in the public-download path.
    import httpx

    spec = SOURCES[source]
    path = output / spec["filename"]
    receipt = path.with_name(path.name + ".source.json")
    if path.exists() and receipt.exists():
        verify_source_receipt(source, path)
        print(f"Verified cached public release: {source}", flush=True)
        return path
    partial = path.with_name(path.name + ".partial")
    offset = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={offset}-"} if offset else {}
    free_guard(output)
    with httpx.stream("GET", spec["url"], headers=headers,
                      timeout=httpx.Timeout(120, read=120), follow_redirects=True) as response:
        response.raise_for_status()
        append = offset > 0 and response.status_code == 206
        if append and not response.headers.get("content-range", "").startswith(f"bytes {offset}-"):
            raise ValueError("Unexpected public download byte range")
        count = offset if append else 0
        expected = response.headers.get("content-length")
        free_guard(output, int(expected) if expected and expected.isdigit() else 0)
        report_at = count + 100_000_000
        with partial.open("ab" if append else "wb") as handle:
            for chunk in response.iter_bytes(1024 * 1024):
                handle.write(chunk)
                count += len(chunk)
                if count >= report_at:
                    free_guard(output)
                    print(f"Downloading public {source}: {count // 1_000_000} MB", flush=True)
                    report_at = count + 100_000_000
        total = response.headers.get("content-range", "").split("/")[-1] if append else response.headers.get("content-length")
        if total and total.isdigit() and count != int(total):
            raise ValueError("Incomplete public reference download")
        partial.replace(path)
        checksum = None
        if "ncbi.nlm.nih.gov" in spec["url"] and spec["url"].endswith(".gz"):
            check = httpx.get(spec["url"] + ".md5", timeout=60)
            check.raise_for_status()
            checksum = re.search(r"\b[0-9a-fA-F]{32}\b", check.text).group().lower()
            with path.open("rb") as handle:
                if hashlib.file_digest(handle, "md5").hexdigest() != checksum:
                    raise ValueError("Public release MD5 checksum mismatch; possible changing weekly snapshot")
        receipt.write_text(json.dumps({"url": spec["url"], "resolved_url": str(response.url),
            "release": spec["release"], "retrieved_utc": now(), "bytes": count,
            "publisher_md5": checksum,
            "sha256": sha256(path), "scope": "complete public release, no patient inputs",
            "headers": {k: response.headers[k] for k in ("etag", "last-modified", "x-goog-hash") if k in response.headers}}, indent=2) + "\n")
    print(f"Downloaded public release: {source}", flush=True)
    return path


def canonical_chrom(chrom):
    chrom = str(chrom).removeprefix("chr")
    return "M" if chrom == "MT" else chrom


def normalize_allele(chrom, pos, ref, alt, fasta=None):
    """Minimal alleles plus reference-verified left alignment of ordinary indels."""
    chrom, pos, ref, alt = canonical_chrom(chrom), int(pos), ref.upper(), alt.upper()
    if chrom not in VALID_CHROMS or pos < 1 or not re.fullmatch("[ACGT]+", ref) or not re.fullmatch("[ACGT]+", alt) or ref == alt:
        raise ValueError("Unsupported public allele representation")
    contig = None
    if fasta is not None:
        contig = chrom if chrom in fasta.references else "chr" + chrom
        if fasta.fetch(contig, pos - 1, pos - 1 + len(ref)).upper() != ref:
            raise ValueError("Public allele REF mismatch")
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if fasta is not None and len(ref) != len(alt):
        while pos > 1 and ref[-1] == alt[-1]:
            before = fasta.fetch(contig, pos - 2, pos - 1).upper()
            ref, alt, pos = before + ref[:-1], before + alt[:-1], pos - 1
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
    return chrom, pos, ref, alt


def connect(path):
    db = sqlite3.connect(path)
    db.execute("PRAGMA journal_mode=DELETE")
    db.execute("PRAGMA synchronous=NORMAL")
    db.execute("PRAGMA cache_size=-100000")
    db.execute("CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL)")
    db.execute("CREATE TABLE IF NOT EXISTS population(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,payload TEXT,PRIMARY KEY(chrom,pos,ref,alt)) WITHOUT ROWID")
    db.execute("INSERT OR IGNORE INTO metadata VALUES('population',?)", (json.dumps({"complete": False, "scope": "not imported; no population frequency can be inferred from absence"}),))
    db.commit()
    return db


def set_metadata(db, key, value):
    db.execute("INSERT OR REPLACE INTO metadata VALUES(?,?)", (key, json.dumps(value, separators=(",", ":"))))
    db.commit()


def mane_genes(path):
    mapping = {}
    if path is None:
        return mapping
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["Ensembl_nuc"].split(".")[0]] = row["symbol"]
    return mapping


def am_rows(path):
    with gzip.open(path, "rt") as handle:
        header = None
        for line in handle:
            if line.startswith("#CHROM") or line.startswith("CHROM\t"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
        if header is None:
            raise ValueError("AlphaMissense hg38 header missing")
        for row in csv.DictReader(handle, fieldnames=header, delimiter="\t"):
            yield row


def build_alphamissense(db, path, mane=None):
    provenance = {"complete": False, "scope": "entire downloaded AlphaMissense hg38 major-transcript release",
                  "source_file": str(path), "source_sha256": sha256(path), "started_utc": now(),
                  "normalization": "single-nucleotide genomic alleles; chr prefix removed", "coverage": "missense SNVs in released transcripts; not all genes or variant classes"}
    set_metadata(db, "alphamissense", provenance)
    db.executescript("""
        DROP VIEW IF EXISTS alphamissense;
        DROP TABLE IF EXISTS am_score;
        DROP TABLE IF EXISTS am_transcript;
        CREATE TABLE am_transcript(id INTEGER PRIMARY KEY,transcript_id TEXT,uniprot_id TEXT,gene TEXT);
        CREATE TABLE am_score(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,transcript INTEGER,
          protein_variant TEXT,score REAL,class TEXT,
          PRIMARY KEY(chrom,pos,ref,alt,transcript)) WITHOUT ROWID;
    """)
    mapping, transcripts, batch, count = mane_genes(mane), {}, [], 0
    for row in am_rows(path):
        key = normalize_allele(row["CHROM"], row["POS"], row["REF"], row["ALT"])
        if row["genome"] != "hg38" or len(key[2]) != 1 or len(key[3]) != 1:
            raise ValueError("AlphaMissense release is not hg38 missense SNVs")
        score = float(row["am_pathogenicity"])
        if not 0 <= score <= 1:
            raise ValueError("AlphaMissense score out of range")
        tx = row["transcript_id"], row["uniprot_id"]
        if tx not in transcripts:
            tx_id = len(transcripts) + 1
            transcripts[tx] = tx_id
            db.execute("INSERT INTO am_transcript VALUES(?,?,?,?)", (tx_id, *tx, mapping.get(tx[0].split(".")[0], "")))
        batch.append((*key, transcripts[tx], row["protein_variant"], score, row["am_class"]))
        count += 1
        if len(batch) >= 50000:
            db.executemany("INSERT INTO am_score VALUES(?,?,?,?,?,?,?,?)", batch)
            db.commit(); batch.clear()
            if count % 1000000 == 0:
                print(f"Indexed public AlphaMissense rows: {count:,}", flush=True)
    if batch:
        db.executemany("INSERT INTO am_score VALUES(?,?,?,?,?,?,?,?)", batch)
    db.executescript("""
        CREATE VIEW alphamissense AS
        SELECT a.chrom,a.pos,a.ref,a.alt,json_group_array(json_object(
          'gene',t.gene,'transcript_id',t.transcript_id,'uniprot_id',t.uniprot_id,
          'protein_variant',a.protein_variant,'am_pathogenicity',a.score,'am_class',a.class,
          'source','AlphaMissense hg38 2023')) AS payload
        FROM am_score a JOIN am_transcript t ON t.id=a.transcript
        GROUP BY a.chrom,a.pos,a.ref,a.alt;
    """)
    provenance.update(complete=True, finished_utc=now(), rows=count, transcripts=len(transcripts))
    set_metadata(db, "alphamissense", provenance)
    return provenance


def review_stars(status):
    value = status.lower()
    if "practice guideline" in value:
        return 4
    if "expert panel" in value:
        return 3
    if "multiple submitters" in value and "no conflicts" in value:
        return 2
    if "criteria provided" in value and "no assertion" not in value:
        return 1
    return 0


def disease_id(db, identifier):
    prefixes = {"OMIM": "OMIM", "MedGen": "MEDGEN", "Orphanet": "ORPHA", "MONDO": "MONDO"}
    if db not in prefixes or not identifier:
        return None
    identifier = identifier.removeprefix(prefixes[db] + ":").removeprefix("Orphanet:")
    return prefixes[db] + ":" + identifier


def parse_rcv(element, fasta=None):
    """Use only the RCV aggregate, never mislabel an SCV/VCV assertion as RCV."""
    assertion = element.find("ReferenceClinVarAssertion")
    if assertion is None or assertion.findtext("RecordStatus") != "current":
        return None, "noncurrent_or_missing_rcv"
    measureset = assertion.find("MeasureSet")
    if measureset is None or measureset.get("Type") != "Variant" or len(measureset.findall("Measure")) != 1:
        return None, "complex_haplotype_or_genotype"
    accession = assertion.find("ClinVarAccession")
    classification = assertion.find("Classifications/GermlineClassification")
    if accession is None or accession.get("Type") != "RCV" or classification is None:
        return None, "no_germline_rcv"
    measure = measureset.find("Measure")
    locations = [loc for loc in measure.findall("SequenceLocation")
                 if loc.get("Assembly") == "GRCh38" and loc.get("positionVCF") and loc.get("referenceAlleleVCF") and loc.get("alternateAlleleVCF")]
    if len(locations) != 1:
        return None, "no_unique_grch38_vcf_location"
    loc = locations[0]
    original = (loc.get("Chr"), int(loc.get("positionVCF")), loc.get("referenceAlleleVCF"), loc.get("alternateAlleleVCF"))
    try:
        key = normalize_allele(*original, fasta=fasta)
    except ValueError as error:
        return None, "ref_mismatch" if "REF mismatch" in str(error) else "unsupported_allele"
    ids, names = set(), set()
    for trait in assertion.findall("TraitSet/Trait"):
        for xref in trait.findall("XRef"):
            identifier = disease_id(xref.get("DB"), xref.get("ID"))
            if identifier:
                ids.add(identifier)
        for name in trait.findall("Name/ElementValue"):
            if name.get("Type") == "Preferred" and name.text:
                names.add(name.text)
    description = classification.find("Description")
    label = "" if description is None else (description.text or "")
    status = classification.findtext("ReviewStatus") or ""
    payload = {"accession": accession.get("Acc") + "." + accession.get("Version", ""),
               "disease_ids": sorted(ids), "disease_names": sorted(names),
               "classification": label, "review_status": status, "review_stars": review_stars(status),
               "conflict": "conflict" in label.lower() or ("conflict" in status.lower() and "no conflicts" not in status.lower()),
               "assertion_type": "RCV", "source": "ClinVarRCVRelease_2026-10",
               "date_last_evaluated": None if description is None else description.get("DateLastEvaluated"),
               "variation_id": measureset.get("ID"), "original_allele": list(original)}
    return (key, payload), None


def build_clinvar_rcv(db, path, fasta=None, fasta_sha256=None):
    provenance = {"complete": False, "scope": "entire downloaded ClinVar RCV release; simple GRCh38 germline allele records",
                  "source_file": str(path), "source_sha256": sha256(path), "started_utc": now(),
                  "reference_fasta_sha256": fasta_sha256,
                  "normalization": "reference-checked left alignment" if fasta is not None else "minimal alleles only; indel equivalent representations may remain"}
    set_metadata(db, "clinvar_rcv", provenance)
    db.executescript("""
        DROP VIEW IF EXISTS clinvar_rcv;
        DROP TABLE IF EXISTS rcv_records;
        CREATE TABLE rcv_records(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,accession TEXT,payload TEXT,
          PRIMARY KEY(chrom,pos,ref,alt,accession)) WITHOUT ROWID;
    """)
    batch, counts = [], Counter()
    with gzip.open(path, "rb") as handle:
        context = ET.iterparse(handle, events=("start", "end"))
        _, root = next(context)
        provenance["release_xml_date"] = root.get("Dated")
        for event, element in context:
            if event != "end" or element.tag != "ClinVarSet":
                continue
            counts["records_seen"] += 1
            parsed, reason = parse_rcv(element, fasta)
            if parsed:
                key, payload = parsed
                batch.append((*key, payload["accession"], json.dumps(payload, separators=(",", ":"))))
                counts["assertions_indexed"] += 1
            else:
                counts[reason] += 1
            element.clear(); root.clear()
            if len(batch) >= 10000:
                db.executemany("INSERT INTO rcv_records VALUES(?,?,?,?,?,?)", batch)
                db.commit(); batch.clear()
            if counts["records_seen"] % 100000 == 0:
                print(f"Indexed public ClinVar RCV records: {counts['records_seen']:,}", flush=True)
    if batch:
        db.executemany("INSERT INTO rcv_records VALUES(?,?,?,?,?,?)", batch)
    db.execute("""CREATE VIEW clinvar_rcv AS SELECT chrom,pos,ref,alt,json_group_array(json(payload)) AS payload
                  FROM rcv_records GROUP BY chrom,pos,ref,alt""")
    provenance.update(complete=True, finished_utc=now(), counts=dict(counts))
    set_metadata(db, "clinvar_rcv", provenance)
    return provenance


def public_tsv(path):
    # ClinVar includes long submitted descriptions and sequence alleles. Python's
    # 128 KiB CSV default rejects valid release rows; retain their full values.
    csv.field_size_limit(max(csv.field_size_limit(), 64 * 1024 * 1024))
    with gzip.open(path, "rt") as handle:
        header = None
        for line in handle:
            if line.startswith(("#VariationID\t", "#AlleleID\t")):
                header = line.rstrip("\r\n").lstrip("#").split("\t")
                break
        if header is None:
            raise ValueError("ClinVar tabular header missing")
        # NCBI's tab-separated release has literal quotation marks in free text,
        # not CSV-quoted cells. CSV quoting can swallow thousands of valid rows.
        for row in csv.DictReader(handle, fieldnames=header, delimiter="\t", quoting=csv.QUOTE_NONE):
            if None in row or any(value is None for value in row.values()):
                raise ValueError("ClinVar tabular row has an unexpected column count")
            yield row


def build_condition_aliases(db, path):
    aliases = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row["type"] != "phenotype" or not re.fullmatch(r"C\d+", row["MedGenCUI"])
                    or not re.fullmatch(r"\d+", row["#MIM number"])):
                continue
            medgen, omim = "MEDGEN:" + row["MedGenCUI"], "OMIM:" + row["#MIM number"]
            aliases.setdefault(medgen, set()).add(omim)
            aliases.setdefault(omim, set()).add(medgen)
    payload = {key: sorted(values) for key, values in aliases.items()}
    set_metadata(db, "condition_aliases", payload)
    set_metadata(db, "condition_aliases_provenance", {"complete": True, "source": SOURCES["disease_aliases"]["url"],
        "source_file": str(path), "source_sha256": sha256(path), "finished_utc": now(),
        "scope": "NCBI mim2gene_medgen rows whose type is phenotype; gene entries excluded", "identifiers": len(payload)})
    return len(payload)


def scv_condition(row):
    """Extract only submitted-condition IDs, never borrow all VCV conditions."""
    text = row.get("ReportedPhenotypeInfo", "")
    ids = {"MEDGEN:" + identifier for identifier in re.findall(r"\b(C\d+)\b", text)}
    submitted = row.get("SubmittedPhenotypeInfo", "")
    for prefix, identifier in re.findall(r"\b(MONDO|OMIM|Orphanet|ORPHA)[:_](\d+)\b", submitted):
        ids.add(("ORPHA" if prefix in {"Orphanet", "ORPHA"} else prefix) + ":" + identifier)
    return sorted(ids)


def build_clinvar_scv(db, variants_path, submissions_path, fasta=None, fasta_sha256=None):
    provenance = {"complete": False, "index_version": SCV_INDEX_VERSION, "assertion_type": "SCV", "scope": "all simple primary-GRCh38 variants in public variant_summary joined to all current public submissions",
                  "source_sha256": sha256(submissions_path), "variant_summary_sha256": sha256(variants_path),
                  "source_files": {"variants": str(variants_path), "submissions": str(submissions_path)},
                  "source_urls": {key: SOURCES[key]["url"] for key in ("clinvar_variants", "clinvar_submissions")},
                  "reference_fasta_sha256": fasta_sha256,
                  "started_utc": now(), "normalization": "reference-verified and left aligned" if fasta is not None else "minimal only",
                  "limitations": ["SCV assertions are individual submissions, not aggregate RCV review status", "Condition IDs come from each SCV only", "Somatic assertions are not treated as germline classifications"]}
    set_metadata(db, "clinvar_scv", provenance)
    set_metadata(db, "clinvar_rcv", {"complete": False, "assertion_type": "SCV", "alias_of": "clinvar_scv"})
    db.executescript("""
        DROP VIEW IF EXISTS clinvar_rcv;
        DROP VIEW IF EXISTS clinvar_scv;
        DROP TABLE IF EXISTS scv_records;
        DROP TABLE IF EXISTS cv_variant;
        DROP TABLE IF EXISTS cv_component;
        DROP TABLE IF EXISTS cv_unresolved;
        CREATE TABLE cv_variant(variation_id TEXT,chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,
          PRIMARY KEY(variation_id,chrom,pos,ref,alt)) WITHOUT ROWID;
        CREATE TABLE scv_records(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,accession TEXT,payload TEXT,
          PRIMARY KEY(chrom,pos,ref,alt,accession,payload)) WITHOUT ROWID;
        CREATE TABLE cv_component(variation_id TEXT,allele_id TEXT,PRIMARY KEY(variation_id,allele_id)) WITHOUT ROWID;
        CREATE TABLE cv_unresolved(variation_id TEXT PRIMARY KEY) WITHOUT ROWID;
    """)
    counts, batch, components, unresolved = Counter(), [], [], []
    for row in public_tsv(variants_path):
        counts["variant_rows_seen"] += 1
        # Count components before assembly/representation filtering: otherwise a
        # haplotype with one unsupported component could become a single allele.
        components.append((row["VariationID"], row["AlleleID"]))
        if len(components) >= 20000:
            db.executemany("INSERT OR IGNORE INTO cv_component VALUES(?,?)", components)
            components.clear()
        if row.get("Assembly") != "GRCh38":
            continue
        try:
            key = normalize_allele(row["Chromosome"], row["PositionVCF"], row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"], fasta)
        except (ValueError, KeyError):
            counts["variant_rows_unsupported_or_ref_mismatch"] += 1
            if canonical_chrom(row.get("Chromosome", "")) in VALID_CHROMS:
                unresolved.append((row["VariationID"],))
            continue
        batch.append((row["VariationID"], *key)); counts["variant_rows_indexed"] += 1
        if len(batch) >= 20000:
            free_guard(Path(variants_path).parent)
            db.executemany("INSERT OR IGNORE INTO cv_variant VALUES(?,?,?,?,?)", batch)
            db.commit(); batch.clear()
    if batch:
        db.executemany("INSERT OR IGNORE INTO cv_variant VALUES(?,?,?,?,?)", batch)
    db.executemany("INSERT OR IGNORE INTO cv_component VALUES(?,?)", components)
    db.executemany("INSERT OR IGNORE INTO cv_unresolved VALUES(?)", unresolved)
    db.commit()
    # A composite variation with multiple distinct primary alleles is not a single-allele assertion.
    counts["composite_variation_ids_excluded"] = db.execute("SELECT COUNT(*) FROM (SELECT variation_id FROM cv_component GROUP BY variation_id HAVING COUNT(*)>1)").fetchone()[0]
    db.execute("DELETE FROM cv_variant WHERE variation_id IN (SELECT variation_id FROM cv_component GROUP BY variation_id HAVING COUNT(*)>1)")
    counts["unresolved_primary_variation_ids"] = db.execute("SELECT COUNT(*) FROM cv_unresolved").fetchone()[0]
    db.execute("DELETE FROM cv_variant WHERE variation_id IN (SELECT variation_id FROM cv_unresolved)")
    counts["multiple_normalized_mapping_ids_excluded"] = db.execute("SELECT COUNT(*) FROM (SELECT variation_id FROM cv_variant GROUP BY variation_id HAVING COUNT(*)>1)").fetchone()[0]
    db.execute("DELETE FROM cv_variant WHERE variation_id IN (SELECT variation_id FROM cv_variant GROUP BY variation_id HAVING COUNT(*)>1)")
    counts["simple_variant_mappings_retained"] = db.execute("SELECT COUNT(*) FROM cv_variant").fetchone()[0]
    db.execute("DROP TABLE cv_component")
    db.execute("DROP TABLE cv_unresolved")
    db.commit()
    print(f"Indexed public ClinVar variant mappings: {counts['variant_rows_indexed']:,}", flush=True)
    batch = []
    cursor = db.cursor()
    for row in public_tsv(submissions_path):
        counts["submission_rows_seen"] += 1
        label = row.get("ClinicalSignificance", "")
        if label.lower() in {"", "na", "not provided", "-"}:
            counts["no_germline_classification"] += 1
            continue
        accession = row.get("SCV", "")
        if not accession.startswith("SCV"):
            counts["no_scv_accession"] += 1
            continue
        key = cursor.execute("SELECT chrom,pos,ref,alt FROM cv_variant WHERE variation_id=?", (row["VariationID"],)).fetchone()
        if key is None:
            counts["no_simple_grch38_mapping"] += 1
            continue
        status = row.get("ReviewStatus", "")
        payload = {"accession": accession, "disease_ids": scv_condition(row),
                   "disease_names": [row.get("ReportedPhenotypeInfo", "")],
                   "classification": label, "review_status": status, "review_stars": review_stars(status),
                   "conflict": "conflict" in label.lower(), "assertion_type": "SCV",
                   "date_last_evaluated": row.get("DateLastEvaluated"), "submitter": row.get("Submitter"),
                   "collection_method": row.get("CollectionMethod"), "origin_counts": row.get("OriginCounts"),
                   "contributes_to_aggregate": {"yes": True, "no": False}.get(row.get("ContributesToAggregateClassification", "").lower()),
                   "variation_id": row["VariationID"], "source": "ClinVar submission_summary snapshot"}
        batch.append((*key, accession, json.dumps(payload, separators=(",", ":"))))
        counts["submission_assertions_indexed"] += 1
        if len(batch) >= 20000:
            free_guard(Path(submissions_path).parent)
            # A public SCV can recur for the same normalized allele with
            # different condition/provenance rows. Preserve every distinct
            # assertion payload; ignore only byte-identical duplicate rows.
            db.executemany("INSERT OR IGNORE INTO scv_records VALUES(?,?,?,?,?,?)", batch)
            db.commit(); batch.clear()
            if counts["submission_assertions_indexed"] % 200000 == 0:
                print(f"Indexed public SCV assertions: {counts['submission_assertions_indexed']:,}", flush=True)
    if batch:
        db.executemany("INSERT OR IGNORE INTO scv_records VALUES(?,?,?,?,?,?)", batch)
    counts["distinct_assertion_records_retained"] = db.execute("SELECT COUNT(*) FROM scv_records").fetchone()[0]
    counts["identical_assertion_rows_deduplicated"] = counts["submission_assertions_indexed"] - counts["distinct_assertion_records_retained"]
    db.executescript("""
        CREATE VIEW clinvar_scv AS SELECT chrom,pos,ref,alt,json_group_array(json(payload)) payload
          FROM scv_records GROUP BY chrom,pos,ref,alt;
        CREATE VIEW clinvar_rcv AS SELECT * FROM clinvar_scv;
    """)
    provenance.update(complete=True, finished_utc=now(), counts=dict(counts),
                      compatibility_view="clinvar_rcv is a legacy view name; every payload explicitly says SCV")
    set_metadata(db, "clinvar_scv", provenance)
    set_metadata(db, "clinvar_rcv", {"complete": True, "assertion_type": "SCV", "alias_of": "clinvar_scv"})
    return provenance


JAVA_SOURCE = r'''
import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.h2.mvstore.*;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;
import org.monarchinitiative.exomiser.core.genome.dao.serialisers.MvStoreUtil;
import com.fasterxml.jackson.databind.ObjectMapper;

public class ExomiserLookup {
  static int chr(String c) { return switch(c) {case "X"->23;case "Y"->24;case "M","MT"->25;default->Integer.parseInt(c);}; }
  public static void main(String[] args) throws Exception {
    if(args.length!=3) throw new IllegalArgumentException("Expected local database, input and output paths");
    ObjectMapper json = new ObjectMapper();
    try(MVStore store = new MVStore.Builder().fileName(args[0]).readOnly().cacheSize(256).open();
        BufferedReader in=Files.newBufferedReader(Path.of(args[1]));
        BufferedWriter out=Files.newBufferedWriter(Path.of(args[2]))) {
      if(!store.hasMap("alleles")) throw new IllegalArgumentException("Missing public alleles map");
      MVMap<AlleleProto.AlleleKey,AlleleProto.AlleleProperties> map=store.openMap("alleles",MvStoreUtil.alleleMapBuilder());
      String line;
      while((line=in.readLine())!=null) {
        String[] a=line.split("\t",-1);
        if(a.length!=4) throw new IllegalArgumentException("Malformed local input row");
        var key=AlleleProto.AlleleKey.newBuilder().setChr(chr(a[0])).setPosition(Integer.parseInt(a[1])).setRef(a[2]).setAlt(a[3]).build();
        var value=map.get(key);
        Map<String,Object> result=new LinkedHashMap<>();
        result.put("key",String.join(":",a)); result.put("database_match",value!=null);
        List<Map<String,Object>> populations=new ArrayList<>();
        Map<String,Object> predictors=new LinkedHashMap<>();
        if(value!=null) {
          for(var f:value.getFrequenciesList()) {
            Map<String,Object> p=new LinkedHashMap<>();
            boolean counts=f.getAn()>0 && f.getAc()>=0 && f.getAc()<=f.getAn();
            boolean percent=f.getAn()==0 && f.getAc()==0 && f.getFrequency()>0 && f.getFrequency()<=100;
            p.put("source","Exomiser:"+f.getFrequencySource().name());
            Double af=counts?Double.valueOf((double)f.getAc()/f.getAn()):percent?Double.valueOf(f.getFrequency()/100.0):null;
            p.put("af",af);
            p.put("af_basis",counts?"AC/AN":percent?"stored_percentage/100":"unknown_or_invalid");
            p.put("an",counts?f.getAn():null);p.put("ac",counts?f.getAc():null);
            p.put("homozygotes",counts && f.getHom()>0?f.getHom():null);
            p.put("homozygotes_zero_or_missing",counts && f.getHom()==0);
            p.put("filters",null);p.put("callable",null);
            Boolean observed=counts?Boolean.valueOf(f.getAc()>0):percent?Boolean.TRUE:null;
            p.put("observed",observed);
            p.put("quality_scope","source preprocessing; per-allele filters not retained in MVStore");
            populations.add(p);
          }
          for(var s:value.getPathogenicityScoresList()) {
            Map<String,Object> p=new LinkedHashMap<>(); p.put("score",s.getScore());
            p.put("source","Exomiser MVStore");p.put("transcript_unresolved",true);
            predictors.put(s.getPathogenicitySource().name(),p);
          }
        }
        result.put("population",populations);result.put("other_predictors",predictors);
        result.put("alphamissense_annotations",List.of());
        out.write(json.writeValueAsString(result));out.newLine();
      }
    }
  }
}
'''


def java_paths():
    java = ROOT / "references/exomiser/jdk-21.0.12.1+1-jre/Contents/Home/bin/java"
    lib = ROOT / "references/exomiser/exomiser-cli-15.1.0/lib"
    return java, lib


def prepare_adapter(output):
    """Compile public-only Java wrapper; invoked separately before any local lookup."""
    directory = output / "mvstore_adapter"
    directory.mkdir(parents=True, exist_ok=True)
    source = directory / "ExomiserLookup.java"
    source.write_text(JAVA_SOURCE)
    compiler = directory / "ecj-3.39.0.jar"
    url = "https://repo.maven.apache.org/maven2/org/eclipse/jdt/ecj/3.39.0/ecj-3.39.0.jar"
    if not compiler.exists():
        import httpx
        response = httpx.get(url, timeout=60); response.raise_for_status()
        compiler.write_bytes(response.content)
        compiler.with_name(compiler.name + ".source.json").write_text(json.dumps({"url": url,"sha256":sha256(compiler),"retrieved_utc":now()},indent=2)+"\n")
    receipt = json.loads(compiler.with_name(compiler.name + ".source.json").read_text())
    if receipt.get("url") != url or receipt.get("sha256") != sha256(compiler):
        raise ValueError("Cached public Java compiler receipt mismatch")
    receipt["bytes"] = compiler.stat().st_size
    compiler.with_name(compiler.name + ".source.json").write_text(json.dumps(receipt, indent=2)+"\n")
    java, lib = java_paths()
    classpath = ":".join(str(path) for path in sorted(lib.glob("*.jar")))
    subprocess.run([str(java), "-jar", str(compiler), "-21", "-proc:none", "-cp", classpath, "-d", str(directory), str(source)], check=True, capture_output=True)
    (directory / "adapter.source.json").write_text(json.dumps({"java_source_sha256":sha256(source),"class_sha256":sha256(directory/"ExomiserLookup.class"),"exomiser_version":"15.1.0","compiler_sha256":sha256(compiler),"compiled_utc":now()},indent=2)+"\n")
    print("Compiled local Exomiser MVStore adapter", flush=True)


def annotate_keys(keys, mvstore_path, output_dir):
    """OFFLINE ONLY. No downloads, logging of alleles, or network APIs.

    The caller controls private output_dir; this routine is not called by the
    public-reference preparation CLI. Database keys are normalized by the caller.
    """
    directory = DEFAULT_OUT / "mvstore_adapter"
    if not (directory / "ExomiserLookup.class").exists():
        raise FileNotFoundError("Compile the public MVStore adapter before offline analysis")
    receipt = json.loads((directory / "adapter.source.json").read_text())
    if (receipt.get("java_source_sha256") != hashlib.sha256(JAVA_SOURCE.encode()).hexdigest()
            or receipt.get("class_sha256") != sha256(directory / "ExomiserLookup.class")):
        raise ValueError("Compiled public adapter is stale or does not match its receipt")
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    inputs, outputs = output_dir / "mvstore_keys.tsv", output_dir / "mvstore_annotations.jsonl"
    with inputs.open("w") as handle:
        for key in keys:
            normalized = normalize_allele(*key)
            handle.write("\t".join(map(str, normalized)) + "\n")
    java, lib = java_paths()
    with (output_dir / "mvstore_adapter.log").open("w") as log:
        subprocess.run([str(java), "-Xmx1G", "-cp", str(lib/"*")+":"+str(directory),
                        "ExomiserLookup", str(mvstore_path), str(inputs), str(outputs)],
                       stdout=log, stderr=log, check=True)
    result = {}
    with outputs.open() as handle:
        for line in handle:
            row = json.loads(line); key = row.pop("key"); result[key] = row
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--sources", nargs="+", choices=SOURCES, default=["clinvar_variants", "clinvar_submissions", "disease_aliases"])
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--prepare-mvstore-adapter", action="store_true")
    parser.add_argument("--mane", type=Path, default=ROOT / "references/track1_v2_functional/MANE.GRCh38.v1.5.summary.txt.gz")
    parser.add_argument("--fasta", type=Path, default=ROOT / "references/track1_v2/hg38.fa")
    args = parser.parse_args()
    if "clinvar_submissions" in args.sources and "clinvar_variants" not in args.sources:
        parser.error("clinvar_submissions requires clinvar_variants for the complete public join")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.prepare_mvstore_adapter:
        prepare_adapter(args.output)
        return
    paths = {}
    for source in args.sources:
        path = args.output / SOURCES[source]["filename"]
        if args.offline:
            verify_source_receipt(source, path)
        else:
            path = download(source, args.output)
        paths[source] = path
    if args.download_only:
        return
    db = connect(args.output / "public_annotations.sqlite")
    fasta = None
    fasta_hash = None
    if "clinvar_rcv" in paths or "clinvar_submissions" in paths:
        import pysam
        fasta = pysam.FastaFile(str(args.fasta))
        fasta_hash = sha256(args.fasta)
    for source, path in paths.items():
        if source == "disease_aliases":
            print(f"Indexed public disease alias identifiers: {build_condition_aliases(db, path):,}", flush=True)
            continue
        if source == "clinvar_variants":
            continue
        if source == "clinvar_submissions":
            existing = db.execute("SELECT value FROM metadata WHERE key='clinvar_scv'").fetchone()
            previous = json.loads(existing[0]) if existing else {}
            if (previous.get("complete") and previous.get("index_version") == SCV_INDEX_VERSION
                    and previous.get("source_sha256") == sha256(path)
                    and previous.get("variant_summary_sha256") == sha256(paths["clinvar_variants"])
                    and previous.get("reference_fasta_sha256") == fasta_hash):
                print("Public SCV index already complete", flush=True); continue
            result = build_clinvar_scv(db, paths["clinvar_variants"], path, fasta, fasta_hash)
            print(json.dumps({"clinvar_scv": result}), flush=True)
            continue
        existing = db.execute("SELECT value FROM metadata WHERE key=?", (source,)).fetchone()
        previous = json.loads(existing[0]) if existing else {}
        if (previous.get("complete") and previous.get("source_sha256") == sha256(path)
                and (source != "clinvar_rcv" or previous.get("reference_fasta_sha256") == fasta_hash)):
            print(f"Public index already complete: {source}", flush=True)
            continue
        result = build_alphamissense(db, path, args.mane) if source == "alphamissense" else build_clinvar_rcv(db, path, fasta, fasta_hash)
        print(json.dumps({source: result}), flush=True)
    db.close()
    if fasta is not None:
        fasta.close()


if __name__ == "__main__":
    main()
