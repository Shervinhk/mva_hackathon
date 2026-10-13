"""Acquire whole-gene/protein public references without opening patient data.

Inputs are public HGNC gene symbols only. All API requests, response bytes and
SHA-256 receipts are retained. The output is a public reference bundle; joining
it to a private callset is a separate, offline operation.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import time

import httpx

ROOT = Path(__file__).resolve().parents[1]
CACHE_ONLY = False
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_QUERY = """query PublicGeneReference($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    gene_id symbol chrom start stop
    variants(dataset: gnomad_r4) {
      variant_id pos ref alt
      exome { ac an af ac_hom filters populations { id ac an ac_hom } }
      genome { ac an af ac_hom filters populations { id ac an ac_hom } }
      joint { ac an homozygote_count hemizygote_count filters
              populations { id ac an ac_hom } }
    }
  }
}"""


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fetch(url: str, path: Path, payload: dict | None = None) -> bytes:
    """Cache exact public responses and verify both request identity and hash."""
    receipt_path = path.with_name(path.name + ".source.json")
    request = {"method": "POST" if payload is not None else "GET", "url": url}
    if payload is not None:
        request["json"] = payload
    if path.exists() and receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt["request"] != request or receipt["sha256"] != digest(path):
            raise ValueError(f"Cached reference request or hash mismatch: {path}")
        return path.read_bytes()
    if CACHE_ONLY:
        raise FileNotFoundError(f"Offline reference or receipt missing: {path}")
    last_error = None
    for attempt in range(3):
        try:
            with httpx.Client(timeout=120, follow_redirects=True) as client:
                response = client.request(request["method"], url, json=payload)
                response.raise_for_status()
            if payload is not None:
                result = response.json()
                if result.get("errors"):
                    raise ValueError(f"Public GraphQL error: {result['errors']}")
                if not result.get("data"):
                    raise ValueError("Public GraphQL response has no data")
            path.write_bytes(response.content)
            receipt_path.write_text(json.dumps({
                "request": request,
                "resolved_url": str(response.url),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "sha256": digest(path), "bytes": path.stat().st_size,
                "response_headers": {key: response.headers[key] for key in
                    ("etag", "last-modified", "content-type") if key in response.headers},
                "scope": "Public whole-gene/protein or whole-release reference; no patient inputs",
            }, indent=2) + "\n")
            print(f"Downloaded {path.name}: {path.stat().st_size} bytes", flush=True)
            return response.content
        except (httpx.HTTPError, ValueError) as error:
            last_error = error
            if attempt < 2:
                time.sleep(2)
    raise RuntimeError(f"Unable to acquire public reference {url}: {last_error}")


def validate_missense_rows(rows: list[dict], sequence: str, protein_scores=None) -> dict:
    observed = set()
    for row in rows:
        variant = row["protein_variant"]
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", variant)
        if not match:
            raise ValueError(f"Unexpected public protein substitution format: {variant}")
        ref, position, alt = match.groups()
        position = int(position)
        if position < 1 or position > len(sequence) or sequence[position - 1] != ref or ref == alt:
            raise ValueError("Public AlphaMissense substitution does not match canonical sequence")
        score = float(row["am_pathogenicity"])
        if not 0 <= score <= 1 or row["am_class"] not in {"likely_benign", "ambiguous", "likely_pathogenic"}:
            raise ValueError("Invalid public AlphaMissense score/class")
        if protein_scores is not None:
            if variant not in protein_scores or abs(protein_scores[variant] - score) > 1e-8:
                raise ValueError("Public genomic/protein AlphaMissense scores disagree")
        observed.add(variant)
    return {"canonical_sequence_matches": True, "unique_substitutions": len(observed)}


def prepare_gene(gene: str, output: Path) -> dict:
    # Queries are based only on a supplied public gene symbol, never variants.
    if not re.fullmatch(r"[A-Z][A-Z0-9-]*", gene):
        raise ValueError("Only public uppercase gene symbols are accepted")
    query = f"(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)"
    url = str(httpx.URL("https://rest.uniprot.org/uniprotkb/search", params={
        "query": query, "format": "json", "fields": "accession,gene_names,sequence", "size": 25}))
    uniprot = json.loads(fetch(url, output / f"{gene}.uniprot.json"))
    entries = uniprot.get("results", [])
    if len(entries) != 1:
        raise ValueError(f"Expected one reviewed human protein for {gene}, got {len(entries)}")
    accession = entries[0]["primaryAccession"]
    alpha = json.loads(fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}",
                             output / f"{gene}.alphafold_metadata.json"))
    canonical = [row for row in alpha if row.get("uniprotAccession") == accession]
    if len(canonical) != 1:
        raise ValueError(f"Ambiguous canonical AlphaFold entry: {gene}/{accession}")
    entry = canonical[0]
    if entry["sequence"] != entries[0]["sequence"]["value"]:
        raise ValueError(f"UniProt/AlphaFold protein sequence mismatch: {gene}")
    files = {}
    protein_scores = None
    for key, label in (("amAnnotationsUrl", "alphamissense_protein"),
                       ("amAnnotationsHg38Url", "alphamissense_hg38")):
        if key not in entry:
            files[label] = {"available": False}
            continue
        path = output / f"{gene}.{label}.csv"
        raw = fetch(entry[key], path)
        text = raw.decode("utf-8-sig")
        rows = list(csv.DictReader(io.StringIO(text)))
        if not rows:
            raise ValueError(f"Empty AlphaMissense reference: {gene}/{label}")
        validation = validate_missense_rows(rows, entry["sequence"], protein_scores)
        if label == "alphamissense_protein":
            protein_scores = {row["protein_variant"]: float(row["am_pathogenicity"]) for row in rows}
            if len(protein_scores) != 19 * len(entry["sequence"]):
                raise ValueError(f"Incomplete whole-protein substitutions: {gene}")
        files[label] = {"available": True, "file": path.name, "rows": len(rows),
                        "columns": list(rows[0]), "validation": validation}
        if label == "alphamissense_hg38":
            files[label]["transcript_ids"] = sorted({row["transcript_id"] for row in rows})
            if {row["genome"] for row in rows} != {"hg38"} or {row["uniprot_id"] for row in rows} != {accession}:
                raise ValueError("Public genomic AlphaMissense build/protein mismatch")
    gnomad_path = output / f"{gene}.gnomad_r4.json"
    gnomad = json.loads(fetch(GNOMAD_API, gnomad_path,
                       {"query": GNOMAD_QUERY, "variables": {"symbol": gene}}))
    annotation = gnomad["data"]["gene"]
    if annotation is None or annotation["symbol"] != gene:
        raise ValueError(f"gnomAD gene mismatch: {gene}")
    for variant in annotation["variants"]:
        for platform in ("exome", "genome", "joint"):
            counts = variant[platform]
            if counts is None:
                continue
            if counts["ac"] is None or counts["an"] is None or not 0 <= counts["ac"] <= counts["an"]:
                raise ValueError("Public gnomAD allele counts inconsistent")
            if counts["an"] and "af" in counts and abs(counts["af"] - counts["ac"] / counts["an"]) > 1e-8:
                raise ValueError("Public gnomAD allele frequency/count mismatch")
    return {"gene": gene, "uniprot_accession": accession,
            "protein_length": len(entry["sequence"]),
            "alphafold_model_version": entry["latestVersion"],
            "alphamissense": files,
            "gnomad": {"api_dataset": "gnomad_r4", "file": gnomad_path.name,
                       "gene_id": annotation["gene_id"],
                       "variants": len(annotation["variants"]),
                       "reference_genome": "GRCh38"}}


def prepare_mane(output: Path) -> dict:
    # Explicit release pin, rather than silently changing a current symlink.
    base = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
    names = ["MANE.GRCh38.v1.5.summary.txt.gz", "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz"]
    for name in names:
        raw = fetch(base + name, output / name)
        if not gzip.decompress(raw):
            raise ValueError(f"Empty MANE reference: {name}")
    return {"release": "1.5", "files": names}


def main():
    global CACHE_ONLY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "references/track1_v2_functional")
    parser.add_argument("--genes", nargs="+", default=["BUB1B"])
    parser.add_argument("--skip-mane", action="store_true")
    parser.add_argument("--offline", action="store_true", help="Verify and reuse cached public references; forbid network requests")
    args = parser.parse_args()
    CACHE_ONLY = args.offline
    args.output.mkdir(parents=True, exist_ok=True)
    report = {"public_reference_only": True, "genes": {}, "errors": {}}
    for gene in args.genes:
        try:
            report["genes"][gene] = prepare_gene(gene, args.output)
            print(f"Public reference ready: {gene}", flush=True)
        except Exception as error:
            report["errors"][gene] = f"{type(error).__name__}: {error}"
            print(f"Public reference incomplete: {gene}: {error}", flush=True)
    if not args.skip_mane:
        try:
            report["mane"] = prepare_mane(args.output)
        except Exception as error:
            report["errors"]["MANE"] = f"{type(error).__name__}: {error}"
    report["created_utc"] = datetime.now(timezone.utc).isoformat()
    report["limitations"] = [
        "gnomad_r4 is the browser API dataset identifier; the API does not supply a minor release per result.",
        "Absence from a gene catalogue is not proof of absence from a fully callable population dataset.",
        "Preserve per-platform and joint filters; do not add exome and genome counts if joint data are available.",
        "AlphaMissense is a pretrained computational score, not a patient diagnosis probability or functional assay.",
        "Protein and genomic AlphaMissense annotations must be joined with explicit transcript/isoform checks.",
    ]
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest = {path.name: digest(path) for path in sorted(args.output.iterdir())
                if path.is_file() and path.name != "manifest.json"}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"genes_ready": list(report["genes"]),
                      "errors": report["errors"], "mane": report.get("mane")}), flush=True)
    if report["errors"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
