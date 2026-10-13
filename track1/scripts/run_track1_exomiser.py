#!/usr/bin/env python3
"""Install public Exomiser data or run a fully local independent comparator.

Patient execution writes all JVM output and results privately. Stdout is only a
generic status. No lookup service, upload or clinical challenge key is used.
"""
from __future__ import annotations

if __package__:
    from .local_inputs import input_path
else:
    from local_inputs import input_path
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import contextlib
import csv
import datetime
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import struct
import stat
import subprocess
import time
import uuid
import zipfile
import zlib

ROOT = Path(__file__).resolve().parents[1]
SOFTWARE_VERSION = '15.1.0'
DATA_VERSION = '2602'
BASE_URL = 'https://g-879a9f.f5dc97.75bc.dn.glob.us/data/'
MIN_FREE = 6 * 1024**3


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def require_local_reference(path):
    """Reject macOS cloud placeholders before a long blocking hydration read."""
    info = Path(path).stat()
    dataless_flag = getattr(stat, 'SF_DATALESS', 0x40000000)
    if getattr(info, 'st_flags', 0) & dataless_flag:
        raise RuntimeError('Public reference is a cloud placeholder; restore it to local storage before analysis.')


def payload_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def input_hashes(vcf, phenotype, preset):
    return {name: sha256(path) for name, path in
            [('vcf', vcf), ('phenotype', phenotype), ('preset', preset)]}


def code_hashes(references):
    software = references / f'exomiser-cli-{SOFTWARE_VERSION}'
    return {'wrapper': sha256(Path(__file__)),
            'phenotype_extractor': sha256(ROOT / 'scripts/track1_baseline.py'),
            'exomiser_cli_jar': sha256(software / f'exomiser-cli-{SOFTWARE_VERSION}.jar'),
            'exomiser_libraries': {path.name: sha256(path) for path in sorted((software / 'lib').glob('*.jar'))},
            'java_executable': sha256(find_java(references))}


def local_path(value):
    if '://' in str(value):
        raise ValueError('Only local filesystem input paths are supported.')
    return Path(value).resolve()


class RemoteZip(io.RawIOBase):
    """Read only public ZIP metadata through byte ranges, never patient queries."""
    def __init__(self, url):
        import httpx
        if not url.startswith(BASE_URL):
            raise ValueError('Only the fixed public Exomiser archive host is allowed.')
        self.client = httpx.Client(timeout=120, follow_redirects=True,
                                   limits=httpx.Limits(max_connections=80, max_keepalive_connections=80))
        self.url = url
        response = self.client.head(url)
        response.raise_for_status()
        self.size = int(response.headers['content-length'])
        self.pos = 0

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = offset if whence == 0 else self.pos + offset if whence == 1 else self.size + offset
        return self.pos

    def read(self, size=-1):
        if size < 0:
            size = self.size - self.pos
        if size == 0:
            return b''
        response = self.client.get(self.url, headers={'Range': f'bytes={self.pos}-{self.pos + size - 1}'})
        response.raise_for_status()
        if response.status_code != 206:
            raise RuntimeError('Public archive does not support required byte ranges.')
        result = response.content
        if len(result) != size:
            raise RuntimeError('Public archive byte-range length mismatch.')
        self.pos += len(result)
        return result


def safe_member(root, name):
    rel = PurePosixPath(name)
    if rel.is_absolute() or '..' in rel.parts:
        raise ValueError('Unsafe archive member path.')
    result = (root / Path(*rel.parts)).resolve()
    if not result.is_relative_to(root.resolve()):
        raise ValueError('Archive member escapes destination.')
    return result


def verify_archive_data(destination, archive, required):
    """Verify locally acquired complete members before any database is opened."""
    manifest_path = destination / (archive + '.manifest.json')
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('url') != BASE_URL + archive:
        raise ValueError('Public reference manifest has unexpected source URL.')
    members = manifest.get('members', {})
    if not set(required).issubset(members):
        raise ValueError('Required public reference manifest entries are missing.')
    verified = {}
    for name, receipt in members.items():
        path = safe_member(destination, name)
        if not path.is_file() or path.stat().st_size != receipt.get('bytes'):
            raise ValueError('Public reference size verification failed.')
        require_local_reference(path)
        digest = sha256(path)
        if digest != receipt.get('sha256'):
            raise ValueError('Public reference SHA-256 verification failed.')
        verified[name] = {'bytes': path.stat().st_size, 'sha256': digest}
    return {'manifest_sha256': sha256(manifest_path), 'members': verified}


def verify_data(references, benchmark=False):
    destination = references / 'data'
    genome = 'test_hg19' if benchmark else DATA_VERSION + '_hg38'
    phenotype = 'test_phenotype' if benchmark else DATA_VERSION + '_phenotype'
    genome_required = [f'{genome}/{genome}_{suffix}' for suffix in
                       ('variants.mv.db', 'genome.mv.db', 'clinvar.mv.db', 'transcripts_ensembl.ser')]
    phenotype_required = [f'{phenotype}/{phenotype}.mv.db', f'{phenotype}/hp.obo',
                          f'{phenotype}/rw_string_10.mv']
    archives = {'exomiser-test-data.zip': genome_required + phenotype_required} if benchmark else {
        f'{DATA_VERSION}_hg38.zip': genome_required,
        f'{DATA_VERSION}_phenotype.zip': phenotype_required,
    }
    result = {'verified_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'complete': True, 'verification': 'exact byte size and SHA-256 for all manifest members',
              'archives': {name: verify_archive_data(destination, name, required)
                           for name, required in archives.items()}}
    dump(destination / ('benchmark_verification.json' if benchmark else 'production_verification.json'), result)
    return result


def download_member(remote, start, length, target, workers=64, chunk_size=4 * 1024**2):
    """Resume parallel public byte ranges; persist receipts before decompression.

    Keeping the compressed member makes an interrupted deflate stream resumable.
    Only one member is retained, then removed after verified extraction.
    """
    compressed_path = target.with_name(target.name + '.compressed.part')
    receipt_path = compressed_path.with_suffix('.receipt.json')
    identity = {'url': remote.url, 'start': start, 'bytes': length, 'chunk_size': chunk_size}
    receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else {**identity, 'chunks': {}}
    if any(receipt.get(key) != value for key, value in identity.items()):
        raise RuntimeError('Public download receipt does not match requested archive member.')
    descriptor = os.open(compressed_path, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        total_chunks = (length + chunk_size - 1) // chunk_size
        pending = []
        for index in range(total_chunks):
            size = min(chunk_size, length - index * chunk_size)
            expected = receipt['chunks'].get(str(index))
            if expected:
                block = os.pread(descriptor, size, index * chunk_size)
                if len(block) == size and hashlib.sha256(block).hexdigest() == expected:
                    continue
            pending.append(index)

        def download(index):
            size = min(chunk_size, length - index * chunk_size)
            begin = start + index * chunk_size
            for attempt in range(4):
                try:
                    response = remote.client.get(remote.url, headers={'Range': f'bytes={begin}-{begin + size - 1}'})
                    response.raise_for_status()
                    data = response.content
                    if response.status_code != 206 or len(data) != size:
                        raise RuntimeError('Public reference chunk length mismatch.')
                    return index, data
                except Exception:
                    if attempt == 3:
                        raise
                    time.sleep(1 + attempt)

        begin_time = time.monotonic()
        received, last_message = 0, begin_time
        # At most two batches are live, so completed futures do not retain GBs.
        with ThreadPoolExecutor(max_workers=workers) as pool:
            iterator = iter(pending)
            running = {pool.submit(download, index) for index in [next(iterator, None) for _ in range(workers)] if index is not None}
            while running:
                future = next(as_completed(running))
                running.remove(future)
                index, block = future.result()
                if shutil.disk_usage(target.parent).free - len(block) < MIN_FREE:
                    raise RuntimeError('Public download stopped at 6 GiB free-space guard.')
                count = os.pwrite(descriptor, block, index * chunk_size)
                if count != len(block):
                    raise RuntimeError('Incomplete local public-reference write.')
                receipt['chunks'][str(index)] = hashlib.sha256(block).hexdigest()
                temporary = receipt_path.with_suffix('.tmp')
                dump(temporary, receipt)
                temporary.replace(receipt_path)
                received += len(block)
                if time.monotonic() - last_message >= 60:
                    speed = received / 1024**2 / (time.monotonic() - begin_time)
                    done = len(receipt['chunks']) / total_chunks
                    print(f'Public member {target.name}: {done:.1%}, {speed:.2f} MiB/s this session', flush=True)
                    last_message = time.monotonic()
                next_index = next(iterator, None)
                if next_index is not None:
                    running.add(pool.submit(download, next_index))
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return compressed_path, receipt_path


def stream_archive(name, destination, skip_phenix=True, only_prefix=None):
    """Resume downloads and CRC-check extraction with a six GiB disk reserve."""
    destination.mkdir(parents=True, exist_ok=True)
    remote = RemoteZip(BASE_URL + name)
    manifest_path = destination / (name + '.manifest.json')
    previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else {'members': {}}
    manifest = {'url': remote.url, 'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'archive_bytes': remote.size, 'members': previous['members'],
                'skipped_optional_phenix': skip_phenix}
    with zipfile.ZipFile(remote) as archive:
        entries = archive.infolist()
    for entry in entries:
        if (entry.is_dir() or (skip_phenix and '/phenix/' in entry.filename)
                or (only_prefix and not entry.filename.startswith(only_prefix))):
            continue
        target = safe_member(destination, entry.filename)
        old = manifest['members'].get(entry.filename)
        if (target.exists() and old and target.stat().st_size == entry.file_size
                and old.get('crc32') == entry.CRC and sha256(target) == old.get('sha256')):
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_name(target.name + '.part')
        # Incomplete decompression cannot be resumed; the compressed chunks can.
        part.unlink(missing_ok=True)
        compressed_path = target.with_name(target.name + '.compressed.part')
        allocated = min(entry.compress_size, compressed_path.stat().st_blocks * 512) if compressed_path.exists() else 0
        needed = entry.file_size + entry.compress_size - allocated
        if shutil.disk_usage(destination).free - needed < MIN_FREE:
            raise RuntimeError('Insufficient free disk for public reference member plus 6 GiB reserve.')
        remote.seek(entry.header_offset)
        header = remote.read(30)
        if header[:4] != b'PK\x03\x04':
            raise RuntimeError('Invalid ZIP local header.')
        filename_length, extra_length = struct.unpack_from('<HH', header, 26)
        start = entry.header_offset + 30 + filename_length + extra_length
        decompressor = zlib.decompressobj(-15) if entry.compress_type == zipfile.ZIP_DEFLATED else None
        if entry.compress_type not in (zipfile.ZIP_DEFLATED, zipfile.ZIP_STORED):
            raise RuntimeError('Unsupported archive compression.')
        crc, written, compressed = 0, 0, 0
        digest = hashlib.sha256()
        print('Downloading public reference member: ' + entry.filename, flush=True)
        compressed_path, receipt_path = download_member(remote, start, entry.compress_size, target)
        with compressed_path.open('rb') as source:
            with part.open('wb') as handle:
                while block := source.read(4 * 1024**2):
                    compressed += len(block)
                    data = decompressor.decompress(block) if decompressor else block
                    handle.write(data)
                    digest.update(data)
                    crc = zlib.crc32(data, crc)
                    written += len(data)
                    if shutil.disk_usage(destination).free < MIN_FREE:
                        raise RuntimeError('Stopped public reference download at 6 GiB free-space guard.')
                if decompressor:
                    data = decompressor.flush()
                    handle.write(data)
                    digest.update(data)
                    crc = zlib.crc32(data, crc)
                    written += len(data)
        if compressed != entry.compress_size or written != entry.file_size or crc != entry.CRC:
            raise RuntimeError('Public reference size or CRC verification failed.')
        part.replace(target)
        manifest['members'][entry.filename] = {'bytes': written, 'crc32': crc, 'sha256': digest.hexdigest()}
        dump(manifest_path, manifest)
        compressed_path.unlink()
        receipt_path.unlink()
        print('Verified public reference member: ' + entry.filename, flush=True)


def build_analysis(vcf, hpo_ids, output, preset, assembly='hg38', proband=None):
    """Use the vendor's exome model without a candidate gene whitelist."""
    import yaml
    if not hpo_ids or any(not re.fullmatch(r'HP:[0-9]{7}', value) for value in hpo_ids):
        raise ValueError('Expected nonempty valid proband HPO identifiers.')
    analysis = yaml.safe_load(Path(preset).read_text())
    analysis.update(genomeAssembly=assembly, vcf=str(local_path(vcf)),
                    hpoIds=sorted(set(hpo_ids)))
    if proband:
        analysis['proband'] = proband
    return {'analysis': analysis, 'outputOptions': {
        'outputContributingVariantsOnly': False, 'numGenes': 0,
        'outputDirectory': str(local_path(output)), 'outputFileName': 'independent_exomiser',
        'outputFormats': ['JSON', 'TSV_GENE', 'TSV_VARIANT', 'VCF'],
    }}


def find_java(root):
    matches = sorted(root.glob('jdk-*/Contents/Home/bin/java'))
    if not matches:
        raise FileNotFoundError('Portable Java runtime missing.')
    return matches[-1]


@contextlib.contextmanager
def private_output(path):
    """Capture Python, native-library and subprocess output at OS descriptors."""
    with Path(path).open('a') as handle:
        saved = os.dup(1), os.dup(2)
        try:
            os.dup2(handle.fileno(), 1)
            os.dup2(handle.fileno(), 2)
            with contextlib.redirect_stdout(handle), contextlib.redirect_stderr(handle):
                yield handle
        finally:
            handle.flush()
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            os.close(saved[0])
            os.close(saved[1])


def normalize_allele(chrom, pos, ref, alt):
    chrom = chrom if chrom.startswith('chr') else 'chr' + chrom
    chrom = 'chrM' if chrom == 'chrMT' else chrom
    pos, ref, alt = int(pos), ref.upper(), alt.upper()
    if not ref or not alt or not re.fullmatch('[ACGTN]+', ref + alt):
        raise ValueError('Explicit sequence alleles required for comparator join.')
    # Minimal representation only. Repeat indels still need FASTA left alignment.
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    return f'{chrom}:{pos}:{ref}:{alt}'


def export_hypotheses(gene_tsv, variant_tsv, output):
    """Join vendor-contributing variants by gene and mode, entirely locally."""
    with Path(gene_tsv).open() as handle:
        genes = {(row['GENE_SYMBOL'], row['MOI']): row for row in csv.DictReader(handle, delimiter='\t')}
    grouped = {}
    skipped_nonsequence = 0
    with Path(variant_tsv).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['CONTRIBUTING_VARIANT'] != '1' or row['FILTER'] != 'PASS':
                continue
            try:
                allele = normalize_allele(row['CONTIG'], row['START'], row['REF'], row['ALT'])
            except ValueError:
                skipped_nonsequence += 1
                continue
            grouped.setdefault((row['GENE_SYMBOL'], row['MOI']), set()).add(allele)
    hypotheses = []
    for key, alleles in grouped.items():
        if key not in genes:
            raise ValueError('Contributing variant gene/mode missing from gene output.')
        gene = genes[key]
        hypotheses.append({'rank': int(gene['#RANK']), 'gene': key[0], 'inheritance': key[1],
                           'score': float(gene['EXOMISER_GENE_COMBINED_SCORE']),
                           'phenotype_score': float(gene['EXOMISER_GENE_PHENO_SCORE']),
                           'variant_ids': sorted(alleles),
                           'phase': 'not_verified_by_comparator_export',
                           'challenge_row_compatible': len(alleles) <= 2,
                           'score_is_calibrated_probability': False})
    hypotheses.sort(key=lambda item: (item['rank'], item['gene'], item['inheritance']))
    dump(output, {'hypotheses': hypotheses, 'skipped_nonsequence_rows': skipped_nonsequence,
                  'normalization': 'chr prefix and minimal allele trimming; repeat indels not left-aligned',
                  'phase_confirmation': False})
    return bool(hypotheses)


def properties_text(references, assembly='hg38', genome_version=DATA_VERSION,
                    phenotype_version=DATA_VERSION):
    return '\n'.join([
        'exomiser.data-directory=' + str(references / 'data'),
        f'exomiser.{assembly}.data-version=' + genome_version,
        f'exomiser.{assembly}.transcript-source=ensembl',
        'exomiser.phenotype.data-version=' + phenotype_version,
        'spring.cache.type=none', 'logging.level.com.zaxxer.hikari=ERROR',
        f'exomiser.{assembly}.use-clinvar-white-list=false',
    ]) + '\n'


def java_command(references, properties, analysis, heap='12g'):
    software = references / f'exomiser-cli-{SOFTWARE_VERSION}'
    return [str(find_java(references)), '-Xms2g', '-Xmx' + heap, '-XX:+UseG1GC',
            '-Dspring.config.location=' + properties.as_uri(),
            '-jar', str(software / f'exomiser-cli-{SOFTWARE_VERSION}.jar'),
            'analyse', '--analysis', str(analysis)]


def run_comparator(args):
    """No patient-derived material is returned from this function or printed."""
    import yaml
    args.output.mkdir(parents=True, exist_ok=True)
    log = args.output / 'private_exomiser.log'
    status_path = args.output / 'comparator_status.json'
    started = time.monotonic()
    status = {'software': 'Exomiser', 'version': SOFTWARE_VERSION, 'data_version': DATA_VERSION,
              'provenance_schema_version': 1, 'run_id': uuid.uuid4().hex,
              'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'completed': False, 'independent_implementation': True,
              'independent_of_prior_case_knowledge': False,
              'clinical_accuracy_established': False, 'official_challenge_score_obtained': False,
              'scope': 'All input variants evaluated by vendor exome preset; no fixed gene panel. Noncoding classes filtered by preset.',
              'submitted': False}
    # Invalidate prior completion before doing anything that can fail or be interrupted.
    dump(status_path, status)
    with private_output(log) as handle:
        try:
            references = verify_data(args.references)
            status['reference_size_and_hash_verified'] = True
            status['reference_provenance'] = references['archives']
            status['reference_provenance_sha256'] = payload_sha256(references['archives'])
            software = args.references / f'exomiser-cli-{SOFTWARE_VERSION}'
            preset = software / 'examples/preset-exome-analysis.yml'
            phenotype_input = args.phenotypes or args.phenotype_docx
            before = input_hashes(args.vcf, phenotype_input, preset)
            status['analysis_input_hashes'] = before
            status['phenotype_input_kind'] = 'json' if args.phenotypes else 'docx'
            status['code_provenance'] = code_hashes(args.references)
            from track1_baseline import extract_phenotypes
            if args.phenotypes:
                terms = json.loads(args.phenotypes.read_text())
            else:
                terms = extract_phenotypes(args.phenotype_docx)
            hpo = [entry['hpo_id'] for entry in terms
                   if entry.get('subject') == 'proband' and entry.get('present', True)]
            import pysam
            with pysam.VariantFile(str(args.vcf)) as variants:
                samples = list(variants.header.samples)
            if len(samples) != 1:
                raise ValueError('Comparator requires one explicitly selected sample.')
            # Fresh output directory prevents an old successful TSV from satisfying this run.
            run_output = args.output / 'runs' / status['run_id']
            run_output.mkdir(parents=True)
            analysis = build_analysis(args.vcf, hpo, run_output, preset, proband=samples[0])
            status['hpo_payload_sha256'] = payload_sha256(analysis['analysis']['hpoIds'])
            status['analysis_payload_sha256'] = payload_sha256(analysis)
            yaml_path = run_output / 'private_analysis.yml'
            yaml_path.write_text(yaml.safe_dump(analysis, sort_keys=False))
            properties = run_output / 'application.properties'
            properties.write_text(properties_text(args.references))
            status['configuration_hashes'] = {'analysis_yaml': sha256(yaml_path),
                                              'properties': sha256(properties)}
            dump(status_path, status)
            command = java_command(args.references, properties, yaml_path)
            process = subprocess.run(command, cwd=run_output, stdout=handle, stderr=handle)
            status['process_exit_code'] = process.returncode
            if input_hashes(args.vcf, phenotype_input, preset) != before:
                raise ValueError('Comparator inputs changed during execution.')
            status['inputs_unchanged_during_run'] = True
            gene_files = list(run_output.glob('independent_exomiser*.genes.tsv'))
            variant_files = list(run_output.glob('independent_exomiser*.variants.tsv'))
            status['completed'] = process.returncode == 0 and len(gene_files) == len(variant_files) == 1
            if status['completed']:
                hypotheses_path = args.output / 'comparator_hypotheses.json'
                status['hypotheses_exported'] = export_hypotheses(gene_files[0], variant_files[0],
                                                                hypotheses_path)
                status['output_hashes'] = {'hypotheses': sha256(hypotheses_path),
                                           'genes_tsv': sha256(gene_files[0]),
                                           'variants_tsv': sha256(variant_files[0])}
            if not status['completed']:
                status['failure'] = 'Local comparator did not complete; inspect private log locally.'
        except Exception:
            import traceback
            traceback.print_exc()
            status['completed'] = False
            status['failure'] = 'Local setup or execution failed; inspect private log locally.'
    status['runtime_seconds'] = round(time.monotonic() - started, 3)
    dump(status_path, status)
    print('Local Exomiser comparator completed.' if status['completed'] else 'Local Exomiser comparator failed; details remain in private local files.')
    return 0 if status['completed'] else 1


def run_benchmark(args):
    """Run the vendor's public synthetic Pfeiffer example, never the child data."""
    import yaml
    verify_data(args.references, benchmark=True)
    output = ROOT / 'research/track1_exomiser_benchmark'
    output.mkdir(parents=True, exist_ok=True)
    software = args.references / f'exomiser-cli-{SOFTWARE_VERSION}'
    examples = software / 'examples'
    source = yaml.safe_load((examples / 'test-analysis-exome.yml').read_text())
    analysis = build_analysis(examples / 'Pfeiffer.vcf', source['analysis']['hpoIds'], output,
                              examples / 'preset-exome-analysis-human-only.yml', assembly='hg19')
    job = output / 'benchmark_analysis.yml'
    job.write_text(yaml.safe_dump(analysis, sort_keys=False))
    properties = output / 'application.properties'
    properties.write_text(properties_text(args.references, 'hg19', 'test', 'test'))
    started = time.monotonic()
    with (output / 'benchmark.log').open('w') as handle:
        process = subprocess.run(java_command(args.references, properties, job, heap='6g'),
                                 cwd=output, stdout=handle, stderr=handle)
    genes = list(output.glob('independent_exomiser*.genes.tsv'))
    rank = None
    if process.returncode == 0 and len(genes) == 1:
        with genes[0].open() as handle:
            ranks = [int(row['#RANK']) for row in csv.DictReader(handle, delimiter='\t')
                     if row['GENE_SYMBOL'] == 'FGFR2']
        rank = min(ranks) if ranks else None
    result = {'benchmark': 'Official Exomiser Pfeiffer public demonstration',
              'dataset_type': 'Synthetic spike-in demonstration, not independent clinical accuracy cohort',
              'completed': process.returncode == 0 and bool(genes), 'exit_code': process.returncode,
              'expected_public_gene': 'FGFR2', 'expected_gene_rank': rank,
              'expected_gene_rank_one': rank == 1, 'clinical_accuracy_established': False,
              'runtime_seconds': round(time.monotonic() - started, 3),
              'limitations': ['Test genomic database contains only demonstration variants',
                             'Human-only phenotype database; no cross-species validation',
                             'One supplied positive example does not measure specificity or generalization'],
              'input_sha256': sha256(examples / 'Pfeiffer.vcf')}
    dump(output / 'benchmark_status.json', result)
    print('Public Exomiser demonstration completed.' if result['completed'] else 'Public Exomiser demonstration failed; see public benchmark log.')
    return 0 if result['completed'] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-data', action='store_true', help='Public downloads only; no patient processing.')
    parser.add_argument('--prepare-benchmark-data', action='store_true', help='Only official public demonstration hg19 reference members.')
    parser.add_argument('--benchmark', action='store_true', help='Run only the public demonstration, not patient data.')
    parser.add_argument('--verify-data', action='store_true', help='Verify public production member sizes and SHA-256 hashes offline.')
    parser.add_argument('--references', type=local_path, default=ROOT / 'references/exomiser')
    parser.add_argument('--vcf', type=local_path, default=input_path('vcf'))
    parser.add_argument('--phenotypes', type=local_path)
    parser.add_argument('--phenotype-docx', type=local_path, default=input_path('phenotype'))
    parser.add_argument('--output', type=local_path, default=ROOT / 'work/private_analysis.nosync/track1_exomiser')
    args = parser.parse_args()
    if args.prepare_benchmark_data:
        stream_archive('exomiser-test-data.zip', args.references / 'data')
        return 0
    if args.verify_data:
        verify_data(args.references)
        print('Public Exomiser production reference sizes and hashes verified.')
        return 0
    if args.benchmark:
        return run_benchmark(args)
    if args.prepare_data:
        # Sequential archives keep disk reservations straightforward.
        for archive in [f'{DATA_VERSION}_hg38.zip', f'{DATA_VERSION}_phenotype.zip']:
            stream_archive(archive, args.references / 'data')
        verify_data(args.references)
        print('Public Exomiser production references prepared.')
        return 0
    return run_comparator(args)


if __name__ == '__main__':
    raise SystemExit(main())
