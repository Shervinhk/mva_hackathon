#!/usr/bin/env python3
"""Restore one fixed public Exomiser DB without retaining a full ZIP beside it.

Downloads resume from independently SHA-256-checked compressed chunk files.
Inflation does NOT resume: a subsequent run discards its incomplete output and
downloads any chunks consumed by the previous attempt before starting again.
The existing acquisition manifest is never rewritten.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import contextlib
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import time
import zipfile
import zlib

if __package__:
    from .run_track1_exomiser import BASE_URL, MIN_FREE, ROOT, RemoteZip, require_local_reference, safe_member
else:
    from run_track1_exomiser import BASE_URL, MIN_FREE, ROOT, RemoteZip, require_local_reference, safe_member

ARCHIVE = '2602_hg38.zip'
MEMBER = '2602_hg38/2602_hg38_variants.mv.db'
CACHE_NAME = 'local_cache.nosync/2602_hg38_variants.mv.db'
CHUNK_SIZE = 4 * 1024**2
MAX_OUTPUT = 16 * 1024**2


def fsync_directory(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def regular_file_or_absent(path):
    if path.is_symlink() or (path.exists() and not path.is_file()):
        raise RuntimeError('Unsafe public cache work-file path.')
    if path.exists():
        require_local_reference(path)


def atomic_json(path, value):
    regular_file_or_absent(path)
    temporary = path.with_name(path.name + '.tmp')
    regular_file_or_absent(temporary)
    with temporary.open('w') as handle:
        json.dump(value, handle, sort_keys=True)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)
    fsync_directory(path.parent)


def check_free(path, additional=0, min_free=MIN_FREE):
    if shutil.disk_usage(path).free - additional < min_free:
        raise RuntimeError('Public cache restoration stopped at the 6 GiB free-space guard.')


def exact_range(remote, start, length):
    """Validate absolute archive ranges before accepting any response body."""
    if start < 0 or length <= 0 or start + length > remote.size:
        raise RuntimeError('Invalid public archive byte range.')
    end = start + length - 1
    headers = {'Range': f'bytes={start}-{end}', 'Accept-Encoding': 'identity'}
    with remote.client.stream('GET', remote.url, headers=headers) as response:
        response.raise_for_status()
        if response.status_code != 206:
            raise RuntimeError('Public archive did not return a partial response.')
        # The official endpoint reports an unknown total ("*") on valid ranges.
        # Its HEAD size is checked against the acquisition manifest by restore().
        match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+|\*)', response.headers.get('content-range', ''))
        if (not match or tuple(map(int, match.groups()[:2])) != (start, end)
                or (match.group(3) != '*' and int(match.group(3)) != remote.size)):
            raise RuntimeError('Public archive Content-Range mismatch.')
        if response.headers.get('content-encoding', 'identity').lower() != 'identity':
            raise RuntimeError('Public archive range has unexpected content encoding.')
        declared = response.headers.get('content-length')
        if declared is not None and (not declared.isdigit() or int(declared) != length):
            raise RuntimeError('Public archive range Content-Length mismatch.')
        result = bytearray()
        for block in response.iter_raw(chunk_size=CHUNK_SIZE):
            if len(result) + len(block) > length:
                raise RuntimeError('Public archive range body exceeds requested length.')
            result.extend(block)
        if len(result) != length:
            raise RuntimeError('Public archive range body is incomplete.')
        return bytes(result)


class ExactRemoteZip(RemoteZip):
    """Reuse the fixed-host ZIP reader with strict checks on every HTTP range."""

    def read(self, size=-1):
        if size < 0:
            size = self.size - self.pos
        if size == 0:
            return b''
        result = exact_range(self, self.pos, size)
        self.pos += len(result)
        return result

    def close(self):
        if not self.closed and hasattr(self, 'client'):
            self.client.close()
        super().close()


def read_manifest(data_root):
    path = safe_member(data_root, ARCHIVE + '.manifest.json')
    require_local_reference(path)
    manifest = json.loads(path.read_text())
    if manifest.get('url') != BASE_URL + ARCHIVE:
        raise RuntimeError('Unexpected public acquisition manifest source.')
    expected = manifest.get('members', {}).get(MEMBER, {})
    if (type(expected.get('bytes')) is not int or expected['bytes'] <= 0
            or type(expected.get('crc32')) is not int or not 0 <= expected['crc32'] <= 0xffffffff
            or not re.fullmatch('[0-9a-f]{64}', expected.get('sha256', ''))
            or type(manifest.get('archive_bytes')) is not int or manifest['archive_bytes'] <= 0):
        raise RuntimeError('Incomplete public acquisition manifest.')
    return manifest, expected


def zip_member(remote, expected):
    with zipfile.ZipFile(remote) as archive:
        matches = [entry for entry in archive.infolist() if entry.filename == MEMBER]
        if len(matches) != 1:
            raise RuntimeError('Public archive must contain exactly one requested member.')
        entry = matches[0]
        next_header = min((item.header_offset for item in archive.infolist()
                           if item.header_offset > entry.header_offset), default=archive.start_dir)
        boundary = min(next_header, archive.start_dir)
    if (entry.is_dir() or entry.compress_type not in (zipfile.ZIP_DEFLATED, zipfile.ZIP_STORED)
            or entry.flag_bits & (1 | 32 | 64)
            or entry.file_size != expected['bytes'] or entry.CRC != expected['crc32']):
        raise RuntimeError('Public ZIP metadata does not match the acquisition manifest.')
    header = exact_range(remote, entry.header_offset, 30)
    fields = struct.unpack('<4s5H3I2H', header)
    signature, _, flags, compression, _, _, crc, compressed, size, name_len, extra_len = fields
    if signature != b'PK\x03\x04' or flags != entry.flag_bits or compression != entry.compress_type:
        raise RuntimeError('Public ZIP local header mismatch.')
    variable = exact_range(remote, entry.header_offset + 30, name_len + extra_len)
    name = variable[:name_len].decode('utf-8' if flags & 0x800 else 'cp437')
    if name != entry.filename:
        raise RuntimeError('Public ZIP local filename mismatch.')
    if not flags & 8:
        extra, position = variable[name_len:], 0
        while position + 4 <= len(extra):
            kind, length = struct.unpack_from('<HH', extra, position)
            position += 4
            value = extra[position:position + length]
            if len(value) != length:
                raise RuntimeError('Invalid public ZIP extra field.')
            position += length
            if kind == 1:
                offset = 0
                if size == 0xffffffff:
                    if len(value) < offset + 8:
                        raise RuntimeError('Incomplete public ZIP64 size.')
                    size = struct.unpack_from('<Q', value, offset)[0]
                    offset += 8
                if compressed == 0xffffffff:
                    if len(value) < offset + 8:
                        raise RuntimeError('Incomplete public ZIP64 compressed size.')
                    compressed = struct.unpack_from('<Q', value, offset)[0]
        if (crc, compressed, size) != (entry.CRC, entry.compress_size, entry.file_size):
            raise RuntimeError('Public ZIP local sizes or CRC mismatch.')
    start = entry.header_offset + 30 + name_len + extra_len
    if entry.compress_size <= 0 or start + entry.compress_size > boundary:
        raise RuntimeError('Public ZIP compressed member range is invalid.')
    return entry, start


def member_identity(remote, entry, start, expected, chunk_size=CHUNK_SIZE):
    return {'format': 1, 'url': remote.url, 'archive_bytes': remote.size,
            'member': entry.filename, 'start': start, 'compressed_bytes': entry.compress_size,
            'compression': entry.compress_type, 'bytes': entry.file_size, 'crc32': entry.CRC,
            'sha256': expected['sha256'], 'chunk_size': chunk_size}


class ChunkStore:
    def __init__(self, directory, identity):
        self.directory, self.identity = directory, identity
        if directory.is_symlink():
            raise RuntimeError('Public chunk directory must not be a symlink.')
        directory.mkdir(parents=True, exist_ok=True)
        identity_path = directory / 'identity.json'
        regular_file_or_absent(identity_path)
        if identity_path.exists():
            if json.loads(identity_path.read_text()) != identity:
                raise RuntimeError('Public chunk receipt belongs to a different archive member.')
        else:
            # Refuse to adopt unrelated preexisting chunk files without provenance.
            if any(path.name != 'identity.json.tmp' for path in directory.iterdir()):
                raise RuntimeError('Public chunk directory has no identity receipt.')
            atomic_json(identity_path, identity)
        self.chunk_size = identity['chunk_size']
        self.length = identity['compressed_bytes']
        self.count = (self.length + self.chunk_size - 1) // self.chunk_size

    def path(self, index):
        return self.directory / f'{index:08d}.chunk'

    def receipt(self, index):
        return self.directory / f'{index:08d}.json'

    def size(self, index):
        return min(self.chunk_size, self.length - index * self.chunk_size)

    def read_verified(self, index):
        path, receipt = self.path(index), self.receipt(index)
        regular_file_or_absent(path)
        regular_file_or_absent(receipt)
        if not path.exists() or not receipt.exists() or path.stat().st_size != self.size(index):
            return None
        try:
            value = json.loads(receipt.read_text())
        except (ValueError, UnicodeError):
            return None
        if not isinstance(value, dict) or value.get('bytes') != self.size(index) or value.get('index') != index:
            return None
        block = path.read_bytes()
        if hashlib.sha256(block).hexdigest() != value.get('sha256'):
            return None
        return block

    def write(self, index, block, min_free=MIN_FREE):
        if len(block) != self.size(index):
            raise RuntimeError('Public compressed chunk size mismatch.')
        path = self.path(index)
        temporary = path.with_suffix('.part')
        regular_file_or_absent(path)
        regular_file_or_absent(temporary)
        check_free(self.directory, len(block), min_free)
        with temporary.open('wb') as handle:
            handle.write(block)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
        fsync_directory(self.directory)
        atomic_json(self.receipt(index), {'index': index, 'bytes': len(block),
                                         'sha256': hashlib.sha256(block).hexdigest()})

    def pending(self):
        pending, retained = [], 0
        for index in range(self.count):
            if self.read_verified(index) is None:
                pending.append(index)
            else:
                retained += self.size(index)
        return pending, retained


def download_chunks(remote, store, workers=32, min_free=MIN_FREE, pending=None):
    if not 1 <= workers <= 64:
        raise ValueError('Public download workers must be between 1 and 64.')
    if pending is None:
        pending, _ = store.pending()

    def fetch(index):
        for attempt in range(4):
            try:
                return index, exact_range(remote, store.identity['start'] + index * store.chunk_size,
                                          store.size(index))
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(attempt + 1)

    done = store.count - len(pending)
    last_report = time.monotonic()
    iterator = iter(pending)
    # Completed futures retain at most one worker batch of 4 MiB blocks.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        running = {pool.submit(fetch, index) for index in
                   [next(iterator, None) for _ in range(workers)] if index is not None}
        try:
            while running:
                future = next(as_completed(running))
                running.remove(future)
                index, block = future.result()
                store.write(index, block, min_free)
                done += 1
                if time.monotonic() - last_report >= 30:
                    print(f'Public cache download: {done / store.count:.1%}', flush=True)
                    last_report = time.monotonic()
                following = next(iterator, None)
                if following is not None:
                    running.add(pool.submit(fetch, following))
                del block, future
        except BaseException:
            for future in running:
                future.cancel()
            raise


def discard_incomplete(part):
    regular_file_or_absent(part)
    if part.exists():
        part.unlink()
        fsync_directory(part.parent)
        print('Restarting public cache inflation; consumed chunks will be downloaded again.', flush=True)


def verify_file(path, expected):
    regular_file_or_absent(path)
    if not path.exists() or path.stat().st_size != expected['bytes']:
        raise RuntimeError('Public cached member size verification failed.')
    digest, crc = hashlib.sha256(), 0
    with path.open('rb') as handle:
        while block := handle.read(MAX_OUTPUT):
            digest.update(block)
            crc = zlib.crc32(block, crc)
    if crc != expected['crc32'] or digest.hexdigest() != expected['sha256']:
        raise RuntimeError('Public cached member SHA-256 or CRC verification failed.')


def inflate_chunks(store, part, target, expected, min_free=MIN_FREE):
    """Consume only durable chunks; always restart the deflate stream at zero."""
    regular_file_or_absent(part)
    regular_file_or_absent(target)
    if part.exists() or target.exists():
        raise RuntimeError('Public cache output already exists; refusing to overwrite it.')
    identity = store.identity
    if (identity['bytes'], identity['crc32'], identity['sha256']) != (
            expected['bytes'], expected['crc32'], expected['sha256']):
        raise RuntimeError('Public compressed chunk identity does not match the manifest.')
    compression = identity['compression']
    if compression not in (zipfile.ZIP_DEFLATED, zipfile.ZIP_STORED):
        raise RuntimeError('Unsupported public member compression.')
    inflater = zlib.decompressobj(-15) if compression == zipfile.ZIP_DEFLATED else None
    digest, crc, written, consumed = hashlib.sha256(), 0, 0, 0
    last_report = time.monotonic()
    with part.open('xb', buffering=0) as handle:
        def write_output(block):
            nonlocal crc, written
            if len(block) > MAX_OUTPUT or written + len(block) > expected['bytes']:
                raise RuntimeError('Public cache inflation exceeds the expected size.')
            if block:
                check_free(part.parent, len(block), min_free)
                view = memoryview(block)
                while view:
                    count = handle.write(view)
                    if not count:
                        raise RuntimeError('Incomplete public cache output write.')
                    view = view[count:]
                digest.update(block)
                crc = zlib.crc32(block, crc)
                written += len(block)

        for index in range(store.count):
            block = store.read_verified(index)
            if block is None:
                raise RuntimeError('Public compressed chunk is missing or has failed its SHA-256 receipt.')
            consumed += len(block)
            if inflater:
                pending = block
                while True:
                    data = inflater.decompress(pending, MAX_OUTPUT)
                    pending = inflater.unconsumed_tail
                    if inflater.unused_data:
                        raise RuntimeError('Public compressed member contains trailing data.')
                    write_output(data)
                    if pending or len(data) == MAX_OUTPUT:
                        continue
                    break
            else:
                write_output(block)
            # Output must reach durable storage BEFORE its source chunk is freed.
            os.fsync(handle.fileno())
            store.path(index).unlink()
            fsync_directory(store.directory)
            if time.monotonic() - last_report >= 30:
                print(f'Public cache inflation: {written / expected["bytes"]:.1%}', flush=True)
                last_report = time.monotonic()
        if (consumed != identity['compressed_bytes'] or written != expected['bytes']
                or crc != expected['crc32'] or digest.hexdigest() != expected['sha256']
                or (inflater is not None and not inflater.eof)):
            raise RuntimeError('Public cache SHA-256, ZIP CRC, size, or deflate-end verification failed.')
        os.fsync(handle.fileno())
    # Hard-link publication is atomic and cannot replace a concurrently created DB.
    os.link(part, target)
    fsync_directory(target.parent)
    part.unlink()
    fsync_directory(target.parent)


def cache_paths(data_root):
    data_root = data_root.resolve()
    target = safe_member(data_root, CACHE_NAME)
    raw_cache = data_root / 'local_cache.nosync'
    if raw_cache.is_symlink():
        raise RuntimeError('Public cache directory must not be a symlink.')
    # Inspect the lexical canonical path so a materialized file is never replaced.
    canonical = data_root / MEMBER
    if safe_member(data_root, MEMBER) != canonical.resolve():
        raise RuntimeError('Unsafe public canonical path.')
    if canonical.is_symlink():
        if canonical.resolve() != target:
            raise RuntimeError('Existing public reference symlink points to a different target.')
    elif canonical.exists():
        raise RuntimeError('Canonical public DB already exists; refusing to replace it, including cloud placeholders.')
    if target.is_symlink() or (data_root / CACHE_NAME).is_symlink():
        raise RuntimeError('Public cache DB must not be a symlink.')
    return canonical, target


def publish_symlink(data_root, canonical, target):
    if safe_member(data_root, CACHE_NAME) != target.resolve():
        raise RuntimeError('Public cache target escaped the data root.')
    canonical.parent.mkdir(parents=True, exist_ok=True)
    relative = os.path.relpath(target, canonical.parent)
    if not (canonical.parent / relative).resolve().is_relative_to(data_root.resolve()):
        raise RuntimeError('Public canonical symlink would escape the data root.')
    if canonical.is_symlink():
        if canonical.resolve() != target.resolve():
            raise RuntimeError('Existing public canonical symlink has a different target.')
    else:
        # symlink() refuses to overwrite any existing file, even a cloud placeholder.
        canonical.symlink_to(relative)
        fsync_directory(canonical.parent)
    if safe_member(data_root, MEMBER) != target.resolve():
        raise RuntimeError('Public canonical symlink verification failed.')


@contextlib.contextmanager
def cache_lock(target):
    path = target.with_name(target.name + '.restore.lock')
    regular_file_or_absent(path)
    with path.open('a') as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Public cache restoration is already running.') from exc
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def restore(data_root, workers=32, download_only=False):
    data_root = Path(data_root).resolve()
    manifest, expected = read_manifest(data_root)
    canonical, target = cache_paths(data_root)
    target.parent.mkdir(parents=True, exist_ok=True)
    with cache_lock(target):
        canonical, target = cache_paths(data_root)
        with ExactRemoteZip(BASE_URL + ARCHIVE) as remote:
            if remote.size != manifest['archive_bytes']:
                raise RuntimeError('Public archive byte size differs from the acquisition manifest.')
            entry, start = zip_member(remote, expected)
            if target.exists():
                verify_file(target, expected)
            else:
                store = ChunkStore(target.with_name(target.name + '.chunks'),
                                   member_identity(remote, entry, start, expected))
                part = target.with_name(target.name + '.part')
                discard_incomplete(part)
                pending, retained = store.pending()
                # This is a preflight estimate; every later write still enforces
                # the reserve. Budget filesystem allocation for receipts and
                # rounded chunk tails separately from retained payload bytes.
                allocation_unit = max(4096, os.statvfs(target.parent).f_frsize)
                metadata_allowance = (store.count + 4) * allocation_unit * 2
                required = (max(entry.compress_size, entry.file_size) - retained
                            + MAX_OUTPUT + CHUNK_SIZE + metadata_allowance)
                check_free(target.parent, max(0, required))
                print('Downloading fixed public Exomiser cache member.', flush=True)
                download_chunks(remote, store, workers, pending=pending)
                if download_only:
                    print('Public compressed chunks verified; inflation remains pending.', flush=True)
                    return
                print('Inflating public cache; verified source chunks are freed after durable writes.', flush=True)
                inflate_chunks(store, part, target, expected)
        publish_symlink(data_root, canonical, target)
        print('Public Exomiser variant cache verified and linked locally.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=ROOT / 'references/exomiser/data')
    parser.add_argument('--workers', type=int, default=32, choices=range(1, 65), metavar='1..64')
    parser.add_argument('--download-only', action='store_true')
    args = parser.parse_args()
    try:
        restore(args.data_dir, args.workers, args.download_only)
    except KeyboardInterrupt:
        parser.exit(130, 'Public cache restoration interrupted; rerun to restart safely.\n')
    except RuntimeError as exc:
        parser.exit(1, f'Public cache restoration failed: {exc}\n')
    except Exception as exc:
        # Do not print filesystem paths, remote response bodies, or tracebacks.
        parser.exit(1, f'Public cache restoration failed ({type(exc).__name__}).\n')


if __name__ == '__main__':
    main()
