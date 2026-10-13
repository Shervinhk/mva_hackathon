"""Synthetic public-cache tests: no remote requests or production DB reads."""
import hashlib
import io
import json
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest import mock
import zipfile
import zlib

from scripts import restore_exomiser_public_cache as cache
from scripts.run_track1_exomiser import safe_member, verify_archive_data


class Response:
    def __init__(self, client, start, end):
        self.client = client
        self.status_code = client.status
        self.body = client.payload[start:end + 1]
        self.headers = {'content-range': f'bytes {start}-{end}/{len(client.payload)}',
                        'content-length': str(end - start + 1)}
        self.headers.update(client.override_headers)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def raise_for_status(self):
        pass

    def iter_raw(self, chunk_size):
        self.client.bodies_read += 1
        body = self.body + self.client.extra_body
        for offset in range(0, len(body), chunk_size):
            yield body[offset:offset + chunk_size]


class Client:
    def __init__(self, payload):
        self.payload = payload
        self.calls = []
        self.status = 206
        self.override_headers = {}
        self.extra_body = b''
        self.bodies_read = 0

    def stream(self, method, url, headers):
        assert method == 'GET'
        assert url == cache.BASE_URL + cache.ARCHIVE
        assert headers['Accept-Encoding'] == 'identity'
        start, end = map(int, headers['Range'].removeprefix('bytes=').split('-'))
        self.calls.append((start, end))
        return Response(self, start, end)

    def close(self):
        pass


class Remote(cache.ExactRemoteZip):
    def __init__(self, payload):
        io.RawIOBase.__init__(self)
        self.client = Client(payload)
        self.url = cache.BASE_URL + cache.ARCHIVE
        self.size = len(payload)
        self.pos = 0


class PublicCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='synthetic-public-cache-')
        self.root = Path(self.temporary.name)
        self.payload = random.Random(2026).randbytes(24_000)

    def tearDown(self):
        self.temporary.cleanup()

    def fixture(self, payload=None, compression=zipfile.ZIP_DEFLATED, chunk_size=1024):
        payload = self.payload if payload is None else payload
        archive = io.BytesIO()
        with zipfile.ZipFile(archive, 'w', compression=compression) as handle:
            handle.writestr('preceding-public-member.txt', b'Public synthetic archive prefix')
            handle.writestr(cache.MEMBER, payload)
        expected = {'bytes': len(payload), 'crc32': zlib.crc32(payload),
                    'sha256': hashlib.sha256(payload).hexdigest()}
        remote = Remote(archive.getvalue())
        entry, start = cache.zip_member(remote, expected)
        identity = cache.member_identity(remote, entry, start, expected, chunk_size)
        store = cache.ChunkStore(self.root / 'chunks', identity)
        remote.client.calls.clear()
        return remote, store, expected

    def download(self, remote, store):
        cache.download_chunks(remote, store, workers=4, min_free=0)

    def test_chunk_resume_redownloads_missing_tampered_and_bad_receipt(self):
        remote, store, _ = self.fixture()
        self.download(remote, store)
        self.assertEqual(len(remote.client.calls), store.count)
        self.assertEqual(store.pending(), ([], store.length))
        remote.client.calls.clear()
        self.download(remote, store)
        self.assertEqual(remote.client.calls, [])
        # A valid chunk is retained; three independently damaged chunks are fetched.
        store.path(0).unlink()
        block = bytearray(store.path(1).read_bytes())
        block[0] ^= 1
        store.path(1).write_bytes(block)
        store.receipt(2).write_text('{invalid JSON')
        self.assertEqual(store.pending()[0], [0, 1, 2])
        self.download(remote, store)
        start = store.identity['start']
        self.assertGreater(start, 0)
        self.assertEqual(set(remote.client.calls), {
            (start + index * store.chunk_size, start + index * store.chunk_size + store.size(index) - 1)
            for index in range(3)})
        self.assertEqual(store.pending(), ([], store.length))

    def test_receipts_cannot_be_reused_for_different_archive_range(self):
        _, store, _ = self.fixture()
        identity = dict(store.identity, start=store.identity['start'] + 1)
        with self.assertRaisesRegex(RuntimeError, 'different archive member'):
            cache.ChunkStore(store.directory, identity)

    def test_absolute_content_range_and_status_checked_before_body(self):
        remote = Remote(b'0123456789abcdef')
        self.assertEqual(cache.exact_range(remote, 5, 4), b'5678')
        # The official server supplies exact offsets with its total reported as *.
        remote.client.override_headers = {'content-range': 'bytes 5-8/*'}
        self.assertEqual(cache.exact_range(remote, 5, 4), b'5678')
        invalid_headers = [
            {'content-range': 'bytes 0-3/16'},  # same length, wrong absolute offset
            {'content-range': 'bytes 5-8/17'},
            {'content-range': 'bytes 0-3/*'},
            {'content-length': '5'},
            {'content-encoding': 'gzip'},
        ]
        for headers in invalid_headers:
            with self.subTest(headers=headers):
                remote.client.override_headers = headers
                before = remote.client.bodies_read
                with self.assertRaises(RuntimeError):
                    cache.exact_range(remote, 5, 4)
                self.assertEqual(remote.client.bodies_read, before)
        remote.client.override_headers = {}
        remote.client.status = 200
        with self.assertRaisesRegex(RuntimeError, 'partial response'):
            cache.exact_range(remote, 5, 4)
        remote.client.status = 206
        remote.client.extra_body = b'x'
        with self.assertRaisesRegex(RuntimeError, 'exceeds requested'):
            cache.exact_range(remote, 5, 4)

    def test_zip_metadata_requires_manifest_crc_and_size(self):
        remote, _, expected = self.fixture()
        for changed in (dict(expected, bytes=expected['bytes'] + 1),
                        dict(expected, crc32=expected['crc32'] ^ 1)):
            with self.assertRaisesRegex(RuntimeError, 'metadata'):
                cache.zip_member(remote, changed)

    def test_ordered_inflate_frees_chunks_only_after_output_fsync(self):
        remote, store, expected = self.fixture()
        self.download(remote, store)
        part, target = self.root / 'output.part', self.root / 'output.db'
        removed, durable_size = [], -1
        actual_fsync, actual_unlink = os.fsync, Path.unlink

        def observe_fsync(descriptor):
            nonlocal durable_size
            actual_fsync(descriptor)
            if part.exists() and os.fstat(descriptor).st_ino == part.stat().st_ino:
                durable_size = part.stat().st_size

        def observe_unlink(path, *args, **kwargs):
            if path.suffix == '.chunk':
                self.assertEqual(durable_size, part.stat().st_size)
                removed.append(int(path.stem))
            return actual_unlink(path, *args, **kwargs)

        with mock.patch.object(cache.os, 'fsync', side_effect=observe_fsync), \
                mock.patch.object(Path, 'unlink', autospec=True, side_effect=observe_unlink):
            cache.inflate_chunks(store, part, target, expected, min_free=0)
        self.assertEqual(removed, list(range(store.count)))
        self.assertEqual(target.read_bytes(), self.payload)
        self.assertFalse(part.exists())
        self.assertFalse(list(store.directory.glob('*.chunk')))
        cache.verify_file(target, expected)

    def test_interruption_requires_output_restart_and_consumed_chunk_redownload(self):
        remote, store, expected = self.fixture()
        self.download(remote, store)
        part, target = self.root / 'output.part', self.root / 'output.db'
        part.write_bytes(b'Interrupted deflate output is not a resumable checkpoint.')
        store.path(0).unlink()
        with self.assertRaisesRegex(RuntimeError, 'already exists'):
            cache.inflate_chunks(store, part, target, expected, min_free=0)
        cache.discard_incomplete(part)
        remote.client.calls.clear()
        self.download(remote, store)
        self.assertEqual(len(remote.client.calls), 1)
        cache.inflate_chunks(store, part, target, expected, min_free=0)
        self.assertEqual(target.read_bytes(), self.payload)

    def test_chunk_corruption_refused_before_consumption(self):
        remote, store, expected = self.fixture()
        self.download(remote, store)
        store.path(0).write_bytes(b'x' * store.size(0))
        part, target = self.root / 'output.part', self.root / 'output.db'
        with self.assertRaisesRegex(RuntimeError, 'SHA-256 receipt'):
            cache.inflate_chunks(store, part, target, expected, min_free=0)
        self.assertTrue(store.path(0).exists())
        self.assertFalse(target.exists())

    def test_corruption_with_updated_chunk_receipt_still_fails_zip_crc(self):
        remote, store, expected = self.fixture(compression=zipfile.ZIP_STORED)
        self.download(remote, store)
        block = bytearray(store.path(0).read_bytes())
        block[0] ^= 1
        store.write(0, bytes(block), min_free=0)
        target = self.root / 'output.db'
        with self.assertRaisesRegex(RuntimeError, 'ZIP CRC'):
            cache.inflate_chunks(store, self.root / 'output.part', target, expected, min_free=0)
        self.assertFalse(target.exists())

    def test_manifest_sha_independently_required_before_promotion(self):
        remote, store, expected = self.fixture()
        self.download(remote, store)
        expected = dict(expected, sha256='0' * 64)
        store.identity = dict(store.identity, sha256=expected['sha256'])
        target = self.root / 'output.db'
        with self.assertRaisesRegex(RuntimeError, 'SHA-256'):
            cache.inflate_chunks(store, self.root / 'output.part', target, expected, min_free=0)
        self.assertFalse(target.exists())

    def test_truncated_deflate_stream_cannot_be_promoted(self):
        remote, store, expected = self.fixture(payload=b'a' * 100_000)
        start = store.identity['start']
        compressed = remote.client.payload[start:start + store.length - 1]
        store.identity = dict(store.identity, compressed_bytes=len(compressed))
        store.length = len(compressed)
        store.count = (store.length + store.chunk_size - 1) // store.chunk_size
        for index in range(store.count):
            store.write(index, compressed[index * store.chunk_size:(index + 1) * store.chunk_size], min_free=0)
        target = self.root / 'output.db'
        with self.assertRaisesRegex(RuntimeError, 'verification failed'):
            cache.inflate_chunks(store, self.root / 'output.part', target, expected, min_free=0)
        self.assertFalse(target.exists())

    def test_complete_restore_uses_retained_chunks_and_preserves_manifest(self):
        remote, _, expected = self.fixture()
        public_archive = remote.client.payload
        data_root = self.root / 'public-data'
        data_root.mkdir()
        manifest = {'url': cache.BASE_URL + cache.ARCHIVE, 'archive_bytes': len(public_archive),
                    'members': {cache.MEMBER: expected}}
        manifest_path = data_root / (cache.ARCHIVE + '.manifest.json')
        manifest_path.write_text(json.dumps(manifest))
        original_manifest = manifest_path.read_bytes()
        check_calls = []
        with mock.patch.object(cache, 'ExactRemoteZip', side_effect=lambda url: Remote(public_archive)), \
                mock.patch.object(cache, 'check_free', side_effect=lambda path, additional=0, min_free=cache.MIN_FREE: check_calls.append(additional)):
            cache.restore(data_root, workers=2, download_only=True)
            initial_preflight = check_calls[0]
            check_calls.clear()
            canonical, target = cache.cache_paths(data_root)
            self.assertFalse(canonical.exists())
            self.assertFalse(target.exists())
            cache.restore(data_root, workers=2)
            retained_preflight = check_calls[0]
            chunks_dir = target.with_name(target.name + '.chunks')
            identity = json.loads((chunks_dir / 'identity.json').read_text())
            self.assertEqual(initial_preflight - retained_preflight, identity['compressed_bytes'])
            count = (identity['compressed_bytes'] + cache.CHUNK_SIZE - 1) // cache.CHUNK_SIZE
            unit = max(4096, os.statvfs(target.parent).f_frsize)
            self.assertEqual(initial_preflight,
                             max(identity['compressed_bytes'], expected['bytes'])
                             + cache.MAX_OUTPUT + cache.CHUNK_SIZE + (count + 4) * unit * 2)
            self.assertEqual(canonical.read_bytes(), self.payload)
            self.assertTrue(canonical.is_symlink())
            self.assertEqual(manifest_path.read_bytes(), original_manifest)
            # A complete rerun validates and reuses the final cache.
            cache.restore(data_root, workers=2)
        self.assertEqual(manifest_path.read_bytes(), original_manifest)

    def test_inflate_bounds_output_and_enforces_free_space_guard(self):
        payload = b'a' * (cache.MAX_OUTPUT * 2 + 123)
        remote, store, expected = self.fixture(payload=payload, chunk_size=cache.CHUNK_SIZE)
        self.download(remote, store)
        sizes = []
        with mock.patch.object(cache, 'check_free', side_effect=lambda path, additional, reserve: sizes.append(additional)):
            cache.inflate_chunks(store, self.root / 'output.part', self.root / 'output.db', expected, min_free=0)
        self.assertIn(cache.MAX_OUTPUT, sizes)
        self.assertTrue(all(0 < size <= cache.MAX_OUTPUT for size in sizes))
        cache.verify_file(self.root / 'output.db', expected)
        with mock.patch.object(cache.shutil, 'disk_usage', return_value=type('Usage', (), {'free': cache.MIN_FREE + 10})()):
            cache.check_free(self.root, additional=10)
            with self.assertRaisesRegex(RuntimeError, '6 GiB'):
                cache.check_free(self.root, additional=11)

    def test_internal_relative_link_works_with_existing_verifier(self):
        data_root = self.root / 'public-data'
        data_root.mkdir()
        canonical, target = cache.cache_paths(data_root)
        target.parent.mkdir()
        target.write_bytes(self.payload)
        expected = {'bytes': len(self.payload), 'crc32': zlib.crc32(self.payload),
                    'sha256': hashlib.sha256(self.payload).hexdigest()}
        manifest = {'url': cache.BASE_URL + cache.ARCHIVE, 'members': {cache.MEMBER: expected}}
        (data_root / (cache.ARCHIVE + '.manifest.json')).write_text(json.dumps(manifest))
        cache.publish_symlink(data_root, canonical, target)
        self.assertTrue(canonical.is_symlink())
        self.assertFalse(Path(os.readlink(canonical)).is_absolute())
        self.assertEqual(safe_member(data_root, cache.MEMBER), target)
        verified = verify_archive_data(data_root, cache.ARCHIVE, [cache.MEMBER])
        self.assertEqual(verified['members'][cache.MEMBER]['sha256'], expected['sha256'])
        self.assertEqual(cache.cache_paths(data_root), (canonical, target))

    def test_canonical_materialized_file_and_external_symlink_are_refused(self):
        data_root = self.root / 'public-data'
        canonical = data_root / cache.MEMBER
        canonical.parent.mkdir(parents=True)
        canonical.write_bytes(b'Existing materialized DB must survive.')
        with self.assertRaisesRegex(RuntimeError, 'refusing to replace'):
            cache.cache_paths(data_root)
        self.assertEqual(canonical.read_bytes(), b'Existing materialized DB must survive.')
        canonical.unlink()  # Synthetic fixture only.
        external = self.root / 'external.db'
        external.write_bytes(b'Outside the data root.')
        canonical.symlink_to(external)
        with self.assertRaises(ValueError):
            cache.cache_paths(data_root)
        with self.assertRaises(ValueError):
            safe_member(data_root, cache.MEMBER)


if __name__ == '__main__':
    unittest.main()
