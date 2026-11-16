"""Checksum memo of files.sha256_file / files.cached_by_stat."""
import json
import os
import time

import pytest

from ftscripts import files


@pytest.fixture(autouse=True)
def fresh_memo(monkeypatch):
    monkeypatch.setattr(files, '_HASH_MEMO', {})
    monkeypatch.setattr(files, '_HASH_STATE',
                        {'path': None, 'pid': None, 'dirty': False, 'last_save': 0.0})
    monkeypatch.setattr(files.atexit, 'register', lambda *args, **kwargs: None)


def _old_file(path, content):
    path.write_bytes(content)
    old = time.time_ns() - 3600 * 10**9
    os.utime(path, ns=(old, old))
    return path


def _forbid_reads(monkeypatch):
    def fail(path):
        raise AssertionError(f'{path} was re-read')
    monkeypatch.setattr(files, '_sha256_uncached', fail)


def test_unchanged_file_is_hashed_once(tmp_path, monkeypatch):
    path = _old_file(tmp_path / 'a.gbk', b'ACGT')
    first = files.sha256_file(path)
    _forbid_reads(monkeypatch)
    assert files.sha256_file(path) == first


def test_modified_file_is_rehashed(tmp_path):
    path = _old_file(tmp_path / 'a.gbk', b'ACGT')
    first = files.sha256_file(path)
    _old_file(path, b'ACGTT')
    assert files.sha256_file(path) != first


def test_same_size_rewrite_with_new_mtime_is_rehashed(tmp_path):
    path = _old_file(tmp_path / 'a.gbk', b'ACGT')
    first = files.sha256_file(path)
    path.write_bytes(b'TTTT')
    assert files.sha256_file(path) != first


def test_recently_modified_file_is_not_memoized(tmp_path, monkeypatch):
    path = tmp_path / 'fresh.gbk'
    path.write_bytes(b'ACGT')
    files.sha256_file(path)
    assert files._HASH_MEMO == {}


def test_memo_persists_across_runs(tmp_path, monkeypatch):
    cache = tmp_path / 'cache.json'
    path = _old_file(tmp_path / 'a.gbk', b'ACGT')
    files.enable_hash_cache(cache)
    first = files.sha256_file(path)
    files.save_hash_cache(prune=True)
    assert json.loads(cache.read_text())['schema'] == 1

    monkeypatch.setattr(files, '_HASH_MEMO', {})
    monkeypatch.setattr(files, '_HASH_STATE',
                        {'path': None, 'pid': None, 'dirty': False, 'last_save': 0.0})
    files.enable_hash_cache(cache)
    _forbid_reads(monkeypatch)
    assert files.sha256_file(path) == first


def test_prune_drops_deleted_files(tmp_path):
    cache = tmp_path / 'cache.json'
    keep = _old_file(tmp_path / 'keep', b'1')
    gone = _old_file(tmp_path / 'gone', b'2')
    files.enable_hash_cache(cache)
    files.sha256_file(keep)
    files.sha256_file(gone)
    gone.unlink()
    files.save_hash_cache(prune=True)
    assert list(json.loads(cache.read_text())['files']) == [str(keep)]


def test_corrupt_cache_file_is_ignored(tmp_path):
    cache = tmp_path / 'cache.json'
    cache.write_text('{not json')
    files.enable_hash_cache(cache)
    path = _old_file(tmp_path / 'a', b'ACGT')
    assert files.sha256_file(path) == files._sha256_uncached(path)


def test_cached_by_stat_reuses_and_invalidates(tmp_path):
    path = _old_file(tmp_path / 'a.gff', b'x')
    calls = []

    def compute(p):
        calls.append(p)
        return ['seq1']

    assert files.cached_by_stat('orphans', path, compute) == ['seq1']
    assert files.cached_by_stat('orphans', path, compute) == ['seq1']
    assert len(calls) == 1
    _old_file(path, b'xy')
    files.cached_by_stat('orphans', path, compute)
    assert len(calls) == 2
