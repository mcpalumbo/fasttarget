"""Database preparation: retries, taxonomy grouping and download integrity.

Most of this module downloads multi-gigabyte catalogues, which no fast test can
exercise. What is testable without network is the part that decides whether a
download is usable and how microbiome hits are normalised, and both are worth
pinning: a partial file presented as complete, or a wrong denominator, changes
results without failing.
"""

import hashlib
import socket

import pandas as pd
import pytest
import requests

import databases


# --------------------------------------------------------------- retrying

def test_a_call_that_succeeds_is_not_retried():
    calls = {'n': 0}

    def once():
        calls['n'] += 1
        return 'done'

    assert databases.retry_on_timeout(once, max_retries=3, delay=0) == 'done'
    assert calls['n'] == 1


def test_a_timeout_is_retried_and_can_succeed():
    calls = {'n': 0}

    def flaky():
        calls['n'] += 1
        if calls['n'] < 3:
            raise requests.exceptions.Timeout('slow')
        return 'done'

    assert databases.retry_on_timeout(flaky, max_retries=5, delay=0) == 'done'
    assert calls['n'] == 3


def test_a_socket_timeout_counts_as_a_timeout():
    calls = {'n': 0}

    def flaky():
        calls['n'] += 1
        if calls['n'] < 2:
            raise socket.timeout('slow')
        return 'done'

    assert databases.retry_on_timeout(flaky, max_retries=3, delay=0) == 'done'


def test_persistent_timeouts_raise_after_the_last_attempt():
    calls = {'n': 0}

    def always():
        calls['n'] += 1
        raise requests.exceptions.Timeout('slow')

    with pytest.raises(Exception, match='Failed after 3 retries'):
        databases.retry_on_timeout(always, max_retries=3, delay=0)
    assert calls['n'] == 3


def test_an_error_that_is_not_a_timeout_is_not_retried():
    """Retrying a permanent failure only delays the report."""
    calls = {'n': 0}

    def broken():
        calls['n'] += 1
        raise ValueError('wrong argument')

    with pytest.raises(ValueError):
        databases.retry_on_timeout(broken, max_retries=3, delay=0)
    assert calls['n'] == 1


def test_the_backoff_decorator_retries_network_errors(monkeypatch):
    monkeypatch.setattr(databases.time, 'sleep', lambda seconds: None)
    calls = {'n': 0}

    @databases.retry_with_backoff(max_retries=4, initial_delay=0)
    def flaky():
        calls['n'] += 1
        if calls['n'] < 3:
            raise requests.exceptions.ConnectionError('reset')
        return 'done'

    assert flaky() == 'done'
    assert calls['n'] == 3


def test_the_backoff_decorator_reraises_the_last_error(monkeypatch):
    monkeypatch.setattr(databases.time, 'sleep', lambda seconds: None)

    @databases.retry_with_backoff(max_retries=2, initial_delay=0)
    def always():
        raise requests.exceptions.Timeout('slow')

    with pytest.raises(requests.exceptions.Timeout):
        always()


def test_the_backoff_delay_grows_and_is_capped(monkeypatch):
    """
    Without a cap the wait between attempts grows without bound; without growth
    a struggling server is hammered.
    """
    waits = []
    monkeypatch.setattr(databases.time, 'sleep', lambda seconds: waits.append(seconds))

    @databases.retry_with_backoff(max_retries=6, initial_delay=2, backoff_factor=2, max_delay=8)
    def always():
        raise requests.exceptions.ConnectionError('reset')

    with pytest.raises(requests.exceptions.ConnectionError):
        always()

    assert waits == [2, 4, 8, 8, 8], 'doubling until the cap, then held there'


def test_a_non_network_error_is_not_retried_by_the_decorator(monkeypatch):
    monkeypatch.setattr(databases.time, 'sleep', lambda seconds: None)
    calls = {'n': 0}

    @databases.retry_with_backoff(max_retries=4, initial_delay=0)
    def broken():
        calls['n'] += 1
        raise KeyError('missing')

    with pytest.raises(KeyError):
        broken()
    assert calls['n'] == 1


# -------------------------------------------------------------- chain ids

@pytest.mark.parametrize('value,expected', [
    ('A', ['A']),
    ('A;B;C', ['A', 'B', 'C']),
    ('A; B ; C', ['A', 'B', 'C']),
    ('A;;B', ['A', 'B']),
    ('', []),
    ('   ', []),
    (None, []),
    (float('nan'), []),
])
def test_chain_identifiers_are_split_and_cleaned(value, expected):
    assert databases.parse_chain_ids(value) == expected


def test_a_numeric_chain_identifier_is_kept_as_text():
    assert databases.parse_chain_ids(1) == ['1']


# ------------------------------------------------- taxonomy normalisation

def _taxonomy(rows):
    columns = ['representative_genome_id', *databases.TAXONOMY_RANKS]
    frame = pd.DataFrame(rows, columns=columns)
    # 'species' is one of the ranks, and the representative id is separate.
    return frame


def test_species_groups_are_the_representative_genomes():
    taxonomy = _taxonomy([
        ('MGYG001', 'Bacteria', 'Firmicutes', 'Clostridia', 'Eubacteriales',
         'Lachnospiraceae', 'Blautia', 'Blautia obeum'),
        ('MGYG002', 'Bacteria', 'Firmicutes', 'Clostridia', 'Eubacteriales',
         'Lachnospiraceae', 'Blautia', 'Blautia wexlerae'),
    ])

    groups = databases.microbiome_taxonomy_group_ids(taxonomy)

    assert groups['species'] == ['MGYG001', 'MGYG002']


def test_a_higher_rank_group_is_its_full_path():
    """
    Two genera of the same name in different families are different groups, so
    the identifier is the path and not the name.
    """
    taxonomy = _taxonomy([
        ('MGYG001', 'Bacteria', 'P1', 'C1', 'O1', 'F1', 'G1', 'S1'),
        ('MGYG002', 'Bacteria', 'P1', 'C1', 'O1', 'F2', 'G1', 'S2'),
    ])

    groups = databases.microbiome_taxonomy_group_ids(taxonomy)

    assert groups['genus'] == ['Bacteria;P1;C1;O1;F1;G1', 'Bacteria;P1;C1;O1;F2;G1']
    assert len(groups['genus']) == 2, 'the same genus name in two families is two groups'


def test_duplicate_paths_are_one_group():
    taxonomy = _taxonomy([
        ('MGYG001', 'Bacteria', 'P1', 'C1', 'O1', 'F1', 'G1', 'S1'),
        ('MGYG002', 'Bacteria', 'P1', 'C1', 'O1', 'F1', 'G1', 'S2'),
    ])

    groups = databases.microbiome_taxonomy_group_ids(taxonomy)

    assert groups['genus'] == ['Bacteria;P1;C1;O1;F1;G1']
    assert len(groups['species']) == 2


def test_an_unclassified_path_is_excluded_from_its_rank():
    """
    An unclassified lineage cannot be a denominator: counting it would make the
    normalised prevalence depend on how much of the catalogue is annotated.
    """
    taxonomy = _taxonomy([
        ('MGYG001', 'Bacteria', 'P1', 'C1', 'O1', 'F1', 'G1', 'S1'),
        ('MGYG002', 'Bacteria', 'P1', 'C1', 'O1', 'unclassified', 'G2', 'S2'),
    ])

    groups = databases.microbiome_taxonomy_group_ids(taxonomy)

    assert groups['family'] == ['Bacteria;P1;C1;O1;F1']
    assert groups['genus'] == ['Bacteria;P1;C1;O1;F1;G1'], 'the unclassified family breaks the path'
    assert len(groups['species']) == 2, 'every representative is still a species group'


def test_an_unclassified_value_at_a_higher_rank_removes_the_deeper_ones():
    taxonomy = _taxonomy([
        ('MGYG001', 'Bacteria', 'unclassified', 'C1', 'O1', 'F1', 'G1', 'S1'),
    ])

    groups = databases.microbiome_taxonomy_group_ids(taxonomy)

    assert groups['domain'] == ['Bacteria']
    for rank in ('phylum', 'class', 'order', 'family', 'genus'):
        assert groups[rank] == [], rank


def test_every_rank_is_reported_even_when_empty():
    groups = databases.microbiome_taxonomy_group_ids(_taxonomy([]))

    assert set(groups) == set(databases.TAXONOMY_RANKS)
    assert all(value == [] for value in groups.values())


def test_group_identifiers_are_sorted_for_stability():
    """The order is a denominator list written to disk, so it must be stable."""
    taxonomy = _taxonomy([
        ('MGYG003', 'Bacteria', 'P2', 'C1', 'O1', 'F1', 'G1', 'S1'),
        ('MGYG001', 'Bacteria', 'P1', 'C1', 'O1', 'F1', 'G1', 'S1'),
    ])

    groups = databases.microbiome_taxonomy_group_ids(taxonomy)

    assert groups['species'] == sorted(groups['species'])
    assert groups['phylum'] == sorted(groups['phylum'])


# ------------------------------------------------------- download integrity

def test_a_checksum_is_computed_over_the_whole_file(tmp_path):
    path = tmp_path / 'payload.bin'
    content = b'x' * (10 * 1024 * 1024)   # larger than the read chunk
    path.write_bytes(content)

    assert databases._file_sha256(path) == hashlib.sha256(content).hexdigest()


def test_an_existing_file_with_the_right_checksum_is_reused(tmp_path):
    path = tmp_path / 'chembl.tar.gz'
    path.write_bytes(b'payload')
    digest = hashlib.sha256(b'payload').hexdigest()

    result = databases._download_ligand_file('https://example.org/x', path,
                                             expected_sha256=digest)

    assert result == path


def test_an_existing_file_with_the_wrong_checksum_is_refused(tmp_path):
    """
    A truncated or tampered archive must not be used: it would fail later, deep
    inside the database preparation, with an unrelated error.
    """
    path = tmp_path / 'chembl.tar.gz'
    path.write_bytes(b'corrupted')

    with pytest.raises(ValueError, match='Checksum mismatch'):
        databases._download_ligand_file('https://example.org/x', path,
                                        expected_sha256=hashlib.sha256(b'payload').hexdigest())


def test_a_partial_download_needing_aria2_is_not_resumed_blindly(tmp_path, monkeypatch):
    """
    An .aria2 control file means the partial file has holes, so appending to it
    with a plain ranged request would produce a corrupt archive.
    """
    destination = tmp_path / 'bindingdb.tar.gz'
    partial = tmp_path / 'bindingdb.tar.gz.part'
    partial.write_bytes(b'incomplete')
    (tmp_path / 'bindingdb.tar.gz.part.aria2').write_bytes(b'control')

    # The accelerator is looked up under sys.prefix/bin, so pointing that at an
    # empty directory makes it unavailable without patching pathlib itself.
    empty_prefix = tmp_path / 'prefix'
    (empty_prefix / 'bin').mkdir(parents=True)
    monkeypatch.setattr(databases.sys, 'prefix', str(empty_prefix))

    with pytest.raises(RuntimeError, match='aria2c'):
        databases._download_ligand_file('https://example.org/x', destination)
