"""Parameter manifests: what counts as the same parameters, and what does not.

Each test here corresponds to a way the resume logic silently reused a result
built with other settings.
"""

import json
import os

import pytest

from ftscripts import files


@pytest.fixture
def output(tmp_path):
    path = os.path.join(str(tmp_path), 'result.tsv')
    with open(path, 'w') as handle:
        handle.write('content\n')
    return path


def test_missing_output_is_a_miss(tmp_path):
    missing = os.path.join(str(tmp_path), 'absent.tsv')
    assert files.manifest_matches(missing, {'a': 1}) is False


def test_empty_output_is_a_miss_even_with_a_manifest(tmp_path):
    empty = os.path.join(str(tmp_path), 'empty.tsv')
    open(empty, 'w').close()
    files.write_manifest(empty, {'a': 1})
    assert files.manifest_matches(empty, {'a': 1}) is False


def test_same_parameters_hit(output):
    files.write_manifest(output, {'min_identity': 95})
    assert files.manifest_matches(output, {'min_identity': 95}) is True


def test_changed_parameter_misses_and_is_named(output):
    files.write_manifest(output, {'min_identity': 95})
    assert files.manifest_matches(output, {'min_identity': 90}) is False
    assert files.manifest_differences(output, {'min_identity': 90}) == ['min_identity 95 -> 90']


@pytest.mark.parametrize('probe', [
    {'min_core_freq': 99, 'min_identity': 95},
    {'min_core_freq': 99.0, 'min_identity': '95'},
    {'min_core_freq': ' 99 ', 'min_identity': 95.0},
])
def test_values_compare_equal_across_yaml_nextflow_and_json(output, probe):
    """The same setting arrives as 99, 99.0 or '99' depending on the caller."""
    files.write_manifest(output, {'min_core_freq': 99, 'min_identity': 95})
    assert files.manifest_matches(output, probe) is True


def test_a_flag_is_not_the_number_one(output):
    """True == 1.0 in Python, so booleans are stored tagged."""
    files.write_manifest(output, {'uniparc_fallback': True})
    assert files.manifest_matches(output, {'uniparc_fallback': True}) is True
    assert files.manifest_matches(output, {'uniparc_fallback': 1}) is False
    assert files.manifest_differences(output, {'uniparc_fallback': False}) == [
        'uniparc_fallback True -> False'
    ]


def test_none_and_text_survive_a_round_trip(output):
    files.write_manifest(output, {'engine': 'docker', 'reference': None})
    assert files.manifest_matches(output, {'engine': 'docker', 'reference': None}) is True
    assert files.manifest_matches(output, {'engine': 'singularity', 'reference': None}) is False


def test_a_field_added_later_is_unknown_not_changed(output):
    """
    A manifest written before a field existed must not stop an expensive run.

    Treating an absent field as a difference made every added parameter break
    existing output directories.
    """
    files.write_manifest(output, {'min_core_freq': 99})
    current = {'min_core_freq': 99, 'genome_count': 19, 'genome_set': 'abc123'}

    assert files.manifest_differences(output, current) == []
    assert files.manifest_unrecorded(output, current) == ['genome_count', 'genome_set']
    # A real change is still a difference, alongside the unknown fields.
    assert files.manifest_differences(output, dict(current, min_core_freq=90)) == [
        'min_core_freq 99 -> 90'
    ]


def test_no_manifest_reports_every_parameter_as_unknown(output):
    assert files.read_manifest(output) is None
    assert files.manifest_differences(output, {'a': 1, 'b': 2}) == []
    assert files.manifest_unrecorded(output, {'a': 1, 'b': 2}) == ['a', 'b']
    # Cheap outputs still rebuild, to replace a result of unknown provenance.
    assert files.manifest_matches(output, {'a': 1}) is False


def test_corrupt_manifest_means_recompute_not_crash(output):
    files.write_manifest(output, {'a': 1})
    with open(files.manifest_path(output), 'w') as handle:
        handle.write('{not json')
    assert files.read_manifest(output) is None
    assert files.manifest_matches(output, {'a': 1}) is False


def test_manifest_records_resolved_values_and_extra(output):
    files.write_manifest(output, {'min_identity': 95}, extra={'genomes_analyzed': 19})
    manifest = files.read_manifest(output)
    assert manifest['parameters'] == {'min_identity': 95.0}
    assert manifest['extra'] == {'genomes_analyzed': 19}
    assert manifest['output'] == 'result.tsv'
    assert 'written' in manifest


def test_fingerprint_notices_additions_removals_and_substitutions():
    base = files.fileset_fingerprint(['/x/a.gff', '/x/b.gff'])

    assert base == files.fileset_fingerprint(['/y/b.gff', '/y/a.gff']), 'order and path must not matter'
    assert base != files.fileset_fingerprint(['/x/a.gff', '/x/b.gff', '/x/c.gff'])
    assert base != files.fileset_fingerprint(['/x/a.gff'])
    # A substitution keeps the count, so the count alone is not enough.
    assert base != files.fileset_fingerprint(['/x/a.gff', '/x/z.gff'])


def test_step_completed_distinguishes_empty_from_crashed(tmp_path):
    """An empty output is a valid result only when the step recorded finishing."""
    empty = os.path.join(str(tmp_path), 'search.tsv')
    open(empty, 'w').close()

    assert files.file_check(empty) is False
    assert files.step_completed(empty) is False, 'empty with no marker may be a crash'

    files.mark_complete(empty)
    assert files.step_completed(empty) is True

    with open(empty, 'w') as handle:
        handle.write('a hit\n')
    assert files.file_check(empty) is True
    assert files.step_completed(empty) is True
