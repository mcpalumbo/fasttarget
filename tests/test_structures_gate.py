"""The structure pipeline's parameter gate.

The whole structures stage is gated on one summary file, so enabling
uniparc_fallback or pocket_full_mode over an existing output directory used to
reuse a summary built without it. The stage downloads structures and runs
FPocket, P2Rank and optionally ColabFold, so a mismatch stops the run and names
what to remove rather than deleting hours of work on its own.
"""

import os

import pytest

from ftscripts import files, structures


BASE = {
    'specie_taxid': 2104,
    'strain_taxid': 272634,
    'full_mode': False,
    'uniparc_fallback': True,
    'colabfold': False,
    'colabfold_all_models': False,
    'amber_option': False,
    'gpu_option': False,
    'resolution_cutoff': 3.5,
    'coverage_cutoff': 40.0,
}


@pytest.fixture
def summary(organism):
    """An existing structure summary with its parameters recorded."""
    output_path, organism_name = organism
    structure_dir = os.path.join(output_path, organism_name, 'structures')
    path = os.path.join(structure_dir, f'{organism_name}_final_structure_summary.tsv')
    with open(path, 'w') as handle:
        handle.write('gene\tstructure\nSYN_0001\tAF_X\n')
    return output_path, organism_name, structure_dir, path


def _check(summary_fixture, parameters):
    output_path, organism_name, structure_dir, path = summary_fixture
    structures._check_structure_pipeline_parameters(
        output_path, organism_name, structure_dir, path, parameters
    )


def test_matching_parameters_pass_silently(summary, capsys):
    files.write_manifest(summary[3], BASE)

    _check(summary, BASE)

    assert capsys.readouterr().out == ''


def test_no_summary_yet_is_not_a_mismatch(organism):
    """Nothing to reuse means nothing to check."""
    output_path, organism_name = organism
    structure_dir = os.path.join(output_path, organism_name, 'structures')
    absent = os.path.join(structure_dir, 'absent.tsv')

    structures._check_structure_pipeline_parameters(
        output_path, organism_name, structure_dir, absent, BASE
    )


@pytest.mark.parametrize('option,value', [
    ('uniparc_fallback', False),
    ('full_mode', True),
    ('colabfold', True),
    ('colabfold_all_models', True),
    ('amber_option', True),
    ('gpu_option', True),
    ('resolution_cutoff', 2.0),
    ('coverage_cutoff', 80.0),
])
def test_every_recorded_option_is_a_cache_key(summary, option, value):
    files.write_manifest(summary[3], BASE)

    with pytest.raises(RuntimeError) as error:
        _check(summary, dict(BASE, **{option: value}))
    assert option in str(error.value)


def test_a_refusal_names_the_files_to_remove(summary):
    files.write_manifest(summary[3], BASE)
    output_path, organism_name, structure_dir, path = summary

    with pytest.raises(RuntimeError) as error:
        _check(summary, dict(BASE, uniparc_fallback=False))

    message = str(error.value)
    assert f'{organism_name}_final_structure_summary.tsv' in message
    assert f'{organism_name}_structure_data.json' in message
    assert '_structure_summary.tsv' in message, 'the per-gene summaries'
    assert 'human_foldseek_dict.json' in message, 'the Foldseek cache outside this module'
    assert 'does not download them again' in message, 'say what is kept'


def test_nothing_is_deleted_by_a_refusal(summary):
    files.write_manifest(summary[3], BASE)

    with pytest.raises(RuntimeError):
        _check(summary, dict(BASE, full_mode=True))

    assert files.file_check(summary[3]), 'the summary must survive the refusal'


def test_a_changed_taxid_also_invalidates_the_downloaded_proteome(summary):
    """
    The taxids select the UniProt proteome and order the mapping tiers, so the
    proteome, its CD-HIT clusters and the BLAST databases are for the wrong taxon.
    """
    files.write_manifest(summary[3], BASE)

    with pytest.raises(RuntimeError) as error:
        _check(summary, dict(BASE, strain_taxid=2097))

    assert 'uniprot_files' in str(error.value)


def test_a_threshold_change_keeps_the_proteome(summary):
    """Re-downloading a proteome that is still valid would be wasteful."""
    files.write_manifest(summary[3], BASE)

    with pytest.raises(RuntimeError) as error:
        _check(summary, dict(BASE, uniparc_fallback=False))

    assert 'uniprot_files' not in str(error.value)


def test_an_option_recorded_by_an_older_version_warns_and_reuses(summary, capsys):
    """
    A summary built before an option was recorded has unknown parameters, not
    wrong ones, so it is reused with a warning naming what is unknown.
    """
    recorded = {key: value for key, value in BASE.items()
                if key not in ('specie_taxid', 'strain_taxid')}
    files.write_manifest(summary[3], recorded)

    _check(summary, BASE)

    output = capsys.readouterr().out
    assert 'does not record how it was built' in output
    assert 'specie_taxid' in output and 'strain_taxid' in output


def test_no_manifest_at_all_warns_and_reuses(summary, capsys):
    _check(summary, BASE)

    output = capsys.readouterr().out
    assert 'does not record how it was built' in output
    assert files.file_check(summary[3])
