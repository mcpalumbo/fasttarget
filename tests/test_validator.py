"""The test validator itself.

A validator that cannot fail is worse than none: it reports success and hides
the thing it was meant to catch. Three of its checks were in that state, so
these tests pin both what passes and what must fail.
"""

import json
import os

import pytest

import validate_test
from validate_test import EXPECTED_GENES


def _empty_validation():
    return {'passed': [], 'warnings': [], 'failed': []}


def _messages(validation, key):
    return ' | '.join(validation[key])


# -------------------------------------------- expectation from the structures

@pytest.fixture
def offtarget_dir(tmp_path):
    """An offtarget directory with Foldseek results for four genes."""
    path = os.path.join(str(tmp_path), 'offtarget')
    foldseek = os.path.join(path, 'foldseek_results')

    searched = EXPECTED_GENES[:4]
    for gene in searched:
        gene_dir = os.path.join(foldseek, f'AF_{gene}_output_foldseek')
        os.makedirs(gene_dir, exist_ok=True)
        with open(os.path.join(gene_dir, f'AF_{gene}_foldseek_results.tsv'), 'w') as handle:
            handle.write('query\ttarget\n')
            handle.write(f'AF_{gene}\tHUMAN\n')

    with open(os.path.join(foldseek, 'human_foldseek_dict.json'), 'w') as handle:
        json.dump({gene: {'target_foldseek': 'HUMAN'} for gene in searched}, handle)

    return path, searched


def _breakdown(alphafold=(), pdb=(), colabfold=(), none=()):
    return {'pdb': list(pdb), 'alphafold': list(alphafold),
            'colabfold': list(colabfold), 'none': list(none)}


def test_the_expected_foldseek_count_comes_from_the_genes_with_a_structure(
        offtarget_dir, monkeypatch):
    """
    A flat 'expect 20' failed on every run where a gene had no structure, which
    is legitimate unless colabfold is enabled to model it.
    """
    path, searched = offtarget_dir
    monkeypatch.setattr(validate_test, 'reference_structure_breakdown',
                        lambda _: _breakdown(alphafold=searched, none=EXPECTED_GENES[4:]))

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=1, organism_dir='ignored')

    foldseek_failures = [m for m in validation['failed'] if 'Foldseek' in m]
    assert foldseek_failures == [], foldseek_failures
    assert '4 non-empty result TSVs for 4 searched genes' in _messages(validation, 'passed')


def test_the_genes_without_a_structure_are_a_warning_naming_colabfold(
        offtarget_dir, monkeypatch):
    path, searched = offtarget_dir
    monkeypatch.setattr(validate_test, 'reference_structure_breakdown',
                        lambda _: _breakdown(alphafold=searched, none=EXPECTED_GENES[4:]))

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=1, organism_dir='ignored')

    warnings = _messages(validation, 'warnings')
    assert '16 of 20 genes have no reference structure' in warnings
    assert 'colabfold' in warnings


def test_the_structure_sources_are_broken_down(offtarget_dir, monkeypatch):
    path, searched = offtarget_dir
    monkeypatch.setattr(validate_test, 'reference_structure_breakdown',
                        lambda _: _breakdown(pdb=searched[:1], alphafold=searched[1:],
                                             none=EXPECTED_GENES[4:]))

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=1, organism_dir='ignored')

    assert 'PDB 1, AlphaFold 3, ColabFold 0, none 16' in _messages(validation, 'passed')


def test_a_missing_result_for_a_searched_gene_fails(offtarget_dir, monkeypatch):
    """One more gene has a structure than there are results: that is a failure."""
    path, searched = offtarget_dir
    monkeypatch.setattr(validate_test, 'reference_structure_breakdown',
                        lambda _: _breakdown(alphafold=EXPECTED_GENES[:5],
                                             none=EXPECTED_GENES[5:]))

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=1, organism_dir='ignored')

    assert any('Expected at least 5 non-empty' in m for m in validation['failed'])


def test_an_empty_foldseek_dictionary_fails(offtarget_dir, monkeypatch):
    """The real defect: three stale files kept every run reporting no hits."""
    path, searched = offtarget_dir
    monkeypatch.setattr(validate_test, 'reference_structure_breakdown',
                        lambda _: _breakdown(alphafold=searched, none=EXPECTED_GENES[4:]))

    with open(os.path.join(path, 'foldseek_results', 'human_foldseek_dict.json'), 'w') as handle:
        json.dump({}, handle)

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=1, organism_dir='ignored')

    assert any('should have 4 keys' in m for m in validation['failed'])


def test_a_dictionary_with_keys_but_no_hits_fails(offtarget_dir, monkeypatch):
    path, searched = offtarget_dir
    monkeypatch.setattr(validate_test, 'reference_structure_breakdown',
                        lambda _: _breakdown(alphafold=searched, none=EXPECTED_GENES[4:]))

    with open(os.path.join(path, 'foldseek_results', 'human_foldseek_dict.json'), 'w') as handle:
        json.dump({gene: {'target_foldseek': None} for gene in searched}, handle)

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=1, organism_dir='ignored')

    assert any('no gene has a hit' in m for m in validation['failed'])


def test_the_colabfold_outputs_are_only_checked_when_present(offtarget_dir, monkeypatch):
    path, searched = offtarget_dir
    monkeypatch.setattr(validate_test, 'reference_structure_breakdown',
                        lambda _: _breakdown(alphafold=searched, none=EXPECTED_GENES[4:]))

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=1, organism_dir='ignored')
    assert 'ColabFold' not in _messages(validation, 'failed')

    # With the file present and empty, it is a failure.
    with open(os.path.join(path, 'foldseek_results', 'human_foldseek_colabfold_dict.json'), 'w') as handle:
        json.dump({}, handle)

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=1, organism_dir='ignored')
    assert any('ColabFold' in m for m in validation['failed'])


# ------------------------------------------------------- microbiome manifest

def _microbiome(path, manifest=None, tsv_count=0):
    species_dir = os.path.join(path, 'microbiomes', 'human-gut', 'species_blast_results')
    os.makedirs(species_dir, exist_ok=True)

    for index in range(tsv_count):
        with open(os.path.join(species_dir, f'MGYG{index:09d}_offtarget_id50_cov70.tsv'), 'w') as handle:
            handle.write('hit\n')

    if manifest is not None:
        with open(os.path.join(species_dir, 'human-gut_offtarget_manifest.json'), 'w') as handle:
            json.dump(manifest, handle)

    for name in ['human_gut_offtarget_counts', 'human_gut_offtarget_norm',
                 'human_gut_genomes_analyzed']:
        for extension in ('tsv', 'csv'):
            with open(os.path.join(species_dir, f'{name}.{extension}'), 'w') as handle:
                handle.write('gene\tvalue\n' if extension == 'tsv' else 'gene,value\n')

    return species_dir


def test_the_searched_genome_count_is_read_from_the_manifest(tmp_path):
    """
    Counting per-genome TSVs reported zero on every completed run, because the
    pipeline deletes them after consolidating into Parquet.
    """
    path = os.path.join(str(tmp_path), 'offtarget')
    _microbiome(path, manifest={
        'indexed_representatives': 4744, 'representatives_with_hits': 4318,
        'hit_rows': 8521, 'genes_with_hits': 5, 'raw_tsv_cleanup': 'completed',
    }, tsv_count=0)

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=4744)

    assert any('4744 representatives searched' in m for m in validation['passed'])
    assert not any('MGYG' in m for m in validation['failed'])


def test_too_few_searched_representatives_fail(tmp_path):
    path = os.path.join(str(tmp_path), 'offtarget')
    _microbiome(path, manifest={
        'indexed_representatives': 12, 'representatives_with_hits': 3,
        'hit_rows': 20, 'genes_with_hits': 1, 'raw_tsv_cleanup': 'completed',
    })

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=4744)

    assert any('12 representatives searched' in m for m in validation['failed'])


def test_an_interrupted_consolidation_is_reported(tmp_path):
    """Cleanup not completed and the TSVs gone means something stopped halfway."""
    path = os.path.join(str(tmp_path), 'offtarget')
    _microbiome(path, manifest={
        'indexed_representatives': 4744, 'representatives_with_hits': 4318,
        'hit_rows': 8521, 'genes_with_hits': 5, 'raw_tsv_cleanup': 'pending',
    }, tsv_count=2)

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=4744)

    assert any('raw_tsv_cleanup' in m for m in validation['warnings'])


def test_without_a_manifest_the_per_genome_files_are_still_counted(tmp_path):
    path = os.path.join(str(tmp_path), 'offtarget')
    _microbiome(path, manifest=None, tsv_count=3)

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=4744)

    assert any('no consolidation manifest' in m for m in validation['failed'])


def test_a_corrupt_manifest_is_reported_not_ignored(tmp_path):
    path = os.path.join(str(tmp_path), 'offtarget')
    species_dir = _microbiome(path, manifest={'indexed_representatives': 4744})
    with open(os.path.join(species_dir, 'human-gut_offtarget_manifest.json'), 'w') as handle:
        handle.write('{not json')

    validation = _empty_validation()
    validate_test.check_offtarget(path, validation, expected_count=len(EXPECTED_GENES),
                                  min_microbiome_genomes=4744)

    assert any('Failed to parse' in m for m in validation['failed'])


# ------------------------------------------- breakdown from the pipeline itself

def test_the_breakdown_returns_none_when_the_organism_cannot_be_read(tmp_path):
    """
    The selection is not reimplemented in the validator, so when the pipeline's
    own function cannot run the validator falls back rather than guessing.
    """
    assert validate_test.reference_structure_breakdown(
        os.path.join(str(tmp_path), 'no_such_organism')
    ) is None
