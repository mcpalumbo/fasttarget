"""Choosing the one structure that represents a gene.

Every downstream step uses this choice: pockets are found on it, Foldseek
searches with it, and the druggability score comes from it. The rules are a
priority order with two cutoffs, so the cases that matter are the boundaries and
the fallbacks.

The first five tests are ported from an earlier local suite; the rest extend
them to the cutoff edges and the tie-breaks.
"""

import pytest

from ftscripts.structures import collect_structures_for_uniprot, select_reference_structure


def pdb(method, resolution, coverage):
    return {'structure_type': 'PDB', 'method': method,
            'resolution': resolution, 'coverage': coverage}


def alphafold(coverage=100):
    return {'structure_type': 'AlphaFold', 'method': 'AlphaFold', 'coverage': coverage}


# ------------------------------------------------------------ priority order

def test_an_eligible_xray_beats_a_better_resolution_em_structure():
    """X-ray is preferred over electron microscopy even at worse resolution."""
    structures = [
        pdb('Electron microscopy', 1.8, 100),
        pdb('X-ray diffraction', 1.9, 60),
        pdb('X-ray diffraction', 2.0, 90),
    ]

    assert select_reference_structure(structures) == 2


def test_coverage_decides_between_structures_of_similar_resolution():
    structures = [
        pdb('X-ray diffraction', 2.6, 60),
        pdb('X-ray diffraction', 3.2, 90),
        pdb('X-ray diffraction', 2.9, 80),
    ]

    assert select_reference_structure(structures) == 2


def test_alphafold_is_used_when_no_experimental_structure_passes():
    structures = [
        pdb('X-ray diffraction', 4.0, 100),   # resolution too poor
        pdb('X-ray diffraction', 2.0, 20),    # coverage too low
        alphafold(),
    ]

    assert select_reference_structure(structures) == 2


def test_nothing_is_selected_when_nothing_qualifies():
    """
    A gene with only unusable structures has no reference, which is a result:
    the pipeline reports it rather than picking a bad structure.
    """
    assert select_reference_structure([
        pdb('X-ray diffraction', 4.0, 100),
        pdb('X-ray diffraction', 2.0, 20),
    ]) is None
    assert select_reference_structure([]) is None


# ----------------------------------------------------------------- cutoffs

@pytest.mark.parametrize('resolution,selected', [
    (3.49, True),
    (3.5, True),     # inclusive: exactly at the limit still qualifies
    (3.51, False),
])
def test_the_resolution_cutoff_is_inclusive(resolution, selected):
    structures = [pdb('X-ray diffraction', resolution, 100)]

    assert (select_reference_structure(structures) == 0) is selected


@pytest.mark.parametrize('coverage,selected', [
    (40.1, True),
    (40.0, True),    # inclusive, as passes_cutoffs implements it
    (39.9, False),
])
def test_the_coverage_cutoff_is_inclusive(coverage, selected):
    structures = [pdb('X-ray diffraction', 2.0, coverage)]

    assert (select_reference_structure(structures) == 0) is selected


def test_the_cutoffs_can_be_relaxed_or_tightened():
    structures = [pdb('X-ray diffraction', 4.0, 50)]

    assert select_reference_structure(structures) is None
    assert select_reference_structure(structures, resolution_cutoff=5.0) == 0
    assert select_reference_structure(structures, resolution_cutoff=5.0,
                                      coverage_cutoff=60) is None


def test_a_structure_without_a_resolution_does_not_pass_as_perfect():
    """
    A missing resolution must not read as zero, which would make an unmeasured
    structure the best candidate.
    """
    structures = [
        {'structure_type': 'PDB', 'method': 'X-ray diffraction',
         'resolution': None, 'coverage': 100},
        pdb('X-ray diffraction', 2.0, 80),
    ]

    assert select_reference_structure(structures) == 1


# ---------------------------------------------------- building the candidates

def test_a_good_pdb_structure_is_enough_on_its_own():
    """No AlphaFold candidate is added when an experimental one qualifies."""
    info = {
        'Sequence': 'A' * 100,
        'AlphaFoldDB': 'AF-P1',
        'PDB_id': [{'ID': '1ABC', 'method': 'X-ray', 'resolution': '2.0 A',
                    'chains': 'A=1-80'}],
    }

    rows = collect_structures_for_uniprot('P1', info)

    assert [row['structure_type'] for row in rows] == ['PDB']
    assert rows[0]['coverage'] == 80.0
    assert rows[0]['residue_range'] == '1-80'


def test_alphafold_is_added_when_the_experimental_structure_is_inadequate():
    info = {
        'Sequence': 'A' * 100,
        'PDB_id': [{'ID': '1ABC', 'method': 'X-ray', 'resolution': '4.0 A',
                    'chains': 'A=1-20'}],
    }

    rows = collect_structures_for_uniprot('P1', info)

    assert [row['structure_type'] for row in rows] == ['PDB', 'AlphaFold']
    assert rows[1]['structure_id'] == 'P1', 'the accession is the model identifier'


def test_coverage_is_computed_from_the_chain_range_and_the_sequence():
    info = {
        'Sequence': 'A' * 200,
        'PDB_id': [{'ID': '1ABC', 'method': 'X-ray', 'resolution': '2.0 A',
                    'chains': 'A=1-100'}],
    }

    rows = collect_structures_for_uniprot('P1', info)

    assert rows[0]['coverage'] == 50.0


def test_a_uniprot_entry_with_no_structures_still_offers_the_model():
    """
    With neither a PDB entry nor an AlphaFoldDB cross-reference, the accession
    itself is tried against AlphaFold DB.
    """
    rows = collect_structures_for_uniprot('P1', {'Sequence': 'A' * 100})

    assert [row['structure_type'] for row in rows] == ['AlphaFold']
    assert rows[0]['structure_id'] == 'P1'
