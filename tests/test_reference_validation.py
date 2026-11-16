"""Reference structures are checked against the residues actually resolved."""
import csv
import os

import pandas as pd

from ftscripts import structures


def _pdb(path, chain, n_residues, ligand=None):
    """Write a minimal PDB with n_residues alanines (N, CA, C) in one chain."""
    lines, serial = [], 1
    for i in range(1, n_residues + 1):
        for atom in ('N', 'CA', 'C'):
            lines.append(f"ATOM  {serial:5d}  {atom:<3} ALA {chain}{i:4d}    "
                         f"{i * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           {atom[0]}")
            serial += 1
    if ligand:
        lines.append(f"HETATM{serial:5d}  C1  {ligand:<3} {chain}{900:4d}    "
                     f"{0.0:8.3f}{5.0:8.3f}{0.0:8.3f}  1.00 20.00           C")
    lines.append('END')
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


def _summary(locus_dir, locus, rows):
    columns = ['uniprot_id', 'structure_type', 'structure_id', 'method', 'resolution', 'chain',
               'residue_range', 'coverage', 'sequence_length', 'is_reference']
    pd.DataFrame(rows, columns=columns).to_csv(
        os.path.join(locus_dir, f'{locus}_structure_summary.tsv'), sep='\t', index=False,
        quoting=csv.QUOTE_NONNUMERIC)


def _reference(locus_dir, locus):
    df = pd.read_csv(os.path.join(locus_dir, f'{locus}_structure_summary.tsv'), sep='\t',
                     keep_default_na=False, dtype=str)
    ref = df[df.is_reference.str.lower() == 'true']
    assert len(ref) == 1
    return ref.iloc[0]


def test_observed_coverage_counts_resolved_residues(tmp_path):
    path = tmp_path / 'PDB_1ABC.pdb'
    _pdb(path, 'A', 10)
    assert structures.observed_chain_coverage(str(path), ['A'], 100) == 10.0
    assert structures.observed_chain_coverage(str(path), ['B'], 100) == 0.0


def test_fragment_reference_falls_back_to_alphafold(tmp_path):
    locus = 'LOC_1'
    locus_dir = tmp_path / locus
    (locus_dir / 'P1').mkdir(parents=True)
    _pdb(locus_dir / 'P1' / 'PDB_1ABC.pdb', 'A', 7)       # declared 1-300, 7 resolved
    _pdb(locus_dir / 'P1' / 'AF_P1.pdb', 'A', 300)
    _summary(str(locus_dir), locus, [
        ['P1', 'PDB', '1ABC', 'X-ray', 2.0, 'A', '1-300', 100.0, 300, True],
        ['P1', 'AlphaFold', 'P1', 'AlphaFold', '', 'A', '1-300', 100.0, 300, False],
    ])
    (locus_dir / 'pockets').mkdir()
    stale = locus_dir / 'pockets' / f'best_DS_fpocket_{locus}.json'
    stale.write_text('{"1ABC": "No_pockets"}')
    assert structures.validate_reference_structure_for_locus(locus, str(locus_dir), 3.5, 40.0)
    ref = _reference(str(locus_dir), locus)
    assert ref.structure_type == 'AlphaFold'
    assert not stale.exists()  # best pockets of the former reference are not reused


def test_fragment_reference_prefers_a_pdb_that_really_passes(tmp_path):
    locus = 'LOC_2'
    locus_dir = tmp_path / locus
    (locus_dir / 'P1').mkdir(parents=True)
    _pdb(locus_dir / 'P1' / 'PDB_1ABC.pdb', 'A', 7)
    _pdb(locus_dir / 'P1' / 'PDB_2DEF.pdb', 'A', 200)
    _summary(str(locus_dir), locus, [
        ['P1', 'PDB', '1ABC', 'X-ray', 1.5, 'A', '1-300', 100.0, 300, True],
        ['P1', 'PDB', '2DEF', 'X-ray', 3.0, 'A', '1-200', 66.7, 300, False],
    ])
    assert structures.validate_reference_structure_for_locus(locus, str(locus_dir), 3.5, 40.0)
    ref = _reference(str(locus_dir), locus)
    assert ref.structure_id == '2DEF'
    assert os.path.isfile(locus_dir / 'P1' / 'PDB_2DEF_A_ref.pdb')


def test_complete_reference_is_kept(tmp_path):
    locus = 'LOC_3'
    locus_dir = tmp_path / locus
    (locus_dir / 'P1').mkdir(parents=True)
    _pdb(locus_dir / 'P1' / 'PDB_1ABC.pdb', 'A', 280)
    _summary(str(locus_dir), locus, [
        ['P1', 'PDB', '1ABC', 'X-ray', 2.0, 'A', '1-300', 100.0, 300, True],
    ])
    assert not structures.validate_reference_structure_for_locus(locus, str(locus_dir), 3.5, 40.0)
    assert _reference(str(locus_dir), locus).structure_id == '1ABC'


def test_missing_reference_row_is_selected_again(tmp_path):
    # The strain AlphaFold row chosen as reference was dropped after a 404.
    locus = 'LOC_4'
    locus_dir = tmp_path / locus
    (locus_dir / 'P1').mkdir(parents=True)
    _pdb(locus_dir / 'P1' / 'PDB_1ABC.pdb', 'A', 60)       # 20%: not a valid reference
    _pdb(locus_dir / 'P1' / 'AF_P1.pdb', 'A', 300)
    _summary(str(locus_dir), locus, [
        ['P1', 'PDB', '1ABC', 'X-ray', 2.0, 'A', '1-60', 20.0, 300, False],
        ['P1', 'AlphaFold', 'P1', 'AlphaFold', '', 'A', '1-300', 100.0, 300, False],
    ])
    assert structures.validate_reference_structure_for_locus(locus, str(locus_dir), 3.5, 40.0)
    assert _reference(str(locus_dir), locus).structure_type == 'AlphaFold'


def test_long_ccd_codes_are_detected_and_filtered(tmp_path):
    from Bio.PDB.Residue import Residue
    good = tmp_path / 'good_ref.pdb'
    _pdb(good, 'A', 5, ligand='HEM')
    assert not structures.pdb_file_has_overflowing_resnames(str(good))
    # A five-character CCD code as Biopython writes it: it fills column 21
    bad = tmp_path / 'bad_ref.pdb'
    bad.write_text(good.read_text() +
                   'HETATM 5213  CAB A1H49 A 803     -73.024 -14.895  -3.443  0.83 32.19           C\n')
    assert structures.pdb_file_has_overflowing_resnames(str(bad))
    assert structures.pdb_format_residue(Residue(('H_HEM', 900, ' '), 'HEM', ' '))
    assert not structures.pdb_format_residue(Residue(('H_A1H49', 803, ' '), 'A1H49', ' '))
