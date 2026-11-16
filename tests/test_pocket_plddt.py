"""Pocket pLDDT: per-pocket metrics and the mean pLDDT >= 70 best-pocket filter."""
import json

from ftscripts import structures


def _atom(serial, name, resnum, bfactor, chain='A', resname='ALA'):
    return (f'ATOM  {serial:5d} {name:<4} {resname} {chain}{resnum:4d}    '
            f'{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{bfactor:6.2f}           C\n')


def _model(path, plddts):
    lines, serial = [], 1
    for resnum, value in enumerate(plddts, 1):
        for name in ('N', 'CA'):
            lines.append(_atom(serial, name, resnum, value))
            serial += 1
    path.write_text(''.join(lines) + 'END\n')


def _fpocket_locus(tmp_path, name, plddts, pockets):
    """Locus with a model or PDB, an fpocket folder and pocketN_atm.pdb files.

    pockets: {number: (druggability, [residue numbers])}
    """
    locus = tmp_path / 'PA14_00001'
    (locus / 'Q1').mkdir(parents=True)
    _model(locus / 'Q1' / f'{name}.pdb', plddts)
    folder = locus / 'pockets' / f'{name}_fpocket'
    (folder / 'pockets').mkdir(parents=True)
    info = []
    for number, (score, residues) in pockets.items():
        info.append(f'Pocket {number} :\n\tScore : \t0.1\n\tDruggability Score : \t{score:.3f}\n')
        # fpocket writes B-factor 0 in the atm files: the model file must be used.
        atoms = ''.join(_atom(i, 'CA', r, 0.0) for i, r in enumerate(residues, 1))
        (folder / 'pockets' / f'pocket{number}_atm.pdb').write_text(atoms)
    (folder / f'{name}_info.txt').write_text('\n'.join(info))
    return folder


def test_residue_plddt_map_reads_ca_bfactor(tmp_path):
    _model(tmp_path / 'AF_Q1.pdb', [90.0, 40.0])
    assert structures.residue_plddt_map(tmp_path / 'AF_Q1.pdb') == {('A', 1, ''): 90.0, ('A', 2, ''): 40.0}


def test_residue_plddt_map_rescales_fractions(tmp_path):
    _model(tmp_path / 'AF_Q1.pdb', [0.9, 0.4])
    assert structures.residue_plddt_map(tmp_path / 'AF_Q1.pdb')[('A', 1, '')] == 90.0


def test_pocket_plddt_stats():
    plddt = {('A', 1, ''): 90.0, ('A', 2, ''): 40.0, ('A', 3, ''): 80.0, ('A', 4, ''): 30.0}
    stats = structures.pocket_plddt_stats({('A', 1, ''), ('A', 2, ''), ('A', 3, ''), ('A', 4, '')}, plddt)
    assert stats == {'mean_plddt': 60.0, 'frac_plddt_lt50': 0.5}
    assert structures.pocket_plddt_stats({('B', 1, '')}, plddt) == {'mean_plddt': None, 'frac_plddt_lt50': None}


def test_p2rank_residue_ids_parsing():
    assert structures.p2rank_pocket_residues('A_12 A_13B B_-1') == {('A', 12, ''), ('A', 13, 'B'), ('B', -1, '')}


def test_fpocket_best_pocket_skips_low_plddt(tmp_path):
    # Residues 1-3 disordered (pLDDT 30), 4-6 confident (pLDDT 95).
    folder = _fpocket_locus(tmp_path, 'AF_Q1', [30, 30, 30, 95, 95, 95],
                            {1: (0.9, [1, 2, 3]), 2: (0.4, [4, 5, 6]), 3: (0.1, [4, 5])})
    pockets = structures.pockets_data_to_dict(str(folder))
    assert pockets['Q1']['Pocket 1']['mean_plddt'] == 30.0
    assert pockets['Q1']['Pocket 1']['frac_plddt_lt50'] == 1.0
    best = structures.pockets_filter(pockets)['Q1']
    assert best['pocket'] == 'Pocket 2'
    assert best['maxDS'] == 0.4
    assert best['mean_plddt'] == 95.0 and best['frac_plddt_lt50'] == 0.0
    assert best['plddt_filter'] == 'passed'


def test_fpocket_best_pocket_falls_back_when_none_passes(tmp_path):
    folder = _fpocket_locus(tmp_path, 'AF_Q1', [30, 30, 60, 60],
                            {1: (0.3, [1, 2]), 2: (0.8, [3, 4])})
    best = structures.pockets_filter(structures.pockets_data_to_dict(str(folder)))['Q1']
    assert best['pocket'] == 'Pocket 2'
    assert best['maxDS'] == 0.8
    assert best['mean_plddt'] == 60.0
    assert best['plddt_filter'] == 'none_passed'


def test_experimental_structure_keeps_top_pocket(tmp_path):
    # B-factors of a crystal structure are not pLDDT: no filter, empty columns.
    folder = _fpocket_locus(tmp_path, 'PDB_1ABC_A_ref', [30, 30, 95, 95],
                            {1: (0.9, [1, 2]), 2: (0.4, [3, 4])})
    best = structures.pockets_filter(structures.pockets_data_to_dict(str(folder)))['1ABC']
    assert best['pocket'] == 'Pocket 1'
    assert best['mean_plddt'] is None and best['frac_plddt_lt50'] is None
    assert best['plddt_filter'] == 'not_model'


def test_cached_fpocket_json_gets_plddt(tmp_path):
    folder = _fpocket_locus(tmp_path, 'AF_Q1', [30, 95], {1: (0.9, [1]), 2: (0.4, [2])})
    old = {'Q1': {'Pocket 1': {'Druggability Score': 0.9}, 'Pocket 2': {'Druggability Score': 0.4}}}
    (folder / 'all_pockets.json').write_text(json.dumps(old))
    pockets = structures.pockets_data_to_dict(str(folder))
    assert pockets['Q1']['Pocket 2']['mean_plddt'] == 95.0
    assert json.loads((folder / 'all_pockets.json').read_text())['Q1']['Pocket 1']['mean_plddt'] == 30.0


def test_p2rank_best_pocket_skips_low_plddt(tmp_path):
    locus = tmp_path / 'PA14_00001'
    (locus / 'Q1').mkdir(parents=True)
    _model(locus / 'Q1' / 'AF_Q1.pdb', [30, 30, 90, 90])
    folder = locus / 'pockets' / 'AF_Q1_p2rank'
    folder.mkdir(parents=True)
    (folder / 'AF_Q1.pdb_predictions.csv').write_text(
        'name     ,  rank,   score, probability, residue_ids\n'
        'pocket1  ,     1,   10.0,       0.800, A_1 A_2\n'
        'pocket2  ,     2,    5.0,       0.500, A_3 A_4\n')
    pockets = structures.p2rank_data_to_dict(str(folder))
    best = structures.p2rank_filter(pockets)['Q1']
    assert best['pocket'].strip() == 'Pocket pocket2'
    assert best['max_probability'] == 0.5
    assert best['mean_plddt'] == 90.0
    assert best['plddt_filter'] == 'passed'


def test_best_pocket_without_filter_marker_is_stale():
    assert structures._best_pockets_missing_plddt({'Q1': {'maxDS': 0.5, 'pocket': 'Pocket 1', 'mean_plddt': 90}})
    current = structures.select_pocket_by_score({'Pocket 1': {'probability': 0.5, 'mean_plddt': 90}}, 'probability')
    assert not structures._best_pockets_missing_plddt({'Q1': current})
    assert not structures._best_pockets_missing_plddt({'Q1': 'No_pockets'})


def test_full_mode_prefers_structures_with_a_passing_pocket():
    failed_model = {'Q1': {'maxDS': 0.9, 'pocket': 'Pocket 1', 'plddt_filter': 'none_passed'}}
    crystal = {'1ABC': {'maxDS': 0.5, 'pocket': 'Pocket 3', 'plddt_filter': 'not_model'}}
    assert structures.select_best_pocket_dict([failed_model, crystal]) == crystal
    assert structures.select_best_pocket_dict([failed_model]) == failed_model


def _agreement_locus(tmp_path, fpocket_residues, p2rank_residues):
    locus = tmp_path / 'PA14_00001'
    (locus / 'Q1').mkdir(parents=True)
    _model(locus / 'Q1' / 'AF_Q1.pdb', [90] * 10)
    fp = locus / 'pockets' / 'AF_Q1_fpocket' / 'pockets'
    fp.mkdir(parents=True)
    (fp / 'pocket3_atm.pdb').write_text(''.join(_atom(i, 'CA', r, 0.0) for i, r in enumerate(fpocket_residues, 1)))
    pr = locus / 'pockets' / 'AF_Q1_p2rank'
    pr.mkdir(parents=True)
    ids = ' '.join(f'A_{r}' for r in p2rank_residues)
    (pr / 'AF_Q1.pdb_predictions.csv').write_text(f'name     , probability, residue_ids\npocket1  ,       0.800, {ids}\n')
    return locus


def test_fpocket_p2rank_same_site(tmp_path):
    locus = _agreement_locus(tmp_path, [1, 2, 3, 4], [3, 4, 5, 6, 7, 8])
    result = structures.pocket_method_agreement(str(locus), 'Pocket 3', 'Pocket pocket1  ')
    assert result == {'fpocket_p2rank_shared_residues': 2, 'fpocket_p2rank_overlap': 50.0,
                      'fpocket_p2rank_same_site': True}


def test_fpocket_p2rank_different_site(tmp_path):
    locus = _agreement_locus(tmp_path, [1, 2, 3, 4], [4, 5, 6, 7, 8])
    result = structures.pocket_method_agreement(str(locus), 'Pocket 3', 'Pocket pocket1  ')
    assert result['fpocket_p2rank_shared_residues'] == 1
    assert result['fpocket_p2rank_same_site'] is False


def test_fpocket_p2rank_agreement_needs_both_pockets(tmp_path):
    locus = _agreement_locus(tmp_path, [1, 2], [1, 2])
    result = structures.pocket_method_agreement(str(locus), 'Pocket 3', 'No_pockets')
    assert set(result.values()) == {None}
