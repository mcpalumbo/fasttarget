import json

import pandas as pd
import pytest

from ftscripts import alphafill, pocket_overlap
from ftscripts.pocket_overlap import SequenceMapper, is_ion, run_pocket_overlap

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLI"  # 30 residues
THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
         "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
         "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


def pdb_line(serial, name, resname, chain, number, x, y=0.0, z=0.0, record="ATOM"):
    return (f"{record:<6}{serial:>5} {name:^4} {resname:>3} {chain}{number:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")


def write_structure(path, offset=100):
    """PDB chain with a 3-residue His tag and numbering shifted by offset."""
    lines, serial = [], 1
    for index, number in enumerate((offset - 5, offset - 4, offset - 3)):
        lines.append(pdb_line(serial, "CA", "HIS", "A", number, -10.0 * (index + 1)))
        serial += 1
    for position, letter in enumerate(SEQ, 1):
        lines.append(pdb_line(serial, "CA", THREE[letter], "A", position + offset, 3.8 * position))
        serial += 1
    path.write_text("".join(lines) + "END\n")


def write_pocket_atm(path, numbers, offset=100):
    path.write_text("".join(pdb_line(i, "CB", THREE[SEQ[n - 1]], "A", n + offset, 0.0)
                            for i, n in enumerate(numbers, 1)))


def write_cif(path, ligands):
    rows = []
    serial = 1
    for position, letter in enumerate(SEQ, 1):
        rows.append(f"ATOM {serial} C CA {THREE[letter]} A 1 {position} {3.8 * position:.3f} 0.000 0.000 1")
        serial += 1
    for asym, comp, coordinates in ligands:
        for x, y in coordinates:
            rows.append(f"HETATM {serial} C C1 {comp} {asym} 2 . {x:.3f} {y:.3f} 0.000 1")
            serial += 1
    header = ["data_test", "loop_"] + [f"_atom_site.{name}" for name in (
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id",
        "label_entity_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num")]
    path.write_text("\n".join(header + rows) + "\n#\n")


@pytest.fixture
def organism(tmp_path):
    base = tmp_path / "ORG"
    (base / "genome").mkdir(parents=True)
    (base / "genome" / "ORG.faa").write_text(f">G1\n{SEQ}\n")
    locus = base / "structures" / "G1"
    (locus / "P1").mkdir(parents=True)
    write_structure(locus / "P1" / "PDB_1ABC_A_ref.pdb")
    fpocket = locus / "pockets" / "PDB_1ABC_A_ref_fpocket"
    (fpocket / "pockets").mkdir(parents=True)
    (fpocket / "all_pockets.json").write_text(json.dumps({"1ABC": {
        "Pocket 1": {"Druggability Score": 0.9}, "Pocket 2": {"Druggability Score": 0.2}}}))
    write_pocket_atm(fpocket / "pockets" / "pocket1_atm.pdb", range(1, 6))
    write_pocket_atm(fpocket / "pockets" / "pocket2_atm.pdb", range(11, 17))
    p2rank = locus / "pockets" / "PDB_1ABC_A_ref_p2rank"
    p2rank.mkdir()
    residues = " ".join(f"A_{n + 100}" for n in range(11, 17))
    (p2rank / "PDB_1ABC_A_ref.pdb_predictions.csv").write_text(
        "name     ,  rank, probability, residue_ids\n"
        f"pocket1  ,     1,       0.700, {residues}\n"
        "pocket2  ,     2,       0.100, A_101 A_102\n")
    alpha = base / "ligands" / "alphafill"
    (alpha / "models" / "G1").mkdir(parents=True)
    # LIG atoms over residues 12-15: CA within 4.5 A of residues 11-16 only.
    write_cif(alpha / "models" / "G1" / "P1_alphafill.cif", [
        ("B", "LIG", [(3.8 * k, 2.0) for k in range(12, 16)]),
        ("C", "ZN", [(3.8 * 2, 2.0)]),
        ("D", "TIN", [(3.8 * 25, 2.0)]),
    ])
    pd.DataFrame([
        {"gene": "G1", "uniprot": "P1", "selection_source": "direct", "ligand_instance_id": f"AF-{asym}",
         "compound_code": comp, "asym_id": asym, "confidence_category": category,
         "cif_path": "models/G1/P1_alphafill.cif"}
        for asym, comp, category in (("B", "LIG", "high"), ("C", "ZN", "high"), ("D", "TIN", "high"))
    ]).to_csv(alpha / "alphafill_ligands_all.tsv", sep="\t", index=False)
    return tmp_path


def test_mapper_ignores_tag_and_offsets():
    mapping, identity = SequenceMapper().map("HHH" + SEQ, SEQ)
    assert identity == 100.0
    assert mapping[3] == 1 and mapping[3 + 29] == 30


def test_ion_detection():
    assert is_ion("ZN", 1) and is_ion("SO4", 5) and is_ion("XYZ", 1)
    assert not is_ion("FMN", 31)


def test_overlap_best_and_match_pockets(organism):
    summary = run_pocket_overlap(organism, "ORG")
    out = organism / "ORG" / "ligands" / "alphafill" / "pocket_overlap"
    rows = pd.read_csv(out / "alphafill_pocket_overlap.tsv", sep="\t", keep_default_na=False)
    lig = rows[rows["compound_code"] == "LIG"].iloc[0]
    assert lig["status"] == "compared"
    assert lig["contact_residues"] == "11;12;13;14;15;16"
    # Best druggability pocket does not hold the ligand; pocket 2 (rank 2) does.
    assert lig["fpocket_best_pocket"] == "PDB_1ABC_A_ref/Pocket 1"
    assert float(lig["fpocket_best_ligand_overlap"]) == 0.0
    assert lig["fpocket_match_pocket"] == "PDB_1ABC_A_ref/Pocket 2"
    assert float(lig["fpocket_match_score"]) == 0.2
    assert int(float(lig["fpocket_match_score_rank"])) == 2
    assert float(lig["fpocket_match_ligand_overlap"]) == 100.0
    assert float(lig["fpocket_match_pocket_overlap"]) == 100.0
    assert lig["fpocket_match_is_best"] == "False"
    # P2Rank residue IDs use PDB numbering (+100) and map to the same positions.
    assert lig["p2rank_best_pocket"] == "PDB_1ABC_A_ref/pocket1"
    assert float(lig["p2rank_best_ligand_overlap"]) == 100.0
    assert rows.set_index("compound_code").loc["ZN", "status"] == "excluded_ion"
    assert rows.set_index("compound_code").loc["TIN", "status"] == "excluded_ion"

    row = summary.set_index("gene").loc["G1"]
    assert row["alphafill_pocket_fpocket_best_pocket"] == "PDB_1ABC_A_ref/Pocket 1"
    assert row["alphafill_pocket_fpocket_best_max_ligand_overlap"] == 0.0
    assert row["alphafill_pocket_fpocket_best_compounds"] == ""
    assert row["alphafill_pocket_fpocket_match_pocket"] == "PDB_1ABC_A_ref/Pocket 2"
    assert row["alphafill_pocket_fpocket_match_compounds"] == "LIG"
    assert row["alphafill_pocket_p2rank_best_compounds"] == "LIG"


def test_few_contacts_are_reported_but_not_summarized(organism):
    summary = run_pocket_overlap(organism, "ORG", min_contact_residues=10)
    rows = pd.read_csv(organism / "ORG" / "ligands" / "alphafill" / "pocket_overlap" /
                       "alphafill_pocket_overlap.tsv", sep="\t")
    lig = rows[rows["compound_code"] == "LIG"].iloc[0]
    assert lig["status"] == "few_contacts"
    assert lig["fpocket_match_ligand_overlap"] == 100.0
    assert summary.set_index("gene").loc["G1", "alphafill_pocket_fpocket_match_pocket"] is None


def test_keep_ions_option(organism):
    run_pocket_overlap(organism, "ORG", exclude_ions=False, min_contact_residues=1)
    rows = pd.read_csv(organism / "ORG" / "ligands" / "alphafill" / "pocket_overlap" /
                       "alphafill_pocket_overlap.tsv", sep="\t")
    assert rows.set_index("compound_code").loc["ZN", "status"] in {"compared", "few_contacts"}


def test_mismatched_model_is_rejected(organism):
    faa = organism / "ORG" / "genome" / "ORG.faa"
    faa.write_text(">G1\nWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n")
    run_pocket_overlap(organism, "ORG")
    rows = pd.read_csv(organism / "ORG" / "ligands" / "alphafill" / "pocket_overlap" /
                       "alphafill_pocket_overlap.tsv", sep="\t")
    assert set(rows["status"]) == {"model_mapping_failed"}


@pytest.mark.parametrize("key, value", [
    ("pocket_overlap", "yes"), ("pocket_contact_cutoff", 0),
    ("pocket_min_overlap", 150), ("pocket_min_contact_residues", 0),
])
def test_invalid_pocket_options(key, value):
    with pytest.raises(ValueError):
        alphafill.validate_options({key: value})


def test_cli(organism, capsys):
    pocket_overlap.main(["--output_path", str(organism), "--organism", "ORG"])
    assert "1 loci written" in capsys.readouterr().out


def test_ligands_are_read_by_gene_in_chunks_with_only_the_needed_columns(tmp_path):
    """Chunk boundaries must not split a gene, and the heavy JSON columns are skipped."""
    rows = [{"gene": f"G{index // 3:02d}", "uniprot": "P1", "selection_source": "direct",
             "ligand_instance_id": f"AF-{index:03d}", "compound_code": "ATP", "asym_id": "B",
             "confidence_category": "high", "cif_path": "", "hit_json": "{" + "x" * 50 + "}"}
            for index in range(20)]
    path = tmp_path / "ligands.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    chunked = list(pocket_overlap._ligands_by_gene(path, chunksize=4))

    assert [gene for gene, _ in chunked] == sorted({row["gene"] for row in rows})
    assert sum(len(group) for _, group in chunked) == len(rows)
    for gene, group in chunked:
        assert set(group["gene"]) == {gene}
        assert "hit_json" not in group.columns


def test_an_unsorted_ligand_table_still_keeps_every_row(tmp_path):
    genes = ["B", "A", "B", "C", "A"]
    rows = [{column: "" for column in pocket_overlap.LIGAND_INPUT_COLUMNS} | {"gene": gene, "ligand_instance_id": str(i)}
            for i, gene in enumerate(genes)]
    path = tmp_path / "ligands.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    chunked = list(pocket_overlap._ligands_by_gene(path, chunksize=2))

    counts = {}
    for gene, group in chunked:
        counts[gene] = counts.get(gene, 0) + len(group)
    assert counts == {"A": 2, "B": 2, "C": 1}


def _forbid_recompute(monkeypatch):
    from ftscripts import pocket_overlap

    def fail(*args, **kwargs):
        raise AssertionError("pocket overlap was recomputed")
    monkeypatch.setattr(pocket_overlap, "load_locus_pockets", fail)


def test_unchanged_inputs_reuse_previous_tables(organism, monkeypatch):
    first = run_pocket_overlap(organism, "ORG")
    _forbid_recompute(monkeypatch)
    again = run_pocket_overlap(organism, "ORG")
    assert list(again.columns) == list(first.columns)
    assert again["gene"].tolist() == first["gene"].tolist()


@pytest.mark.parametrize("change", ["pocket_file", "parameter", "missing_output"])
def test_changed_inputs_recompute(organism, monkeypatch, change):
    import os
    run_pocket_overlap(organism, "ORG")
    out = organism / "ORG" / "ligands" / "alphafill" / "pocket_overlap"
    kwargs = {}
    if change == "pocket_file":
        atm = (organism / "ORG" / "structures" / "G1" / "pockets" / "PDB_1ABC_A_ref_fpocket"
               / "pockets" / "pocket1_atm.pdb")
        st = atm.stat()
        os.utime(atm, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000))
    elif change == "parameter":
        kwargs = {"min_overlap": 10.0}
    else:
        (out / "alphafill_pocket_overlap.tsv").unlink()
    calls = []
    from ftscripts import pocket_overlap
    original = pocket_overlap.load_locus_pockets
    monkeypatch.setattr(pocket_overlap, "load_locus_pockets",
                        lambda *a, **k: calls.append(1) or original(*a, **k))
    run_pocket_overlap(organism, "ORG", **kwargs)
    assert calls
    assert (out / "alphafill_pocket_overlap.tsv").is_file()


def test_best_pocket_follows_structures_module(tmp_path):
    """pocket_overlap's best is the pocket the structures module chose (pLDDT filter)."""
    pockets_dir = tmp_path / "PA14_00001" / "pockets"
    pockets_dir.mkdir(parents=True)
    (pockets_dir / "best_DS_fpocket_PA14_00001.json").write_text(
        json.dumps({"Q1": {"maxDS": 0.4, "pocket": "Pocket 2", "plddt_filter": "passed"}}))
    (pockets_dir / "best_prob_p2rank_PA14_00001.json").write_text(
        json.dumps({"Q1": {"max_probability": 0.5, "pocket": "Pocket pocket2  "}}))
    assert pocket_overlap._pipeline_best(pockets_dir, "fpocket") == ("Q1", "2")
    assert pocket_overlap._pipeline_best(pockets_dir, "p2rank") == ("Q1", "pocket2")
    assert pocket_overlap._pipeline_structure_key("AF_Q1") == "Q1"
    assert pocket_overlap._pipeline_structure_key("PDB_1ABC_A_ref") == "1ABC"

    top = {"id": "AF_Q1/Pocket 1", "score": 0.9, "rank": 1, "residues": {1, 2, 3}}
    chosen = {"id": "AF_Q1/Pocket 2", "score": 0.4, "rank": 2, "residues": {7, 8, 9}}
    row = pocket_overlap.compare({7, 8, 9}, {"status": "ok", "pockets": [top, chosen], "best": chosen}, "fpocket")
    assert row["fpocket_best_pocket"] == "AF_Q1/Pocket 2"
    assert row["fpocket_match_is_best"] is True


def test_ligand_classes():
    from ftscripts.pocket_overlap import ligand_class
    assert ligand_class("ZN", 1) == "ion"
    assert ligand_class("BF2", 3) == "ion"
    assert ligand_class("HOH", 1) == "ion"          # one heavy atom
    for code in ("GOL", "EDO", "PEG", "1PG", "C8E", "CO2", "9JE", "DMS"):
        assert ligand_class(code, 6) == "artifact", code
    for code in ("HEM", "SF4", "FES", "MGD", "B12"):
        assert ligand_class(code, 20) == "prosthetic", code
    for code in ("ATP", "NAD", "FAD", "FMN", "COA", "SAM", "PLP"):
        assert ligand_class(code, 20) == "cofactor", code
    for code in ("NOV", "ICB", "TCA", "CB3", "LIG"):
        assert ligand_class(code, 12) == "other", code


def _add_ligands(organism, ligands):
    """Add (asym, comp, x positions) ligands to the fixture model and ligand table."""
    alpha = organism / "ORG" / "ligands" / "alphafill"
    base = [("B", "LIG", [(3.8 * k, 2.0) for k in range(12, 16)]),
            ("C", "ZN", [(3.8 * 2, 2.0)]), ("D", "TIN", [(3.8 * 25, 2.0)])]
    extra = [(asym, comp, [(3.8 * k, 2.0) for k in range(12, 16)]) for asym, comp in ligands]
    write_cif(alpha / "models" / "G1" / "P1_alphafill.cif", base + extra)
    table = pd.read_csv(alpha / "alphafill_ligands_all.tsv", sep="\t")
    rows = [{"gene": "G1", "uniprot": "P1", "selection_source": "direct",
             "ligand_instance_id": f"AF-{asym}", "compound_code": comp, "asym_id": asym,
             "confidence_category": "high", "cif_path": "models/G1/P1_alphafill.cif"}
            for asym, comp in ligands]
    pd.concat([table, pd.DataFrame(rows)]).to_csv(alpha / "alphafill_ligands_all.tsv", sep="\t", index=False)


def test_artifacts_and_prosthetic_groups_are_not_compared(organism):
    _add_ligands(organism, [("E", "GOL"), ("F", "HEM"), ("G", "ATP")])
    summary = run_pocket_overlap(organism, "ORG")
    rows = pd.read_csv(organism / "ORG" / "ligands" / "alphafill" / "pocket_overlap" /
                       "alphafill_pocket_overlap.tsv", sep="\t", keep_default_na=False).set_index("compound_code")
    assert rows.loc["GOL", "status"] == "excluded_artifact"
    assert rows.loc["HEM", "status"] == "excluded_prosthetic"
    assert rows.loc["ATP", "status"] == "compared" and rows.loc["ATP", "ligand_class"] == "cofactor"
    assert rows.loc["LIG", "ligand_class"] == "other"
    row = summary.set_index("gene").loc["G1"]
    assert row["alphafill_pocket_p2rank_best_compounds"] == "ATP;LIG"
    assert row["alphafill_pocket_p2rank_best_non_cofactor_compounds"] == "LIG"
    assert row["alphafill_pocket_fpocket_match_non_cofactor_compounds"] == "LIG"


def test_keep_artifacts_and_prosthetic_options(organism):
    _add_ligands(organism, [("E", "GOL"), ("F", "HEM")])
    run_pocket_overlap(organism, "ORG", exclude_artifacts=False, exclude_prosthetic=False)
    rows = pd.read_csv(organism / "ORG" / "ligands" / "alphafill" / "pocket_overlap" /
                       "alphafill_pocket_overlap.tsv", sep="\t").set_index("compound_code")
    assert rows.loc["GOL", "status"] == "compared"
    assert rows.loc["HEM", "status"] == "compared"
