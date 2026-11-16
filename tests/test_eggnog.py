"""eggNOG-mapper annotations: parsing, COG resolution and the confidence filter."""
import pandas as pd
import pytest

import configuration
from ftscripts import eggnog

HEADER = ("#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\ttax_ceiling\tfarthest_donor_lineage\t"
          "COG_category\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\t"
          "KEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\tannotation_confidence\n")
LEGEND = ("## emapper-3.0.0-beta6\n"
          "## confidence field order: Preferred_name GOs EC KEGG_ko KEGG_Pathway KEGG_Module "
          "KEGG_Reaction KEGG_rclass BRITE KEGG_TC CAZy BiGG_Reaction PFAMs\n")


def _row(query, cog, name, ec, ko, pathway, pfam, confidence):
    return (f"{query}\t1.X\t1e-50\t200\tOG1@2|A-1,OG2@1236|B-1\tProkaryota\t1;2\t{cog}\t{name}\t-\t{ec}\t"
            f"{ko}\t{pathway}\t-\t-\t-\t-\t-\t-\t-\t{pfam}\t{confidence}\n")


def _file(tmp_path, rows, legend=LEGEND):
    path = tmp_path / "query.emapper.annotations"
    path.write_text(legend + HEADER + "".join(rows) + "## 3 queries scanned\n")
    return path


def test_cog_id_resolves_category_from_ncbi_tables(tmp_path):
    path = _file(tmp_path, [_row("G1", "COG0043", "ubiD", "ec:4.1.1.98", "K03182", "00130,01100",
                                 "UbiD_10_431", "h-hhh---h---h")])
    table = eggnog.build_table(path, ["G1", "G2"])
    g1 = table.set_index("gene").loc["G1"]
    assert g1.eggnog_cog == "COG0043"
    assert g1.eggnog_cog_category == "H"
    assert g1.eggnog_cog_category_desc == "Coenzyme transport and metabolism"
    assert g1.eggnog_cog_pathway == "Ubiquinone biosynthesis"
    assert "eggnog_ec" not in table.columns
    assert g1.eggnog_kegg_pathway == "00130"  # 01100 is a generic overview map
    assert g1.eggnog_og == "OG1@2|A-1; OG2@1236|B-1"
    # A locus without an eggNOG hit keeps a row with empty columns.
    assert table.set_index("gene").loc["G2"].isna().all()
    assert "eggnog_description" not in table.columns


def test_category_letters_from_older_emapper(tmp_path):
    path = _file(tmp_path, [_row("G1", "EH", "x", "-", "ko:K00001", "map00010", "-", "")], legend="")
    g1 = eggnog.build_table(path, ["G1"]).iloc[0]
    assert pd.isna(g1.eggnog_cog) and g1.eggnog_cog_category == "EH"
    assert g1.eggnog_cog_category_desc == "Amino acid transport and metabolism; Coenzyme transport and metabolism"
    assert g1.eggnog_kegg_ko == "K00001" and g1.eggnog_kegg_pathway == "00010"


def test_low_confidence_values_are_dropped_and_reported(tmp_path):
    # Preferred_name low, KEGG_ko medium, PFAMs high.
    path = _file(tmp_path, [_row("G1", "S", "abcD", "1.1.1.1", "K00001", "-", "PF1", "l-hm--------h")])
    g1 = eggnog.build_table(path, ["G1"]).iloc[0]
    assert pd.isna(g1.eggnog_preferred_name)
    assert g1.eggnog_low_confidence_dropped == "preferred_name"
    assert g1.eggnog_kegg_ko == "K00001" and g1.eggnog_medium_confidence_fields == "kegg_ko"
    assert g1.eggnog_pfams == "PF1"
    assert g1.eggnog_cog_category_desc == "Function unknown"


@pytest.mark.parametrize("level, kept", [("high", None), ("low", "abcD")])
def test_min_confidence_option(tmp_path, level, kept):
    path = _file(tmp_path, [_row("G1", "-", "abcD", "-", "K00001", "-", "-", "l--m---------")])
    g1 = eggnog.build_table(path, ["G1"], min_confidence=level).iloc[0]
    assert (pd.isna(g1.eggnog_preferred_name) if kept is None else g1.eggnog_preferred_name == kept)
    assert pd.isna(g1.eggnog_kegg_ko) is (level == "high")


def test_file_from_another_organism_is_rejected(tmp_path):
    path = _file(tmp_path, [_row(f"X{i}", "-", "-", "-", "-", "-", "-", "-" * 13) for i in range(3)])
    with pytest.raises(ValueError, match="this organism"):
        eggnog.build_table(path, ["G1"])


def test_config_validation(tmp_path):
    from tests.test_configuration import _minimal
    config = _minimal(tmp_path)
    config["metabolism-eggnog"] = {"enabled": True, "annotations_file": str(tmp_path / "missing")}
    with pytest.raises(ValueError, match="annotations_file not found"):
        configuration.validate_config(config)
    path = _file(tmp_path, [])
    config["metabolism-eggnog"] = {"enabled": True, "annotations_file": str(path), "min_confidence": "x"}
    with pytest.raises(ValueError, match="min_confidence"):
        configuration.validate_config(config)
    config["metabolism-eggnog"]["min_confidence"] = "high"
    configuration.validate_config(config)
    assert configuration.Config(config).metabolism_eggnog["min_confidence"] == "high"
    del config["metabolism-eggnog"]
    assert configuration.Config(config).metabolism_eggnog is False
