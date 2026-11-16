"""Compare AlphaFill ligand placements with fpocket/P2Rank pockets.

Pockets and AlphaFill models are often different structures (a PDB chain vs an
AlphaFold model, or a >=95% identical UniProt fallback). Coordinates are
therefore never compared across files. Instead, both sides are mapped onto the
locus protein sequence (genome/<organism>.faa) by sequence alignment:

* ligand contact residues: protein residues of the AlphaFill model within
  ``contact_cutoff`` Angstrom of any ligand heavy atom (same file, same frame);
* pocket residues: residues lining each fpocket pocket (pocketN_atm.pdb) or
  listed by P2Rank (residue_ids in *_predictions.csv).

For every ligand and method two pockets are reported:

* best: the best pocket the structures module selected (best_DS_fpocket_<locus>.json,
  best_prob_p2rank_<locus>.json): the highest score (fpocket Druggability Score,
  P2Rank probability) among pockets of predicted models with mean pLDDT >= 70.
  Without those files, the highest score across all structures of the locus;
* match: the pocket that best coincides with the ligand (most shared residues),
  whatever its score. It is not the best pocket unless match_is_best is True.

ligand_overlap = shared / ligand contact residues (%), pocket_overlap =
shared / pocket residues (%). Ranking by shared residues is the same as
ranking by ligand_overlap, since the ligand's contact set is fixed.

Ligands are classified (ligand_class) and some are never compared (status
excluded_<class>): ions (known codes or single heavy atom), tiny and easily false
positive matches; crystallization artifacts and noise (water, CO2, PEG, glycerol,
detergents, buffers: BioLiP list plus additions); prosthetic groups that stay bound
in the structure (heme, Fe-S clusters...), whose cavity is no drug site. Mobile
cofactors (ATP, NAD, FAD, CoA...) are compared and flagged. Ligands with
fewer than min_contact_residues contacts are compared (status few_contacts)
but left out of the locus summary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1


METHODS = {"fpocket": "Druggability Score", "p2rank": "probability"}
# Monoatomic and small inorganic ions. Any ligand with a single heavy atom is
# also treated as an ion, so unlisted metal/halide codes are excluded too.
ION_CODES = frozenset({
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "MN3", "FE", "FE2",
    "CO", "3CO", "NI", "3NI", "CU", "CU1", "ZN", "CD", "HG", "PB", "AL", "GA", "AG",
    "AU", "PT", "TL", "YB", "SM", "GD", "LU", "TB", "EU", "HO", "ER", "LA", "CE",
    "PR", "ND", "Y1", "SB", "CR", "V", "W", "MO", "CL", "BR", "IOD", "F", "OH",
    "SO4", "SUL", "PO4", "PI", "2HP", "NO3", "NO2", "CO3", "BCT", "NH4", "SCN",
    "AZI", "CYN", "WO4", "MOO", "VO4", "SE4", "PER", "IUM", "OXL",
})
# Ligand classes. Only "cofactor" and "other" ligands are compared with pockets:
# ions are too small (false positive matches), crystallization artifacts and
# noise (water, CO2, PEG, glycerol, detergents, buffers...) say nothing about a
# binding site, and prosthetic groups (heme, Fe-S clusters...) stay bound in the
# structure, so their cavity is not a site for a drug. Mobile cofactors (ATP,
# NAD, FAD, CoA...) are compared and flagged: their sites can be competed for.
ARTIFACTS_FILE = Path(__file__).resolve().parent / "data" / "biolip_artifacts.tsv"
# Not in the BioLiP list but artifacts or noise in AlphaFill transplants.
EXTRA_ARTIFACT_CODES = frozenset({
    "C8E", "9JE", "DXC", "OXY", "O2", "H2S", "NO", "CMO", "PEO", "UNL", "UNX", "DMU",
    "P33", "PE3", "PE8", "PG0", "PG5", "PG6", "PGF", "XPE", "7PE", "BU3", "CPS",
})
# In the BioLiP list but genuine ligands: inhibitors and drugs, and metabolites that
# are substrates or products of bacterial enzymes (they mark a real binding site).
# Buffers and crystallization salts with a metabolic role (citrate, tartrate,
# malonate, urea, pyrophosphate) stay excluded: in crystals they are mostly additives.
ARTIFACT_EXCEPTIONS = frozenset({
    # inhibitors, drugs
    "CB3", "HAE", "OXM",
    # amino acids, polyamines, osmolytes
    "HCS", "PUT", "SPD", "SPM", "SPK", "BET", "CHT", "TAU", "KAI",
    # central metabolism
    "SIN", "MLT", "LMR", "ICT", "HCA", "LAC", "2OP", "GLV", "3PG", "13P", "MAE",
    # nucleotides, bases
    "ORO", "GUN", "DUD",
    # aromatic acids (substrates of decarboxylases, hydroxylases, CoA ligases)
    "BEZ", "PHB", "DHB", "SAL", "HC4", "TCA", "IPH",
})
# Anion mimics and small polyatomic ions, treated as ions.
EXTRA_ION_CODES = frozenset({"TSU", "SO3", "BF2", "BEF", "ALF", "AF3", "MGF", "VO3", "FLL"})
PROSTHETIC_CODES = frozenset({
    # hemes and porphyrins
    "HEM", "HEC", "HEA", "HEB", "HDD", "HAS", "HEO", "HEV", "HDM", "DHE", "SRM", "HNI", "CYC",
    "COH", "ZNH", "MHM", "1FH", "HCO", "CLA", "CHL", "BCL", "BCB", "BPH", "PHO", "F43",
    # iron-sulfur and other metal clusters
    "SF4", "FES", "F3S", "FS4", "SF3", "CLF", "ICS", "CFM", "CFN", "HCN", "NFS", "XCC", "WCC",
    # molybdopterin, cobalamin
    "MGD", "MTE", "PCD", "MCN", "MTV", "MSS", "B12", "CBL", "COB", "CNC", "5AD",
    # covalently attached groups
    "LPA", "LPB", "BTN", "BTI", "PNS", "RET", "TPQ", "TOQ", "LTQ", "CTQ",
})
COFACTOR_CODES = frozenset({
    "ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "APC", "GTP", "GDP", "GNP", "GCP", "GSP",
    "GMP", "CTP", "CDP", "CMP", "UTP", "UDP", "U5P", "TTP", "DTP", "DGT", "DCP", "TMP",
    "COA", "ACO", "SCA", "MLC", "CAA", "HMG", "MCA", "CO8", "MYA", "LCV",
    "NAD", "NAI", "NAP", "NDP", "NMN", "NAJ", "NAE", "FAD", "FDA", "FAE", "FMN", "FNR",
    "BYN", "4LU", "4LS", "RBF", "PLP", "PMP", "PLR", "SAM", "SAH", "SFG", "TPP", "TDP",
    "TPW", "THF", "THG", "FOL", "DHF", "MEF", "FON", "H4B", "BH4", "GSH", "GDS", "PQQ",
    "UPG", "UD1", "A3P", "PAP", "PPS", "MTA", "ADN", "CO7", "NA7",
})


def _load_artifact_codes() -> frozenset[str]:
    codes = set()
    try:
        for line in ARTIFACTS_FILE.read_text(encoding="utf-8").splitlines():
            if line and not line.startswith("#"):
                codes.add(line.split("\t", 1)[0].strip().upper())
    except OSError as error:
        logging.warning("Pocket overlap: cannot read %s (%s); only the built-in artifact codes are used",
                        ARTIFACTS_FILE, error)
    return frozenset((codes | EXTRA_ARTIFACT_CODES) - ARTIFACT_EXCEPTIONS)


ARTIFACT_CODES = _load_artifact_codes()

NON_STANDARD = {"MSE": "M", "SEC": "U", "PYL": "O", "CSO": "C", "CME": "C",
                "KCX": "K", "LLP": "K", "MLY": "K", "TPO": "T", "SEP": "S", "PTR": "Y"}
OUTPUT_DIR = "pocket_overlap"
LIGAND_FILE = "alphafill_pocket_overlap.tsv"
SUMMARY_FILE = "alphafill_pocket_overlap_summary.tsv"
MANIFEST_FILE = "pocket_overlap_manifest.json"


def _method_columns(method: str) -> list[str]:
    return [f"{method}_status", f"{method}_best_pocket", f"{method}_best_score",
            f"{method}_best_shared_residues", f"{method}_best_ligand_overlap",
            f"{method}_best_pocket_overlap", f"{method}_match_pocket", f"{method}_match_score",
            f"{method}_match_score_rank", f"{method}_match_shared_residues",
            f"{method}_match_ligand_overlap", f"{method}_match_pocket_overlap",
            f"{method}_match_is_best"]


LIGAND_COLUMNS = [
    "gene", "uniprot", "selection_source", "ligand_instance_id", "compound_code",
    "asym_id", "confidence_category", "ligand_class", "status", "model_mapping_identity",
    "contact_residue_count", "contact_residues",
    *(column for method in METHODS for column in _method_columns(method)),
]


def _summary_columns(method: str) -> list[str]:
    prefix = f"alphafill_pocket_{method}"
    return [f"{prefix}_best_pocket", f"{prefix}_best_score",
            f"{prefix}_best_max_ligand_overlap", f"{prefix}_best_compounds",
            f"{prefix}_best_non_cofactor_compounds",
            f"{prefix}_match_pocket", f"{prefix}_match_score", f"{prefix}_match_ligand_overlap",
            f"{prefix}_match_compounds", f"{prefix}_match_non_cofactor_compounds"]


SUMMARY_COLUMNS = ["gene", "alphafill_pocket_status",
                   *(column for method in METHODS for column in _summary_columns(method))]


def residue_letter(name: str) -> str:
    name = name.strip().upper()
    if name in NON_STANDARD:
        return NON_STANDARD[name]
    letter = seq1(name)
    return letter if letter and letter != "X" else "X"


def is_ion(compound: str, heavy_atoms: int) -> bool:
    code = compound.strip().upper()
    return code in ION_CODES or code in EXTRA_ION_CODES or heavy_atoms <= 1


def ligand_class(compound: str, heavy_atoms: int) -> str:
    """'ion', 'prosthetic', 'cofactor', 'artifact' or 'other' (see the code lists above)."""
    code = compound.strip().upper()
    if is_ion(code, heavy_atoms):
        return "ion"
    if code in PROSTHETIC_CODES:
        return "prosthetic"
    if code in COFACTOR_CODES:
        return "cofactor"
    if code in ARTIFACT_CODES:
        return "artifact"
    return "other"


class SequenceMapper:
    """Map residue indices of one sequence onto the locus sequence (1-based)."""

    def __init__(self) -> None:
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        self.aligner = aligner

    def map(self, source: str, target: str) -> tuple[dict[int, int], float]:
        """Return {source index: target position (1-based)} and identity (%)."""
        source = source.replace("X", "A") if source else source
        if not source or not target:
            return {}, 0.0
        alignment = self.aligner.align(source, target)[0]
        mapping, identical = {}, 0
        for (s_start, s_end), (t_start, t_end) in zip(*alignment.aligned):
            for offset in range(s_end - s_start):
                mapping[s_start + offset] = t_start + offset + 1
                identical += source[s_start + offset] == target[t_start + offset]
        identity = 100.0 * identical / len(mapping) if mapping else 0.0
        return mapping, identity


# ---------------------------------------------------------------- structures

def _residue_key(chain: str, number: Any, icode: str = "") -> tuple[str, int, str]:
    return chain.strip(), int(number), (icode or "").strip()


def structure_residue_map(pdb_path: Path, locus_seq: str, mapper: SequenceMapper,
                          min_identity: float) -> dict[tuple[str, int, str], int]:
    """Map (chain, resseq, icode) of a pocket structure to locus positions.

    Each chain is aligned separately; chains of other proteins (complex
    partners) fall below min_identity and are ignored.
    """
    structure = PDBParser(QUIET=True).get_structure("s", str(pdb_path))
    result: dict[tuple[str, int, str], int] = {}
    model = next(iter(structure), None)
    if model is None:
        return result
    for chain in model:
        keys, letters = [], []
        for residue in chain:
            if not is_aa(residue, standard=False):
                continue
            hetflag, number, icode = residue.id
            keys.append(_residue_key(chain.id, number, icode))
            letters.append(residue_letter(residue.get_resname()))
        if len(keys) < 10:
            continue
        mapping, identity = mapper.map("".join(letters), locus_seq)
        if identity < min_identity:
            continue
        for index, position in mapping.items():
            result[keys[index]] = position
    return result


def _pocket_atm_residues(path: Path) -> set[tuple[str, int, str]]:
    residues = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 27:
                try:
                    residues.add(_residue_key(line[21], line[22:26], line[26]))
                except ValueError:
                    continue
    return residues


def _fpocket_scores(folder: Path, name: str) -> dict[str, float]:
    data = None
    json_path = folder / "all_pockets.json"
    if json_path.is_file():
        try:
            loaded = json.loads(json_path.read_text(encoding="utf-8"))
            data = next(iter(loaded.values())) if isinstance(loaded, dict) and loaded else None
        except (OSError, ValueError):
            data = None
    if not isinstance(data, dict):
        data, current = {}, None
        info = folder / f"{name}_info.txt"
        if info.is_file():
            for line in info.read_text(encoding="utf-8", errors="replace").splitlines():
                match = re.match(r"^Pocket (\d+) :", line)
                if match:
                    current = f"Pocket {match.group(1)}"
                    data[current] = {}
                elif current and "Druggability Score" in line:
                    data[current]["Druggability Score"] = float(line.split(":")[1])
    scores = {}
    for label, props in data.items():
        if isinstance(props, dict) and props.get("Druggability Score") is not None:
            scores[label] = float(props["Druggability Score"])
    return scores


def _parse_p2rank_residue(token: str) -> tuple[str, int, str] | None:
    chain, _, number = token.strip().rpartition("_")
    match = re.fullmatch(r"(-?\d+)([A-Za-z]?)", number)
    if not chain or not match:
        return None
    return _residue_key(chain, match.group(1), match.group(2))


def _find_input_structure(locus_dir: Path, name: str, folder: Path) -> Path | None:
    for path in sorted(locus_dir.rglob(f"{name}.pdb")):
        if "pockets" not in path.relative_to(locus_dir).parts:
            return path
    fallback = folder / f"{name}_out.pdb"
    return fallback if fallback.is_file() else None


BEST_POCKET_FILES = {"fpocket": "best_DS_fpocket_{gene}.json", "p2rank": "best_prob_p2rank_{gene}.json"}


def _pipeline_structure_key(name: str) -> str:
    """Key the structures module gives a pocket folder (see pockets_data_to_dict)."""
    parts = name.split("_")
    if parts[0] == "PDB" and len(parts) >= 2:
        return parts[1]
    if parts[0] == "AF" and len(parts) >= 2:
        return "_".join(parts[1:])
    if parts[0] == "CB" and len(parts) >= 2:
        return name.rsplit("_", 1)[0]
    return name


def _pocket_label_key(label: str) -> str:
    """'Pocket 3' -> '3', 'Pocket pocket1  ' / 'pocket1' -> 'pocket1'."""
    return re.sub(r"^Pocket\s+", "", str(label).strip()).strip()


def _pipeline_best(pockets_dir: Path, method: str) -> tuple[str, str] | None:
    """(structure key, pocket label key) of the best pocket the structures module chose.

    That choice applies the pLDDT filter of predicted models, so it is not always
    the top-scoring pocket. None when the file is missing or names no pocket.
    """
    path = pockets_dir / BEST_POCKET_FILES[method].format(gene=pockets_dir.parent.name)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict):
        return None
    for key, props in data.items():
        if isinstance(props, dict) and props.get("pocket"):
            return str(key), _pocket_label_key(props["pocket"])
    return None


def load_locus_pockets(locus_dir: Path, locus_seq: str, mapper: SequenceMapper,
                       min_identity: float) -> dict[str, dict[str, Any]]:
    """Return, per method, a status, the pockets in locus numbering and the best one.

    The best pocket is the one the structures module selected (best_*.json, which
    applies the pLDDT filter); without that file, the top-scoring pocket.
    """
    pockets_dir = locus_dir / "pockets"
    result = {method: {"status": "no_pocket_data", "pockets": []} for method in METHODS}
    if not pockets_dir.is_dir():
        return result
    structure_maps: dict[Path, dict] = {}

    def residue_map(path: Path) -> dict:
        if path not in structure_maps:
            try:
                structure_maps[path] = structure_residue_map(path, locus_seq, mapper, min_identity)
            except Exception as error:  # malformed structures must not stop the run
                logging.warning("Pocket overlap: cannot read %s: %s", path, error)
                structure_maps[path] = {}
        return structure_maps[path]

    for method, suffix in (("fpocket", "_fpocket"), ("p2rank", "_p2rank")):
        folders = sorted(path for path in pockets_dir.iterdir()
                         if path.is_dir() and path.name.endswith(suffix))
        if not folders:
            continue
        pockets, unmapped = [], 0
        for folder in folders:
            name = folder.name[: -len(suffix)]
            structure = _find_input_structure(locus_dir, name, folder)
            mapping = residue_map(structure) if structure else {}
            if not mapping:
                unmapped += 1
                continue
            raw: list[tuple[str, float, set]] = []
            if method == "fpocket":
                scores = _fpocket_scores(folder, name)
                for label, score in scores.items():
                    number = label.split()[-1]
                    atm = folder / "pockets" / f"pocket{number}_atm.pdb"
                    if atm.is_file():
                        raw.append((label, score, _pocket_atm_residues(atm)))
            else:
                for csv_path in folder.glob("*_predictions.csv"):
                    frame = pd.read_csv(csv_path, skipinitialspace=True)
                    frame.columns = frame.columns.str.strip()
                    for _, row in frame.iterrows():
                        tokens = str(row.get("residue_ids", "")).split()
                        keys = {key for key in map(_parse_p2rank_residue, tokens) if key}
                        raw.append((str(row["name"]).strip(), float(row["probability"]), keys))
            for label, score, keys in raw:
                pockets.append({"id": f"{name}/{label}", "score": score,
                                "structure_key": _pipeline_structure_key(name),
                                "label_key": _pocket_label_key(label),
                                "residues": {mapping[key] for key in keys if key in mapping}})
        if pockets:
            ordered = sorted(pockets, key=lambda item: -item["score"])
            for rank, pocket in enumerate(ordered, 1):
                pocket["rank"] = rank
            chosen = _pipeline_best(pockets_dir, method)
            best = next((pocket for pocket in ordered
                         if chosen and (pocket["structure_key"], pocket["label_key"]) == chosen), ordered[0])
            result[method] = {"status": "ok", "pockets": ordered, "best": best}
        elif unmapped:
            result[method] = {"status": "structure_mapping_failed", "pockets": []}
        else:
            result[method] = {"status": "no_pockets", "pockets": []}
    return result


# ---------------------------------------------------------------- AlphaFill

def read_alphafill_model(cif_path: Path) -> dict[str, Any]:
    """Return protein residues/coordinates and heavy-atom ligands by label_asym_id."""
    data = MMCIF2Dict(str(cif_path))
    columns = ("_atom_site.group_PDB", "_atom_site.type_symbol", "_atom_site.label_comp_id",
               "_atom_site.label_asym_id", "_atom_site.label_seq_id",
               "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z")
    group, element, comp, asym, seq_id, x, y, z = (data[column] for column in columns)
    models = data.get("_atom_site.pdbx_PDB_model_num")
    first_model = models[0] if models else None
    protein_xyz, protein_seq_ids, residues = [], [], {}
    ligands: dict[str, dict[str, Any]] = {}
    for index in range(len(group)):
        if models and models[index] != first_model:
            continue
        if element[index].upper() in {"H", "D"}:
            continue
        xyz = (float(x[index]), float(y[index]), float(z[index]))
        if group[index] == "ATOM" and seq_id[index] not in {".", "?"}:
            number = int(seq_id[index])
            protein_xyz.append(xyz)
            protein_seq_ids.append(number)
            residues.setdefault(number, residue_letter(comp[index]))
        else:
            ligand = ligands.setdefault(asym[index], {"compound": comp[index], "xyz": []})
            ligand["xyz"].append(xyz)
    numbers = sorted(residues)
    return {"sequence": "".join(residues[number] for number in numbers),
            "seq_ids": numbers,
            "protein_xyz": np.asarray(protein_xyz, dtype=float).reshape(-1, 3),
            "protein_seq_ids": np.asarray(protein_seq_ids, dtype=int),
            "ligands": ligands}


def contact_residues(model: dict[str, Any], ligand_xyz: list, cutoff: float) -> set[int]:
    """label_seq_id of protein residues with any heavy atom within cutoff."""
    if not len(ligand_xyz) or not len(model["protein_xyz"]):
        return set()
    ligand = np.asarray(ligand_xyz, dtype=float)
    protein = model["protein_xyz"]
    # Bounding box prefilter keeps the distance matrix small for big models.
    low, high = ligand.min(axis=0) - cutoff, ligand.max(axis=0) + cutoff
    inside = np.all((protein >= low) & (protein <= high), axis=1)
    if not inside.any():
        return set()
    distances = np.linalg.norm(protein[inside][:, None, :] - ligand[None, :, :], axis=2)
    close = distances.min(axis=1) <= cutoff
    return {int(number) for number in model["protein_seq_ids"][inside][close]}


# ---------------------------------------------------------------- overlap

def _percent(numerator: int, denominator: int) -> float | None:
    return round(100.0 * numerator / denominator, 1) if denominator else None


def compare(contacts: set[int], method_data: dict[str, Any], method: str) -> dict[str, Any]:
    row: dict[str, Any] = {column: None for column in _method_columns(method)}
    row[f"{method}_status"] = method_data["status"]
    pockets = method_data["pockets"]
    if method_data["status"] != "ok" or not pockets:
        return row
    best = method_data.get("best") or pockets[0]

    def describe(prefix: str, pocket: dict[str, Any]) -> None:
        shared = len(contacts & pocket["residues"])
        row[f"{method}_{prefix}_pocket"] = pocket["id"]
        row[f"{method}_{prefix}_score"] = pocket["score"]
        row[f"{method}_{prefix}_shared_residues"] = shared
        row[f"{method}_{prefix}_ligand_overlap"] = _percent(shared, len(contacts))
        row[f"{method}_{prefix}_pocket_overlap"] = _percent(shared, len(pocket["residues"]))

    describe("best", best)
    match = max(pockets, key=lambda pocket: (
        len(contacts & pocket["residues"]),
        len(contacts & pocket["residues"]) / max(1, len(pocket["residues"])),
        pocket["score"], -pocket["rank"]))
    if contacts & match["residues"]:
        describe("match", match)
        row[f"{method}_match_score_rank"] = match["rank"]
        row[f"{method}_match_is_best"] = match is best
        row[f"{method}_status"] = "ok"
    else:
        row[f"{method}_status"] = "no_overlap"
    return row


def _compound_list(rows: pd.DataFrame, non_cofactor: bool = False) -> str:
    if non_cofactor and "ligand_class" in rows:
        rows = rows[rows["ligand_class"] != "cofactor"]
    return ";".join(sorted(set(rows["compound_code"].astype(str).str.upper())))


def summarize(ligand_rows: pd.DataFrame, min_overlap: float) -> pd.DataFrame:
    """One row per locus from accepted (high/intermediate) ligands that were compared.

    Ions, artifacts and prosthetic groups are never compared. *_compounds lists
    every compound in the pocket, mobile cofactors included; *_non_cofactor_compounds
    leaves cofactors out (substrates, inhibitors, other ligands).
    """
    rows = []
    for gene, group in ligand_rows.groupby("gene", sort=True):
        compared = group[group["status"] == "compared"]
        accepted = compared[compared["confidence_category"].isin(["high", "intermediate"])]
        row: dict[str, Any] = {"gene": gene,
                               "alphafill_pocket_status": "compared" if len(accepted)
                               else "no_accepted_ligands" if len(compared)
                               else ";".join(sorted(set(group["status"])))}
        for method in METHODS:
            prefix = f"alphafill_pocket_{method}"
            row.update({column: None for column in _summary_columns(method)})
            known = group[group[f"{method}_best_pocket"].notna()]
            if len(known):
                row[f"{prefix}_best_pocket"] = known[f"{method}_best_pocket"].iloc[0]
                row[f"{prefix}_best_score"] = known[f"{method}_best_score"].iloc[0]
            if not len(accepted):
                continue
            best_overlap = accepted[f"{method}_best_ligand_overlap"].astype(float)
            if best_overlap.notna().any():
                row[f"{prefix}_best_max_ligand_overlap"] = best_overlap.max()
                hits = accepted[best_overlap >= min_overlap]
                row[f"{prefix}_best_compounds"] = _compound_list(hits)
                row[f"{prefix}_best_non_cofactor_compounds"] = _compound_list(hits, non_cofactor=True)
            match_overlap = accepted[f"{method}_match_ligand_overlap"].astype(float)
            if match_overlap.notna().any():
                leader = accepted.loc[match_overlap.idxmax()]
                pocket = leader[f"{method}_match_pocket"]
                row[f"{prefix}_match_pocket"] = pocket
                row[f"{prefix}_match_score"] = leader[f"{method}_match_score"]
                row[f"{prefix}_match_ligand_overlap"] = leader[f"{method}_match_ligand_overlap"]
                in_match = accepted[(accepted[f"{method}_match_pocket"] == pocket)
                                  & (match_overlap >= min_overlap)]
                row[f"{prefix}_match_compounds"] = _compound_list(in_match)
                row[f"{prefix}_match_non_cofactor_compounds"] = _compound_list(in_match, non_cofactor=True)
        rows.append(row)
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


# Only these columns of alphafill_ligands_all.tsv are needed. The donor/hit/transplant
# JSON columns make up most of the file (hundreds of MB) and are never read here.
LIGAND_INPUT_COLUMNS = ["gene", "uniprot", "selection_source", "ligand_instance_id",
                        "compound_code", "asym_id", "confidence_category", "cif_path"]


def _ligands_by_gene(ligands_path: Path, chunksize: int = 20000):
    """
    Yield (gene, ligands) from the AlphaFill ligand table, reading it in chunks.

    Only LIGAND_INPUT_COLUMNS are loaded. The table is written sorted by gene, so a
    gene is complete once the next one starts; the last gene of a chunk is carried
    over to the next chunk. Genes are yielded in sorted order.
    """
    header = pd.read_csv(ligands_path, sep="\t", nrows=0).columns
    missing = [column for column in LIGAND_INPUT_COLUMNS if column not in header]
    if missing:
        raise ValueError(f"{ligands_path} lacks columns: {', '.join(missing)}")
    groups: dict[str, list[pd.DataFrame]] = {}
    pending = None
    for chunk in pd.read_csv(ligands_path, sep="\t", dtype=str, keep_default_na=False,
                             usecols=LIGAND_INPUT_COLUMNS, chunksize=chunksize):
        if pending is not None:
            chunk = pd.concat([pending, chunk], ignore_index=True)
        last_gene = chunk["gene"].iloc[-1] if len(chunk) else None
        pending = chunk[chunk["gene"] == last_gene]
        for gene, group in chunk[chunk["gene"] != last_gene].groupby("gene", sort=False):
            groups.setdefault(gene, []).append(group)
        # Genes already complete (an unsorted file is still handled correctly, at
        # the cost of keeping its groups until the end).
        if groups and _is_sorted_prefix(groups, last_gene):
            for gene in sorted(groups):
                yield gene, pd.concat(groups[gene], ignore_index=True)
            groups = {}
    if pending is not None and len(pending):
        groups.setdefault(pending["gene"].iloc[0], []).append(pending)
    for gene in sorted(groups):
        yield gene, pd.concat(groups[gene], ignore_index=True)


def _is_sorted_prefix(groups: dict[str, list[pd.DataFrame]], next_gene: str | None) -> bool:
    """True if every buffered gene sorts before the gene still being read."""
    return next_gene is None or max(groups) < next_gene


def _input_fingerprint(base: Path, root: Path, ligands_path: Path, faa: Path,
                       params: dict[str, Any]) -> str:
    """Digest of everything the overlap tables depend on.

    The comparison takes minutes, while its inputs rarely change between runs.
    The ligand table, protein FASTA and this module's code are hashed (memoized
    by files.sha256_file); the pocket/structure files of the loci with ligands
    and the AlphaFill model CIFs are fingerprinted by size and mtime_ns, which
    change whenever fpocket, P2Rank or AlphaFill rewrite them.
    """
    from ftscripts import files
    table = pd.read_csv(ligands_path, sep="\t", usecols=["gene", "cif_path"], dtype=str,
                        keep_default_na=False)
    digest = hashlib.sha256(json.dumps({
        "params": params,
        "code": files.sha256_file(__file__),
        "ligands": files.sha256_file(ligands_path),
        "faa": files.sha256_file(faa),
        "artifacts": files.sha256_file(ARTIFACTS_FILE) if ARTIFACTS_FILE.is_file() else None,
    }, sort_keys=True).encode())

    def add(path: Path, label: str) -> None:
        try:
            st = os.stat(path)
            digest.update(f"{label}\t{st.st_size}\t{st.st_mtime_ns}\n".encode())
        except FileNotFoundError:
            digest.update(f"{label}\tmissing\n".encode())

    structures_dir = base / "structures"
    for gene in sorted(set(table["gene"]) - {""}):
        locus_dir = structures_dir / gene
        add(locus_dir, gene)
        for folder, dirs, names in os.walk(locus_dir):
            dirs.sort()
            for name in sorted(names):
                path = Path(folder) / name
                add(path, str(path.relative_to(structures_dir)))
    for cif in sorted(set(table["cif_path"]) - {""}):
        add(root / cif, cif)
    return digest.hexdigest()


def _cached_summary(out: Path, fingerprint: str) -> pd.DataFrame | None:
    try:
        manifest = json.loads((out / MANIFEST_FILE).read_text())
    except (OSError, ValueError):
        return None
    if not isinstance(manifest, dict) or manifest.get("fingerprint") != fingerprint:
        return None
    if not (out / LIGAND_FILE).is_file() or not (out / SUMMARY_FILE).is_file():
        return None
    try:
        # "NA" is the CCD code of sodium: only empty cells are missing values.
        summary = pd.read_csv(out / SUMMARY_FILE, sep="\t", keep_default_na=False, na_values=[""])
    except (OSError, ValueError, pd.errors.ParserError):
        return None
    if list(summary.columns) != list(SUMMARY_COLUMNS):
        return None
    return summary


def run_pocket_overlap(output_path: str | Path, organism_name: str,
                       contact_cutoff: float = 4.5, min_overlap: float = 50.0,
                       min_mapping_identity: float = 90.0,
                       exclude_ions: bool = True,
                       min_contact_residues: int = 4,
                       exclude_artifacts: bool = True,
                       exclude_prosthetic: bool = True) -> pd.DataFrame:
    """Write per-ligand and per-locus overlap tables; return the locus table."""
    base = Path(output_path).resolve() / organism_name
    root = base / "ligands" / "alphafill"
    ligands_path = root / "alphafill_ligands_all.tsv"
    if not ligands_path.is_file():
        raise FileNotFoundError(f"Run AlphaFill first: {ligands_path} not found")
    faa = base / "genome" / f"{organism_name}.faa"
    out = root / OUTPUT_DIR
    params = {"contact_cutoff": contact_cutoff, "min_overlap": min_overlap,
              "min_mapping_identity": min_mapping_identity, "exclude_ions": exclude_ions,
              "min_contact_residues": min_contact_residues,
              "exclude_artifacts": exclude_artifacts, "exclude_prosthetic": exclude_prosthetic}
    excluded_classes = {name for name, flag in (("ion", exclude_ions), ("artifact", exclude_artifacts),
                                                ("prosthetic", exclude_prosthetic)) if flag}
    fingerprint = _input_fingerprint(base, root, ligands_path, faa, params)
    cached = _cached_summary(out, fingerprint)
    if cached is not None:
        logging.info("Pocket overlap reused from %s: inputs and parameters unchanged", out)
        return cached
    locus_seqs = {record.id: str(record.seq).rstrip("*").upper()
                  for record in SeqIO.parse(str(faa), "fasta")}
    mapper = SequenceMapper()
    structures_dir = base / "structures"
    rows = []
    for gene, gene_ligands in _ligands_by_gene(ligands_path):
        locus_seq = locus_seqs.get(gene, "")
        pockets = (load_locus_pockets(structures_dir / gene, locus_seq, mapper, min_mapping_identity)
                   if locus_seq else {m: {"status": "no_locus_sequence", "pockets": []} for m in METHODS})
        for cif_rel, group in gene_ligands.groupby("cif_path", sort=True):
            model, model_map, identity, model_error = None, {}, None, ""
            cif = root / cif_rel if cif_rel else None
            if cif is None or not cif.is_file():
                model_error = "no_model_cif"
            elif not locus_seq:
                model_error = "no_locus_sequence"
            else:
                try:
                    model = read_alphafill_model(cif)
                    index_map, identity = mapper.map(model["sequence"], locus_seq)
                    model_map = {model["seq_ids"][index]: position
                                 for index, position in index_map.items()}
                    if identity < min_mapping_identity:
                        model_error = "model_mapping_failed"
                except Exception as error:
                    logging.warning("Pocket overlap: cannot read %s: %s", cif, error)
                    model_error = "model_read_error"
            for _, ligand in group.iterrows():
                row = {key: ligand.get(key, "") for key in (
                    "gene", "uniprot", "selection_source", "ligand_instance_id",
                    "compound_code", "asym_id", "confidence_category")}
                row["model_mapping_identity"] = round(identity, 1) if identity is not None else None
                for method in METHODS:
                    row.update({column: None for column in _method_columns(method)})
                    row[f"{method}_status"] = pockets[method]["status"]
                placed = model["ligands"].get(ligand["asym_id"]) if model else None
                if placed is not None:
                    row["ligand_class"] = ligand_class(ligand["compound_code"] or placed["compound"],
                                                       len(placed["xyz"]))
                if model_error:
                    row["status"] = model_error
                elif placed is None:
                    row["status"] = "ligand_not_in_model"
                elif row["ligand_class"] in excluded_classes:
                    row["status"] = f"excluded_{row['ligand_class']}"
                else:
                    contacts = {model_map[number] for number in
                                contact_residues(model, placed["xyz"], contact_cutoff)
                                if number in model_map}
                    row["contact_residue_count"] = len(contacts)
                    row["contact_residues"] = ";".join(map(str, sorted(contacts)))
                    if not contacts:
                        row["status"] = "no_contacts"
                    else:
                        # Overlap is still reported for few-contact ligands, but
                        # 1-3 residues give meaningless percentages, so they are
                        # kept out of the locus summary.
                        row["status"] = ("compared" if len(contacts) >= min_contact_residues
                                         else "few_contacts")
                        for method in METHODS:
                            row.update(compare(contacts, pockets[method], method))
                rows.append(row)
    ligand_df = pd.DataFrame(rows, columns=LIGAND_COLUMNS)
    summary = summarize(ligand_df, min_overlap) if len(ligand_df) else pd.DataFrame(columns=SUMMARY_COLUMNS)
    out.mkdir(parents=True, exist_ok=True)
    (out / MANIFEST_FILE).unlink(missing_ok=True)
    ligand_df.to_csv(out / LIGAND_FILE, sep="\t", index=False)
    summary.to_csv(out / SUMMARY_FILE, sep="\t", index=False)
    (out / "README.txt").write_text(
        "AlphaFill ligand vs pocket overlap (FastTarget)\n\n"
        "Residues are compared in locus sequence numbering (genome .faa), so the\n"
        "pocket structure (PDB chain or AlphaFold model) and the AlphaFill model need\n"
        "not be the same file. Ligand contacts: model residues with a heavy atom within\n"
        f"{contact_cutoff} A of the ligand. Pocket residues: fpocket pocketN_atm.pdb or\n"
        "P2Rank residue_ids. Sequence mappings below "
        f"{min_mapping_identity}% identity are rejected.\n\n"
        "best_*: best pocket chosen by the structures module (best_DS_fpocket_<locus>.json,\n"
        "best_prob_p2rank_<locus>.json): highest fpocket Druggability Score or P2Rank\n"
        "probability, among pockets with mean pLDDT >= 70 for predicted models (the\n"
        "highest-scoring pocket when none passes). Without those files, the highest\n"
        "scoring pocket of the locus across all its structures.\n"
        "match_*: pocket that coincides most with the ligand (most shared residues),\n"
        "whatever its score: it is not the best pocket unless match_is_best is True.\n"
        "match_score_rank is its rank by score (1 = best pocket).\n"
        "ligand_overlap: % of ligand contact residues inside the pocket.\n"
        "pocket_overlap: % of pocket residues contacted by the ligand.\n\n"
        f"Summary: accepted (high/intermediate) ligands only; *_compounds lists those\n"
        f"with ligand_overlap >= {min_overlap}%; *_non_cofactor_compounds leaves out mobile\n"
        "cofactors (ATP, NAD, FAD, CoA...). Not compared (column ligand_class, status\n"
        "excluded_<class>): ions (known codes or one heavy atom), crystallization\n"
        "artifacts and noise (water, CO2, PEG, glycerol, detergents, buffers; BioLiP list\n"
        "in ftscripts/data/biolip_artifacts.tsv) and prosthetic groups that stay bound\n"
        "(heme, Fe-S clusters, molybdopterin, cobalamin). Ligands with fewer than "
        f"{min_contact_residues} contact\nresidues get status few_contacts: overlap is "
        "reported but not summarized.\n", encoding="utf-8")
    # Written last: interrupted output is never reused.
    temporary = out / f".{MANIFEST_FILE}.tmp"
    temporary.write_text(json.dumps({"fingerprint": fingerprint, "params": params}, indent=2) + "\n")
    os.replace(temporary, out / MANIFEST_FILE)
    logging.info("Pocket overlap: %s ligands compared, results in %s",
                 int((ligand_df["status"] == "compared").sum()), out)
    return summary


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output_path", required=True, help="Directory containing <organism>/")
    parser.add_argument("--organism", required=True)
    parser.add_argument("--contact_cutoff", type=float, default=4.5)
    parser.add_argument("--min_overlap", type=float, default=50.0)
    parser.add_argument("--min_mapping_identity", type=float, default=90.0)
    parser.add_argument("--min_contact_residues", type=int, default=4)
    parser.add_argument("--keep_ions", action="store_true")
    parser.add_argument("--keep_artifacts", action="store_true")
    parser.add_argument("--keep_prosthetic", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
    summary = run_pocket_overlap(args.output_path, args.organism, args.contact_cutoff,
                                 args.min_overlap, args.min_mapping_identity,
                                 not args.keep_ions, args.min_contact_residues,
                                 not args.keep_artifacts, not args.keep_prosthetic)
    print(f"{len(summary)} loci written to "
          f"{Path(args.output_path) / args.organism / 'ligands' / 'alphafill' / OUTPUT_DIR}")


if __name__ == "__main__":
    main()
