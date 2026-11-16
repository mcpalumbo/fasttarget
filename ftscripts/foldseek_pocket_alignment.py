#!/usr/bin/env python3
"""Screen PA14 pockets against the residue-level Foldseek human alignments.

One row is written for every pocket intersecting a significant Foldseek hit.
This tests structural correspondence, not the presence of a human binding pocket
or binding of a particular compound. Run after the Foldseek TSVs have finished.

Usage:
    python -m ftscripts.foldseek_pocket_alignment \
      organism/PaPA14 databases/human_structures \
      organism/PaPA14/offtarget/pocket_human_alignment.tsv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from statistics import mean

from Bio.Align import substitution_matrices


AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCEQGHILKMFPSTWYV",
))
AA3.update(MSE="M", SEC="U", PYL="O")
# Deliberately overlapping classes: histidine is aromatic and may be basic;
# methionine is aliphatic and sulfur-containing. These are descriptive labels,
# not evidence that two residues make the same ligand contacts.
CHEMICAL_GROUPS = {
    "basic": set("KRH"), "acidic": set("DE"), "aromatic": set("FWYH"),
    "aliphatic": set("AVILM"), "polar_uncharged": set("STNQ"),
    "hydroxyl": set("STY"), "amide": set("NQ"), "sulfur": set("CM"),
    "small": set("AGST"),
}
BLOSUM62 = substitution_matrices.load("BLOSUM62")
PREDICTED_PREFIXES = ("AF_", "AF3_", "CB_")
SCREEN_VERSION = 1
POCKET_ID = re.compile(r"([A-Za-z0-9])_(-?\d+)([A-Za-z]?)$")
FIELDS = [
    "gene", "query", "target", "method", "pocket", "pocket_score", "score_kind",
    "pocket_residues", "mapped_residues", "mapped_fraction", "mapped_pairs",
    "pocket_identity", "pocket_blosum62_positives", "pocket_blosum62_mean",
    "pocket_substitutions",
    "pocket_lddt_mean", "pocket_lddt_low_fraction",
    "query_pocket_plddt", "query_mapped_plddt", "target_mapped_plddt",
    "query_model_plddt", "target_model_plddt", "query_low_conf_fraction",
    "target_low_conf_fraction", "query_is_model", "target_is_model",
    "evalue", "prob", "alnlen", "qcov", "tcov", "lddt", "qtmscore",
    "ttmscore", "alntmscore", "assessment", "reason",
]


def number(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def rounded(values, digits=3):
    values = [v for v in values if v is not None]
    return round(mean(values), digits) if values else None


def key_label(key: tuple[str, int, str]) -> str:
    return f"{key[0]}:{key[1]}{key[2]}"


def shared_chemistry(first: str, second: str) -> list[str]:
    """Chemical properties shared by two different amino acids."""
    return [name for name, group in CHEMICAL_GROUPS.items()
            if first in group and second in group]


def blosum62_score(first: str, second: str) -> float | None:
    try:
        return float(BLOSUM62[first, second])
    except (KeyError, IndexError):
        return None


@dataclass(frozen=True)
class Residue:
    key: tuple[str, int, str]
    aa: str
    confidence: float | None


@lru_cache(maxsize=256)
def pdb_residues(path: str) -> tuple[Residue, ...]:
    """Residues in Foldseek sequence order, retaining original PDB numbering."""
    if Path(path).suffix.lower() in (".cif", ".mmcif"):
        from Bio.PDB import MMCIFParser

        structure = MMCIFParser(QUIET=True).get_structure("human", path)
        model = next(structure.get_models(), None)
        residues = []
        if model is not None:
            for chain in model:
                for residue in chain:
                    if "CA" not in residue or residue.resname.strip() not in AA3:
                        continue
                    key = (chain.id.strip(), residue.id[1], residue.id[2].strip())
                    residues.append(Residue(key, AA3[residue.resname.strip()],
                                            float(residue["CA"].bfactor)))
        return tuple(residues)
    residues = []
    seen = set()
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM  ", "HETATM")) or line[12:16].strip() != "CA":
                continue
            if line.startswith("HETATM") and line[17:20].strip() not in AA3:
                continue
            try:
                key = (line[21].strip(), int(line[22:26]), line[26].strip())
                if key in seen:
                    continue
                confidence = float(line[60:66])
            except (ValueError, IndexError):
                continue
            seen.add(key)
            residues.append(Residue(key, AA3.get(line[17:20].strip(), "X"), confidence))
    if residues and max(r.confidence for r in residues) <= 1:
        residues = [Residue(r.key, r.aa, r.confidence * 100) for r in residues]
    return tuple(residues)


def pocket_key(token: str) -> tuple[str, int, str] | None:
    match = POCKET_ID.fullmatch(token.strip())
    return (match[1], int(match[2]), match[3]) if match else None


def load_pockets(pocket_root: Path, query: str):
    """Yield method, pocket name, score, and PDB residue keys."""
    p2rank = pocket_root / f"{query}_p2rank" / f"{query}.pdb_predictions.csv"
    if p2rank.is_file():
        with p2rank.open(newline="") as handle:
            for row in csv.DictReader(handle, skipinitialspace=True):
                row = {k.strip(): v.strip() for k, v in row.items() if k}
                keys = {key for token in row.get("residue_ids", "").split()
                        if (key := pocket_key(token)) is not None}
                if keys:
                    yield "p2rank", row["name"], number(row.get("probability")), "probability", keys

    fp_dir = pocket_root / f"{query}_fpocket"
    properties_path = fp_dir / "all_pockets.json"
    if properties_path.is_file():
        properties = json.loads(properties_path.read_text())
        properties = properties.get(query.removeprefix("AF_").removeprefix("PDB_"),
                                    next(iter(properties.values()), {}))
        if not isinstance(properties, dict):
            properties = {}
        for name, values in properties.items():
            match = re.search(r"(\d+)$", name)
            if not match or not isinstance(values, dict):
                continue
            atom_file = fp_dir / "pockets" / f"pocket{match[1]}_atm.pdb"
            if not atom_file.is_file():
                continue
            keys = set()
            with atom_file.open(errors="replace") as handle:
                for line in handle:
                    if line.startswith("ATOM"):
                        try:
                            keys.add((line[21].strip(), int(line[22:26]), line[26].strip()))
                        except (ValueError, IndexError):
                            pass
            if keys:
                yield "fpocket", name, number(values.get("Druggability Score")), "druggability_score", keys


def query_index(structures_dir: Path):
    """Index only primary query PDBs, never the numerous pocket PDB files."""
    result = defaultdict(list)
    for gene_dir in structures_dir.glob("PA14_*"):
        if not gene_dir.is_dir():
            continue
        for subdir in gene_dir.iterdir():
            if not subdir.is_dir() or subdir.name == "pockets":
                continue
            for pdb in subdir.glob("*.pdb"):
                pocket_root = subdir / "pockets" if subdir.name == "AF_model" else gene_dir / "pockets"
                result[pdb.stem].append((gene_dir.name, pdb, pocket_root))
    return result


def alignment_pairs(row: dict, query: tuple[Residue, ...], target: tuple[Residue, ...]):
    """Return (query residue, target residue, lDDT) for each aligned pair.

    Foldseek starts are 1-based sequence indices; verify letters to catch an
    incompatible PDB or shifted numbering before using a residue map.
    """
    qaln, taln = row.get("qaln", ""), row.get("taln", "")
    if not qaln or len(qaln) != len(taln):
        return None, "missing_or_unequal_alignments"
    try:
        qpos, tpos = int(row["qstart"]) - 1, int(row["tstart"]) - 1
        local_lddt = [float(v) for v in row["lddtfull"].split(",") if v]
    except (KeyError, ValueError):
        return None, "invalid_alignment_coordinates_or_lddt"
    pairs = []
    bad_letters = 0
    compared = 0
    li = 0
    for qa, ta in zip(qaln, taln):
        qr = query[qpos] if qa != "-" and 0 <= qpos < len(query) else None
        tr = target[tpos] if ta != "-" and 0 <= tpos < len(target) else None
        if qa != "-" and qr is None or ta != "-" and tr is None:
            return None, "alignment_outside_structure"
        for letter, residue in ((qa, qr), (ta, tr)):
            if residue is not None and residue.aa != "X" and letter != "X":
                compared += 1
                bad_letters += residue.aa != letter
        if qr is not None and tr is not None:
            if li >= len(local_lddt):
                return None, "lddtfull_shorter_than_aligned_pairs"
            pairs.append((qr, tr, local_lddt[li]))
            li += 1
        qpos += qa != "-"
        tpos += ta != "-"
    if li != len(local_lddt):
        return None, "lddtfull_longer_than_aligned_pairs"
    if compared and bad_letters / compared > 0.05:
        return None, "structure_sequence_mismatch"
    return pairs, None


def assess(row, method, pocket_name, score, score_kind, keys, query, target, pairs, error,
           min_pocket_coverage, min_local_lddt, min_qtm, min_plddt):
    qmodel = row["query"].startswith(PREDICTED_PREFIXES)
    tmodel = row["target"].startswith(PREDICTED_PREFIXES)
    q_by_key = {r.key: r for r in query}
    pocket_pairs = [(q, t, l) for q, t, l in (pairs or []) if q.key in keys]
    mapped = {q.key for q, _, _ in pocket_pairs}
    q_all_conf = [q_by_key[key].confidence for key in keys if key in q_by_key] if qmodel else []
    q_conf = [q.confidence for q, _, _ in pocket_pairs] if qmodel else []
    t_conf = [t.confidence for _, t, _ in pocket_pairs] if tmodel else []
    local_lddt = [l for _, _, l in pocket_pairs]
    fraction = len(mapped) / len(keys)
    identity = sum(q.aa == t.aa for q, t, _ in pocket_pairs) / len(pocket_pairs) if pocket_pairs else None
    substitution_scores = [blosum62_score(q.aa, t.aa) for q, t, _ in pocket_pairs]
    known_scores = [v for v in substitution_scores if v is not None]
    positives = sum(v > 0 for v in known_scores) / len(known_scores) if known_scores else None
    substitutions = ";".join(
        f"{key_label(q.key)}{q.aa}>{key_label(t.key)}{t.aa}"
        f"[{','.join(shared_chemistry(q.aa, t.aa)) or 'different'}]"
        for q, t, _ in pocket_pairs if q.aa != t.aa
    )
    qtm = number(row.get("qtmscore"))
    reasons = []
    if error:
        assessment = "invalid_mapping"
        reasons.append(error)
    else:
        if fraction < min_pocket_coverage:
            reasons.append("low_pocket_coverage")
        if rounded(local_lddt) is None or rounded(local_lddt) < min_local_lddt:
            reasons.append("low_local_lddt")
        if qmodel and (rounded(q_conf) is None or rounded(q_conf) < min_plddt):
            reasons.append("low_or_missing_query_plddt")
        if tmodel and (rounded(t_conf) is None or rounded(t_conf) < min_plddt):
            reasons.append("low_or_missing_target_plddt")
        if reasons:
            assessment = "weak_or_uncertain"
        elif qtm is None or qtm < min_qtm:
            assessment = "local_match_review"
            reasons.append("low_query_tm_global")
        else:
            assessment = "structural_candidate"
    return {
        "method": method, "pocket": pocket_name, "pocket_score": score, "score_kind": score_kind,
        "pocket_residues": len(keys), "mapped_residues": len(mapped),
        "mapped_fraction": round(fraction, 3),
        "mapped_pairs": ";".join(f"{key_label(q.key)}>{key_label(t.key)}" for q, t, _ in pocket_pairs),
        "pocket_identity": round(identity, 3) if identity is not None else None,
        "pocket_blosum62_positives": round(positives, 3) if positives is not None else None,
        "pocket_blosum62_mean": rounded(known_scores),
        "pocket_substitutions": substitutions,
        "pocket_lddt_mean": rounded(local_lddt),
        "pocket_lddt_low_fraction": round(sum(x < min_local_lddt for x in local_lddt) / len(local_lddt), 3)
        if local_lddt else None,
        "query_pocket_plddt": rounded(q_all_conf) if qmodel else None,
        "query_mapped_plddt": rounded(q_conf) if qmodel else None,
        "target_mapped_plddt": rounded(t_conf) if tmodel else None,
        "query_model_plddt": rounded(r.confidence for r in query) if qmodel else None,
        "target_model_plddt": rounded(r.confidence for r in target) if tmodel else None,
        "query_low_conf_fraction": round(sum(v < min_plddt for v in q_conf) / len(q_conf), 3)
        if q_conf else None,
        "target_low_conf_fraction": round(sum(v < min_plddt for v in t_conf) / len(t_conf), 3)
        if t_conf else None,
        "query_is_model": qmodel, "target_is_model": tmodel,
        "assessment": assessment, "reason": ";".join(reasons),
    }


def run(args, writer):
    org = Path(args.organism)
    queries = query_index(org / "structures")
    human = Path(args.human_structures)
    files = sorted((org / "offtarget" / "foldseek_results").glob("*_output_foldseek/*results.tsv"))
    if args.query:
        files = [p for p in files if p.parent.name == f"{args.query}_output_foldseek"]
        if not files:
            raise FileNotFoundError(f"Foldseek TSV not found for query {args.query}")
    elif not files:
        raise FileNotFoundError(f"No Foldseek TSVs found under {org / 'offtarget' / 'foldseek_results'}")
    stats = defaultdict(int)
    for result_file in files:
        query_name = result_file.parent.name.removesuffix("_output_foldseek")
        if query_name not in queries:
            blank = dict.fromkeys(FIELDS)
            blank.update(query=query_name, assessment="query_structure_missing",
                         reason="query_pdb_not_found_under_structures")
            writer.writerow(blank)
            stats["query_structure_missing"] += 1
            continue
        for gene, query_path, pocket_root in queries.get(query_name, []):
            pockets = [p for p in load_pockets(pocket_root, query_name)
                       if p[2] is not None and p[2] >=
                       (args.min_p2rank_probability if p[0] == "p2rank" else args.min_fpocket_ds)]
            if not pockets:
                blank = dict.fromkeys(FIELDS)
                blank.update(gene=gene, query=query_name, assessment="no_selected_pocket",
                             reason="no_pocket_above_score_threshold_or_pocket_files_missing")
                writer.writerow(blank)
                stats["queries_without_pockets"] += 1
                continue
            query = pdb_residues(str(query_path))
            key_to_index = {r.key: i + 1 for i, r in enumerate(query)}
            indexed_pockets = [(m, n, s, kind, keys, {key_to_index[k] for k in keys if k in key_to_index})
                               for m, n, s, kind, keys in pockets]
            seen = set()
            with result_file.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for hit in reader:
                    evalue = number(hit.get("evalue"))
                    if evalue is None or evalue > args.max_evalue:
                        continue
                    try:
                        qstart, qend = int(hit["qstart"]), int(hit["qend"])
                    except (KeyError, ValueError):
                        continue
                    intersecting = [p for p in indexed_pockets if any(qstart <= i <= qend for i in p[5])]
                    if not intersecting:
                        continue
                    target_path = human / f'{hit["target"]}.pdb'
                    if not target_path.is_file():
                        target_path = human / f'{hit["target"]}.cif'
                    if not target_path.is_file():
                        stats["missing_human_structures"] += 1
                        continue
                    target = pdb_residues(str(target_path))
                    pairs, error = alignment_pairs(hit, query, target)
                    for method, name, score, kind, keys, _ in intersecting:
                        seen.add((method, name))
                        result = assess(hit, method, name, score, kind, keys, query, target,
                                        pairs, error, args.min_pocket_coverage, args.min_local_lddt,
                                        args.min_qtm, args.min_plddt)
                        result.update(gene=gene, query=query_name, target=hit["target"])
                        result.update({k: hit.get(k) for k in
                                       ("evalue", "prob", "alnlen", "qcov", "tcov", "lddt",
                                        "qtmscore", "ttmscore", "alntmscore")})
                        writer.writerow(result)
                        stats[result["assessment"]] += 1
            q_by_key = {r.key: r for r in query}
            for method, name, score, kind, keys, _ in indexed_pockets:
                if (method, name) in seen:
                    continue
                qmodel = query_name.startswith(PREDICTED_PREFIXES)
                blank = dict.fromkeys(FIELDS)
                blank.update(gene=gene, query=query_name, method=method, pocket=name,
                             pocket_score=score, score_kind=kind, pocket_residues=len(keys),
                             mapped_residues=0, mapped_fraction=0.0,
                             query_is_model=qmodel,
                             query_pocket_plddt=rounded(q_by_key[k].confidence for k in keys if k in q_by_key)
                             if qmodel else None,
                             query_model_plddt=rounded(r.confidence for r in query) if qmodel else None,
                             assessment="no_evaluable_hit",
                             reason="no_significant_hit_intersects_pocket_or_human_structure_missing")
                writer.writerow(blank)
                stats["no_evaluable_hit"] += 1
            stats["queries_processed"] += 1
    return dict(stats)


def screen_to_tsv(organism: Path, human_structures: Path, output: Path):
    """Run from FastTarget, refreshing the TSV when Foldseek/pocket inputs change."""
    from . import files as cache_files

    organism = Path(organism)
    human_structures = Path(human_structures)
    output = Path(output)
    result_files = sorted((organism / "offtarget" / "foldseek_results").glob(
        "*_output_foldseek/*results.tsv"))
    if not result_files:
        raise FileNotFoundError("Foldseek pocket screen needs completed Foldseek TSVs")
    indexed = query_index(organism / "structures")
    inputs = list(result_files)
    for result_file in result_files:
        query_name = result_file.parent.name.removesuffix("_output_foldseek")
        for _, query_path, pocket_root in indexed.get(query_name, []):
            inputs.extend((query_path,
                           pocket_root / f"{query_name}_p2rank" / f"{query_name}.pdb_predictions.csv",
                           pocket_root / f"{query_name}_fpocket" / "all_pockets.json"))
    inputs.append(human_structures / "DB_foldseek" / "DB_human_reference")
    digest = hashlib.sha256()
    for path in sorted(set(inputs)):
        try:
            stat = path.stat()
        except FileNotFoundError:
            continue
        digest.update(f"{path}:{stat.st_size}:{stat.st_mtime_ns}\n".encode())
    parameters = {"screen_version": SCREEN_VERSION, "input_signature": digest.hexdigest()[:16],
                  "foldseek_files": len(result_files)}
    if cache_files.manifest_matches(str(output), parameters, label=output.name):
        print(f"Pocket-human alignment already current: {output}")
        return None
    args = argparse.Namespace(organism=organism, human_structures=human_structures,
                              query=None, max_evalue=0.01, min_pocket_coverage=0.7,
                              min_local_lddt=0.5, min_qtm=0.5, min_plddt=70.0,
                              min_p2rank_probability=0.5, min_fpocket_ds=0.5)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        stats = run(args, writer)
    os.replace(temporary, output)
    cache_files.write_manifest(str(output), parameters, extra={"summary": stats}, hash_output=False)
    print(f"Pocket-human alignment saved: {output} ({stats})")
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("organism", type=Path)
    parser.add_argument("human_structures", type=Path)
    parser.add_argument("output", type=Path, help="TSV destination, or '-' for stdout")
    parser.add_argument("--query", help="Only process one Foldseek query identifier")
    parser.add_argument("--max-evalue", type=float, default=0.01)
    parser.add_argument("--min-pocket-coverage", type=float, default=0.7)
    parser.add_argument("--min-local-lddt", type=float, default=0.5)
    parser.add_argument("--min-qtm", type=float, default=0.5)
    parser.add_argument("--min-plddt", type=float, default=70.0)
    parser.add_argument("--min-p2rank-probability", type=float, default=0.5)
    parser.add_argument("--min-fpocket-ds", type=float, default=0.5)
    args = parser.parse_args()
    if str(args.output) == "-":
        import sys
        writer = csv.DictWriter(sys.stdout, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        stats = run(args, writer)
        print(json.dumps(stats), file=sys.stderr)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            writer.writeheader()
            stats = run(args, writer)
        print(json.dumps(stats))


if __name__ == "__main__":
    main()
