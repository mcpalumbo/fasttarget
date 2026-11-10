"""Local BindingDB extraction and protein mapping, independent of FastTarget execution.

The original MySQL dump is retained. ``prepare`` creates an indexed SQLite copy
of the scientific tables, without executing SQL from the dump. No network APIs
or DuckDB are required. Activity categories are deliberately not inferred.

Run with ``conda activate fasttarget``. Prepare once per downloaded release::

    python -m ftscripts.bindingdb prepare --dump databases/bindingdb/202611/BDB-mySQL_All_202611.dmp

Extract a species and map its proteins (FASTA IDs must be locus tags)::

    python -m ftscripts.bindingdb run --database databases/bindingdb/202611/BDB-mySQL_All_202611.sqlite \
        --taxid 287 --taxonomy-dir /path/to/ncbi_taxonomy --proteins /path/to/proteins.faa \
        --output-dir /path/to/organism/ligands/bindingdb --include-human

Taxonomy files are separate data inputs: nodes.dmp is required, merged.dmp and
names.dmp are recommended. Python expands descendants if TaxonKit is absent
from this environment. An explicit --taxids-file is supported and recorded as
not verified complete. BLAST tools are resolved only in the Python environment.

Outputs are TSV and CSV: source proteins, raw/normalized activities, unique
BindingDB monomer IDs, sequence mapping, locus-to-ligand activities, locus
summary, references and optional human evidence. A measurement can have several
endpoints; measurement counts deduplicate its ID. Counts indicate records, not
proven activity. No compound-parent collapse or activity thresholds are applied.
ITC values retain their source fields without an inferred unit conversion.
Human off-target evidence includes individual proteins and complexes containing
human components; mixed/unknown component taxids are reported. Complex activity
is not attributed to any individual subunit. MIC/organism-level assays are deferred. run_metadata.json records completion,
source fingerprint, taxonomy scope and mapping thresholds. Default BLAST filters
are 90% identity and 70% coverage on BOTH genome and BindingDB sequences.
Each BindingDB protein is queried against the genome. Among passing hits, keep
the highest bit score, then lowest E-value. Exact best-hit ties map to all tied
loci and are reported in mapping/bindingdb_best_hit_ties.tsv. Missing evidence is not biological inactivity. This module is standalone;
it does not modify the FastTarget configuration or final pipeline table.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import html
import json
import logging
import math
from pathlib import Path
import re
import sqlite3
import sys
import subprocess
from typing import Iterator

import pandas as pd
from Bio import SeqIO

LOG = logging.getLogger(__name__)
PREPARATION_VERSION = 1
TABLES = {
    "polymer", "monomer", "enzyme_reactant_set", "ki_result", "assay", "entry",
    "entry_citation", "article", "itc_result_a_b_ab", "complex", "complex_component",
}
REQUIRED = {"polymer", "monomer", "enzyme_reactant_set", "ki_result", "assay"}
# Bounded to literal values emitted by mysqldump; unsupported syntax is an error.
_LITERAL = re.compile(r"'(?:[^'\\]|\\.|'')*'|NULL|[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", re.DOTALL)
_ESCAPES = {"0": "\0", "b": "\b", "n": "\n", "r": "\r", "t": "\t", "Z": "\x1a",
            "\\": "\\", "'": "'", '"': '"', "%": "\\%", "_": "\\_"}
ENDPOINTS = ("ki", "kd", "ic50", "ec50")
PROTEIN_COLUMNS = ["polymerid", "unpid1", "unpid2", "display_name", "scientific_name",
                   "source_organism", "taxid", "type", "sequence", "chembl_id"]
MAPPING_COLUMNS = ["gene", "polymerid", "bindingdb_taxid", "identity", "genome_coverage", "bindingdb_coverage",
                   "evalue", "bitscore", "mapping_method", "best_hit_tie_count"]
ACTIVITY_COLUMNS = ["measurement_id", "measurement_source", "target_id", "target_type", "polymerid", "complexid", "monomerid", "compound_id",
                    "entryid", "assayid", "endpoint", "value_raw", "relation_raw", "relation",
                    "value", "units", "pactivity", "value_status", "assay_name", "assay_description",
                    "raw_record_json", "reactant_set_json"]
LIGAND_COLUMNS = ["monomerid", "compound_id", "display_name", "chembl_id", "inchi_key", "inchi",
                  "smiles_string", "bacterial_target_count", "bacterial_measurement_count",
                  "bacterial_quantitative_measurement_count", "human_target_count", "human_protein_target_count",
                  "human_complex_target_count", "human_measurement_count"]
SUMMARY_COLUMNS = ["gene", "bindingdb_status", "bindingdb_target_count", "bindingdb_compounds",
                   "bindingdb_compound_count"]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".part")
    partial.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    partial.replace(path)


def _environment_tool(name: str, required: bool = True) -> str | None:
    """Resolve tools only within the Python environment executing this module."""
    candidate = Path(sys.prefix) / "bin" / name
    if candidate.is_file():
        return str(candidate)
    if required:
        raise FileNotFoundError(f"{name} is missing from environment {sys.prefix}; add it to requirements.yml")
    return None


def _quote(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def iter_mysql_rows(values: str) -> Iterator[list]:
    """Decode mysqldump literals, preserving Unicode, NULL and backslash escapes."""
    pos, size = 0, len(values)
    while pos < size:
        while pos < size and values[pos].isspace():
            pos += 1
        if pos >= size or (values[pos] == ";" and not values[pos + 1:].strip()):
            return
        if pos >= size or values[pos] != "(":
            raise ValueError(f"Expected a row at SQL offset {pos}")
        pos += 1
        row = []
        while True:
            while pos < size and values[pos].isspace():
                pos += 1
            match = _LITERAL.match(values, pos)
            if not match:
                raise ValueError(f"Unsupported MySQL literal at SQL offset {pos}")
            token = match.group()
            if token.startswith("'"):
                # Quoted values stay strings, including numeric-looking identifiers.
                value = re.sub(r"\\(.)|'(')", lambda m: _ESCAPES.get(m[1], m[1]) if m[1] else "'", token[1:-1])
            elif token == "NULL":
                value = None
            else:
                # SQLite column affinity converts numeric types; DECIMAL remains text.
                value = token
            row.append(value)
            pos = match.end()
            while pos < size and values[pos].isspace():
                pos += 1
            if pos >= size:
                raise ValueError("Unterminated INSERT row")
            delimiter = values[pos]
            pos += 1
            if delimiter == ")":
                break
            if delimiter != ",":
                raise ValueError(f"Expected a value separator at SQL offset {pos - 1}")
        yield row
        while pos < size and values[pos].isspace():
            pos += 1
        if pos < size and values[pos] == ",":
            pos += 1
            if not values[pos:].strip() or values[pos:].strip() == ";":
                raise ValueError("Trailing comma after INSERT row")
        elif values[pos:].strip() in ("", ";"):
            return
        else:
            raise ValueError(f"Unexpected SQL after INSERT row at offset {pos}")


def prepare_bindingdb(dump_path: str | Path, database_path: str | Path | None = None) -> Path:
    """Build a reusable SQLite database atomically from the official MySQL dump."""
    source = Path(dump_path).expanduser().resolve()
    destination = Path(database_path).expanduser().resolve() if database_path else source.with_suffix(".sqlite")
    if source == destination:
        raise ValueError("SQLite output must differ from the original dump")
    if not source.is_file():
        raise FileNotFoundError(source)
    manifest_path = destination.with_suffix(destination.suffix + ".metadata.json")
    fingerprint = _sha256(source)
    if destination.exists() and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("source_sha256") == fingerprint and manifest.get("preparation_version") == PREPARATION_VERSION:
            LOG.info("Reusing prepared BindingDB database: %s", destination)
            return destination
        raise ValueError("Existing SQLite metadata does not match this dump; use another output path")
    if destination.exists():
        raise ValueError(f"Existing SQLite has no preparation metadata: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    if partial.exists():
        partial.unlink()
    connection = sqlite3.connect(partial)
    counts, schema, current, columns, primary = {}, {}, None, [], []
    completed_dump = False
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        with source.open(encoding="utf-8") as incoming:
            for line_number, line in enumerate(incoming, 1):
                create = re.match(r"CREATE TABLE `([^`]+)` \(", line)
                if create:
                    current = create[1] if create[1] in TABLES else None
                    columns, primary = [], []
                    if current:
                        schema[current] = []
                    continue
                if current:
                    schema[current].append(line.rstrip())
                    column = re.match(r"\s+`([^`]+)`\s+([A-Za-z]+)", line)
                    if column:
                        name, mysql_type = column.groups()
                        affinity = "INTEGER" if mysql_type.lower() in {"int", "bigint", "smallint", "tinyint", "mediumint"} else "TEXT"
                        columns.append((name, affinity))
                    elif line.lstrip().startswith("PRIMARY KEY"):
                        primary = re.findall(r"`([^`]+)`", line)
                    if line.startswith(")"):
                        definitions = [f"{_quote(name)} {kind}" for name, kind in columns]
                        if primary:
                            definitions.append("PRIMARY KEY (" + ",".join(map(_quote, primary)) + ")")
                        connection.execute(f"CREATE TABLE {_quote(current)} ({','.join(definitions)})")
                        counts[current] = 0
                        LOG.info("Importing table %s", current)
                        current = None
                    continue
                insert = re.match(r"INSERT INTO `([^`]+)` VALUES ", line)
                if insert and insert[1] in TABLES:
                    table = insert[1]
                    width = len(connection.execute(f"PRAGMA table_info({_quote(table)})").fetchall())
                    statement = f"INSERT INTO {_quote(table)} VALUES ({','.join('?' for _ in range(width))})"
                    batch = []
                    for row in iter_mysql_rows(line[insert.end():]):
                        if len(row) != width:
                            raise ValueError(f"Column mismatch in {table} at dump line {line_number}")
                        batch.append(row)
                        if len(batch) >= 2000:
                            connection.executemany(statement, batch)
                            counts[table] += len(batch)
                            batch.clear()
                    if batch:
                        connection.executemany(statement, batch)
                        counts[table] += len(batch)
                elif re.match(r"INSERT INTO `([^`]+)`", line) and re.match(r"INSERT INTO `([^`]+)`", line)[1] in TABLES:
                    raise ValueError(f"Unsupported INSERT syntax at dump line {line_number}")
                elif line.startswith("-- Dump completed on"):
                    completed_dump = True
        if not completed_dump:
            raise ValueError("Missing mysqldump completion marker; dump may be truncated")
        if not REQUIRED.issubset(counts):
            raise ValueError(f"Missing BindingDB tables: {REQUIRED - set(counts)}")
        for table, column_names in {
            "polymer": [("taxid",)], "enzyme_reactant_set": [("enzyme_polymerid",), ("inhibitor_monomerid",)],
            "ki_result": [("reactant_set_id",)], "entry_citation": [("entryid",)],
            "itc_result_a_b_ab": [("cell_polymerid",), ("syr_polymerid",), ("cell_monomerid",), ("syr_monomerid",)],
        }.items():
            if table in counts:
                for names in column_names:
                    name = "idx_" + table + "_" + "_".join(names)
                    connection.execute(f"CREATE INDEX {_quote(name)} ON {_quote(table)} ({','.join(map(_quote, names))})")
        connection.commit()
        if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("SQLite integrity check failed")
        actual = {table: connection.execute(f"SELECT COUNT(*) FROM {_quote(table)}").fetchone()[0] for table in counts}
        if actual != counts:
            raise ValueError("Imported row counts do not match the parsed dump")
    except Exception:
        connection.close()
        partial.unlink(missing_ok=True)
        raise
    connection.close()
    partial.replace(destination)
    _json(manifest_path, {"preparation_version": PREPARATION_VERSION, "source": source.name,
                         "source_sha256": fingerprint, "database": destination.name,
                         "prepared_at": datetime.now(timezone.utc).isoformat(), "row_counts": counts,
                         "numeric_policy": "INTEGER for integer columns; all other source values preserved as TEXT",
                         "tables": sorted(counts), "integrity_check": "ok"})
    LOG.info("Prepared BindingDB: %s; rows per table: %s", destination, counts)
    return destination


def expand_taxids(taxid: int, taxonomy_dir: str | Path | None = None,
                  taxids: list[int] | None = None) -> tuple[list[int], dict]:
    """Expand NCBI descendants; an explicit set is recorded as such, never implied complete."""
    if isinstance(taxid, bool) or int(taxid) <= 0:
        raise ValueError("taxid must be a positive integer")
    taxid = int(taxid)
    if taxids is not None:
        if any(isinstance(t, bool) or int(t) <= 0 for t in taxids):
            raise ValueError("Expanded taxids must be positive integers")
        return sorted({taxid, *map(int, taxids)}), {"method": "explicit_taxids", "descendants_verified": False}
    directory = Path(taxonomy_dir or Path.home() / ".taxonkit").expanduser().resolve()
    nodes = directory / "nodes.dmp"
    if not nodes.is_file():
        raise FileNotFoundError("Provide --taxonomy-dir with nodes.dmp (and preferably names.dmp, merged.dmp), "
                                "or an explicit --taxids-file. Species expansion is never silently skipped.")
    merged = directory / "merged.dmp"
    redirects = {}
    if merged.is_file():
        for line in merged.open():
            parts = line.split("|")
            redirects[int(parts[0])] = int(parts[1])
    requested, seen = taxid, set()
    while taxid in redirects:
        if taxid in seen:
            raise ValueError("Cycle in merged taxids")
        seen.add(taxid)
        taxid = redirects[taxid]
    taxonkit = _environment_tool("taxonkit", required=False)
    if taxonkit:
        result = subprocess.run([taxonkit, "list", "--ids", str(taxid), "--indent", "", "--data-dir", str(directory)],
                                check=True, capture_output=True, text=True)
        expanded = {int(line.strip()) for line in result.stdout.splitlines() if line.strip().isdigit()}
        method = "taxonkit"
    else:
        children, present = defaultdict(list), False
        with nodes.open() as handle:
            for line in handle:
                parts = line.split("|")
                child, parent = int(parts[0]), int(parts[1])
                if child == taxid:
                    present = True
                if child != parent:
                    children[parent].append(child)
        if not present:
            raise ValueError(f"Unknown or deleted NCBI taxid: {taxid}")
        expanded, pending = set(), [taxid]
        while pending:
            node = pending.pop()
            if node not in expanded:
                expanded.add(node)
                pending.extend(children.get(node, []))
        method = "ncbi_nodes"
    if taxid not in expanded:
        raise ValueError(f"Taxonomy expansion did not include the root taxid {taxid}")
    metadata = {"method": method, "descendants_verified": True, "requested_taxid": requested,
                "resolved_taxid": taxid, "files": {}}
    for name in ("nodes.dmp", "merged.dmp", "names.dmp", "delnodes.dmp"):
        path = directory / name
        if path.is_file():
            metadata["files"][name] = {"sha256": _sha256(path), "bytes": path.stat().st_size}
    names = directory / "names.dmp"
    if names.is_file():
        with names.open() as handle:
            for line in handle:
                parts = [part.strip() for part in line.split("|")]
                if parts[0] == str(taxid) and parts[3] == "scientific name":
                    metadata["scientific_name"] = parts[1]
                    break
    return sorted(expanded), metadata


def _frame(connection, sql: str, params=()) -> pd.DataFrame:
    return pd.read_sql_query(sql, connection, params=params)


def _temporary_ids(connection, table: str, column: str, values) -> None:
    connection.execute(f"DROP TABLE IF EXISTS temp.{table}")
    connection.execute(f"CREATE TEMP TABLE {table} ({column} INTEGER PRIMARY KEY)")
    connection.executemany(f"INSERT INTO {table} VALUES (?)", [(int(value),) for value in sorted(set(values))])


def parse_measurement(raw, uncertainty=None) -> dict:
    """Keep censored values as bounds; only exact positive nM measurements get pactivity."""
    result = {"value_raw": raw, "relation_raw": uncertainty, "relation": None, "value": None,
              "units": "nM", "pactivity": None, "value_status": "unparsed"}
    if raw is None or not str(raw).strip():
        result["value_status"] = "missing"
        return result
    text = html.unescape(str(raw)).strip().replace("≤", "<=").replace("≥", ">=")
    match = re.fullmatch(r"\s*(<=|>=|<|>|=|~)?\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*", text)
    if not match:
        return result
    relation, value = match.groups()
    qualifier = html.unescape(str(uncertainty or "")).strip().replace("≤", "<=").replace("≥", ">=")
    if qualifier and qualifier not in {"=", "<", ">", "<=", ">=", "~"}:
        result.update(value=float(value), value_status="unsupported_uncertainty")
        return result
    if relation and qualifier and relation != qualifier:
        result["value_status"] = "conflicting_relations"
        return result
    relation = relation or qualifier or "="
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        result["value_status"] = "nonpositive_or_nonfinite"
        return result
    result.update(relation=relation, value=value, value_status="exact" if relation == "=" else "bound_or_approximate")
    if relation == "=":
        result["pactivity"] = 9 - math.log10(value)
    return result


def _normalized_sequence(value) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper().rstrip("*")


def _protein_evidence(row: dict) -> str | None:
    """Use sequence only. The free-text type and accessions remain descriptive."""
    sequence = _normalized_sequence(row.get("sequence"))
    if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+", sequence) and set(sequence) - set("ACGTUN"):
        return "amino_acid_sequence"
    return None


def _protein_rows(connection, taxids: list[int]) -> pd.DataFrame:
    _temporary_ids(connection, "selected_taxids", "taxid", taxids)
    frame = _frame(connection, """SELECT p.* FROM polymer p JOIN selected_taxids t
                    ON trim(p.taxid) = CAST(t.taxid AS TEXT) ORDER BY p.polymerid""")
    frame["protein_selection_evidence"] = [_protein_evidence(row) for row in frame.to_dict("records")]
    selected = frame[frame.protein_selection_evidence.notna()].copy()
    LOG.info("Protein selection: %s/%s polymers retained for %s taxids", len(selected), len(frame), len(taxids))
    return selected


def _activities(connection, proteins: pd.DataFrame, monomer_ids=None, *, target_kind="protein") -> pd.DataFrame:
    if target_kind not in {"protein", "complex"}:
        raise ValueError("target_kind must be protein or complex")
    column = "polymerid" if target_kind == "protein" else "complexid"
    table = "selected_proteins" if target_kind == "protein" else "selected_complexes"
    prefix = "BDBP" if target_kind == "protein" else "BDBC"
    target_type = "SINGLE_PROTEIN" if target_kind == "protein" else "PROTEIN_COMPLEX"
    guard = "AND coalesce(r.enzyme_complexid, 0) = 0" if target_kind == "protein" else ""
    _temporary_ids(connection, table, column, proteins[column].tolist())
    ligand_join = ""
    if monomer_ids is not None:
        _temporary_ids(connection, "selected_ligands", "monomerid", monomer_ids)
        ligand_join = "JOIN selected_ligands lig ON lig.monomerid = r.inhibitor_monomerid"
    measurements = _frame(connection, f"""SELECT k.*, r.enzyme_{column} AS target_entity_id,
                        r.inhibitor_monomerid AS ligand_monomerid,
                        a.assay_name, a.description AS assay_description
                    FROM enzyme_reactant_set r JOIN {table} p ON p.{column} = r.enzyme_{column}
                    JOIN ki_result k ON k.reactant_set_id = r.reactant_set_id AND k.entryid = r.entryid
                    {ligand_join}
                    LEFT JOIN assay a ON a.assayid = k.assayid AND a.entryid = k.entryid
                    WHERE r.inhibitor_monomerid > 0 {guard}
                    ORDER BY k.ki_result_id""")
    rows = []
    for raw in measurements.astype(object).where(pd.notna(measurements), None).to_dict("records"):
        base = {"measurement_id": f"ki:{raw['ki_result_id']}", "measurement_source": "ki_result",
                "target_id": f"{prefix}{int(raw['target_entity_id'])}", "target_type": target_type,
                "polymerid": int(raw["target_entity_id"]) if target_kind == "protein" else None,
                "complexid": int(raw["target_entity_id"]) if target_kind == "complex" else None,
                "monomerid": int(raw["ligand_monomerid"]),
                "compound_id": f"BDBM{int(raw['ligand_monomerid'])}", "entryid": raw["entryid"],
                "assayid": raw["assayid"], "assay_name": raw["assay_name"],
                "assay_description": raw["assay_description"], "raw_record_json": json.dumps(raw, default=str),
                "reactant_set_json": None}
        found = False
        for endpoint in ENDPOINTS:
            if raw.get(endpoint) is not None and str(raw[endpoint]).strip():
                rows.append({**base, "endpoint": endpoint.upper(), **parse_measurement(raw[endpoint], raw.get(endpoint + "_uncert"))})
                found = True
        if not found:
            rows.append({**base, "endpoint": "OTHER", **parse_measurement(None), "units": None})
    # ITC entries have their own identifiers and retain their source units/fields.
    tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if "itc_result_a_b_ab" in tables:
        itc_guard = "AND coalesce(i.cell_complexid,0)=0 AND coalesce(i.syr_complexid,0)=0" if target_kind == "protein" else ""
        itc = _frame(connection, f"""SELECT i.*, CASE WHEN i.cell_{column} = p.{column} THEN i.syr_monomerid
                               ELSE i.cell_monomerid END AS ligand_monomerid,
                               p.{column} AS target_entity_id
                           FROM itc_result_a_b_ab i JOIN {table} p
                           ON ((i.cell_{column} = p.{column} AND i.syr_monomerid > 0)
                           OR (i.syr_{column} = p.{column} AND i.cell_monomerid > 0)) {itc_guard}""")
        allowed = set(monomer_ids) if monomer_ids is not None else None
        for raw in itc.astype(object).where(pd.notna(itc), None).to_dict("records"):
            mid = int(raw["ligand_monomerid"])
            if allowed is not None and mid not in allowed:
                continue
            rows.append({"measurement_id": f"itc:{raw['entryid']}:{raw['itc_result_a_b_ab_id']}",
                         "measurement_source": "itc_result_a_b_ab",
                         "target_id": f"{prefix}{int(raw['target_entity_id'])}", "target_type": target_type,
                         "polymerid": int(raw["target_entity_id"]) if target_kind == "protein" else None,
                         "complexid": int(raw["target_entity_id"]) if target_kind == "complex" else None,
                         "monomerid": mid, "compound_id": f"BDBM{mid}", "entryid": raw["entryid"],
                         "assayid": None, "endpoint": "ITC_K_RAW", "value_raw": raw.get("k"),
                         "relation_raw": raw.get("k_uncert"), "relation": None, "value": None,
                         "units": None, "pactivity": None, "value_status": "source_units_not_normalized",
                         "assay_name": None, "assay_description": raw.get("comments"),
                         "raw_record_json": json.dumps(raw, default=str), "reactant_set_json": None})
    frame = pd.DataFrame(rows, columns=ACTIVITY_COLUMNS)
    for identifier in ("polymerid", "complexid"):
        frame[identifier] = pd.array(frame[identifier], dtype="Int64")
    if not measurements.empty:
        _temporary_ids(connection, "selected_reactants", "reactant_set_id", measurements.reactant_set_id.tolist())
        sets = _frame(connection, """SELECT r.* FROM enzyme_reactant_set r
                         JOIN selected_reactants s ON s.reactant_set_id=r.reactant_set_id""")
        lookup = {int(row["reactant_set_id"]): json.dumps(row, default=str)
                  for row in sets.astype(object).where(pd.notna(sets), None).to_dict("records")}
        for index, row in frame[frame.measurement_source == "ki_result"].iterrows():
            source = json.loads(row.raw_record_json)
            frame.at[index, "reactant_set_json"] = lookup.get(int(source["reactant_set_id"]))
    return frame


def _human_complexes(connection) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Find complexes with at least one explicitly human polymer component.

    Mixed-species and missing-taxonomy components remain visible. Measurements
    are joined to the complex once, never multiplied by its component count.
    """
    complexes = _frame(connection, """SELECT x.* FROM complex x WHERE EXISTS (
        SELECT 1 FROM complex_component cc JOIN polymer p ON p.polymerid=cc.polymerid
        WHERE cc.complexid=x.complexid AND trim(p.taxid)='9606') ORDER BY x.complexid""")
    _temporary_ids(connection, "human_complex_ids", "complexid", complexes.complexid.tolist())
    components = _frame(connection, """SELECT cc.*, p.display_name AS protein_name,
        p.unpid1, p.unpid2, p.sequence, p.taxid AS polymer_taxid, p.scientific_name,
        p.type AS polymer_type, m.display_name AS monomer_name
        FROM complex_component cc JOIN human_complex_ids s ON s.complexid=cc.complexid
        LEFT JOIN polymer p ON p.polymerid=cc.polymerid
        LEFT JOIN monomer m ON m.monomerid=cc.monomerid ORDER BY cc.complexid, cc.componentid""")
    annotations = []
    for row in complexes.to_dict("records"):
        polymers = components[(components.complexid == row["complexid"]) & (components.polymerid.fillna(0) > 0)]
        taxids = polymers.polymer_taxid.fillna("").astype(str).str.strip()
        known = taxids.str.fullmatch(r"[1-9][0-9]*")
        human_count = int((taxids == "9606").sum())
        nonhuman_count = int((known & (taxids != "9606")).sum())
        unknown_count = int((~known).sum())
        scope = "mixed_species" if nonhuman_count else ("partially_unknown" if unknown_count else "all_annotated_polymers_human")
        annotations.append({**row, "human_component_count": human_count,
                            "nonhuman_component_count": nonhuman_count, "unknown_taxid_component_count": unknown_count,
                            "component_taxids": ";".join(sorted(set(taxids[known]))), "human_scope": scope})
    columns = list(complexes.columns) + ["human_component_count", "nonhuman_component_count", "unknown_taxid_component_count",
                                        "component_taxids", "human_scope"]
    return pd.DataFrame(annotations, columns=columns), components


def select_best_hits(candidates: pd.DataFrame, tie_policy: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pick the best passing genome hit per BindingDB polymer; report exact ties."""
    if tie_policy not in {"unassigned", "all"}:
        raise ValueError("tie_policy must be 'unassigned' or 'all'")
    if candidates.empty:
        return candidates.copy(), candidates.copy()
    candidates = candidates.drop_duplicates(["polymerid", "gene"]).copy()
    best_score = candidates.groupby("polymerid").bitscore.transform("max")
    best = candidates[candidates.bitscore == best_score].copy()
    best_evalue = best.groupby("polymerid").evalue.transform("min")
    best = best[best.evalue == best_evalue].copy()
    ambiguous = best[best.duplicated("polymerid", keep=False)].sort_values(["polymerid", "gene"])
    if tie_policy == "unassigned":
        best = best[~best.polymerid.isin(ambiguous.polymerid)]
    return best.sort_values(["polymerid", "gene"]).reset_index(drop=True), ambiguous.reset_index(drop=True)


def map_proteins(proteins_fasta: str | Path, targets: pd.DataFrame, workdir: Path,
                 min_identity: float = 90, min_coverage: float = 70, threads: int = 4) -> tuple[pd.DataFrame, list[str]]:
    if not 0 < min_identity <= 100 or not 0 < min_coverage <= 100 or threads < 1:
        raise ValueError("Identity/coverage must be in (0,100]; threads must be positive")
    records = list(SeqIO.parse(str(proteins_fasta), "fasta"))
    genes = [record.id for record in records]
    if not records or len(set(genes)) != len(genes):
        raise ValueError("Protein FASTA must contain unique, nonempty locus_tag identifiers")
    workdir.mkdir(parents=True, exist_ok=True)
    target_fasta = workdir / "bindingdb_targets.faa"
    usable = targets.copy()
    usable["sequence"] = usable.sequence.map(_normalized_sequence)
    usable = usable[usable.sequence.str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+")].copy()
    with target_fasta.open("w") as handle:
        for row in usable.itertuples():
            handle.write(f">BDBP{row.polymerid}\n{_normalized_sequence(row.sequence)}\n")
    if usable.empty:
        return pd.DataFrame(columns=MAPPING_COLUMNS), genes
    makeblastdb = _environment_tool("makeblastdb")
    blastp = _environment_tool("blastp")
    prefix = workdir / "genome_proteins"
    subprocess.run([makeblastdb, "-in", str(Path(proteins_fasta).resolve()), "-dbtype", "prot", "-out", str(prefix)],
                   check=True, capture_output=True, text=True)
    raw_path = workdir / "blastp.tsv"
    # Query each BindingDB donor against the genome. Retrieve candidates before
    # filtering; max_target_seqs=1 could discard the best hit passing coverage.
    subprocess.run([blastp, "-query", str(target_fasta), "-db", str(prefix),
                    "-out", str(raw_path), "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore",
                    "-max_hsps", "1", "-max_target_seqs", str(max(5, len(genes))), "-evalue", "1e-5",
                    "-num_threads", str(threads)], check=True, capture_output=True, text=True)
    rows = []
    donor_taxids = dict(zip(usable.polymerid, usable.taxid))
    with raw_path.open() as handle:
        for line in handle:
            q, s, identity, length, qs, qe, ss, se, qlen, slen, evalue, bitscore = line.rstrip().split("\t")
            qcov = 100 * (abs(int(qe) - int(qs)) + 1) / int(qlen)
            scov = 100 * (abs(int(se) - int(ss)) + 1) / int(slen)
            if float(identity) >= min_identity and min(qcov, scov) >= min_coverage:
                rows.append({"gene": s, "polymerid": int(q.removeprefix("BDBP")), "identity": float(identity),
                             "bindingdb_taxid": donor_taxids[int(q.removeprefix("BDBP"))],
                             "genome_coverage": scov, "bindingdb_coverage": qcov, "evalue": float(evalue),
                             "bitscore": float(bitscore), "mapping_method": "blastp", "best_hit_tie_count": 1})
    candidates = pd.DataFrame(rows, columns=MAPPING_COLUMNS).drop_duplicates(["gene", "polymerid"])
    best, ties = select_best_hits(candidates, tie_policy="all")
    if not best.empty:
        best["best_hit_tie_count"] = best.groupby("polymerid").gene.transform("size")
    if not ties.empty:
        ties["best_hit_tie_count"] = ties.groupby("polymerid").gene.transform("size")
    candidates.to_csv(workdir / "bindingdb_passing_hits.tsv", sep="\t", index=False)
    ties.to_csv(workdir / "bindingdb_best_hit_ties.tsv", sep="\t", index=False)
    LOG.info("Mapping: %s passing pairs, %s selected pairs, %s donors with exact best-hit ties",
             len(candidates), len(best), ties.polymerid.nunique())
    return best, genes


def _export(root: Path, name: str, frame: pd.DataFrame) -> None:
    for suffix, separator in (("tsv", "\t"), ("csv", ",")):
        destination = root / f"{name}.{suffix}"
        partial = destination.with_name(destination.name + ".part")
        frame.to_csv(partial, sep=separator, index=False)
        partial.replace(destination)


def run_bindingdb(database_path: str | Path, taxid: int, output_dir: str | Path, *,
                  proteins_fasta: str | Path | None = None, taxonomy_dir: str | Path | None = None,
                  taxids: list[int] | None = None,
                  include_human: bool = False,
                  min_identity: float = 90,
                  min_coverage: float = 70, threads: int = 4) -> dict:
    """Extract species evidence and optionally map it to locus tags; no API calls."""
    source = Path(database_path).expanduser().resolve()
    if source.suffix.lower() in {".dmp", ".sql"}:
        source = prepare_bindingdb(source)
    if not source.is_file():
        raise FileNotFoundError(source)
    expanded, taxonomy = expand_taxids(taxid, taxonomy_dir, taxids)
    root = Path(output_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    metadata = {"status": "running", "taxid": int(taxid), "taxids": expanded, "taxonomy": taxonomy,
                "database": source.name, "include_human": include_human,
                "mapping": {"min_identity": min_identity, "min_query_and_target_coverage": min_coverage,
                            "threads": threads, "all_passing_hits": False,
                            "direction": "bindingdb_to_genome", "selection": "highest_bitscore_then_lowest_evalue_per_bindingdb_protein",
                            "ties": "assign_all_equal_best", "coverage_scope": ["genome", "bindingdb"]},
                "started_at": datetime.now(timezone.utc).isoformat(),
                "target_selection": "Individual polymers with valid amino-acid sequences; sequences restricted to A/C/G/T/U/N excluded as ambiguous; type is descriptive only; complex measurements excluded",
                "activity_categories": "not_assigned", "compound_identity": "BindingDB monomerid; no parent collapse",
                "human_complex_policy": "Include complexes with at least one component explicitly annotated with taxid 9606; preserve mixed/unknown component taxonomy",
                "organism_assays": {"status": "deferred", "reason": "No reliable assay organism taxid in the inspected dump"}}
    prep_metadata = source.with_suffix(source.suffix + ".metadata.json")
    if prep_metadata.is_file():
        metadata["database_preparation"] = json.loads(prep_metadata.read_text())
    _json(root / "run_metadata.json", metadata)
    connection = sqlite3.connect(source.as_uri() + "?mode=ro", uri=True)
    try:
        proteins = _protein_rows(connection, expanded)
        LOG.info("Taxids: %s; single protein targets: %s", len(expanded), len(proteins))
        activities = _activities(connection, proteins)
        monomers = sorted(set(activities.monomerid.tolist()))
        _temporary_ids(connection, "bacterial_ligands", "monomerid", monomers)
        ligand_source = _frame(connection, "SELECT m.* FROM monomer m JOIN bacterial_ligands b ON b.monomerid=m.monomerid ORDER BY m.monomerid")
        ligand_rows = []
        for raw in ligand_source.to_dict("records"):
            evidence = activities[activities.monomerid == raw["monomerid"]]
            ligand_rows.append({**raw, "compound_id": f"BDBM{raw['monomerid']}",
                                "bacterial_target_count": evidence.polymerid.nunique(),
                                "bacterial_measurement_count": evidence.measurement_id.nunique(),
                                "bacterial_quantitative_measurement_count": evidence[
                                    evidence.value_status.isin(["exact", "bound_or_approximate"])].measurement_id.nunique()})
        ligands = pd.DataFrame(ligand_rows, columns=LIGAND_COLUMNS)
        if len(ligands) != len(monomers):
            raise ValueError("Activity records reference compounds missing from the monomer table")
        if proteins_fasta:
            mapping, genes = map_proteins(proteins_fasta, proteins, root / "mapping", min_identity, min_coverage, threads)
            metadata["mapping"]["proteins_sha256"] = _sha256(Path(proteins_fasta))
        else:
            mapping, genes = pd.DataFrame(columns=MAPPING_COLUMNS), []
        pairs = mapping.merge(activities, on="polymerid", how="inner")
        human_proteins = pd.DataFrame(columns=PROTEIN_COLUMNS)
        human = pd.DataFrame(columns=ACTIVITY_COLUMNS)
        human_complexes = pd.DataFrame(columns=["complexid", "display_name", "human_scope"])
        human_components = pd.DataFrame(columns=["complexid", "componentid", "polymerid", "monomerid", "polymer_taxid"])
        cross = pd.DataFrame(columns=["compound_id", "bacterial_polymerid", "human_target_id", "human_target_type", "human_polymerid", "human_complexid"])
        if include_human:
            human_proteins = _protein_rows(connection, [9606])
            human = _activities(connection, human_proteins, monomers)
            human_proteins = human_proteins[human_proteins.polymerid.isin(human.polymerid)]
            human_complexes, human_components = _human_complexes(connection)
            complex_activities = _activities(connection, human_complexes, monomers, target_kind="complex")
            if not complex_activities.empty:
                human = pd.concat([human, complex_activities], ignore_index=True)
            human_complexes = human_complexes[human_complexes.complexid.isin(human.complexid)]
            human_components = human_components[human_components.complexid.isin(human_complexes.complexid)]
            bacterial_pairs = activities[["compound_id", "polymerid"]].drop_duplicates().rename(columns={"polymerid": "bacterial_polymerid"})
            human_pairs = human[["compound_id", "target_id", "target_type", "polymerid", "complexid"]].drop_duplicates().rename(
                columns={"target_id": "human_target_id", "target_type": "human_target_type",
                         "polymerid": "human_polymerid", "complexid": "human_complexid"})
            cross = bacterial_pairs.merge(human_pairs, on="compound_id", how="inner")
        for column, field, scope in (("human_target_count", "target_id", None),
                                     ("human_protein_target_count", "target_id", "SINGLE_PROTEIN"),
                                     ("human_complex_target_count", "target_id", "PROTEIN_COMPLEX"),
                                     ("human_measurement_count", "measurement_id", None)):
            if include_human:
                evidence = human if scope is None else human[human.target_type == scope]
                counts = evidence.groupby("monomerid")[field].nunique()
                ligands[column] = ligands.monomerid.map(counts).fillna(0).astype("Int64")
            else:
                ligands[column] = pd.Series(pd.NA, index=ligands.index, dtype="Int64")
        summaries = []
        for gene in genes:
            gene_mapping = mapping[mapping.gene == gene]
            evidence = pairs[pairs.gene == gene]
            compounds = set(evidence.monomerid)
            summaries.append({"gene": gene, "bindingdb_status": "no_target_mapping" if gene_mapping.empty else ("no_activities" if evidence.empty else "found"),
                              "bindingdb_target_count": gene_mapping.polymerid.nunique(),
                              "bindingdb_compounds": ";".join(f"BDBM{mid}" for mid in sorted(compounds)),
                              "bindingdb_compound_count": len(compounds)})
        summary = pd.DataFrame(summaries, columns=SUMMARY_COLUMNS)
        entry_ids = set(activities.entryid.dropna()) | set(human.entryid.dropna())
        _temporary_ids(connection, "selected_entries", "entryid", entry_ids)
        references = _frame(connection, """SELECT c.*, a.* FROM entry_citation c JOIN selected_entries s ON s.entryid=c.entryid
                              JOIN article a ON a.articleid=c.articleid""")
        references = references.loc[:, ~references.columns.duplicated()]
        for name, frame in (("bindingdb_proteins", proteins), ("bindingdb_activities", activities),
                            ("bindingdb_ligands", ligands), ("bindingdb_mapping", mapping),
                            ("bindingdb_protein_ligand_activities", pairs), ("bindingdb_summary", summary),
                            ("bindingdb_references", references)):
            _export(root, name, frame)
        # Always replace human outputs: disabling the option cannot leave stale evidence.
        for name, frame in (("bindingdb_human_proteins", human_proteins), ("bindingdb_human_activities", human),
                            ("bindingdb_human_crossactivity", cross), ("bindingdb_human_complexes", human_complexes),
                            ("bindingdb_human_complex_components", human_components)):
            _export(root, name, frame)
        if not proteins.empty:
            with (root / "bindingdb_proteins.faa").open("w") as handle:
                for row in proteins.itertuples():
                    if row.sequence:
                        handle.write(f">BDBP{row.polymerid}\n{_normalized_sequence(row.sequence)}\n")
        else:
            (root / "bindingdb_proteins.faa").write_text("")
        metadata.update(status="completed", finished_at=datetime.now(timezone.utc).isoformat(),
                        counts={"taxids": len(expanded), "targets": len(proteins),
                                "activity_rows": len(activities), "measurements": activities.measurement_id.nunique(),
                                "ligands": len(ligands), "input_proteins": len(genes), "mapped_proteins": mapping.gene.nunique(),
                                "mapping_pairs": len(mapping),
                                "donors_with_best_hit_ties": mapping[mapping.best_hit_tie_count > 1].polymerid.nunique(),
                                "human_activity_rows": len(human), "human_complex_targets": len(human_complexes)})
        if proteins.empty:
            LOG.warning("No single protein BindingDB targets for taxid %s and its selected descendants", taxid)
        LOG.info("BindingDB extraction completed: %s", metadata["counts"])
    except Exception as error:
        metadata.update(status="failed", error=str(error))
        _json(root / "run_metadata.json", metadata)
        raise
    finally:
        connection.close()
    _json(root / "run_metadata.json", metadata)
    return {"summary": summary, "metadata": metadata, "output_dir": root}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="Import the official MySQL dump once into local SQLite")
    prepare.add_argument("--dump", required=True)
    prepare.add_argument("--database")
    run = sub.add_parser("run", help="Extract species evidence and optionally map a locus-tag FASTA")
    run.add_argument("--database", required=True, help="Prepared SQLite or original .dmp/.sql")
    run.add_argument("--taxid", required=True, type=int)
    run.add_argument("--output-dir", required=True, help="Usually <organism>/ligands/bindingdb")
    run.add_argument("--proteins", help="Protein FASTA with locus_tag as each sequence ID")
    taxonomy = run.add_mutually_exclusive_group()
    taxonomy.add_argument("--taxonomy-dir", help="NCBI nodes.dmp, names.dmp and merged.dmp directory")
    taxonomy.add_argument("--taxids-file", help="Explicit taxids, one per line; completeness is not assumed")
    run.add_argument("--include-human", action="store_true")
    run.add_argument("--min-identity", type=float, default=90)
    run.add_argument("--min-coverage", type=float, default=70)
    run.add_argument("--threads", type=int, default=4)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    if args.command == "prepare":
        prepare_bindingdb(args.dump, args.database)
    else:
        taxids = [int(line.strip()) for line in Path(args.taxids_file).read_text().splitlines() if line.strip()] if args.taxids_file else None
        run_bindingdb(args.database, args.taxid, args.output_dir, proteins_fasta=args.proteins,
                      taxonomy_dir=args.taxonomy_dir, taxids=taxids,
                      include_human=args.include_human,
                      min_identity=args.min_identity,
                      min_coverage=args.min_coverage, threads=args.threads)


if __name__ == "__main__":
    main()
