
from ftscripts import programs, metadata, files, structures
from ftscripts.microbiome_catalogues import (
    catalogue_column_prefix,
    catalogue_species_path,
    get_catalogue,
)
import os
import json
import csv
import pandas as pd
import multiprocessing
import glob
import datetime
from tqdm import tqdm
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
import shutil
import subprocess
import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

MICROBIOME_BLAST_COLUMNS = 13
MICROBIOME_HIT_COLUMNS = (
    "gene",
    "representative_genome_id",
    "subject_protein_id",
    "pident",
    "qcovhsp",
    "evalue",
    "bitscore",
)
MICROBIOME_HIT_SCHEMA = pa.schema([
    pa.field("gene", pa.string()),
    pa.field("representative_genome_id", pa.string()),
    pa.field("subject_protein_id", pa.string()),
    pa.field("pident", pa.float64()),
    pa.field("qcovhsp", pa.float64()),
    pa.field("evalue", pa.float64()),
    pa.field("bitscore", pa.float64()),
])
MICROBIOME_EVALUE = "1e-5"
MICROBIOME_MAX_TARGET_SEQS = 1
MICROBIOME_MAX_HSPS = 1
MICROBIOME_PARQUET_BATCH_ROWS = 75000
MICROBIOME_PARQUET_ROW_GROUP_ROWS = 100000
MICROBIOME_TAXONOMY_RANKS = (
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


def _duckdb_identifier(identifier):
    """
    Escapes a column identifier for DuckDB SQL.

    :param identifier: Column identifier.
    :return: Double-quoted DuckDB identifier.
    """
    escaped_identifier = identifier.replace('"', '""')
    return f'"{escaped_identifier}"'


def _microbiome_rank_path(rank, table_alias=None):
    """
    Returns a DuckDB expression for the hierarchical taxonomy path of a rank.

    :param rank: Taxonomy rank to represent.
    :param table_alias: Optional SQL table alias to prepend to column names.
    :return: SQL expression joining the taxonomy path up to the selected rank.
    """
    rank_index = MICROBIOME_TAXONOMY_RANKS.index(rank)
    prefix = f"{table_alias}." if table_alias else ""
    return " || ';' || ".join(
        f"{prefix}{_duckdb_identifier(taxonomy_rank)}"
        for taxonomy_rank in MICROBIOME_TAXONOMY_RANKS[:rank_index + 1]
    )


def _microbiome_classified_path_filter(rank, table_alias=None):
    """
    Returns a DuckDB filter that keeps fully classified paths for a rank.

    :param rank: Taxonomy rank to filter.
    :param table_alias: Optional SQL table alias to prepend to column names.
    :return: SQL boolean expression excluding unclassified values in the path.
    """
    rank_index = MICROBIOME_TAXONOMY_RANKS.index(rank)
    prefix = f"{table_alias}." if table_alias else ""
    return " AND ".join(
        f"{prefix}{_duckdb_identifier(taxonomy_rank)} <> 'unclassified'"
        for taxonomy_rank in MICROBIOME_TAXONOMY_RANKS[:rank_index + 1]
    )


@dataclass(frozen=True)
class GenomeSearchResult:
    genome_id: str
    status: str
    output_path: str
    error: str = None


def _format_filter_value(value):
    """
    Formats a numeric filter value for use in result file names.

    :param value: Numeric filter value.
    :return: Compact string representation of the value.
    """
    return f"{float(value):g}"


def _validate_microbiome_blast_output(output_path):
    """
    Validates the column count of a DIAMOND microbiome result file.

    :param output_path: Path to the DIAMOND output TSV.
    :raises ValueError: If a row does not contain the expected number of columns.
    """
    with open(output_path, "r", encoding="utf-8") as output_file:
        for line_number, line in enumerate(output_file, start=1):
            if len(line.rstrip("\n").split("\t")) != MICROBIOME_BLAST_COLUMNS:
                raise ValueError(
                    f"Invalid DIAMOND output in {output_path} at line {line_number}: "
                    f"expected {MICROBIOME_BLAST_COLUMNS} columns."
                )


def _microbiome_catalogue_status(species_path):
    """
    Returns the available and indexed representatives of a catalogue.

    :param species_path: Path to the species catalogue directory.
    :return: Sets containing genome directories and indexed genome IDs.
    """
    genome_dirs = {
        entry for entry in os.listdir(species_path)
        if os.path.isdir(os.path.join(species_path, entry))
    }
    indexed_genomes = {
        genome for genome in genome_dirs
        if os.path.isfile(os.path.join(species_path, genome, f"{genome}_DB.dmnd"))
    }
    return genome_dirs, indexed_genomes


def _microbiome_results_path(output_path, organism_name, catalogue_name):
    """
    Returns the directory containing per-genome microbiome results.

    :param output_path: Path of the organism output.
    :param organism_name: Name of the organism.
    :param catalogue_name: Name of the microbiome catalogue.
    :return: Path to the species BLAST result directory.
    """
    return os.path.join(
        output_path,
        organism_name,
        "offtarget",
        "microbiomes",
        catalogue_name,
        "species_blast_results",
    )


def _microbiome_consolidated_paths(output_path, organism_name, catalogue_name):
    """
    Returns the paths used by consolidated microbiome results.

    :param output_path: Path of the organism output.
    :param organism_name: Name of the organism.
    :param catalogue_name: Name of the microbiome catalogue.
    :return: Result directory, consolidated Parquet path, and manifest path.
    """
    results_path = _microbiome_results_path(
        output_path,
        organism_name,
        catalogue_name,
    )
    return (
        results_path,
        os.path.join(results_path, f"{catalogue_name}_offtarget_hits.parquet"),
        os.path.join(results_path, f"{catalogue_name}_offtarget_manifest.json"),
    )

def _read_microbiome_manifest(manifest_path):
    """
    Reads a microbiome consolidation manifest.

    :param manifest_path: Path to the manifest JSON file.
    :return: Manifest dictionary, or None if the file is missing or invalid.
    """
    try:
        with open(manifest_path, "r", encoding="utf-8") as manifest_file:
            manifest = json.load(manifest_file)
    except (OSError, json.JSONDecodeError):
        return None
    return manifest if isinstance(manifest, dict) else None


def is_microbiome_consolidation_compatible(
    output_path,
    organism_name,
    catalogue_name,
    identity_filter,
    coverage_filter,
    verify_checksum=True,
):
    """
    Checks whether consolidated microbiome hits match a catalogue and thresholds.

    :param output_path: Path of the organism output.
    :param organism_name: Name of the organism.
    :param catalogue_name: Name of a supported MGnify catalogue.
    :param identity_filter: Identity threshold associated with the results.
    :param coverage_filter: Query coverage threshold associated with the results.
    :param verify_checksum: Whether to recompute and verify the Parquet SHA-256.
    :return: True if the Parquet and manifest are present, valid, and compatible.
    """
    _, parquet_path, manifest_path = _microbiome_consolidated_paths(
        output_path,
        organism_name,
        catalogue_name,
    )
    if not os.path.isfile(parquet_path) or not os.path.isfile(manifest_path):
        return False

    manifest = _read_microbiome_manifest(manifest_path)
    if manifest is None:
        return False

    expected_suffix = microbiome_result_suffix(
        identity_filter,
        coverage_filter,
    )
    catalogue = get_catalogue(catalogue_name)
    try:
        compatible = (
            manifest.get("catalogue") == catalogue_name
            and manifest.get("source_result_suffix") == expected_suffix
            and manifest.get("expected_representatives")
            == catalogue["number_of_species"]
            and _format_filter_value(manifest.get("identity_filter"))
            == _format_filter_value(identity_filter)
            and _format_filter_value(manifest.get("coverage_filter"))
            == _format_filter_value(coverage_filter)
            and manifest.get("hits_parquet") == os.path.basename(parquet_path)
        )
    except (TypeError, ValueError):
        return False
    if not compatible:
        return False

    try:
        parquet_rows = pq.ParquetFile(parquet_path).metadata.num_rows
        if parquet_rows != manifest.get("hit_rows"):
            return False
        return (
            not verify_checksum
            or files.sha256_file(parquet_path)
            == manifest.get("hits_parquet_sha256")
        )
    except (OSError, pa.ArrowException):
        return False


def _validate_microbiome_source_tsv(tsv_path):
    """
    Validates one per-genome DIAMOND result and counts its hit rows.

    Empty files are accepted as valid searches without hits.

    :param tsv_path: Path to the per-genome result TSV.
    :return: Number of non-empty result rows.
    :raises ValueError: If columns, numeric values, or query IDs are invalid.
    """
    row_count = 0
    query_ids = set()
    with open(tsv_path, "r", encoding="utf-8", newline="") as input_file:
        reader = csv.reader(input_file, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) != MICROBIOME_BLAST_COLUMNS:
                raise ValueError(
                    f"Invalid DIAMOND output in {tsv_path} at line {line_number}: "
                    f"expected {MICROBIOME_BLAST_COLUMNS} columns, found {len(row)}."
                )
            if not row[0]:
                raise ValueError(
                    f"Empty qseqid in {tsv_path} at line {line_number}."
                )
            if row[0] in query_ids:
                raise ValueError(
                    f"Duplicate qseqid {row[0]!r} in {tsv_path}; "
                    "--max-target-seqs 1 requires at most one row per query."
                )
            query_ids.add(row[0])
            try:
                for position in (2, 10, 11, 12):
                    float(row[position])
            except ValueError as error:
                raise ValueError(
                    f"Invalid numeric value in {tsv_path} at line {line_number}."
                ) from error
            row_count += 1
    return row_count


def _write_unsorted_microbiome_parquet(
    source_files,
    unsorted_path,
    batch_rows=MICROBIOME_PARQUET_BATCH_ROWS,
):
    """
    Writes per-genome microbiome hits incrementally to an unsorted Parquet file.

    :param source_files: Sequence of representative IDs and result TSV paths.
    :param unsorted_path: Path to the temporary unsorted Parquet file.
    :param batch_rows: Maximum number of rows accumulated before writing a batch.
    """
    columns = {column: [] for column in MICROBIOME_HIT_COLUMNS}

    def flush_batch(writer):
        """
        Writes the current in-memory hit batch and clears its columns.

        :param writer: Open PyArrow Parquet writer.
        """
        if not columns["gene"]:
            return
        writer.write_table(
            pa.Table.from_pydict(columns, schema=MICROBIOME_HIT_SCHEMA),
            row_group_size=MICROBIOME_PARQUET_ROW_GROUP_ROWS,
        )
        for values in columns.values():
            values.clear()

    with pq.ParquetWriter(
        unsorted_path,
        MICROBIOME_HIT_SCHEMA,
        compression="zstd",
    ) as writer:
        for genome_id, tsv_path in source_files:
            with open(tsv_path, "r", encoding="utf-8", newline="") as input_file:
                for row in csv.reader(input_file, delimiter="\t"):
                    if not row:
                        continue
                    columns["gene"].append(row[0])
                    columns["representative_genome_id"].append(genome_id)
                    columns["subject_protein_id"].append(row[1])
                    columns["pident"].append(float(row[2]))
                    columns["qcovhsp"].append(float(row[12]))
                    columns["evalue"].append(float(row[10]))
                    columns["bitscore"].append(float(row[11]))
                    if len(columns["gene"]) >= batch_rows:
                        flush_batch(writer)
        flush_batch(writer)


def _duckdb_path(file_path):
    """
    Escapes a file path for use in a DuckDB SQL string.

    :param file_path: File path to escape.
    :return: Absolute SQL-safe file path.
    """
    return os.path.abspath(file_path).replace("'", "''")


def _sort_microbiome_parquet(unsorted_path, sorted_path, temporary_directory):
    """
    Sorts consolidated microbiome hits using DuckDB external storage.

    :param unsorted_path: Path to the unsorted input Parquet file.
    :param sorted_path: Path to the sorted output Parquet file.
    :param temporary_directory: Directory available for DuckDB temporary data.
    """
    connection = duckdb.connect()
    try:
        connection.execute(
            f"SET temp_directory='{_duckdb_path(temporary_directory)}'"
        )
        connection.execute(
            f"""
            COPY (
                SELECT * FROM read_parquet('{_duckdb_path(unsorted_path)}')
                ORDER BY gene, representative_genome_id
            )
            TO '{_duckdb_path(sorted_path)}'
            (
                FORMAT PARQUET,
                COMPRESSION ZSTD,
                ROW_GROUP_SIZE {MICROBIOME_PARQUET_ROW_GROUP_ROWS}
            )
            """
        )
    finally:
        connection.close()


def _validate_consolidated_microbiome_parquet(
    parquet_path,
    expected_rows,
    expected_representative_ids,
):
    """
    Validates the schema and contents of a consolidated microbiome Parquet file.

    :param parquet_path: Path to the consolidated Parquet file.
    :param expected_rows: Expected number of hit rows.
    :param expected_representative_ids: Valid representative IDs for the catalogue.
    :return: Number of genes and representatives containing hits.
    :raises ValueError: If any schema or content validation fails.
    """
    parquet_file = pq.ParquetFile(parquet_path)
    if parquet_file.schema_arrow != MICROBIOME_HIT_SCHEMA:
        raise ValueError(
            f"Unexpected consolidated hit schema in {parquet_path}: "
            f"{parquet_file.schema_arrow}."
        )
    if parquet_file.metadata.num_rows != expected_rows:
        raise ValueError(
            f"Consolidated Parquet contains {parquet_file.metadata.num_rows} rows; "
            f"{expected_rows} source rows were expected."
        )

    previous_key = None
    out_of_order = 0
    for batch in parquet_file.iter_batches(
        columns=["gene", "representative_genome_id"],
        batch_size=MICROBIOME_PARQUET_BATCH_ROWS,
    ):
        genes = batch.column(0).to_pylist()
        representatives_in_batch = batch.column(1).to_pylist()
        for key in zip(genes, representatives_in_batch):
            if previous_key is not None and key < previous_key:
                out_of_order += 1
            previous_key = key

    connection = duckdb.connect()
    try:
        source = f"read_parquet('{_duckdb_path(parquet_path)}')"
        null_rows = connection.execute(
            f"""
            SELECT count(*) FROM {source}
            WHERE gene IS NULL OR representative_genome_id IS NULL
            """
        ).fetchone()[0]
        duplicate_pairs = connection.execute(
            f"""
            SELECT count(*) FROM (
                SELECT gene, representative_genome_id
                FROM {source}
                GROUP BY gene, representative_genome_id
                HAVING count(*) > 1
            )
            """
        ).fetchone()[0]
        representatives = {
            row[0] for row in connection.execute(
                f"SELECT DISTINCT representative_genome_id FROM {source}"
            ).fetchall()
        }
        genes_with_hits, representatives_with_hits = connection.execute(
            f"""
            SELECT
                count(DISTINCT gene),
                count(DISTINCT representative_genome_id)
            FROM {source}
            """
        ).fetchone()
    finally:
        connection.close()

    if null_rows:
        raise ValueError(
            f"Consolidated Parquet contains {null_rows} rows with null identifiers."
        )
    if duplicate_pairs:
        raise ValueError(
            f"Consolidated Parquet contains {duplicate_pairs} duplicate "
            "gene and representative_genome_id pairs."
        )
    unexpected_representatives = representatives - expected_representative_ids
    if unexpected_representatives:
        raise ValueError(
            f"Consolidated Parquet contains {len(unexpected_representatives)} "
            "representative IDs outside the catalogue."
        )
    if out_of_order:
        raise ValueError(
            f"Consolidated Parquet contains {out_of_order} out-of-order rows."
        )
    return genes_with_hits, representatives_with_hits


def microbiome_result_suffix(identity_filter, coverage_filter):
    """
    Creates the suffix used by per-genome microbiome result files.

    :param identity_filter: Identity threshold associated with the results.
    :param coverage_filter: Query coverage threshold associated with the results.
    :return: Result file suffix containing both thresholds.
    """
    identity = _format_filter_value(identity_filter)
    coverage = _format_filter_value(coverage_filter)
    return f"_offtarget_id{identity}_cov{coverage}.tsv"

def consolidate_microbiome_hits(
    databases_path,
    output_path,
    organism_name,
    catalogue_name,
    identity_filter,
    coverage_filter,
    delete_source_tsvs=False,
):
    """
    Consolidates validated per-genome DIAMOND results into a sorted Parquet file.

    :param databases_path: Path where microbiome catalogues are stored.
    :param output_path: Path of the organism output.
    :param organism_name: Name of the organism.
    :param catalogue_name: Name of a supported MGnify catalogue.
    :param identity_filter: Identity threshold associated with the results.
    :param coverage_filter: Query coverage threshold associated with the results.
    :param delete_source_tsvs: Whether to remove source TSVs after validation.
    :return: Path to the consolidated Parquet file.
    """
    _, existing_parquet_path, _ = _microbiome_consolidated_paths(
        output_path,
        organism_name,
        catalogue_name,
    )
    if is_microbiome_consolidation_compatible(
        output_path,
        organism_name,
        catalogue_name,
        identity_filter,
        coverage_filter,
    ):
        if delete_source_tsvs:
            cleanup_consolidated_microbiome_tsvs(
                databases_path,
                output_path,
                organism_name,
                catalogue_name,
                identity_filter,
                coverage_filter,
                delete_source_tsvs=True,
            )
        return existing_parquet_path

    catalogue = get_catalogue(catalogue_name)
    expected_count = catalogue["number_of_species"]
    species_path = catalogue_species_path(databases_path, catalogue_name)
    _, indexed_genomes = _microbiome_catalogue_status(species_path)
    if len(indexed_genomes) != expected_count:
        raise RuntimeError(
            f"{catalogue_name} has {len(indexed_genomes)} indexed representatives; "
            f"{expected_count} were expected."
        )

    taxonomy_path = os.path.join(
        species_path,
        "representative_taxonomy.parquet",
    )
    if not os.path.isfile(taxonomy_path):
        raise FileNotFoundError(
            f"Representative taxonomy not found: {taxonomy_path}"
        )
    taxonomy_table = pq.read_table(
        taxonomy_path,
        columns=["representative_genome_id"],
    )
    taxonomy_ids = set(
        taxonomy_table.column("representative_genome_id").to_pylist()
    )
    if len(taxonomy_ids) != expected_count or taxonomy_ids != indexed_genomes:
        raise RuntimeError(
            f"{catalogue_name} taxonomy and indexed representatives do not match."
        )

    results_path, parquet_path, manifest_path = _microbiome_consolidated_paths(
        output_path,
        organism_name,
        catalogue_name,
    )
    if not os.path.isdir(results_path):
        raise FileNotFoundError(
            f"Microbiome result directory not found: {results_path}"
        )

    temporary_files = [
        name for name in os.listdir(results_path)
        if name.endswith(".tmp")
    ]
    if temporary_files:
        raise RuntimeError(
            f"Temporary files found in {results_path}: "
            f"{', '.join(sorted(temporary_files)[:5])}."
        )

    result_suffix = microbiome_result_suffix(
        identity_filter,
        coverage_filter,
    )
    source_files = [
        (
            genome_id,
            os.path.join(results_path, f"{genome_id}{result_suffix}"),
        )
        for genome_id in sorted(indexed_genomes)
    ]
    missing_files = [
        tsv_path for _, tsv_path in source_files
        if not os.path.isfile(tsv_path)
    ]
    if missing_files:
        raise RuntimeError(
            f"{catalogue_name} is missing {len(missing_files)} expected result "
            f"files for suffix {result_suffix}. Examples: "
            f"{', '.join(missing_files[:5])}."
        )

    rows_by_representative = {}
    for genome_id, tsv_path in source_files:
        rows_by_representative[genome_id] = _validate_microbiome_source_tsv(
            tsv_path
        )
    total_rows = sum(rows_by_representative.values())
    empty_files = sum(count == 0 for count in rows_by_representative.values())

    unsorted_path = os.path.join(
        results_path,
        f"{catalogue_name}_offtarget_hits.unsorted.parquet.tmp",
    )
    sorted_temporary_path = f"{parquet_path}.tmp"
    duckdb_temporary_path = os.path.join(
        results_path,
        f".{catalogue_name}_duckdb.tmp",
    )
    try:
        _write_unsorted_microbiome_parquet(source_files, unsorted_path)
        os.makedirs(duckdb_temporary_path)
        _sort_microbiome_parquet(
            unsorted_path,
            sorted_temporary_path,
            duckdb_temporary_path,
        )
        genes_with_hits, representatives_with_hits = (
            _validate_consolidated_microbiome_parquet(
                sorted_temporary_path,
                total_rows,
                taxonomy_ids,
            )
        )
        parquet_sha256 = files.sha256_file(sorted_temporary_path)
        os.replace(sorted_temporary_path, parquet_path)
    except Exception:
        for temporary_path in (unsorted_path, sorted_temporary_path):
            if os.path.exists(temporary_path):
                os.remove(temporary_path)
        raise
    finally:
        if os.path.exists(unsorted_path):
            os.remove(unsorted_path)
        if os.path.isdir(duckdb_temporary_path):
            shutil.rmtree(duckdb_temporary_path)

    manifest = {
        "catalogue": catalogue_name,
        "catalogue_url": catalogue["ftp_site"],
        "created_at_utc": datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat(),
        "identity_filter": float(identity_filter),
        "coverage_filter": float(coverage_filter),
        "evalue": MICROBIOME_EVALUE,
        "max_target_seqs": MICROBIOME_MAX_TARGET_SEQS,
        "max_hsps": MICROBIOME_MAX_HSPS,
        "diamond_version": programs.diamond_version(),
        "expected_representatives": expected_count,
        "indexed_representatives": len(indexed_genomes),
        "searched_representatives": len(source_files),
        "empty_result_files": empty_files,
        "hit_rows": total_rows,
        "genes_with_hits": genes_with_hits,
        "representatives_with_hits": representatives_with_hits,
        "hits_parquet": os.path.basename(parquet_path),
        "hits_parquet_sha256": parquet_sha256,
        "source_result_suffix": result_suffix,
        "rows_by_representative": rows_by_representative,
        "raw_tsv_cleanup": (
            "pending" if delete_source_tsvs else "disabled"
        ),
    }
    files.atomic_write_json(manifest, manifest_path)

    if delete_source_tsvs:
        cleanup_consolidated_microbiome_tsvs(
            databases_path,
            output_path,
            organism_name,
            catalogue_name,
            identity_filter,
            coverage_filter,
            delete_source_tsvs=True,
        )

    print(
        f"Consolidated {total_rows} {catalogue_name} hits into {parquet_path}."
    )
    return parquet_path


def cleanup_consolidated_microbiome_tsvs(
    databases_path,
    output_path,
    organism_name,
    catalogue_name,
    identity_filter,
    coverage_filter,
    delete_source_tsvs=False,
):
    """
    Safely removes TSVs represented by a validated consolidated Parquet file.

    :param databases_path: Path where microbiome catalogues are stored.
    :param output_path: Path of the organism output.
    :param organism_name: Name of the organism.
    :param catalogue_name: Name of a supported MGnify catalogue.
    :param identity_filter: Identity threshold associated with the results.
    :param coverage_filter: Query coverage threshold associated with the results.
    :param delete_source_tsvs: Explicit authorization to remove source TSV files.
    :return: Number of TSV files removed in this invocation.
    """
    if not delete_source_tsvs:
        return 0

    results_path, parquet_path, manifest_path = _microbiome_consolidated_paths(
        output_path,
        organism_name,
        catalogue_name,
    )
    if not os.path.isfile(parquet_path) or not os.path.isfile(manifest_path):
        raise RuntimeError(
            "Cannot clean microbiome TSVs without both Parquet and manifest."
        )
    temporary_files = [
        name for name in os.listdir(results_path)
        if name.endswith(".tmp")
    ]
    if temporary_files:
        raise RuntimeError(
            f"Cannot clean microbiome TSVs while temporary files exist in "
            f"{results_path}."
        )

    manifest = _read_microbiome_manifest(manifest_path)
    if manifest is None:
        raise RuntimeError(f"Invalid microbiome manifest: {manifest_path}")
    if not is_microbiome_consolidation_compatible(
        output_path,
        organism_name,
        catalogue_name,
        identity_filter,
        coverage_filter,
    ):
        raise RuntimeError(
            "Consolidated microbiome hits are not compatible with the requested "
            "catalogue and thresholds."
        )

    rows_by_representative = manifest.get("rows_by_representative")
    if not isinstance(rows_by_representative, dict):
        raise RuntimeError(
            "Microbiome manifest does not contain per-representative row counts."
        )
    if any(
        isinstance(count, bool)
        or not isinstance(count, int)
        or count < 0
        for count in rows_by_representative.values()
    ):
        raise RuntimeError("Manifest contains invalid per-representative row counts.")
    if sum(rows_by_representative.values()) != manifest.get("hit_rows"):
        raise RuntimeError(
            "Manifest hit count does not match per-representative row counts."
        )

    species_path = catalogue_species_path(databases_path, catalogue_name)
    _, indexed_genomes = _microbiome_catalogue_status(species_path)
    expected_count = get_catalogue(catalogue_name)["number_of_species"]
    if (
        len(indexed_genomes) != expected_count
        or set(rows_by_representative) != indexed_genomes
    ):
        raise RuntimeError(
            "Manifest representatives do not match the indexed catalogue."
        )

    result_suffix = microbiome_result_suffix(
        identity_filter,
        coverage_filter,
    )
    if manifest.get("source_result_suffix") != result_suffix:
        raise RuntimeError("Manifest result suffix does not match thresholds.")

    removed = 0
    for genome_id in rows_by_representative:
        tsv_path = os.path.join(
            results_path,
            f"{genome_id}{result_suffix}",
        )
        if os.path.isfile(tsv_path):
            os.remove(tsv_path)
            removed += 1

    manifest["raw_tsv_cleanup"] = "completed"
    files.atomic_write_json(manifest, manifest_path)
    return removed


def _available_cpus():
    """
    Returns the number of CPUs available to the current process.

    :return: Number of available CPUs.
    """
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return multiprocessing.cpu_count()


def validate_microbiome_search_environment(query_faa, output_dir):
    """
    Validates resources shared by every representative-genome search.
    """

    if not files.file_check(query_faa):
        raise RuntimeError(f"Query protein FASTA not found or empty: {query_faa}")
    if not os.access(query_faa, os.R_OK):
        raise RuntimeError(f"Query protein FASTA is not readable: {query_faa}")
    if shutil.which("diamond") is None:
        raise RuntimeError("DIAMOND executable was not found in PATH.")

    os.makedirs(output_dir, exist_ok=True)
    if not os.access(output_dir, os.W_OK):
        raise RuntimeError(f"Output directory is not writable: {output_dir}")


def _is_systemic_search_error(error):
    if isinstance(error, OSError):
        return True
    if isinstance(error, subprocess.CalledProcessError):
        message = f"{error.stderr or ''} {error.stdout or ''}".lower()
        systemic_messages = (
            "no space left on device",
            "permission denied",
            "read-only file system",
            "cannot allocate memory",
        )
        return any(text in message for text in systemic_messages)
    return False


def create_microbiome_shards(species_path, output_dir, shard_size):
    """
    Creates size-balanced shard files containing indexed representative genome IDs.

    :return: List of shard file paths.
    """

    if isinstance(shard_size, bool) or int(shard_size) < 1:
        raise ValueError("shard_size must be a positive integer.")

    _, indexed_genomes = _microbiome_catalogue_status(species_path)
    if not indexed_genomes:
        raise RuntimeError(f"No indexed representative genomes found in {species_path}.")

    genomes = []
    for genome_id in indexed_genomes:
        faa_path = os.path.join(species_path, genome_id, f"{genome_id}.faa")
        if not files.file_check(faa_path):
            raise FileNotFoundError(
                f"Representative protein FASTA not found: {faa_path}"
            )
        genomes.append((genome_id, os.path.getsize(faa_path)))

    shard_count = max(1, (len(genomes) + int(shard_size) - 1) // int(shard_size))
    shards = [{"size": 0, "genomes": []} for _ in range(shard_count)]
    for genome_id, faa_size in sorted(genomes, key=lambda item: item[1], reverse=True):
        shard = min(shards, key=lambda item: (item["size"], len(item["genomes"])))
        shard["genomes"].append(genome_id)
        shard["size"] += faa_size

    os.makedirs(output_dir, exist_ok=True)
    shard_paths = []
    for index, shard in enumerate(shards, start=1):
        shard_path = os.path.join(output_dir, f"shard_{index:04d}.txt")
        with open(shard_path, "w", encoding="utf-8") as shard_file:
            for genome_id in sorted(shard["genomes"]):
                shard_file.write(f"{genome_id}\n")
        shard_paths.append(shard_path)

    return shard_paths


def search_one_genome(
    genome_id,
    genome_db,
    query_faa,
    output_path,
    identity_filter,
    coverage_filter,
    threads=4,
):
    """
    Runs and validates one DIAMOND search against a representative genome.

    :return: GenomeSearchResult with success, skipped, error, or system_error status.
    """

    temporary_output_path = f"{output_path}.tmp"

    try:
        if os.path.exists(output_path):
            try:
                _validate_microbiome_blast_output(output_path)
                return GenomeSearchResult(genome_id, "skipped", output_path)
            except ValueError:
                os.remove(output_path)

        if os.path.exists(temporary_output_path):
            os.remove(temporary_output_path)

        programs.run_diamond_blastp(
            blastdb=genome_db,
            query=query_faa,
            output=temporary_output_path,
            evalue=MICROBIOME_EVALUE,
            max_hsps=MICROBIOME_MAX_HSPS,
            outfmt=(
                "6 qseqid sseqid pident length mismatch gapopen qstart qend "
                "sstart send evalue bitscore qcovhsp"
            ),
            cpus=threads,
            identity=identity_filter,
            query_cover=coverage_filter,
            max_target_seqs=MICROBIOME_MAX_TARGET_SEQS,
        )
        _validate_microbiome_blast_output(temporary_output_path)
        os.replace(temporary_output_path, output_path)
        return GenomeSearchResult(genome_id, "success", output_path)
    except Exception as error:
        if os.path.exists(temporary_output_path):
            os.remove(temporary_output_path)
        return GenomeSearchResult(
            genome_id,
            "system_error" if _is_systemic_search_error(error) else "error",
            output_path,
            str(error),
        )


def _warn_incomplete_microbiome_catalogue(
    catalogue_name,
    expected_genomes,
    genome_dirs,
    indexed_genomes,
):
    if len(genome_dirs) != expected_genomes:
        logging.warning(
            "The %s catalogue contains %d genome directories; %d were expected.",
            catalogue_name,
            len(genome_dirs),
            expected_genomes,
        )
        print(
            f"Warning: The {catalogue_name} catalogue contains {len(genome_dirs)} "
            f"genome directories; {expected_genomes} were expected."
        )
    if len(indexed_genomes) != expected_genomes:
        logging.warning(
            "Only %d of %d expected %s genomes have a DIAMOND index.",
            len(indexed_genomes),
            expected_genomes,
            catalogue_name,
        )
        print(
            f"Warning: Only {len(indexed_genomes)} of {expected_genomes} expected "
            f"{catalogue_name} genomes have a DIAMOND index."
        )


def human_offtarget_blast (databases_path, output_path, organism_name, cpus=multiprocessing.cpu_count()):

    """
    Runs NCBI BLASTP against the human proteome.
    This function uses the `run_blastp` function from the `programs` module.
    It uses the HUMAN_DB database created by the `index_db_blast_human` function from the `databases` module.
    The blast output is saved in the 'offtarget' folder of the organism directory.

    :param databases_path: Directory where HUMAN databases are stored.
    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param cpus: Number of threads (CPUs) to use in the blast search.
    
    """

    #Database files
    humanprot_index_path = os.path.join(databases_path, 'HUMAN_DB')

    #Organism files
    organism_path = os.path.join(output_path, organism_name)
    organism_prot_seq_path = os.path.join(organism_path, f'genome/{organism_name}.faa')

    offtarget_path = os.path.join(organism_path, 'offtarget')
    blast_output_path = os.path.join(offtarget_path, 'human_offtarget_blast.tsv')

    programs.run_blastp(
        blastdb= humanprot_index_path,
        query= organism_prot_seq_path,
        output=blast_output_path,
        evalue= '1e-5',
        outfmt= '6 std qcovhsp qcovs',
        cpus=cpus
    )

def microbiome_offtarget_blast_species(
    databases_path,
    output_path,
    organism_name,
    catalogue_name,
    identity_filter,
    coverage_filter,
    cpus=multiprocessing.cpu_count(),
    threads_per_genome=4,
):
    """
    Runs Diamond BLASTP of the organism proteome against each genome in the microbiome species catalogue.
    Each representative genome is stored under
    `databases/microbiomes/<catalogue>/species_catalogue`, with its DIAMOND index.

    For each genome, the BLAST output is stored in the 'offtarget' folder of the organism directory,
    with one result file per genome. The process can be resumed if interrupted.

    :param databases_path: Path where the species catalogue databases are stored.
    :param output_path: Path of the organism output.
    :param organism_name: Name of the organism (folder name under 'organism').
    :param catalogue_name: Name of a supported MGnify catalogue.
    :param identity_filter: Identity threshold associated with the result files.
    :param coverage_filter: Query coverage threshold associated with the result files.
    :param cpus: Total CPU budget for concurrent searches.
    :param threads_per_genome: Maximum DIAMOND threads assigned to each genome.
    """

    catalogue = get_catalogue(catalogue_name)
    expected_genomes = catalogue["number_of_species"]
    species_databases_path = catalogue_species_path(databases_path, catalogue_name)

    if is_microbiome_consolidation_compatible(
        output_path,
        organism_name,
        catalogue_name,
        identity_filter,
        coverage_filter,
    ):
        print(
            f"Compatible consolidated {catalogue_name} results already exist; "
            "skipping DIAMOND searches."
        )
        return []

    # Path to organism proteome (.faa file)
    organism_path = os.path.join(output_path, organism_name)
    organism_prot_seq_path = os.path.join(organism_path, "genome", f"{organism_name}.faa")

    # Output folder
    offtarget_path = os.path.join(
        organism_path,
        "offtarget",
        "microbiomes",
        catalogue_name,
        "species_blast_results",
    )
    os.makedirs(offtarget_path, exist_ok=True)
    validate_microbiome_search_environment(
        organism_prot_seq_path,
        offtarget_path,
    )

    genome_dirs, indexed_genomes = _microbiome_catalogue_status(species_databases_path)
    _warn_incomplete_microbiome_catalogue(
        catalogue_name,
        expected_genomes,
        genome_dirs,
        indexed_genomes,
    )

    if not indexed_genomes:
        raise RuntimeError(
            f"No indexed microbiome genomes were found in {species_databases_path}."
        )

    result_suffix = microbiome_result_suffix(identity_filter, coverage_filter)
    available_cpus = _available_cpus()
    cpu_budget = max(1, min(int(cpus), available_cpus))
    genome_threads = max(1, min(int(threads_per_genome), cpu_budget))
    parallel_genomes = max(1, cpu_budget // genome_threads)

    print(
        f"{catalogue_name}: {parallel_genomes} concurrent searches, "
        f"{genome_threads} DIAMOND threads per genome "
        f"({cpu_budget} CPUs available to this module)."
    )

    search_arguments = []
    for genome_dir in sorted(indexed_genomes):
        genome_path = os.path.join(species_databases_path, genome_dir)
        blast_output_path = os.path.join(offtarget_path, f"{genome_dir}{result_suffix}")
        genome_db = os.path.join(genome_path, f'{genome_dir}_DB')
        search_arguments.append((genome_dir, genome_db, blast_output_path))

    results = []
    systemic_failure = None
    executor = ThreadPoolExecutor(max_workers=parallel_genomes)
    try:
        futures = {
            executor.submit(
                search_one_genome,
                genome_id,
                genome_db,
                organism_prot_seq_path,
                blast_output_path,
                identity_filter,
                coverage_filter,
                genome_threads,
            ): genome_id
            for genome_id, genome_db, blast_output_path in search_arguments
        }
        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc=f"Searching {catalogue_name}",
        ):
            result = future.result()
            results.append(result)
            if result.status == "system_error":
                systemic_failure = result
                for pending_future in futures:
                    if pending_future is not future:
                        pending_future.cancel()
                break
    finally:
        executor.shutdown(wait=True, cancel_futures=True)

    if systemic_failure is not None:
        raise RuntimeError(
            f"Systemic DIAMOND search failure for {systemic_failure.genome_id}: "
            f"{systemic_failure.error}"
        )

    status_counts = {
        status: sum(result.status == status for result in results)
        for status in ("success", "skipped", "error")
    }
    print(
        f"{catalogue_name} search summary: {status_counts['success']} completed, "
        f"{status_counts['skipped']} skipped, {status_counts['error']} failed."
    )

    failed = [result for result in results if result.status == "error"]
    if failed:
        examples = "; ".join(
            f"{result.genome_id}: {result.error}"
            for result in failed[:5]
        )
        raise RuntimeError(
            f"{len(failed)} {catalogue_name} DIAMOND searches failed. {examples}"
        )

    if len(results) != len(indexed_genomes):
        raise RuntimeError(
            f"{catalogue_name} search accounting mismatch: received "
            f"{len(results)} results for {len(indexed_genomes)} indexed genomes."
        )

    consolidate_microbiome_hits(
        databases_path,
        output_path,
        organism_name,
        catalogue_name,
        identity_filter,
        coverage_filter,
    )
    return results

def microbiome_offtarget_blast_allproteins (databases_path, output_path, organism_name, cpus=multiprocessing.cpu_count()):

    """
    Runs ncbi blastp against microbiome proteome.
    This function uses the `run_blastp` function from the `programs` module.
    It uses the MICROBIOME_DB database created by the `index_db_blast_microbiome` function from the `databases` module.
    The blast output is saved in the 'offtarget' folder of the organism directory.

    :param databases_path: Directory where MICROBIOME database is stored.
    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param cpus: Number of threads (CPUs) to use in the blast search.
    
    """
    
    #Database files
    species_databases_path = os.path.join(databases_path,  'species_catalogue')
    microbiome_index_path = os.path.join(species_databases_path, 'MICROBIOME_DB')

    #Organism files
    organism_path = os.path.join(output_path, organism_name)
    organism_prot_seq_path = os.path.join(organism_path, f'genome/{organism_name}.faa')

    offtarget_path = os.path.join(organism_path, 'offtarget')
    blast_output_path = os.path.join(offtarget_path, 'microbiome_offtarget_blast.tsv')

    programs.run_diamond_blastp(
        blastdb= microbiome_index_path,
        query= organism_prot_seq_path,
        output=blast_output_path,
        evalue= '1e-5',
        outfmt= '6 std qcovhsp qcovs',
        cpus=cpus
    )

def human_offtarget_parse (output_path, organism_name):

    """
    Parse NCBI BLASTP results against human proteome, stored in the file 'human_offtarget_blast.tsv'.
    Obtains the hit with the highest percentage of identity for each locus_tag.
    Returns a dictionary with locus_tag as key and highest percentage of identity as value, 
    and a DataFrame with all locus_tags from the genome and their respective values.
    The DataFrame is created using the `metadata_table_with_values` function from the `metadata` module.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    
    :return: Dictionary with locus_tag as key and highest percentage of identity value.
    :return: DataFrame with all locus_tags from the genome and their respective values.
    """

    offtarget_path = os.path.join(output_path, organism_name, 'offtarget')
    human_blast_output = os.path.join(offtarget_path, 'human_offtarget_blast.tsv')
    human_results = os.path.join(offtarget_path, 'human_offtarget.tsv')

    if not files.file_check(human_results):
        blast_output_df = files.read_blast_output(human_blast_output)

        highest_pident_values = {}

        for index,row in blast_output_df.iterrows():
            qseqid = row['qseqid']
            pident = row['pident']

            if qseqid not in highest_pident_values or pident > highest_pident_values[qseqid]:
                highest_pident_values[qseqid] = pident

        df_human = metadata.metadata_table_with_values(output_path, organism_name, highest_pident_values, 
                                            'human_offtarget', offtarget_path, 'no_hit')
    else:
        print('Human offtarget analysis already done, output file found')
        print(human_results)
        df_human = pd.read_csv(human_results, sep='\t', header=0)

    return df_human

def microbiome_species_parse(
    databases_path,
    output_path,
    organism_name,
    catalogue_name,
    identity_filter,
    coverage_filter,
    genome_output_path=None,
):
    """
    Parses consolidated microbiome hits at every supported taxonomy rank.

    Counts distinct classified taxa per protein and normalizes them by the total
    classified taxa available at each rank. Unclassified values are excluded from
    ranks above species. Species are represented by their representative genome ID.

    :param databases_path: Path where microbiome catalogues are stored.
    :param output_path: Path of the organism output.
    :param organism_name: Name of the organism.
    :param catalogue_name: Name of a supported MGnify catalogue.
    :param identity_filter: Identity threshold associated with the consolidated hits.
    :param coverage_filter: Coverage threshold associated with the consolidated hits.
    :param genome_output_path: Optional root containing the staged organism genome.
    :return: Tuple containing normalized and count metadata tables for every rank,
        followed by the species normalization denominator table.
    """
    species_path = catalogue_species_path(databases_path, catalogue_name)
    taxonomy_path = os.path.join(
        species_path,
        "representative_taxonomy.parquet",
    )
    offtarget_path, hits_path, _ = _microbiome_consolidated_paths(
        output_path,
        organism_name,
        catalogue_name,
    )
    if not is_microbiome_consolidation_compatible(
        output_path,
        organism_name,
        catalogue_name,
        identity_filter,
        coverage_filter,
    ):
        raise RuntimeError(
            f"No compatible consolidated microbiome hits were found for "
            f"{catalogue_name} with identity={identity_filter} and "
            f"coverage={coverage_filter}."
        )
    if not os.path.isfile(taxonomy_path):
        raise FileNotFoundError(
            f"Representative taxonomy not found: {taxonomy_path}"
        )

    hits_source = f"read_parquet('{_duckdb_path(hits_path)}')"
    taxonomy_source = f"read_parquet('{_duckdb_path(taxonomy_path)}')"
    counts_by_rank = {}
    denominators = {}
    connection = duckdb.connect()
    try:
        invalid_hits = connection.execute(
            f"""
            SELECT count(*)
            FROM {hits_source}
            WHERE pident IS NULL
                OR qcovhsp IS NULL
                OR pident < {float(identity_filter)}
                OR qcovhsp < {float(coverage_filter)}
            """
        ).fetchone()[0]
        if invalid_hits:
            raise RuntimeError(
                f"Consolidated microbiome hits contain {invalid_hits} rows below "
                f"the requested filters: identity={identity_filter}, "
                f"coverage={coverage_filter}."
            )

        for rank in MICROBIOME_TAXONOMY_RANKS:
            if rank == "species":
                denominator = connection.execute(
                    f"SELECT count(*) FROM {taxonomy_source}"
                ).fetchone()[0]
                rows = connection.execute(
                    f"""
                    SELECT gene, count(DISTINCT representative_genome_id)
                    FROM {hits_source}
                    GROUP BY gene
                    """
                ).fetchall()
            else:
                taxonomy_path_expression = _microbiome_rank_path(rank)
                joined_path_expression = _microbiome_rank_path(rank, "t")
                taxonomy_path_filter = _microbiome_classified_path_filter(rank)
                joined_path_filter = _microbiome_classified_path_filter(rank, "t")
                denominator = connection.execute(
                    f"""
                    SELECT count(DISTINCT {taxonomy_path_expression})
                    FROM {taxonomy_source}
                    WHERE {taxonomy_path_filter}
                    """
                ).fetchone()[0]
                rows = connection.execute(
                    f"""
                    SELECT h.gene, count(DISTINCT {joined_path_expression})
                    FROM {hits_source} h
                    JOIN {taxonomy_source} t
                    USING (representative_genome_id)
                    WHERE {joined_path_filter}
                    GROUP BY h.gene
                    """
                ).fetchall()
            if denominator == 0:
                raise RuntimeError(
                    f"No classified {rank} taxa are available in {taxonomy_path}."
                )
            denominators[rank] = denominator
            counts_by_rank[rank] = dict(rows)
    finally:
        connection.close()

    print(f"Parsing consolidated microbiome hits for {catalogue_name}...")
    column_prefix = catalogue_column_prefix(catalogue_name)
    genome_output_path = genome_output_path or output_path
    result_tables = []
    for rank in MICROBIOME_TAXONOMY_RANKS:
        counts = counts_by_rank[rank]
        normalized = {
            gene: count / denominators[rank]
            for gene, count in counts.items()
        }
        property_prefix = (
            column_prefix
            if rank == "species"
            else f"{column_prefix}_{rank}"
        )
        result_tables.append(
            metadata.metadata_table_with_values(
                genome_output_path,
                organism_name,
                normalized,
                f"{property_prefix}_offtarget_norm",
                offtarget_path,
                0,
            )
        )
        result_tables.append(
            metadata.metadata_table_with_values(
                genome_output_path,
                organism_name,
                counts,
                f"{property_prefix}_offtarget_counts",
                offtarget_path,
                0,
            )
        )

    species_denominator = denominators["species"]
    df_total_genomes = metadata.metadata_table_with_values(
        genome_output_path,
        organism_name,
        {},
        f'{column_prefix}_genomes_analyzed',
        offtarget_path,
        species_denominator,
    )
    result_tables.append(df_total_genomes)
    return tuple(result_tables)


def microbiome_protein_clusters_parse (output_path, organism_name, identity_filter, coverage_filter):

    """
    Parse NCBI BLASTP results against microbiome proteome, stored in the file 'microbiome_offtarget_blast.tsv'.
    Filters results based on identity and coverage thresholds, then counts occurrences of each locus_tag.
    Returns a dictionary with normalized counts of each locus_tag that has at least one hit with UHGP90, 
    and a DataFrame with all locus_tags from the genome and their respective values.
    The DataFrame is created using the `metadata_table_with_values` function from the `metadata` module.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param identity_filter: Percentage identity filter value. Keeps results above this value in the pident column.
    :param coverage_filter: Query coverage filter value. Keeps results above this value in the qcovs column.
    
    :return: Dictionary with normalized counts of locus_tags that have at least one hit with UHGP90.
    :return: DataFrame with all locus_tags from the genome and their respective normalized values.
    """

    offtarget_path = os.path.join(output_path, organism_name, 'offtarget')
    microbiome_blast_output = os.path.join(offtarget_path, 'microbiome_offtarget_blast.tsv')
    microbiome_results = os.path.join(offtarget_path, 'gut_microbiome_offtarget.tsv')

    if not files.file_check(microbiome_results):

        blast_output_df = files.read_blast_output(microbiome_blast_output)

        #Filter % identity and coverage
        filtered_df = blast_output_df[(blast_output_df['pident'] > identity_filter) 
                                    & (blast_output_df['qcovs'] > coverage_filter)]

        value_counts = filtered_df['qseqid'].value_counts()
        max_count = value_counts.max()
        norm_counts = value_counts / max_count

        normalized_counts_dict = norm_counts.to_dict()

        df_microbiome = metadata.metadata_table_with_values(output_path, organism_name, normalized_counts_dict, 
                                            'gut_microbiome_offtarget', offtarget_path, 'no_hit')

    else:
        print('Microbiome offtarget analysis already done, output file found')
        print(microbiome_results)
        df_microbiome = pd.read_csv(microbiome_results, sep='\t', header=0)
        
    return df_microbiome


def run_foldseek_human_structures (databases_path, output_path, organism_name, container_engine='docker'):

    """
    Runs Foldseek easy-search against unified human reference structures database.
    Uses ONLY the reference structures for each locus_tag.
    
    Reference structures are obtained via structures.get_all_reference_structures():
    - PDB reference structures: PDB_{uniprot}_{pdb_id}_{chain}.pdb (extracted chains)
    - AlphaFold models: AF_{uniprot}.pdb (full predictions)
    
    Each reference structure is searched against the unified human reference database
    containing both PDB and AlphaFold structures.

    :param output_path: Directory of the organism output.
    :param databases_path: Directory where FOLDSEEK structure databases are stored.
    :param organism_name: Name of the organism.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    :return: Dictionary with locus_tag as key and path to foldseek results file as value.
    """
    foldseek_results_mapping = {}

    print(f'\n{"="*80}')
    print('FOLDSEEK HUMAN OFFTARGET ANALYSIS')
    print(f'{"="*80}\n')

    # Human unified reference database (PDB + AlphaFold structures)
    # DB is created in DB_foldseek subdirectory by programs.run_foldseek_create_index_db
    db_human_path = os.path.join(databases_path, 'human_structures', 'DB_foldseek')

    # Get all reference structures using helper function
    print('Getting reference structures for all locus_tags...')
    reference_dict = structures.get_all_reference_structures(output_path, organism_name, path_mode=True)
    
    # Filter out None values (locus_tags without structures)
    reference_dict = {k: v for k, v in reference_dict.items() if v is not None}
    
    print(f'Found {len(reference_dict)} locus_tags with reference structures')
    
    # Offtarget path - create directory before any early returns
    offtarget_path = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtarget_path, 'foldseek_results')

    if not os.path.exists(foldseek_results_path):
        os.makedirs(foldseek_results_path, exist_ok=True)
    
    if not reference_dict:
        print('ERROR: No reference structures found. Make sure to run structures pipeline first.')
        return {}

    # Run Foldseek for each reference structure
    success_count = 0
    error_count = 0
    skipped_count = 0
    
    for locus_tag, struct_path in reference_dict.items():
        struct_name = os.path.basename(struct_path)
        struct_dir = os.path.dirname(struct_path)
        
        print(f'\n[{locus_tag}] Processing {struct_name}')
        
        # Search against unified human reference database (PDB + AlphaFold)
        try:
            programs.run_foldseek_search(struct_dir, db_human_path, 'DB_human_reference', struct_name, foldseek_results_path, container_engine=container_engine)
            
            # Validate output file exists before reporting success
            struct_basename = struct_name.split('.')[0]
            result_file = os.path.join(foldseek_results_path, f'{struct_basename}_output_foldseek', f'{struct_basename}_vs_DB_human_reference_foldseek_results.tsv')
            
            if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                print(f'  ✓ Foldseek search completed')
                success_count += 1
                foldseek_results_mapping[locus_tag] = [result_file]
            else:
                print(f'  ✗ Foldseek search failed: output file not found or empty')
                logging.error(f'Foldseek output file missing or empty: {result_file}')
                error_count += 1
                foldseek_results_mapping[locus_tag] = []

        except Exception as e:
            print(f'  ✗ Foldseek search failed')
            logging.exception(f'Error running Foldseek search: {e}')
            error_count += 1
            foldseek_results_mapping[locus_tag] = []
    
    print(f'\n{"="*80}')
    print('FOLDSEEK SUMMARY')
    print(f'{"="*80}')
    print(f'Total locus_tags with structures: {len(reference_dict)}')
    print(f'Successfully processed: {success_count}')
    print(f'Errors: {error_count}')
    print(f'{"="*80}\n')

    return foldseek_results_mapping

def run_foldseek_human_colabfold_structures(databases_path, output_path, organism_name, container_engine='docker'):

    """
    Runs Foldseek easy-search against unified human reference structures database
    using ColabFold structures (CB_*.pdb) for each locus_tag.

    Raw Foldseek results are stored in the shared 'foldseek_results' folder so
    existing CB_* queries can be reused if they were already generated by the
    default analysis.

    :param output_path: Directory of the organism output.
    :param databases_path: Directory where FOLDSEEK structure databases are stored.
    :param organism_name: Name of the organism.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    :return: Dictionary with locus_tag as key and path to foldseek results file as value.
    """
    foldseek_results_mapping = {}

    print(f'\n{"="*80}')
    print('FOLDSEEK HUMAN OFFTARGET ANALYSIS (COLABFOLD)')
    print(f'{"="*80}\n')

    db_human_path = os.path.join(databases_path, 'human_structures', 'DB_foldseek')
    organism_structures_path = os.path.join(output_path, organism_name, 'structures')
    all_locus_tags = metadata.ref_gbk_locus(output_path, organism_name)

    offtarget_path = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtarget_path, 'foldseek_results')

    if not os.path.exists(foldseek_results_path):
        os.makedirs(foldseek_results_path, exist_ok=True)

    success_count = 0
    error_count = 0
    skipped_count = 0

    for locus_tag in all_locus_tags:
        locus_dir = os.path.join(organism_structures_path, locus_tag)
        cb_structures = structures.find_colabfold_for_locus(locus_dir)

        if not cb_structures:
            print(f'\n[{locus_tag}] No ColabFold structure found, skipping')
            skipped_count += 1
            foldseek_results_mapping[locus_tag] = []
            continue

        struct_path = cb_structures[0]
        struct_name = os.path.basename(struct_path)
        struct_dir = os.path.dirname(struct_path)

        print(f'\n[{locus_tag}] Processing {struct_name}')

        try:
            programs.run_foldseek_search(
                struct_dir,
                db_human_path,
                'DB_human_reference',
                struct_name,
                foldseek_results_path,
                container_engine=container_engine
            )

            struct_basename = struct_name.split('.')[0]
            result_file = os.path.join(
                foldseek_results_path,
                f'{struct_basename}_output_foldseek',
                f'{struct_basename}_vs_DB_human_reference_foldseek_results.tsv'
            )

            if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                print('  ✓ Foldseek search completed')
                success_count += 1
                foldseek_results_mapping[locus_tag] = [result_file]
            else:
                print('  ✗ Foldseek search failed: output file not found or empty')
                logging.error(f'Foldseek output file missing or empty: {result_file}')
                error_count += 1
                foldseek_results_mapping[locus_tag] = []

        except Exception as e:
            print('  ✗ Foldseek search failed')
            logging.exception(f'Error running Foldseek search with ColabFold structure: {e}')
            error_count += 1
            foldseek_results_mapping[locus_tag] = []

    print(f'\n{"="*80}')
    print('FOLDSEEK COLABFOLD SUMMARY')
    print(f'{"="*80}')
    print(f'Total locus_tags in genome: {len(all_locus_tags)}')
    print(f'Successfully processed: {success_count}')
    print(f'Skipped (no ColabFold structure): {skipped_count}')
    print(f'Errors: {error_count}')
    print(f'{"="*80}\n')

    return foldseek_results_mapping

def foldseek_human_parser (output_path, organism_name, map_foldseek):

    """
    Parses Foldseek results for reference structures and maps them to locus_tags.
    
    For each locus_tag, selects the best match (highest TM-score) from the unified 
    human reference database search.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param map_foldseek: Dictionary mapping locus_tag to foldseek result file (single file per locus_tag).
    
    :return: Dictionary with locus_tag as key and best foldseek match.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtargets_dir, 'foldseek_results')

    foldseek_dict_file = os.path.join(foldseek_results_path, 'human_foldseek_dict.json')

    if not files.file_check(foldseek_dict_file):
        
        print('\nParsing Foldseek results...')
        results_foldseek_dict = {}
        
        # Guard against empty mapping (no structures found)
        if not map_foldseek:
            print('No Foldseek results to parse (empty mapping)')
            with open(foldseek_dict_file, 'w') as f:
                json.dump({}, f)
            return {}
        
        # Parse results for each locus_tag
        for locus_tag, result_files in map_foldseek.items():
            
            # Skip if no result files (e.g., foldseek search failed)
            if not result_files:
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
                print(f'  {locus_tag}: Foldseek search failed, no results')
                continue
            
            # Process the single result file for this locus_tag
            result_file = result_files[0]
            
            if not files.file_check(result_file):
                print(f'  {locus_tag}: Result file not found: {result_file}')
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
                continue
            
            try:
                df = pd.read_csv(result_file, sep='\t', usecols=['query', 'target', 'alnlen', 'qcov', 'tcov', 'lddt', 'qtmscore', 'ttmscore', 'alntmscore', 'rmsd', 'prob', 'pident', 'evalue'])
                
                if not df.empty:
                    # Sort by TM score (max of query and target TM-scores)
                    df['max_tmscore'] = df[['qtmscore', 'ttmscore']].max(axis=1)
                    best_row = df.sort_values(by='max_tmscore', ascending=False).iloc[0]
                    
                    results_foldseek_dict[locus_tag] = {
                        'query_structure': best_row['query'],
                        'target_foldseek': best_row['target'],
                        'alnlen_foldseek': best_row['alnlen'],
                        'qcov_foldseek': best_row['qcov'],
                        'tcov_foldseek': best_row['tcov'],
                        'lddt_foldseek': best_row['lddt'],
                        'qtmscore_foldseek': best_row['qtmscore'],
                        'ttmscore_foldseek': best_row['ttmscore'],
                        'alntmscore_foldseek': best_row['alntmscore'],
                        'rmsd_foldseek': best_row['rmsd'],
                        'prob_foldseek': best_row['prob'],
                        'pident_foldseek': best_row['pident'],
                        'evalue_foldseek': best_row['evalue']
                    }
                    print(f'  {locus_tag}: Best match = {best_row["target"]} (TM-score={best_row["max_tmscore"]:.3f})')
                else:
                    # Empty result file
                    results_foldseek_dict[locus_tag] = {
                        'query_structure': None,
                        'target_foldseek': None,
                        'alnlen_foldseek': None,
                        'qcov_foldseek': None,
                        'tcov_foldseek': None,
                        'lddt_foldseek': None,
                        'qtmscore_foldseek': None,
                        'ttmscore_foldseek': None,
                        'alntmscore_foldseek': None,
                        'rmsd_foldseek': None,
                        'prob_foldseek': None,
                        'pident_foldseek': None,
                        'evalue_foldseek': None
                    }
                    print(f'  {locus_tag}: No hits found in foldseek results')
                    
            except Exception as e:
                logging.exception(f'Could not read {result_file}: {e}')
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
        
        # Save results
        files.dict_to_json(foldseek_results_path, 'human_foldseek_dict.json', results_foldseek_dict)
        print(f'\nFoldseek results saved to {foldseek_dict_file}')
        print(f'Total locus_tags with matches: {len([v for v in results_foldseek_dict.values() if v["target_foldseek"] is not None])}')
    else:
        print(f'Loading existing foldseek results from {foldseek_dict_file}')
        results_foldseek_dict = files.json_to_dict(foldseek_dict_file)
    
    return results_foldseek_dict

def foldseek_human_colabfold_parser(output_path, organism_name, map_foldseek):

    """
    Parses Foldseek results for ColabFold structures and maps them to locus_tags.

    For each locus_tag, selects the best match (highest TM-score) from the unified
    human reference database search using the ColabFold structure as query.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param map_foldseek: Dictionary mapping locus_tag to foldseek result file.

    :return: Dictionary with locus_tag as key and best foldseek match.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtargets_dir, 'foldseek_results')

    foldseek_dict_file = os.path.join(foldseek_results_path, 'human_foldseek_colabfold_dict.json')

    if not files.file_check(foldseek_dict_file):

        print('\nParsing Foldseek ColabFold results...')
        results_foldseek_dict = {}

        if not map_foldseek:
            print('No Foldseek ColabFold results to parse (empty mapping)')
            with open(foldseek_dict_file, 'w') as f:
                json.dump({}, f)
            return {}

        for locus_tag, result_files in map_foldseek.items():

            if not result_files:
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
                print(f'  {locus_tag}: Foldseek ColabFold search failed, no results')
                continue

            result_file = result_files[0]

            if not files.file_check(result_file):
                print(f'  {locus_tag}: Result file not found: {result_file}')
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
                continue

            try:
                df = pd.read_csv(result_file, sep='\t', usecols=['query', 'target', 'alnlen', 'qcov', 'tcov', 'lddt', 'qtmscore', 'ttmscore', 'alntmscore', 'rmsd', 'prob', 'pident', 'evalue'])

                if not df.empty:
                    df['max_tmscore'] = df[['qtmscore', 'ttmscore']].max(axis=1)
                    best_row = df.sort_values(by='max_tmscore', ascending=False).iloc[0]

                    results_foldseek_dict[locus_tag] = {
                        'query_structure': best_row['query'],
                        'target_foldseek': best_row['target'],
                        'alnlen_foldseek': best_row['alnlen'],
                        'qcov_foldseek': best_row['qcov'],
                        'tcov_foldseek': best_row['tcov'],
                        'lddt_foldseek': best_row['lddt'],
                        'qtmscore_foldseek': best_row['qtmscore'],
                        'ttmscore_foldseek': best_row['ttmscore'],
                        'alntmscore_foldseek': best_row['alntmscore'],
                        'rmsd_foldseek': best_row['rmsd'],
                        'prob_foldseek': best_row['prob'],
                        'pident_foldseek': best_row['pident'],
                        'evalue_foldseek': best_row['evalue']
                    }
                    print(f'  {locus_tag}: Best ColabFold match = {best_row["target"]} (TM-score={best_row["max_tmscore"]:.3f})')
                else:
                    results_foldseek_dict[locus_tag] = {
                        'query_structure': None,
                        'target_foldseek': None,
                        'alnlen_foldseek': None,
                        'qcov_foldseek': None,
                        'tcov_foldseek': None,
                        'lddt_foldseek': None,
                        'qtmscore_foldseek': None,
                        'ttmscore_foldseek': None,
                        'alntmscore_foldseek': None,
                        'rmsd_foldseek': None,
                        'prob_foldseek': None,
                        'pident_foldseek': None,
                        'evalue_foldseek': None
                    }
                    print(f'  {locus_tag}: No hits found in Foldseek ColabFold results')

            except Exception as e:
                logging.exception(f'Could not read {result_file}: {e}')
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }

        files.dict_to_json(foldseek_results_path, 'human_foldseek_colabfold_dict.json', results_foldseek_dict)
        print(f'\nFoldseek ColabFold results saved to {foldseek_dict_file}')
        print(f'Total locus_tags with ColabFold matches: {len([v for v in results_foldseek_dict.values() if v["target_foldseek"] is not None])}')
    else:
        print(f'Loading existing Foldseek ColabFold results from {foldseek_dict_file}')
        results_foldseek_dict = files.json_to_dict(foldseek_dict_file)

    return results_foldseek_dict

def merge_foldseek_data (output_path, organism_name):
    """
    Format foldseek results for final output table.
    
    With the new structure organization, foldseek results are already keyed by locus_tag,
    so this function primarily reformats the data for the final table.
    
    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param id_equivalences: Dictionary with locus_tag and uniprot_id (for compatibility).
    :param uniprot_proteome_annotations: Dictionary with annotations (for compatibility).

    :return: Dictionary with the merged data formatted for final table.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtargets_dir, 'foldseek_results')
    foldseek_res_file = os.path.join(foldseek_results_path, 'human_foldseek_dict.json')
    foldseek_mapped_file = os.path.join(offtargets_dir, f'{organism_name}_final_foldseek_results.json')

    if not files.file_check(foldseek_mapped_file):
        if files.file_check(foldseek_res_file):  

            results_foldseek_dict = files.json_to_dict(foldseek_res_file)

            mapped_dict = {}

            # Results are already organized by locus_tag, just reformat
            for locus_tag, foldseek_data in results_foldseek_dict.items():
                mapped_dict[locus_tag] = {
                    'gene': locus_tag,
                    'query_structure': foldseek_data.get('query_structure'),
                    'structure': foldseek_data.get('query_structure'),  # The reference structure used
                    'target': foldseek_data.get('target_foldseek'),
                    'alnlen': foldseek_data.get('alnlen_foldseek'),
                    'qcov': foldseek_data.get('qcov_foldseek'),
                    'tcov': foldseek_data.get('tcov_foldseek'),
                    'lddt': foldseek_data.get('lddt_foldseek'),
                    'qtmscore': foldseek_data.get('qtmscore_foldseek'),
                    'ttmscore': foldseek_data.get('ttmscore_foldseek'),
                    'alntmscore': foldseek_data.get('alntmscore_foldseek'),
                    'rmsd': foldseek_data.get('rmsd_foldseek'),
                    'prob': foldseek_data.get('prob_foldseek'),
                    'pident': foldseek_data.get('pident_foldseek'),
                    'evalue': foldseek_data.get('evalue_foldseek')
                }

            files.dict_to_json(offtargets_dir, f'{organism_name}_final_foldseek_results.json', mapped_dict)
            print(f'\nFoldseek data merged and saved to {foldseek_mapped_file}')
            print(f'Total genes with foldseek results: {len([v for v in mapped_dict.values() if v["target"] is not None])}')
        else:
            print(f'File {foldseek_res_file} not found.')
            mapped_dict = {}
    else:
        mapped_dict = files.json_to_dict(foldseek_mapped_file)
        print(f'Foldseek results in {foldseek_mapped_file}.')

    return mapped_dict

def merge_foldseek_colabfold_data(output_path, organism_name):
    """
    Format Foldseek ColabFold results for final output table.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.

    :return: Dictionary with the merged data formatted for the final table.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtargets_dir, 'foldseek_results')
    foldseek_res_file = os.path.join(foldseek_results_path, 'human_foldseek_colabfold_dict.json')
    foldseek_mapped_file = os.path.join(offtargets_dir, f'{organism_name}_final_foldseek_colabfold_results.json')

    if not files.file_check(foldseek_mapped_file):
        if files.file_check(foldseek_res_file):

            results_foldseek_dict = files.json_to_dict(foldseek_res_file)

            mapped_dict = {}

            for locus_tag, foldseek_data in results_foldseek_dict.items():
                mapped_dict[locus_tag] = {
                    'gene': locus_tag,
                    'query_structure': foldseek_data.get('query_structure'),
                    'structure': foldseek_data.get('query_structure'),
                    'target': foldseek_data.get('target_foldseek'),
                    'alnlen': foldseek_data.get('alnlen_foldseek'),
                    'qcov': foldseek_data.get('qcov_foldseek'),
                    'tcov': foldseek_data.get('tcov_foldseek'),
                    'lddt': foldseek_data.get('lddt_foldseek'),
                    'qtmscore': foldseek_data.get('qtmscore_foldseek'),
                    'ttmscore': foldseek_data.get('ttmscore_foldseek'),
                    'alntmscore': foldseek_data.get('alntmscore_foldseek'),
                    'rmsd': foldseek_data.get('rmsd_foldseek'),
                    'prob': foldseek_data.get('prob_foldseek'),
                    'pident': foldseek_data.get('pident_foldseek'),
                    'evalue': foldseek_data.get('evalue_foldseek')
                }

            files.dict_to_json(offtargets_dir, f'{organism_name}_final_foldseek_colabfold_results.json', mapped_dict)
            print(f'\nFoldseek ColabFold data merged and saved to {foldseek_mapped_file}')
            print(f'Total genes with Foldseek ColabFold results: {len([v for v in mapped_dict.values() if v["target"] is not None])}')
        else:
            print(f'File {foldseek_res_file} not found.')
            mapped_dict = {}
    else:
        mapped_dict = files.json_to_dict(foldseek_mapped_file)
        print(f'Foldseek ColabFold results in {foldseek_mapped_file}.')

    return mapped_dict

def final_foldseek_structure_table (output_path, organism_name, mapped_dict):
    """
    Create a final table with the merge dictionary of the Foldseek results. 
    Saves the table in a .tsv file named using the organism name followed by '_final_foldseek_results.tsv'in the 'structures' directory.
    Returns a DataFrame with the final results.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param mapped_dict: Dictionary with the merged data.

    :return: DataFrame with the final results.
    """
    
    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')

    all_locus_tags = metadata.ref_gbk_locus(output_path, organism_name)
    
    # Create a DataFrame with final results
    final_foldseek_file = os.path.join(offtargets_dir, f'{organism_name}_final_foldseek_results.tsv')

    if not files.file_check(final_foldseek_file):

        rows = []

        for locus_tag in all_locus_tags:
            if locus_tag in mapped_dict:
                rows.append(mapped_dict[locus_tag])
            else:
                rows.append({'gene': locus_tag, 'query_structure': None, 'structure': 'No hit', 'target': None, 'alnlen': None, 'qcov': None, 'tcov': None, 'lddt': None, 'qtmscore': None, 'ttmscore': None, 'alntmscore': None, 'rmsd': None, 'prob': None, 'pident': None, 'evalue': None })

        final_foldseek_df = pd.DataFrame(rows).rename(columns={
            'query_structure': 'FS_query_structure',
            'structure': 'FS_organism_structure_query',
            'target': 'FS_human_structure_hit',
            'alnlen': 'FS_alnlen',
            'qcov': 'FS_qcov',
            'tcov': 'FS_tcov',
            'lddt': 'FS_lddt',
            'qtmscore': 'FS_qtmscore',
            'ttmscore': 'FS_ttmscore',
            'alntmscore': 'FS_alntmscore',
            'rmsd': 'FS_rmsd',
            'prob': 'FS_prob',
            'pident': 'FS_pident',
            'evalue': 'FS_evalue'
        })

        final_foldseek_df.to_csv(final_foldseek_file, sep='\t', index=False)
   
        print(f'Foldseek final results saved to {final_foldseek_file}.')
    
    else:
        print(f'Foldseek final results in {final_foldseek_file}.')
        final_foldseek_df = pd.read_csv(final_foldseek_file, sep='\t')

    return final_foldseek_df

def final_foldseek_colabfold_structure_table(output_path, organism_name, mapped_dict):
    """
    Create a final table with the merged dictionary of the Foldseek ColabFold results.
    Saves the table in a .tsv file named using the organism name followed by
    '_final_foldseek_colabfold_results.tsv'.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param mapped_dict: Dictionary with the merged data.

    :return: DataFrame with the final results.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    all_locus_tags = metadata.ref_gbk_locus(output_path, organism_name)
    final_foldseek_file = os.path.join(offtargets_dir, f'{organism_name}_final_foldseek_colabfold_results.tsv')

    if not files.file_check(final_foldseek_file):

        rows = []

        for locus_tag in all_locus_tags:
            if locus_tag in mapped_dict:
                rows.append(mapped_dict[locus_tag])
            else:
                rows.append({'gene': locus_tag, 'query_structure': None, 'structure': 'No hit', 'target': None, 'alnlen': None, 'qcov': None, 'tcov': None, 'lddt': None, 'qtmscore': None, 'ttmscore': None, 'alntmscore': None, 'rmsd': None, 'prob': None, 'pident': None, 'evalue': None})

        final_foldseek_df = pd.DataFrame(rows).rename(columns={
            'query_structure': 'FS_CB_query_structure',
            'structure': 'FS_CB_organism_structure_query',
            'target': 'FS_CB_human_structure_hit',
            'alnlen': 'FS_CB_alnlen',
            'qcov': 'FS_CB_qcov',
            'tcov': 'FS_CB_tcov',
            'lddt': 'FS_CB_lddt',
            'qtmscore': 'FS_CB_qtmscore',
            'ttmscore': 'FS_CB_ttmscore',
            'alntmscore': 'FS_CB_alntmscore',
            'rmsd': 'FS_CB_rmsd',
            'prob': 'FS_CB_prob',
            'pident': 'FS_CB_pident',
            'evalue': 'FS_CB_evalue'
        })

        final_foldseek_df.to_csv(final_foldseek_file, sep='\t', index=False)

        print(f'Foldseek ColabFold final results saved to {final_foldseek_file}.')

    else:
        print(f'Foldseek ColabFold final results in {final_foldseek_file}.')
        final_foldseek_df = pd.read_csv(final_foldseek_file, sep='\t')

    return final_foldseek_df
