#!/usr/bin/env python3
"""
FastTarget Pipeline - Test Validation Script

This script validates the output of the test run to ensure
all modules executed correctly and produced expected results.

Usage:
    python validate_test.py --test_dir organism/test
"""

import os
import sys
import contextlib
import argparse
import pandas as pd
import glob
from datetime import datetime

# Expected genes in the test dataset
EXPECTED_GENES = [
    'MPN_RS00140', 'MPN_RS00300', 'MPN_RS00445', 'MPN_RS00890', 'MPN_RS01250',
    'MPN_RS01270', 'MPN_RS01385', 'MPN_RS01390', 'MPN_RS01460', 'MPN_RS01665',
    'MPN_RS01670', 'MPN_RS02025', 'MPN_RS02430', 'MPN_RS03100', 'MPN_RS02775',
    'MPN_RS03450', 'MPN_RS03550', 'MPN_RS03565', 'MPN_RS03570', 'MPN_RS04200'
]

import os
import json
import glob
import pandas as pd

def _is_nonempty_file(path: str) -> bool:
    """Exists and size > 0 bytes."""
    return os.path.isfile(path) and os.path.getsize(path) > 0

def _count_data_rows_tsv(path: str, sep: str = "\t") -> int:
    """
    Count data rows in a TSV/CSV (excluding header).
    Returns 0 if unreadable/empty.
    """
    if not _is_nonempty_file(path):
        return 0
    try:
        df = pd.read_csv(path, sep=sep)
        return int(df.shape[0])
    except Exception:
        return 0

def _read_df_any(path_tsv: str, path_csv: str):
    """Try read TSV first, then CSV. Returns (df, used_path, used_sep) or (None, None, None)."""
    if _is_nonempty_file(path_tsv):
        try:
            return pd.read_csv(path_tsv, sep="\t"), path_tsv, "\t"
        except Exception:
            pass
    if _is_nonempty_file(path_csv):
        try:
            return pd.read_csv(path_csv), path_csv, ","
        except Exception:
            pass
    return None, None, None

def reference_structure_breakdown(organism_dir: str):
    """
    Classify each protein by the reference structure Foldseek searched with.

    Foldseek runs on one reference structure per gene, so the number of results to
    expect is the number of genes that have one, not the number of genes in the
    genome. Enabling colabfold is what models the remainder and brings the two
    numbers together.

    The selection is not reimplemented here: it depends on pocket_full_mode, on
    which PDB chains were extracted and on a priority order, and a copy of it in
    the validator would drift from the pipeline and report expectations the run
    never had. The pipeline's own get_all_reference_structures is used instead, so
    the expectation is by construction the set the run searched.

    :param organism_dir: Path of the organism output directory.
    :return: Dictionary of source name to the sorted locus tags using it, with the
        genes that have no structure under 'none', or None when the pipeline's
        function cannot be used.
    """
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from ftscripts import structures as ft_structures
    except Exception:
        return None

    output_path = os.path.dirname(os.path.abspath(organism_dir))
    organism_name = os.path.basename(os.path.abspath(organism_dir))

    try:
        # The function reports its search per gene, which is useful in a pipeline
        # run and noise in a validation summary.
        with open(os.devnull, 'w') as quiet, contextlib.redirect_stdout(quiet):
            reference_structures = ft_structures.get_all_reference_structures(
                output_path, organism_name, path_mode=True
            )
    except Exception:
        return None

    breakdown = {'pdb': [], 'alphafold': [], 'colabfold': [], 'none': []}
    for locus_tag, path in sorted((reference_structures or {}).items()):
        if not path:
            breakdown['none'].append(locus_tag)
            continue
        name = os.path.basename(path if isinstance(path, str) else str(path))
        if name.startswith('AF_'):
            breakdown['alphafold'].append(locus_tag)
        elif name.startswith('CB_'):
            breakdown['colabfold'].append(locus_tag)
        else:
            breakdown['pdb'].append(locus_tag)

    return breakdown


def check_offtarget(offtarget_dir: str, validation: dict, expected_count: int = 20, min_microbiome_genomes: int = 4744,
                    organism_dir: str = None):
    """
    Checks 3 independent parts:
      1) Foldseek: one result per gene that has a reference structure, and a
         human_foldseek_dict.json with one key per searched gene
      2) Human offtarget files + table shape (2 cols, 20 rows) + blast nonempty
      3) Human gut off-target: the consolidated hits and their manifest
    Only the part that fails appends errors to validation["failed"].
    """

    # Ensure keys exist
    validation.setdefault("passed", [])
    validation.setdefault("warnings", [])
    validation.setdefault("failed", [])

    # --- Base dir ---
    if not os.path.isdir(offtarget_dir):
        validation["failed"].append("❌ Offtarget directory not found")
        return  # nothing else to do

    # ============================================================
    # 1) FOLDSEEK
    # ============================================================
    foldseek_dir = os.path.join(offtarget_dir, "foldseek_results")
    foldseek_errors = []
    foldseek_warnings = []
    foldseek_passed = []

    # How many Foldseek results to expect, broken down by where the structure came
    # from. A flat 'expect 20' failed every run in which one gene had no structure,
    # which is a legitimate outcome unless ColabFold is enabled to model it.
    breakdown = reference_structure_breakdown(organism_dir) if organism_dir else None
    if breakdown is not None and sum(len(v) for v in breakdown.values()):
        searchable = breakdown['pdb'] + breakdown['alphafold'] + breakdown['colabfold']
        expected_foldseek = len(searchable)
        foldseek_passed.append(
            f"✅ Foldseek: reference structures by source - "
            f"PDB {len(breakdown['pdb'])}, AlphaFold {len(breakdown['alphafold'])}, "
            f"ColabFold {len(breakdown['colabfold'])}, none {len(breakdown['none'])}"
        )
        if breakdown['none']:
            foldseek_warnings.append(
                f"⚠️  Foldseek: {len(breakdown['none'])} of {expected_count} genes have no reference "
                f"structure and are not searched ({', '.join(breakdown['none'][:5])}). "
                f"Enable colabfold to model them and reach {expected_count}."
            )
    else:
        # Without the structures directory the old genome-wide expectation is all
        # there is to go on.
        expected_foldseek = expected_count

    if not os.path.isdir(foldseek_dir):
        foldseek_errors.append("❌ Foldseek: foldseek_results directory not found")
    else:
        # Expect subfolders (one per searched protein) and inside each a *_foldseek_results.tsv
        subdirs = [d for d in os.listdir(foldseek_dir) if os.path.isdir(os.path.join(foldseek_dir, d))]
        if len(subdirs) != expected_foldseek:
            foldseek_warnings.append(
                f"⚠️  Foldseek: Expected {expected_foldseek} result folders, found {len(subdirs)}"
            )

        # Find all foldseek result TSVs recursively
        foldseek_tsvs = glob.glob(os.path.join(foldseek_dir, "**", "*_foldseek_results.tsv"), recursive=True)
        nonempty_foldseek_tsvs = [p for p in foldseek_tsvs if _is_nonempty_file(p)]

        if len(nonempty_foldseek_tsvs) < expected_foldseek:
            foldseek_errors.append(
                f"❌ Foldseek: Expected at least {expected_foldseek} non-empty *_foldseek_results.tsv "
                f"(one per gene with a reference structure), found {len(nonempty_foldseek_tsvs)}"
            )
        else:
            foldseek_passed.append(
                f"✅ Foldseek: Found {len(nonempty_foldseek_tsvs)} non-empty result TSVs "
                f"for {expected_foldseek} searched genes"
            )

        # human_foldseek_dict.json: one key per searched gene
        dict_path = os.path.join(foldseek_dir, "human_foldseek_dict.json")
        if not _is_nonempty_file(dict_path):
            foldseek_errors.append("❌ Foldseek: human_foldseek_dict.json not found or empty")
        else:
            try:
                with open(dict_path, "r", encoding="utf-8") as fh:
                    d = json.load(fh)
                if not isinstance(d, dict):
                    foldseek_errors.append("❌ Foldseek: human_foldseek_dict.json is not a JSON dict")
                else:
                    nkeys = len(d)
                    # A gene whose search produced no hit is stored with None fields,
                    # which is a result. An empty dict is not.
                    with_hit = sum(
                        1 for v in d.values()
                        if isinstance(v, dict) and v.get('target_foldseek') is not None
                    )
                    if nkeys != expected_foldseek:
                        foldseek_errors.append(
                            f"❌ Foldseek: human_foldseek_dict.json should have {expected_foldseek} keys "
                            f"(one per searched gene), found {nkeys}"
                        )
                    elif with_hit == 0 and expected_foldseek:
                        foldseek_errors.append(
                            f"❌ Foldseek: human_foldseek_dict.json has {nkeys} keys but no gene has a hit"
                        )
                    else:
                        foldseek_passed.append(
                            f"✅ Foldseek: human_foldseek_dict.json has {nkeys} keys, "
                            f"{with_hit} with a human structural hit"
                        )
            except Exception as e:
                foldseek_errors.append(f"❌ Foldseek: Failed to parse human_foldseek_dict.json ({e})")

        # ColabFold Foldseek outputs are a separate run with its own outputs; check
        # them only when they exist, since colabfold is optional.
        colabfold_dict = os.path.join(foldseek_dir, "human_foldseek_colabfold_dict.json")
        if os.path.exists(colabfold_dict):
            try:
                with open(colabfold_dict, "r", encoding="utf-8") as fh:
                    cb = json.load(fh)
                if isinstance(cb, dict) and cb:
                    foldseek_passed.append(
                        f"✅ Foldseek (ColabFold): human_foldseek_colabfold_dict.json has {len(cb)} keys"
                    )
                else:
                    foldseek_errors.append(
                        "❌ Foldseek (ColabFold): human_foldseek_colabfold_dict.json is empty"
                    )
            except Exception as e:
                foldseek_errors.append(
                    f"❌ Foldseek (ColabFold): Failed to parse human_foldseek_colabfold_dict.json ({e})"
                )

    # Commit foldseek results
    validation["passed"].extend(foldseek_passed)
    validation["warnings"].extend(foldseek_warnings)
    validation["failed"].extend(foldseek_errors)

    # ============================================================
    # 2) HUMAN OFFTARGET
    # ============================================================
    human_errors = []
    human_warnings = []
    human_passed = []

    human_blast = os.path.join(offtarget_dir, "human_offtarget_blast.tsv")
    if not _is_nonempty_file(human_blast):
        human_errors.append("❌ Human offtarget: human_offtarget_blast.tsv not found or empty")
    else:
        nrows = _count_data_rows_tsv(human_blast, sep="\t")
        if nrows < 1:
            human_errors.append("❌ Human offtarget: human_offtarget_blast.tsv has no data rows")
        else:
            human_passed.append(f"✅ Human offtarget: BLAST output has {nrows} rows")

    human_tsv = os.path.join(offtarget_dir, "human_offtarget.tsv")
    human_csv = os.path.join(offtarget_dir, "human_offtarget.csv")

    df_human, used_path, _ = _read_df_any(human_tsv, human_csv)
    if df_human is None:
        # vos dijiste que deben existir 3 archivos: blast.tsv, csv, tsv.
        # acá marcamos error si faltan los resúmenes (tsv/csv)
        missing = []
        if not os.path.exists(human_tsv):
            missing.append("human_offtarget.tsv")
        if not os.path.exists(human_csv):
            missing.append("human_offtarget.csv")
        human_errors.append(f"❌ Human offtarget: Missing/unreadable summary file(s): {', '.join(missing) if missing else 'human_offtarget.[tsv/csv]'}")
    else:
        # checks: 2 columns, 20 rows
        if df_human.shape[0] != expected_count:
            human_errors.append(f"❌ Human offtarget: {os.path.basename(used_path)} should have {expected_count} rows, found {df_human.shape[0]}")
        if df_human.shape[1] != 2:
            human_errors.append(f"❌ Human offtarget: {os.path.basename(used_path)} should have 2 columns, found {df_human.shape[1]}")
        if df_human.shape[0] == expected_count and df_human.shape[1] == 2:
            human_passed.append(f"✅ Human offtarget: {os.path.basename(used_path)} has 2 columns and {expected_count} rows")

        # (opcional pero útil) si esperás columna human_offtarget:
        if "human_offtarget" in df_human.columns:
            non_null = int(df_human["human_offtarget"].notna().sum())
            if non_null == expected_count:
                human_passed.append("✅ human_offtarget: All proteins have values")
            elif non_null > 0:
                human_warnings.append(f"⚠️  human_offtarget: Only {non_null}/{expected_count} have values")
            else:
                human_errors.append("❌ human_offtarget: Column has no data")
        else:
            human_warnings.append("⚠️  Human offtarget: 'human_offtarget' column not present (skipping value completeness check)")

    # Also enforce presence of both files (even if we read one):
    if not os.path.exists(human_tsv):
        human_errors.append("❌ Human offtarget: human_offtarget.tsv not found")
    elif not _is_nonempty_file(human_tsv):
        human_errors.append("❌ Human offtarget: human_offtarget.tsv is empty")

    if not os.path.exists(human_csv):
        human_errors.append("❌ Human offtarget: human_offtarget.csv not found")
    elif not _is_nonempty_file(human_csv):
        human_errors.append("❌ Human offtarget: human_offtarget.csv is empty")

    # Commit human results
    validation["passed"].extend(human_passed)
    validation["warnings"].extend(human_warnings)
    validation["failed"].extend(human_errors)

    # ============================================================
    # 3) MICROBIOME OFFTARGET
    # ============================================================
    micro_errors = []
    micro_warnings = []
    micro_passed = []

    species_dir = os.path.join(
        offtarget_dir,
        "microbiomes",
        "human-gut",
        "species_blast_results",
    )
    if not os.path.isdir(species_dir):
        micro_errors.append("❌ Microbiome offtarget: species_blast_results directory not found")
    else:
        # Per-genome DIAMOND results are consolidated into a Parquet file and then
        # deleted, which the manifest records as raw_tsv_cleanup. Counting the TSVs
        # therefore reported zero on every completed run; the manifest is what says
        # how many genomes were actually searched.
        genome_tsvs = glob.glob(os.path.join(species_dir, "MGYG*_offtarget_id*_cov*.tsv"))
        manifest_paths = glob.glob(os.path.join(species_dir, "*_offtarget_manifest.json"))

        if manifest_paths:
            try:
                with open(manifest_paths[0], "r", encoding="utf-8") as fh:
                    micro_manifest = json.load(fh)
                indexed = micro_manifest.get("indexed_representatives", 0)
                with_hits = micro_manifest.get("representatives_with_hits", 0)
                hit_rows = micro_manifest.get("hit_rows", 0)
                genes_with_hits = micro_manifest.get("genes_with_hits", 0)
                cleanup = micro_manifest.get("raw_tsv_cleanup", "not recorded")

                if indexed < min_microbiome_genomes:
                    micro_errors.append(
                        f"❌ Microbiome offtarget: manifest reports {indexed} representatives searched, "
                        f"expected ≥{min_microbiome_genomes}"
                    )
                else:
                    micro_passed.append(
                        f"✅ Microbiome offtarget: {indexed} representatives searched, "
                        f"{with_hits} with hits, {hit_rows} hit rows, {genes_with_hits} genes with hits"
                    )

                if cleanup != "completed" and len(genome_tsvs) < min_microbiome_genomes:
                    micro_warnings.append(
                        f"⚠️  Microbiome offtarget: raw_tsv_cleanup is '{cleanup}' and only "
                        f"{len(genome_tsvs)} per-genome TSVs remain"
                    )
            except Exception as e:
                micro_errors.append(
                    f"❌ Microbiome offtarget: Failed to parse {os.path.basename(manifest_paths[0])} ({e})"
                )
        elif len(genome_tsvs) < min_microbiome_genomes:
            micro_errors.append(
                f"❌ Microbiome offtarget: no consolidation manifest, and only {len(genome_tsvs)} "
                f"per-genome result files found (expected ≥{min_microbiome_genomes})"
            )
        else:
            micro_passed.append(f"✅ Microbiome offtarget: Found {len(genome_tsvs)} completed MGYG genome TSVs")

        # Human gut summaries: must exist, not empty, 2 columns, 20 rows
        gut_files = sorted(glob.glob(os.path.join(species_dir, "human_gut_offtarget*")))
        # Filter to tsv/csv only (your tree shows both)
        gut_files = [p for p in gut_files if p.endswith(".tsv") or p.endswith(".csv")]

        if not gut_files:
            micro_errors.append("❌ Microbiome offtarget: No human_gut_offtarget* summary files found")
        else:
            bad = []
            ok = 0
            for p in gut_files:
                if not _is_nonempty_file(p):
                    bad.append(f"{os.path.basename(p)} (empty)")
                    continue
                try:
                    if p.endswith(".tsv"):
                        df_gut = pd.read_csv(p, sep="\t")
                    else:
                        df_gut = pd.read_csv(p)
                    if df_gut.shape[0] != expected_count or df_gut.shape[1] != 2:
                        bad.append(f"{os.path.basename(p)} ({df_gut.shape[1]} cols, {df_gut.shape[0]} rows)")
                    else:
                        ok += 1
                except Exception as e:
                    bad.append(f"{os.path.basename(p)} (unreadable: {e})")

            if bad:
                micro_errors.append(
                    "❌ Microbiome offtarget: Some human_gut_offtarget* files are invalid (need 2 cols, 20 rows): "
                    + "; ".join(bad[:8]) + (" ..." if len(bad) > 8 else "")
                )
            if ok > 0 and not bad:
                micro_passed.append(f"✅ Microbiome offtarget: All {ok} human_gut_offtarget* summary files have 2 cols and {expected_count} rows")
            elif ok > 0 and bad:
                micro_warnings.append(f"⚠️  Microbiome offtarget: {ok}/{len(gut_files)} human_gut_offtarget* files look OK")

        analyzed_csv = os.path.join(species_dir, "human_gut_genomes_analyzed.csv")
        analyzed_tsv = os.path.join(species_dir, "human_gut_genomes_analyzed.tsv")
        for p in (analyzed_csv, analyzed_tsv):
            if not os.path.exists(p):
                micro_warnings.append(f"⚠️  Microbiome offtarget: {os.path.basename(p)} not found")
            elif not _is_nonempty_file(p):
                micro_warnings.append(f"⚠️  Microbiome offtarget: {os.path.basename(p)} is empty")
            else:
                micro_passed.append(f"✅ Microbiome offtarget: {os.path.basename(p)} found")

    # Commit microbiome results
    validation["passed"].extend(micro_passed)
    validation["warnings"].extend(micro_warnings)
    validation["failed"].extend(micro_errors)

def validate_test_results(results_dir, check_all_modules=True, organism_dir=None):
    """
    Validate that expected output files were created and data integrity is maintained.
    
    :param results_dir: Path to results directory (test_results_TIMESTAMP)
    :param check_all_modules: Whether to check all module outputs (False for stub mode)
    :return: Dictionary with validation results
    """
    print("\n" + "="*80)
    print("VALIDATING TEST RESULTS")
    print("="*80)
    print(f"📁 Results directory: {results_dir}")
    print(f"📅 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    validation = {
        "passed": [],
        "failed": [],
        "warnings": []
    }
    
    if not os.path.exists(results_dir):
        validation["failed"].append(f"❌ Results directory not found: {results_dir}")
        return validation
    
    # Check main results table - CRITICAL
    main_results = os.path.join(results_dir, "test_results_table.tsv")
    if not os.path.exists(main_results):
        validation["failed"].append("❌ Main results table not found")
        return validation
    
    validation["passed"].append("✅ Main results table found")
    
    df = None
    expected_count = len(EXPECTED_GENES)
    try:
        df = pd.read_csv(main_results, sep='\t')
        
        # CRITICAL: Validate gene count and IDs
        print("\n" + "-"*80)
        print("DATA INTEGRITY CHECKS")
        print("-"*80)
        
        actual_count = len(df)
        
        if actual_count != expected_count:
            validation["failed"].append(
                f"❌ Gene count mismatch: Expected {expected_count}, got {actual_count}"
            )
        else:
            validation["passed"].append(
                f"✅ Gene count correct: {actual_count} proteins"
            )
        
        # Check gene IDs
        if 'gene' in df.columns:
            actual_genes = set(df['gene'].tolist())
            expected_genes_set = set(EXPECTED_GENES)
            
            missing_genes = expected_genes_set - actual_genes
            extra_genes = actual_genes - expected_genes_set
            
            if missing_genes:
                validation["failed"].append(
                    f"❌ Missing genes: {sorted(list(missing_genes))}"
                )
            if extra_genes:
                validation["warnings"].append(
                    f"⚠️  Unexpected genes: {sorted(list(extra_genes))}"
                )
            if not missing_genes and not extra_genes:
                validation["passed"].append(
                    f"✅ All expected gene IDs present"
                )
        else:
            validation["failed"].append("❌ 'gene' column missing from results table")
        
        # Check column count
        col_count = len(df.columns)
        validation["passed"].append(f"✅ Results table has {col_count} columns")
        
        print(f"\n📊 Results table shape: {df.shape}")
        print(f"📊 Columns: {col_count}")

        print("\n" + "-"*80)
        print("RESULTS TABLE MODULE CHECKS")
        print("-"*80)

        # Metabolism columns
        ptools_bc = 'PTOOLS_betweenness_centrality'
        ptools_edges = 'PTOOLS_edges'
        if ptools_bc in df.columns:
            numeric_vals = df[ptools_bc].dropna()
            if len(numeric_vals) > 0 and (numeric_vals >= 0).all() and (numeric_vals <= 1).all():
                validation["passed"].append(f"✅ {ptools_bc}: Valid range 0-1")
            elif len(numeric_vals) == 0:
                validation["warnings"].append(f"⚠️  {ptools_bc}: All values are empty")
            else:
                validation["failed"].append(f"❌ {ptools_bc}: Values outside 0-1 range")

        if ptools_edges in df.columns:
            edge_vals = df[ptools_edges].dropna()
            if len(edge_vals) > 0 and (edge_vals == edge_vals.astype(int)).all():
                validation["passed"].append(f"✅ {ptools_edges}: Valid integers")
            elif len(edge_vals) == 0:
                validation["warnings"].append(f"⚠️  {ptools_edges}: All values are empty")
            else:
                validation["failed"].append(f"❌ {ptools_edges}: Contains non-integer values")

        mgt_bc = 'MGT_betweenness_centrality'
        mgt_edges = 'MGT_edges'
        if mgt_bc in df.columns:
            numeric_vals = df[mgt_bc].dropna()
            if len(numeric_vals) > 0 and (numeric_vals >= 0).all() and (numeric_vals <= 1).all():
                validation["passed"].append(f"✅ {mgt_bc}: Valid range 0-1")
            elif len(numeric_vals) == 0:
                validation["warnings"].append(f"⚠️  {mgt_bc}: All values are empty")
            else:
                validation["failed"].append(f"❌ {mgt_bc}: Values outside 0-1 range")

        if mgt_edges in df.columns:
            edge_vals = df[mgt_edges].dropna()
            if len(edge_vals) > 0 and (edge_vals == edge_vals.astype(int)).all():
                validation["passed"].append(f"✅ {mgt_edges}: Valid integers")
            elif len(edge_vals) == 0:
                validation["warnings"].append(f"⚠️  {mgt_edges}: All values are empty")
            else:
                validation["failed"].append(f"❌ {mgt_edges}: Contains non-integer values")

        # Structures columns
        structure_cols = [
            'uniprot', 'structure', 'druggability_score',
            'fpocket_pocket', 'p2rank_probability', 'p2rank_pocket'
        ]
        for col in structure_cols:
            if col not in df.columns:
                validation["failed"].append(f"❌ Structure column missing: {col}")

        if 'druggability_score' in df.columns:
            drug_scores = df['druggability_score'].dropna()
            if len(drug_scores) > 0 and (drug_scores >= 0).all() and (drug_scores <= 1).all():
                validation["passed"].append(
                    f"✅ druggability_score: Valid range 0-1 ({len(drug_scores)} values)"
                )
            elif len(drug_scores) == 0:
                validation["warnings"].append("⚠️  druggability_score: No values found")
            else:
                validation["failed"].append("❌ druggability_score: Values outside 0-1 range")

        # Conservation columns
        if 'core_corecruncher' in df.columns:
            true_count = (df['core_corecruncher'] == True).sum()
            if true_count > 0:
                validation["passed"].append(f"✅ core_corecruncher: {true_count} core genes")
            else:
                validation["warnings"].append("⚠️  core_corecruncher: No core genes found")
        if 'core_roary' in df.columns:
            true_count = (df['core_roary'] == True).sum()
            if true_count > 0:
                validation["passed"].append(f"✅ core_roary: {true_count} core genes")
            else:
                validation["warnings"].append("⚠️  core_roary: No core genes found")

        # Offtarget column
        if 'human_offtarget' in df.columns:
            non_null = df['human_offtarget'].notna().sum()
            if non_null == expected_count:
                validation["passed"].append("✅ human_offtarget: All proteins have values")
            elif non_null > 0:
                validation["warnings"].append(
                    f"⚠️  human_offtarget: Only {non_null}/{expected_count} have values"
                )
            else:
                validation["failed"].append("❌ human_offtarget: Column has no data")
        else:
            validation["failed"].append("❌ human_offtarget: Column missing")

        # Essentiality column
        if 'hit_in_deg' in df.columns:
            true_count = (df['hit_in_deg'] == True).sum()
            total_not_null = df['hit_in_deg'].notna().sum()
            if true_count > 0:
                validation["passed"].append(f"✅ hit_in_deg: {true_count} essential genes found")
            elif total_not_null > 0:
                validation["warnings"].append(
                    f"⚠️  hit_in_deg: No essential genes found ({total_not_null} values)"
                )
            else:
                validation["failed"].append("❌ hit_in_deg: Column has no data")
        else:
            validation["failed"].append("❌ hit_in_deg: Column missing")

        # Localization column
        psortb_col_found = False
        for col in df.columns:
            if 'psort' in col.lower():
                psortb_col_found = True
                non_empty = df[col].notna().sum()
                if non_empty == expected_count:
                    validation["passed"].append(f"✅ {col}: All {non_empty} proteins have values")
                elif non_empty > 0:
                    validation["warnings"].append(
                        f"⚠️  {col}: Only {non_empty}/{expected_count} proteins have values"
                    )
                else:
                    validation["failed"].append(f"❌ {col}: Column has no data")
                break

        if not psortb_col_found:
            validation["warnings"].append("⚠️  PSORTb column not found in results table")
        
    except Exception as e:
        validation["failed"].append(f"❌ Error reading results table: {e}")
        return validation
    
    # Check tables_for_TP directory
    print("\n" + "-"*80)
    print("TARGET PATHOGEN TABLES")
    print("-"*80)
    
    tp_dir = os.path.join(results_dir, "tables_for_TP")
    if not os.path.exists(tp_dir):
        validation["failed"].append("❌ tables_for_TP directory not found")
    else:
        tsv_files = [f for f in os.listdir(tp_dir) if f.endswith('.tsv')]
        if len(tsv_files) < 22:
            validation["warnings"].append(
                f"⚠️  Only {len(tsv_files)} tables in tables_for_TP (expected more)"
            )
        else:
            validation["passed"].append(
                f"✅ tables_for_TP has {len(tsv_files)} tables"
            )
        print(f"📊 TP tables found: {len(tsv_files)}")
    
    return validation


def validate_module_outputs(organism_dir):
    """
    Validate module outputs once using filesystem artifacts only.
    """
    print("\n" + "-"*80)
    print("MODULE OUTPUT VALIDATION")
    print("-"*80)

    validation = {
        "passed": [],
        "failed": [],
        "warnings": []
    }


    expected_count = len(EXPECTED_GENES)

    # Check genome files
    genome_dir = os.path.join(organism_dir, "genome")
    if os.path.exists(genome_dir):
        required_genome_files = ['test.faa', 'test.gff', 'test.fna']
        found = sum(1 for f in required_genome_files if os.path.exists(os.path.join(genome_dir, f)))
        if found == len(required_genome_files):
            validation["passed"].append(f"✅ Genome files complete ({found}/{len(required_genome_files)})")
        else:
            validation["failed"].append(f"❌ Genome files incomplete ({found}/{len(required_genome_files)})")
    else:
        validation["failed"].append("❌ Genome directory not found")

    # Check metabolism
    metabolism_dir = os.path.join(organism_dir, "metabolism")
    if os.path.exists(metabolism_dir):
        metabolism_files = [f for f in os.listdir(metabolism_dir) if f.endswith('.tsv')]
        if len(metabolism_files) >= 5:
            validation["passed"].append(f"✅ Metabolism analysis complete ({len(metabolism_files)} files)")
        else:
            validation["warnings"].append(f"⚠️  Metabolism incomplete ({len(metabolism_files)} files)")

    else:
        validation["failed"].append("❌ Metabolism directory not found")

    # Check structures
    structures_dir = os.path.join(organism_dir, "structures")
    # No fixed 'one gene lacks a structure' assumption: which genes have one is
    # read from the run itself, below.

    if os.path.exists(structures_dir):
        final_structure = os.path.join(structures_dir, "test_final_structure_summary.tsv")
        if os.path.exists(final_structure):
            validation["passed"].append("✅ Structure analysis complete")
        else:
            validation["warnings"].append("⚠️  Structure analysis incomplete")

        # Every gene gets a directory, including one with no structure: the gene
        # that lacks a structure still has its summary and pockets folder. The
        # expectation used to be len(EXPECTED_GENES) - 1, which confused 'no
        # structure' with 'no directory' and reported a warning on a correct run.
        # A missing directory means the per-gene setup did not run, so it is an error.
        protein_dirs = [d for d in os.listdir(structures_dir)
                       if os.path.isdir(os.path.join(structures_dir, d)) and d.startswith('MPN_')]
        if len(protein_dirs) == expected_count:
            validation["passed"].append(f"✅ Protein structure directories: {len(protein_dirs)} found")
        else:
            missing_dirs = sorted(set(EXPECTED_GENES) - set(protein_dirs))
            validation["failed"].append(
                f"❌ Protein directories: {len(protein_dirs)}/{expected_count} found"
                + (f", missing {missing_dirs[:5]}" if missing_dirs else "")
            )

        missing_pdb = []
        missing_fpocket = []
        missing_p2rank = []
        for gene in EXPECTED_GENES:
            gene_dir = os.path.join(structures_dir, gene)
            if os.path.exists(gene_dir):
                if not glob.glob(os.path.join(gene_dir, "*", "*.pdb")):
                    missing_pdb.append(gene)
                if not glob.glob(os.path.join(gene_dir, "pockets", "*_fpocket")):
                    missing_fpocket.append(gene)
                if not glob.glob(os.path.join(gene_dir, "pockets", "*_p2rank")):
                    missing_p2rank.append(gene)

        # These used to pass whenever fewer than 19 genes were missing a result, so
        # 18 of 20 genes could have no structure, no fpocket and no p2rank and all
        # three still reported success. What is actually expected is that the genes
        # missing a result are exactly the genes with no reference structure.
        breakdown = reference_structure_breakdown(organism_dir)
        without_structure = set(breakdown['none']) if breakdown else set()

        for label, missing in (
            ('PDB structures', missing_pdb),
            ('fpocket results', missing_fpocket),
            ('p2rank results', missing_p2rank),
        ):
            unexpected = sorted(set(missing) - without_structure)
            if unexpected:
                validation["failed"].append(
                    f"❌ {len(unexpected)} protein(s) with a reference structure are missing "
                    f"{label}: {unexpected[:5]}"
                )
            elif missing:
                validation["passed"].append(
                    f"✅ Expected proteins have {label} "
                    f"({len(missing)} without, all lacking a reference structure)"
                )
            else:
                validation["passed"].append(f"✅ All {expected_count} proteins have {label}")

    else:
        validation["failed"].append("❌ Structures directory not found")

    # Check conservation
    conservation_dir = os.path.join(organism_dir, "conservation")
    if os.path.exists(conservation_dir):
        core_files = [f for f in os.listdir(conservation_dir) if f.startswith('core_') and f.endswith('.tsv')]
        if len(core_files) >= 1:
            validation["passed"].append(f"✅ Conservation analysis complete ({len(core_files)} core tables)")
        else:
            validation["warnings"].append("⚠️  Conservation analysis incomplete")

        corecruncher_dir = os.path.join(conservation_dir, "corecruncher_output")
        if os.path.exists(corecruncher_dir):
            families_core = os.path.join(corecruncher_dir, "families_core.txt")
            if os.path.exists(families_core):
                validation["passed"].append("✅ CoreCruncher families_core.txt found")
            else:
                validation["failed"].append("❌ CoreCruncher families_core.txt not found")

        roary_dir = os.path.join(conservation_dir, "roary_output")
        if os.path.exists(roary_dir):
            gene_presence = os.path.join(roary_dir, "results", "gene_presence_absence.csv")
            if os.path.exists(gene_presence):
                validation["passed"].append("✅ Roary gene_presence_absence.csv found")
            else:
                validation["failed"].append("❌ Roary gene_presence_absence.csv not found")

    else:
        validation["warnings"].append("⚠️  Conservation directory not found")

    # Check offtarget
    offtarget_dir = os.path.join(organism_dir, "offtarget")
    check_offtarget(offtarget_dir, validation, expected_count=expected_count, min_microbiome_genomes=4744,
                    organism_dir=organism_dir)


    # Check essentiality
    essentiality_dir = os.path.join(organism_dir, "essentiality")
    if os.path.exists(essentiality_dir):
        deg_file = os.path.join(essentiality_dir, "hit_in_deg.tsv")
        if os.path.exists(deg_file):
            validation["passed"].append("✅ Essentiality analysis complete")
        else:
            validation["failed"].append("❌ Essentiality analysis incomplete")

        deg_blast = os.path.join(essentiality_dir, "deg_blast.tsv")
        if os.path.exists(deg_blast):
            validation["passed"].append("✅ DEG BLAST results found")
        else:
            validation["failed"].append("❌ deg_blast.tsv not found")

    else:
        validation["failed"].append("❌ Essentiality directory not found")

    # Check localization
    localization_dir = os.path.join(organism_dir, "localization")
    if os.path.exists(localization_dir):
        psortb_file = os.path.join(localization_dir, "psortb_localization.tsv")
        if os.path.exists(psortb_file):
            validation["passed"].append("✅ Localization analysis complete")
        else:
            validation["warnings"].append("⚠️  Localization analysis incomplete")

    else:
        validation["failed"].append("❌ Localization directory not found")

    return validation


def print_validation_summary(validation):
    """Print a summary of validation results"""
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    
    if validation["passed"]:
        print("\n✅ PASSED CHECKS:")
        for check in validation["passed"]:
            print(f"  {check}")
    
    if validation["warnings"]:
        print("\n⚠️  WARNINGS:")
        for warning in validation["warnings"]:
            print(f"  {warning}")
    
    if validation["failed"]:
        print("\n❌ FAILED CHECKS:")
        for failure in validation["failed"]:
            print(f"  {failure}")
    
    print("\n" + "="*80)
    total_passed = len(validation["passed"])
    total_warnings = len(validation["warnings"])
    total_failed = len(validation["failed"])
    
    print(f"Total: {total_passed} passed, {total_warnings} warnings, {total_failed} failed")
    print("="*80)
    
    # Determine exit code
    if total_failed > 0:
        print("\n❌ TEST VALIDATION FAILED")
        return 1
    elif total_warnings > 0:
        print("\n✅ TEST VALIDATION PASSED WITH WARNINGS")
        print("\n⚠️  NOTE: Some warnings are expected if:")
        print("   • Certain databases are not downloaded (human offtarget, microbiome)")
        print("   • Some optional analyses were skipped")
        print("   • Running in stub mode")
        return 0
    else:
        print("\n✅ ALL VALIDATION CHECKS PASSED!")
        print("\nCongratulations! FastTarget pipeline is working correctly.")
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Validate FastTarget test results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate all results in the test directory
  python validate_test.py
  
  # Validate a specific test directory
  python validate_test.py --test_dir organism/test
  
  # Validate stub mode results (skip module checks)
  python validate_test.py --stub
        """
    )
    parser.add_argument(
        '--stub',
        action='store_true',
        help='Stub mode validation (only check main results, skip module outputs)'
    )
    parser.add_argument(
        '--test_dir',
        type=str,
        default='organism/test',
        help='Path to test directory (default: organism/test)'
    )
    
    args = parser.parse_args()

    test_dir = args.test_dir
    if not os.path.isdir(test_dir):
        print(f"❌ Test directory not found: {test_dir}")
        sys.exit(2)

    # If a results directory is provided by mistake, validate it directly.
    if os.path.basename(test_dir).startswith("test_results_"):
        results_dirs = [test_dir]
        organism_dir = os.path.dirname(test_dir)
    else:
        pattern = os.path.join(test_dir, "test_results_*")
        results_dirs = sorted(glob.glob(pattern), key=os.path.getctime)
        organism_dir = test_dir
        if not results_dirs:
            print(f"❌ No results directories found matching: {pattern}")
            print("\nMake sure you have run the test pipeline first.")
            sys.exit(2)

    check_all = not args.stub
    any_failed = False
    any_warnings = False

    if check_all:
        module_validation = validate_module_outputs(organism_dir)
        exit_code = print_validation_summary(module_validation)
        if exit_code == 1:
            any_failed = True
        if module_validation["warnings"]:
            any_warnings = True

    for results_dir in results_dirs:
        validation = validate_test_results(
            results_dir,
            check_all_modules=False,
            organism_dir=organism_dir
        )
        exit_code = print_validation_summary(validation)
        if exit_code == 1:
            any_failed = True
        if validation["warnings"]:
            any_warnings = True

    if any_failed:
        sys.exit(1)
    if any_warnings:
        sys.exit(0)
    sys.exit(0)
