import configuration
import fasttarget
import os
import pandas as pd
import argparse
import multiprocessing
from ftscripts import files, structures, pathways, offtargets, genome, essentiality, metadata
from datetime import datetime
import logging
import sys
from ftscripts.logger import logger 
import shutil
import time
import glob

# Expected genes in the test dataset
EXPECTED_GENES = [
    'MPN_RS00140', 'MPN_RS00300', 'MPN_RS00445', 'MPN_RS00890', 'MPN_RS01250',
    'MPN_RS01270', 'MPN_RS01385', 'MPN_RS01390', 'MPN_RS01460', 'MPN_RS01665',
    'MPN_RS01670', 'MPN_RS02025', 'MPN_RS02430', 'MPN_RS03100', 'MPN_RS02775',
    'MPN_RS03450', 'MPN_RS03550', 'MPN_RS03565', 'MPN_RS03570', 'MPN_RS04200'
]

def validate_test_results(test_path, organism_name, test_config):
    """
    Validate that expected output files were created and data integrity is maintained.
    
    :param test_path: Path to test organism directory
    :param organism_name: Name of the test organism
    :param test_config: Configuration object to check enabled modules
    :return: Dictionary with validation results
    """
    print("\n" + "="*80)
    print("VALIDATING TEST RESULTS")
    print("="*80)
    
    validation = {
        "passed": [],
        "failed": [],
        "warnings": []
    }
    
    # Find the latest results directory
    results_pattern = os.path.join(test_path, f"{organism_name}_results_*")
    results_dirs = glob.glob(results_pattern)
    if not results_dirs:
        validation["failed"].append("❌ No results directory found")
        return validation
    
    results_dir = max(results_dirs, key=os.path.getctime)  # Get most recent
    print(f"📁 Checking results in: {os.path.basename(results_dir)}")
    
    # Check main results table first - this is critical
    main_results = os.path.join(results_dir, f"{organism_name}_results_table.tsv")
    if not os.path.exists(main_results):
        validation["failed"].append("❌ Main results table not found")
        return validation
    
    try:
        df = pd.read_csv(main_results, sep='\t')
        
        # CRITICAL: Validate gene count and IDs
        print("\n" + "-"*80)
        print("DATA INTEGRITY CHECKS")
        print("-"*80)
        
        expected_count = len(EXPECTED_GENES)
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
        
        # GENERAL: Check genome files
        print("\n" + "-"*80)
        print("GENOME FILES")
        print("-"*80)
        genome_dir = os.path.join(test_path, "genome")
        genome_files = {
            'FAA': f"{organism_name}.faa",
            'GFF': f"{organism_name}.gff",
            'FNA': f"{organism_name}.fna"
        }
        for file_type, filename in genome_files.items():
            filepath = os.path.join(genome_dir, filename)
            if os.path.exists(filepath):
                size = os.path.getsize(filepath)
                validation["passed"].append(f"✅ Genome {file_type}: {size:,} bytes")
            else:
                validation["failed"].append(f"❌ Genome {file_type} not found: {filename}")
        
        # MODULE-BASED VALIDATION
        print("\n" + "-"*80)
        print("MODULE VALIDATION")
        print("-"*80)
        
        # METABOLIC MODULE
        if test_config.metabolism_pathwaytools['enabled'] or test_config.metabolism_sbml['enabled']:
            print("\n--- Metabolism Module ---")
            metabolism_dir = os.path.join(test_path, "metabolism")
            
            # Check PathwayTools files
            if test_config.metabolism_pathwaytools['enabled']:
                ptools_files = glob.glob(os.path.join(metabolism_dir, "PTOOLS_*"))
                if ptools_files:
                    validation["passed"].append(f"✅ PathwayTools files: {len(ptools_files)} found")
                else:
                    validation["failed"].append("❌ No PathwayTools output files found")
                
                # Validate PTOOLS columns
                ptools_bc = 'PTOOLS_betweenness_centrality'
                ptools_edges = 'PTOOLS_edges'
                if ptools_bc in df.columns:
                    # Check numeric 0-1
                    numeric_vals = df[ptools_bc].dropna()
                    if len(numeric_vals) > 0:
                        if (numeric_vals >= 0).all() and (numeric_vals <= 1).all():
                            validation["passed"].append(f"✅ {ptools_bc}: Valid range 0-1")
                        else:
                            validation["failed"].append(f"❌ {ptools_bc}: Values outside 0-1 range")
                    else:
                        validation["warnings"].append(f"⚠️  {ptools_bc}: All values are empty")
                
                if ptools_edges in df.columns:
                    # Check integer
                    edge_vals = df[ptools_edges].dropna()
                    if len(edge_vals) > 0:
                        if (edge_vals == edge_vals.astype(int)).all():
                            validation["passed"].append(f"✅ {ptools_edges}: Valid integers")
                        else:
                            validation["failed"].append(f"❌ {ptools_edges}: Contains non-integer values")
                    else:
                        validation["warnings"].append(f"⚠️  {ptools_edges}: All values are empty")
                
                # Check chokepoint columns (None or reaction IDs)
                ptools_choke_cols = ['PTOOLS_producing_chokepoints', 'PTOOLS_consuming_chokepoints', 
                                     'PTOOLS_both_chokepoints']
                for col in ptools_choke_cols:
                    if col in df.columns:
                        has_values = df[col].notna().any()
                        if has_values:
                            validation["passed"].append(f"✅ {col}: Contains reaction IDs")
                        # Empty is OK for chokepoints
            
            # Check MetaGraphTools files
            if test_config.metabolism_sbml['enabled']:
                mgt_dirs = glob.glob(os.path.join(metabolism_dir, "MGT_*"))
                if mgt_dirs:
                    validation["passed"].append(f"✅ MetaGraphTools output: {len(mgt_dirs)} directory(ies)")
                else:
                    validation["failed"].append("❌ No MetaGraphTools output directories found")
                
                # Validate MGT columns
                mgt_bc = 'MGT_betweenness_centrality'
                mgt_edges = 'MGT_edges'
                if mgt_bc in df.columns:
                    # Check numeric 0-1
                    numeric_vals = df[mgt_bc].dropna()
                    if len(numeric_vals) > 0:
                        if (numeric_vals >= 0).all() and (numeric_vals <= 1).all():
                            validation["passed"].append(f"✅ {mgt_bc}: Valid range 0-1")
                        else:
                            validation["failed"].append(f"❌ {mgt_bc}: Values outside 0-1 range")
                    else:
                        validation["warnings"].append(f"⚠️  {mgt_bc}: All values are empty")
                
                if mgt_edges in df.columns:
                    # Check integer
                    edge_vals = df[mgt_edges].dropna()
                    if len(edge_vals) > 0:
                        if (edge_vals == edge_vals.astype(int)).all():
                            validation["passed"].append(f"✅ {mgt_edges}: Valid integers")
                        else:
                            validation["failed"].append(f"❌ {mgt_edges}: Contains non-integer values")
                    else:
                        validation["warnings"].append(f"⚠️  {mgt_edges}: All values are empty")
                
                # Check MGT chokepoint columns (TRUE/FALSE, at least one TRUE)
                mgt_choke_cols = ['MGT_consuming_chokepoints', 'MGT_producing_chokepoints']
                for col in mgt_choke_cols:
                    if col in df.columns:
                        true_count = (df[col] == True).sum()
                        if true_count > 0:
                            validation["passed"].append(f"✅ {col}: {true_count} TRUE values found")
                        else:
                            validation["warnings"].append(f"⚠️  {col}: No TRUE values found")
        
        # LOCALIZATION MODULE
        if test_config.psortb['enabled']:
            print("\n--- Localization Module ---")
            localization_dir = os.path.join(test_path, "localization")
            
            # Check PSORTb CSV file
            psortb_csv = glob.glob(os.path.join(localization_dir, "*_psortb_*.txt"))
            if psortb_csv:
                validation["passed"].append(f"✅ PSORTb output: {len(psortb_csv)} file(s)")
            else:
                validation["failed"].append("❌ PSORTb output file not found")
            
            # Check localization column (might not be in main table if not merged yet)
            psortb_col_found = False
            for col in df.columns:
                if 'psort' in col.lower():
                    psortb_col_found = True
                    non_empty = df[col].notna().sum()
                    if non_empty == expected_count:
                        validation["passed"].append(f"✅ {col}: All {non_empty} proteins have values")
                    elif non_empty > 0:
                        validation["warnings"].append(f"⚠️  {col}: Only {non_empty}/{expected_count} proteins have values")
                    else:
                        validation["failed"].append(f"❌ {col}: Column has no data")
                    break
            
            if not psortb_col_found:
                validation["warnings"].append("⚠️  PSORTb column not found in results table (may need manual integration)")
        
        # ESSENTIALITY MODULE
        if test_config.deg['enabled']:
            print("\n--- Essentiality Module ---")
            essentiality_dir = os.path.join(test_path, "essentiality")
            
            # Check DEG blast file
            deg_blast = os.path.join(essentiality_dir, "deg_blast.tsv")
            if os.path.exists(deg_blast):
                size = os.path.getsize(deg_blast)
                validation["passed"].append(f"✅ DEG BLAST results: {size:,} bytes")
            else:
                validation["failed"].append("❌ deg_blast.tsv not found")
            
            # Check hit_in_deg column (at least one TRUE)
            if 'hit_in_deg' in df.columns:
                true_count = (df['hit_in_deg'] == True).sum()
                total_not_null = df['hit_in_deg'].notna().sum()
                if true_count > 0:
                    validation["passed"].append(f"✅ hit_in_deg: {true_count} essential genes found")
                else:
                    validation["warnings"].append(f"⚠️  hit_in_deg: No essential genes found ({total_not_null} values)")
            else:
                validation["failed"].append("❌ hit_in_deg: Column missing")
        
        # CONSERVATION MODULE
        if test_config.core['enabled']:
            print("\n--- Conservation Module ---")
            core_dir = os.path.join(test_path, "conservation")
            
            # Check CoreCruncher
            if test_config.core.get('corecruncher', False):
                corecruncher_dir = os.path.join(core_dir, "corecruncher_output")
                if os.path.exists(corecruncher_dir):
                    faa_files = glob.glob(os.path.join(corecruncher_dir, "faa", "*.faa"))
                    families_core = os.path.join(corecruncher_dir, "families_core.txt")
                    
                    if len(faa_files) > 70:
                        validation["passed"].append(f"✅ CoreCruncher FAA files: {len(faa_files)} found")
                    else:
                        validation["warnings"].append(f"⚠️  CoreCruncher FAA files: Only {len(faa_files)} found (expected >70)")
                    
                    if os.path.exists(families_core):
                        validation["passed"].append("✅ CoreCruncher families_core.txt found")
                    else:
                        validation["failed"].append("❌ CoreCruncher families_core.txt not found")
                else:
                    validation["failed"].append("❌ CoreCruncher output directory not found")
                
                # Check core_corecruncher column
                if 'core_corecruncher' in df.columns:
                    true_count = (df['core_corecruncher'] == True).sum()
                    if true_count > 0:
                        validation["passed"].append(f"✅ core_corecruncher: {true_count} core genes")
                    else:
                        validation["warnings"].append("⚠️  core_corecruncher: No core genes found")
            
            # Check Roary
            if test_config.core.get('roary', False):
                roary_dir = os.path.join(core_dir, "roary_output")
                if os.path.exists(roary_dir):
                    gff_dir = os.path.join(roary_dir, "gff")
                    results_dir = os.path.join(roary_dir, "results")
                    
                    if os.path.exists(gff_dir):
                        gff_files = glob.glob(os.path.join(gff_dir, "*.gff"))
                        if len(gff_files) > 70:
                            validation["passed"].append(f"✅ Roary GFF files: {len(gff_files)} found")
                        else:
                            validation["warnings"].append(f"⚠️  Roary GFF files: Only {len(gff_files)} found (expected >70)")
                    else:
                        validation["failed"].append("❌ Roary gff directory not found")
                    
                    gene_presence = os.path.join(results_dir, "gene_presence_absence.csv")
                    if os.path.exists(gene_presence):
                        validation["passed"].append("✅ Roary gene_presence_absence.csv found")
                    else:
                        validation["failed"].append("❌ Roary gene_presence_absence.csv not found")
                else:
                    validation["failed"].append("❌ Roary output directory not found")
                
                # Check core_roary column
                if 'core_roary' in df.columns:
                    true_count = (df['core_roary'] == True).sum()
                    if true_count > 0:
                        validation["passed"].append(f"✅ core_roary: {true_count} core genes")
                    else:
                        validation["warnings"].append("⚠️  core_roary: No core genes found")
        
        # OFFTARGET MODULE
        if test_config.offtarget['enabled']:
            print("\n--- Offtarget Module ---")
            offtarget_dir = os.path.join(test_path, "offtarget")
            
            # Check human offtarget files
            if test_config.offtarget.get('human', False):
                human_blast = os.path.join(offtarget_dir, "human_offtarget_blast.tsv")
                human_csv = os.path.join(offtarget_dir, "human_offtarget.csv")
                
                if os.path.exists(human_blast):
                    size = os.path.getsize(human_blast)
                    validation["passed"].append(f"✅ Human offtarget BLAST: {size:,} bytes")
                else:
                    validation["failed"].append("❌ human_offtarget_blast.tsv not found")
                
                if os.path.exists(human_csv):
                    validation["passed"].append("✅ Human offtarget CSV found")
                else:
                    validation["warnings"].append("⚠️  human_offtarget.csv not found")
                
                # Check human_offtarget column (numbers or "no_hit")
                if 'human_offtarget' in df.columns:
                    non_null = df['human_offtarget'].notna().sum()
                    if non_null == expected_count:
                        # Check values are numeric or "no_hit"
                        has_no_hit = (df['human_offtarget'] == 'no_hit').any()
                        has_numeric = df['human_offtarget'][df['human_offtarget'] != 'no_hit'].apply(
                            lambda x: str(x).replace('.', '').replace('-', '').isdigit() if pd.notna(x) else False
                        ).any()
                        validation["passed"].append(f"✅ human_offtarget: All proteins have values")
                    else:
                        validation["warnings"].append(f"⚠️  human_offtarget: Only {non_null}/{expected_count} have values")
                else:
                    validation["failed"].append("❌ human_offtarget: Column missing")
            
            # Check FoldSeek results
            if test_config.offtarget.get('foldseek_human', False):
                foldseek_dir = os.path.join(offtarget_dir, "foldseek_results")
                if os.path.exists(foldseek_dir):
                    validation["passed"].append("✅ FoldSeek results directory found")
                else:
                    validation["warnings"].append("⚠️  FoldSeek results directory not found")
            
            # Check microbiome offtarget
            if test_config.offtarget.get('microbiome', False):
                microbiome_files = glob.glob(os.path.join(offtarget_dir, "*microbiome*.tsv")) + \
                                  glob.glob(os.path.join(offtarget_dir, "species_blast_results", "*microbiome*.tsv"))
                if microbiome_files:
                    validation["passed"].append(f"✅ Microbiome offtarget: {len(microbiome_files)} file(s)")
                else:
                    validation["warnings"].append("⚠️  No microbiome offtarget files found")
        
        # STRUCTURES MODULE
        if test_config.structures['enabled']:
            print("\n--- Structures Module ---")
            structures_dir = os.path.join(test_path, "structures")
            
            if os.path.exists(structures_dir):
                # Check for 20 protein directories
                protein_dirs = [d for d in os.listdir(structures_dir) 
                               if os.path.isdir(os.path.join(structures_dir, d)) and d.startswith('MPN_')]
                
                if len(protein_dirs) == expected_count:
                    validation["passed"].append(f"✅ Protein structure directories: {len(protein_dirs)} found")
                else:
                    validation["warnings"].append(f"⚠️  Protein directories: {len(protein_dirs)}/{expected_count} found")
                
                # Check each directory has PDB, fpocket, and p2rank
                missing_pdb = []
                missing_fpocket = []
                missing_p2rank = []
                
                for gene in EXPECTED_GENES:
                    gene_dir = os.path.join(structures_dir, gene)
                    if os.path.exists(gene_dir):
                        # Check for PDB files
                        pdb_files = glob.glob(os.path.join(gene_dir, "*", "*.pdb"))
                        if not pdb_files:
                            missing_pdb.append(gene)
                        
                        # Check for fpocket results
                        fpocket_dirs = glob.glob(os.path.join(gene_dir, "pockets", "*_fpocket"))
                        if not fpocket_dirs:
                            missing_fpocket.append(gene)
                        
                        # Check for p2rank results
                        p2rank_dirs = glob.glob(os.path.join(gene_dir, "pockets", "*_p2rank"))
                        if not p2rank_dirs:
                            missing_p2rank.append(gene)
                
                if not missing_pdb:
                    validation["passed"].append(f"✅ All proteins have PDB structures")
                else:
                    validation["warnings"].append(f"⚠️  {len(missing_pdb)} proteins missing PDB: {missing_pdb[:3]}")
                
                if not missing_fpocket:
                    validation["passed"].append(f"✅ All proteins have fpocket results")
                else:
                    validation["warnings"].append(f"⚠️  {len(missing_fpocket)} proteins missing fpocket: {missing_fpocket[:3]}")
                
                if not missing_p2rank:
                    validation["passed"].append(f"✅ All proteins have p2rank results")
                else:
                    validation["warnings"].append(f"⚠️  {len(missing_p2rank)} proteins missing p2rank: {missing_p2rank[:3]}")
                
                # Check structure columns in dataframe
                structure_cols = ['uniprot', 'structure', 'druggability_score', 
                                 'fpocket_pocket', 'p2rank_probability', 'p2rank_pocket']
                
                for col in structure_cols:
                    if col not in df.columns:
                        validation["failed"].append(f"❌ Structure column missing: {col}")
                
                # Check druggability_score is 0-1
                if 'druggability_score' in df.columns:
                    drug_scores = df['druggability_score'].dropna()
                    if len(drug_scores) > 0:
                        if (drug_scores >= 0).all() and (drug_scores <= 1).all():
                            validation["passed"].append(f"✅ druggability_score: Valid range 0-1 ({len(drug_scores)} values)")
                        else:
                            validation["failed"].append(f"❌ druggability_score: Values outside 0-1 range")
                
                # Check pocket data columns have values
                pocket_cols = ['fpocket_pocket', 'p2rank_probability', 'p2rank_pocket']
                for col in pocket_cols:
                    if col in df.columns:
                        non_null = df[col].notna().sum()
                        if non_null > 0:
                            validation["passed"].append(f"✅ {col}: {non_null} proteins with data")
                        else:
                            validation["warnings"].append(f"⚠️  {col}: No pocket data found")
            else:
                validation["failed"].append("❌ Structures directory not found")
        
    except Exception as e:
        validation["failed"].append(f"❌ Error validating results table: {str(e)}")
        logging.exception(e)
    
    # Print summary
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    for item in validation["passed"]:
        print(item)
    for item in validation["warnings"]:
        print(item)
    for item in validation["failed"]:
        print(item)
    print("="*80)
    
    total_checks = len(validation["passed"]) + len(validation["failed"])
    success_rate = len(validation["passed"]) / total_checks * 100 if total_checks > 0 else 0
    print(f"\n📊 Result: {len(validation['passed'])}/{total_checks} checks passed ({success_rate:.1f}%)")
    print(f"   Warnings: {len(validation['warnings'])}")
    
    return validation

def run_test(test_config, databases_path, output_path):
    """
    Run the fasttarget pipeline with the test configuration.
    
    :param test_config: Configuration object
    :param databases_path: Path to the databases directory
    :param output_path: Path to the output directory
    :return: Tuple of (result, validation_dict, elapsed_time)
    """
    try:
        print("\n" + "="*80)
        print("FASTTARGET INSTALLATION TEST")
        print("="*80)
        print(f"This test runs the complete pipeline with a 20-protein dataset")
        print(f"to verify that all dependencies and databases are installed correctly.")
        print(f"\nTest organism: {test_config.organism['name']}")
        print(f"Expected proteins: {len(EXPECTED_GENES)}")
        print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*80)
        
        start_time = time.time()
        result = fasttarget.main(test_config, databases_path, output_path)
        elapsed_time = time.time() - start_time
        
        print("\n" + "="*80)
        print(f"Pipeline completed in {elapsed_time/60:.2f} minutes ({elapsed_time:.1f} seconds)")
        print("="*80)
        
        # Validate results
        test_path = os.path.join(output_path, test_config.organism['name'])
        validation = validate_test_results(test_path, test_config.organism['name'], test_config)
        
        return result, validation, elapsed_time
        
    except Exception as e:
        logging.exception(f"An error occurred during the test run: {e}")
        raise


if __name__ == "__main__":

    # Define paths
    base_path = os.path.dirname(os.path.abspath(__file__))
    test_path = os.path.join(base_path, 'organism', 'test')


    databases_default_path = os.path.join(base_path, 'databases')
    output_default_path = os.path.join(base_path, 'organism')

    parser = argparse.ArgumentParser(description='Test script')
    parser.add_argument('--databases_path', type=str, default=databases_default_path, help='Path to the databases directory')
    parser.add_argument('--output_path', type=str, default=output_default_path, help='Path to the output directory')
    parser.add_argument('--container_engine', type=str, default='docker', help='Container engine to use (docker or singularity)')
    args = parser.parse_args()
    
    # COMPLETE TEST CONFIGURATION - All modules enabled
    config_dict = {
        "organism": {
            "name": "test",
            "tax_id": 2104,
            "strain_taxid": 272634,
            "gbk_file": f'{test_path}/test.gbk'
        },
        "cpus": None,
        "container_engine": args.container_engine,
        "metabolism-PathwayTools": {
            "enabled": True,
            "sbml_file": f'{test_path}/test.sbml',
            "chokepoint_file": f'{test_path}/test_chokepoints.txt',
            "smarttable_file": f'{test_path}/test_genes_smarttable.tsv'
        },
        "metabolism-SBML": {
            "enabled": True,
            "sbml_file": f'{test_path}/iEG158_mpn.xml',
            "filter_file": None
        },
        "structures": {
            "enabled": True
        },
        "core": {
            "enabled": True,
            "roary": True,
            "corecruncher": True,
            "min_identity": 95,
            "min_core_freq": 99
        },
        "offtarget": {
            "enabled": True,
            "human": True,
            "microbiome": True,
            "microbiome_identity_filter": 50,
            "microbiome_coverage_filter": 70,
            "foldseek_human": True
        },
        "deg": {
            "enabled": True,
            "deg_identity_filter": 50,
            "deg_coverage_filter": 70
        },
        "psortb": {
            "enabled": True,
            "gram_type": "n"
        },
        "metadata": {
            "enabled": True,
            "meta_tables": [f'{test_path}/test_meta_table1.txt']
        }
    }

    test_config = configuration.Config(config_dict)
    
    # Run test and validate
    try:
        result, validation, elapsed_time = run_test(test_config, args.databases_path, args.output_path)
        
        # Exit with appropriate code
        if validation["failed"]:
            print("\n" + "="*80)
            print("❌ TEST FAILED - Some Checks Did Not Pass")
            print("="*80)
            print("\nCommon reasons for test failures:")
            print("\n1. MISSING DATABASES (Expected if not downloaded):")
            print("   • Human proteome database for offtarget analysis for example")
            print("\n2. MISSING DEPENDENCIES:")
            print("\n3. CONFIGURATION ISSUES:")
            print("   • Invalid parameters")
            print("\n⚠️  NOTE: If you intentionally didn't download certain databases,")
            print("   failures in those modules are EXPECTED and NOT a problem.")

            sys.exit(1)
        elif validation["warnings"]:
            print("\n" + "="*80)
            print("✅ TEST PASSED WITH WARNINGS")
            print("="*80)
            print("The pipeline completed successfully but some optional features")
            print("Check warnings above.")
            print("="*80)
            sys.exit(0)
        else:
            print("\n" + "="*80)
            print("✅ ALL TESTS PASSED!")
            print("="*80)
            print("Congratulations! FastTarget is fully installed and working.")
            print("All modules completed successfully.")
            print("="*80)
            sys.exit(0)
            
    except Exception as e:
        print("\n" + "="*80)
        print(f"❌ TEST CRASHED")
        print("="*80)
        print(f"Error: {str(e)}")
        print("\nThe pipeline encountered a fatal error.")
        print("Check the log file for details:")
        log_files = glob.glob(os.path.join(base_path, "logs", "fasttarget_*.log"))
        if log_files:
            latest_log = max(log_files, key=os.path.getctime)
            print(f"  {latest_log}")
        print("="*80)
        sys.exit(2)
