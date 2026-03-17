#!/usr/bin/env nextflow

/*
 * Module: Structures - ColabFold (Parallelized)
 * ==============================================
 * Stage 2.5: Generates ColabFold models (optional) - ONE PROTEIN PER TASK
 *
 * Runs either:
 * - generate_colabfold_missing_single (only missing structures)
 * - generate_colabfold_all_single (all proteins)
 *
 * Each task processes a single locus_tag independently in its isolated working directory.
 * No race conditions, no shared state, full parallelization across GPU nodes.
 */

process COLABFOLD_SINGLE {
    tag "${locus_tag}"
    label 'gpu_process'

    input:
    tuple val(locus_tag), path(locus_structure_dir)
    val organism_name
    val output_path
    path gbk_file
    val amber_option
    val gpu_option
    val colabfold_all_models

    output:
    tuple val(locus_tag), path("colabfold_output"), emit: colabfold_results, optional: true

    script:
    def base_path = workflow.projectDir.parent
    def amber_option_py = amber_option ? 'True' : 'False'
    def gpu_option_py = gpu_option ? 'True' : 'False'
    def colabfold_all_models_py = colabfold_all_models ? 'True' : 'False'
    """#!/usr/bin/env python3

import sys
import os
import glob
import shutil
import json

# Add parent directory to path to import ftscripts
sys.path.insert(0, '${base_path}')

from ftscripts import structures, files, metadata

# Define locus_tag first
locus_tag = '${locus_tag}'

print('=' * 80)
print(f'COLABFOLD SINGLE: Processing {locus_tag}'.center(80))
print('─' * 80)

# Setup work directory for this protein
work_dir = os.getcwd()
locus_structure_dir = '${locus_structure_dir}'
organism_dir = os.path.join(work_dir, '${organism_name}')
structures_dir = os.path.join(organism_dir, 'structures')
genome_dir = os.path.join(organism_dir, 'genome')
os.makedirs(structures_dir, exist_ok=True)
os.makedirs(genome_dir, exist_ok=True)

# Load UniProt mapping
proteome_ids_file = os.path.join(os.path.dirname(locus_structure_dir), 'uniprot_files', f'uniprot_${organism_name}_id_mapping.json')
if files.file_check(proteome_ids_file):
    map_results = files.json_to_dict(proteome_ids_file)
else:
    map_results = {}
    print(f'WARNING: UniProt ID mapping file not found: {proteome_ids_file}')

# Copy input structure to work directory
task_locus_dir = os.path.join(structures_dir, locus_tag)
if os.path.exists(locus_structure_dir):
    shutil.copytree(locus_structure_dir, task_locus_dir, dirs_exist_ok=True)
    print(f'Copied locus structure: {locus_tag}')
else:
    print(f'ERROR: Input structure directory not found: {locus_structure_dir}')
    sys.exit(1)

# Copy GBK file for sequence extraction
gbk_source = '${gbk_file}'
gbk_dest = os.path.join(genome_dir, '${organism_name}.gbk')
shutil.copy2(gbk_source, gbk_dest)
print(f'Copied GBK file: {os.path.basename(gbk_dest)}')

# Run ColabFold
print(f'[2.5] Running ColabFold for {locus_tag}...')

if ${colabfold_all_models_py}:
    structures.generate_colabfold_all_single(
        locus_tag,
        work_dir,
        '${organism_name}',
        structures_dir,
        map_results,
        amber_option=${amber_option_py},
        gpu_option=${gpu_option_py}
    )
else:
    structures.generate_colabfold_missing_single(
        locus_tag,
        work_dir,
        '${organism_name}',
        structures_dir,
        map_results,
        amber_option=${amber_option_py},
        gpu_option=${gpu_option_py}
    )

# Prepare output - copy both the relaxed/unrelaxed model AND the colabfold_models directory
os.makedirs('colabfold_output', exist_ok=True)
colabfold_models_src = os.path.join(task_locus_dir, 'colabfold_models')

# Copy the colabfold_models directory
if os.path.exists(colabfold_models_src):
    shutil.copytree(colabfold_models_src, 'colabfold_output/models', dirs_exist_ok=True)
    print(f'Copied colabfold_models directory for {locus_tag}')
else:
    print(f'WARNING: No colabfold_models directory found for {locus_tag}')

# Copy the relaxed/unrelaxed model that was extracted to the root of task_locus_dir
relaxed_file = os.path.join(task_locus_dir, f"CB_{locus_tag}_relaxed1.pdb")
unrelaxed_file = os.path.join(task_locus_dir, f"CB_{locus_tag}_unrelaxed1.pdb")

if os.path.exists(relaxed_file):
    shutil.copy2(relaxed_file, 'colabfold_output')
    print(f'Copied relaxed model for {locus_tag}')
elif os.path.exists(unrelaxed_file):
    shutil.copy2(unrelaxed_file, 'colabfold_output')
    print(f'Copied unrelaxed model for {locus_tag}')
else:
    print(f'WARNING: No relaxed/unrelaxed model found for {locus_tag}')

# Copy any other PDB files from task_locus_dir root (extracted chains, etc.)
for pdb_file in glob.glob(os.path.join(task_locus_dir, '*.pdb')):
    shutil.copy2(pdb_file, 'colabfold_output')

print(f'COLABFOLD_SINGLE completed successfully for {locus_tag}')
"""

    stub:
    """
    mkdir -p colabfold_output/models
    touch colabfold_output/models/${locus_tag}_unrelaxed1.pdb
    echo "STUB: ColabFold model for ${locus_tag}"
    """
}


process COLABFOLD_COLLECT {
    tag "${organism_name}"
    label 'low_resources'
    publishDir "${output_path}/${organism_name}", mode: 'copy', pattern: "structures/**"

    input:
    val organism_name
    val output_path
    path structure_dir
    val colabfold_results

    output:
    path "${organism_name}/structures", emit: structure_dir
    path "${organism_name}/structures/**", emit: all_structures
    val organism_name, emit: organism_name

    script:
    """#!/usr/bin/env python3

import sys
import os
import glob
import shutil
from pathlib import Path

print('=' * 80)
print('STAGE 2.5 COLLECT: Merging ColabFold results'.center(80))
print('─' * 80)

organism_name = '${organism_name}'
base_path = os.path.join(organism_name, 'structures')

# Ensure output directory exists
os.makedirs(base_path, exist_ok=True)

# First, copy ALL existing structures from upstream (STRUCTURES_EXTRACT_CHAINS_COLLECT output)
staged_structures = '${structure_dir}'
if os.path.exists(staged_structures):
    if os.path.islink(staged_structures):
        staged_structures = os.path.realpath(staged_structures)
    print(f'Copying existing structures from: {staged_structures}')
    for item in os.listdir(staged_structures):
        src_path = os.path.join(staged_structures, item)
        dst_path = os.path.join(base_path, item)
        if os.path.isdir(src_path):
            shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
        else:
            shutil.copy2(src_path, dst_path)
    print(f'Copied all existing structures to {base_path}')
else:
    print(f'WARNING: input structures directory not found: {staged_structures}')

colabfold_data_flat = ${groovy.json.JsonOutput.toJson(colabfold_results)}

# Convert flat list into pairs: [tag1, path1, tag2, path2, ...] -> [(tag1, path1), (tag2, path2), ...]
colabfold_data = [(colabfold_data_flat[i], colabfold_data_flat[i+1]) for i in range(0, len(colabfold_data_flat), 2)]

print(f'Processing {len(colabfold_data)} ColabFold results...')

successful_count = 0
failed_count = 0
for locus_tag, result_dir in colabfold_data:
    colabfold_result_dir = Path(result_dir)
    
    if not colabfold_result_dir.exists():
        print(f'WARNING: ColabFold output not found for {locus_tag} at {result_dir}')
        failed_count += 1
        continue

    # Merge ColabFold models into existing locus structure directory
    locus_structure_dir = os.path.join(base_path, locus_tag)
    os.makedirs(locus_structure_dir, exist_ok=True)
    
    # Copy colabfold_models directory
    models_src = colabfold_result_dir / 'models'
    models_dest = os.path.join(locus_structure_dir, 'colabfold_models')

    if models_src.exists():
        shutil.copytree(models_src, models_dest, dirs_exist_ok=True)
        print(f'Added ColabFold models directory for {locus_tag}')
    else:
        print(f'WARNING: No models directory for {locus_tag} at {models_src}')
        failed_count += 1
        continue
    
    # Copy relaxed/unrelaxed model files from colabfold_output root
    for pdb_file in colabfold_result_dir.glob('*.pdb'):
        dest_pdb = os.path.join(locus_structure_dir, pdb_file.name)
        shutil.copy2(pdb_file, dest_pdb)
        print(f'Copied model file {pdb_file.name} for {locus_tag}')
    
    successful_count += 1

print(f'COLLECT: Successfully merged {successful_count}/{len(colabfold_data)} ColabFold results')
if failed_count > 0:
    print(f'WARNING: {failed_count} results were not merged')
print('Stage 2.5 COLLECT completed successfully')
"""

    stub:
    """
    mkdir -p ${organism_name}/structures/gene_example
    touch ${organism_name}/structures/gene_example/CB_gene_example_unrelaxed1.pdb
    echo "STUB: ColabFold collect for ${organism_name}"
    """
}
