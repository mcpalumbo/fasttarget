#!/usr/bin/env nextflow

/*
 * Module: Structures - ColabFold
 * ===============================
 * Stage 2.5: Generates ColabFold models (optional)
 *
 * Runs either:
 * - make_models_colabfold (only missing structures)
 * - make_models_colabfold_all_proteins (all proteins)
 */

process STRUCTURES_COLABFOLD {
    tag "${organism_name}"
    label 'high_resources'
    publishDir "${output_path}", mode: 'copy', pattern: "${organism_name}/structures/**"

    input:
    val organism_name
    val output_path
    path genome_files
    path structure_dir
    val amber_option
    val gpu_option
    val colabfold_all_models

    output:
    path "${organism_name}/structures", emit: structure_dir
    path "${organism_name}/structures/**", emit: all_structures
    val organism_name, emit: organism_name

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

# Add parent directory to path to import ftscripts
sys.path.insert(0, '${base_path}')

from ftscripts import structures

print('=' * 80)
print('STAGE 2.5: COLABFOLD MODEL GENERATION'.center(80))
print('─' * 80)

# Setup work directory
work_dir = os.getcwd()
organism_dir = os.path.join(work_dir, '${organism_name}')
structures_dir = os.path.join(organism_dir, 'structures')
os.makedirs(structures_dir, exist_ok=True)

# Copy genome files to expected location
genome_dir = os.path.join(organism_dir, 'genome')
os.makedirs(genome_dir, exist_ok=True)
for genome_file in glob.glob('${organism_name}.*'):
    target_file = os.path.join(genome_dir, os.path.basename(genome_file))
    shutil.copy2(genome_file, target_file)
    print(f'Copied genome file: {os.path.basename(genome_file)}')

# Copy staged structures from upstream process (Nextflow channel input)
staged_structures = '${structure_dir}'
if os.path.exists(staged_structures):
    if os.path.islink(staged_structures):
        staged_structures = os.path.realpath(staged_structures)
    print(f'Copying staged structures from: {staged_structures}')
    shutil.copytree(staged_structures, structures_dir, dirs_exist_ok=True)
else:
    print(f'WARNING: staged structures directory not found: {staged_structures}')

print('[2.5] Running ColabFold models...')
if ${colabfold_all_models_py}:
    structures.make_models_colabfold_all_proteins(
        work_dir,
        '${organism_name}',
        amber_option=${amber_option_py},
        gpu_option=${gpu_option_py}
    )
else:
    structures.make_models_colabfold(
        work_dir,
        '${organism_name}',
        amber_option=${amber_option_py},
        gpu_option=${gpu_option_py}
    )

print('Stage 2.5 completed successfully')
"""

    stub:
    """
    mkdir -p ${organism_name}/structures/gene_example
    touch ${organism_name}/structures/gene_example/CB_gene_example_unrelaxed1.pdb
    echo "STUB: ColabFold models for ${organism_name}"
    """
}
