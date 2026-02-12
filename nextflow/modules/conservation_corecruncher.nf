#!/usr/bin/env nextflow

/*
 * Module: Conservation - CoreCruncher Analysis
 * ============================================
 * Runs CoreCruncher core genome analysis on downloaded genomes
 * 
 * Steps:
 * 1. Prepare FAA files in corecruncher_output/faa/
 * 2. Run CoreCruncher with specified parameters
 * 3. Parse CoreCruncher output to generate core genes table
 * 
 * Can run in parallel with Roary
 */

process CONSERVATION_CORECRUNCHER {
    tag "${organism_name}"
    label 'high_resources'
    publishDir "${output_path}", mode: 'copy', pattern: "${organism_name}/conservation/**"
    
    input:
    val organism_name
    val output_path
    path faa_dir
    path genome_files
    val min_core_freq
    val min_identity
    val cpus
    val container_engine
    
    output:
    path "${organism_name}/conservation/core_corecruncher.tsv", emit: core_table
    path "${organism_name}/conservation/corecruncher_output/", emit: corecruncher_dir
    path "${organism_name}/conservation/**", emit: all_corecruncher
    val organism_name, emit: organism_name
    
    script:
    def base_path = workflow.projectDir.parent
    """#!/usr/bin/env python3
    
import sys
import os
import shutil
import glob

# Add parent directory to path to import ftscripts
sys.path.insert(0, '${base_path}')

import shutil
import glob
from ftscripts import genome

print('=' * 80)
print('CORECRUNCHER CORE GENOME ANALYSIS'.center(80))
print('=' * 80)

# Setup work directory structure
work_dir = os.getcwd()
conservation_dir = os.path.join(work_dir, '${organism_name}', 'conservation')
corecruncher_output_dir = os.path.join(conservation_dir, 'corecruncher_output')
os.makedirs(os.path.join(corecruncher_output_dir, 'faa'), exist_ok=True)

print(f'Working in: {work_dir}')
print(f'Conservation directory: {conservation_dir}')
print(f'CoreCruncher output directory: {corecruncher_output_dir}')

# Copy genome files to expected location
genome_dir = os.path.join(work_dir, '${organism_name}', 'genome')
os.makedirs(genome_dir, exist_ok=True)
for genome_file in glob.glob('${organism_name}.*'):
    target_file = os.path.join(genome_dir, os.path.basename(genome_file))
    shutil.copy2(genome_file, target_file)
    print(f'Copied genome file: {os.path.basename(genome_file)}')

# Copy FAA files from staged faa directory to corecruncher faa directory
staged_faa_dir = '${faa_dir}'
faa_target_dir = os.path.join(corecruncher_output_dir, 'faa')

print(f'Copying FAA files from: {staged_faa_dir}')
faa_files = glob.glob(os.path.join(staged_faa_dir, '*.faa'))

for faa_file in faa_files:
    target_file = os.path.join(faa_target_dir, os.path.basename(faa_file))
    shutil.copy2(faa_file, target_file)
    print(f'  - Copied: {os.path.basename(faa_file)}')

print(f'Total FAA files copied: {len(faa_files)}')
print('')
print(f'Parameters:')
print(f'  - Min core frequency: ${min_core_freq}%')
print(f'  - Min identity: ${min_identity}%')
print(f'  - CPUs: ${cpus}')

print('[1] Running CoreCruncher...')
print('')
genome.core_genome_programs(
    work_dir,
    '${organism_name}',
    ${min_core_freq},
    ${min_identity},
    ${cpus},
    program_list=['corecruncher'],
    container_engine='${container_engine}'
)

print('[2] Parsing CoreCruncher results...')
print('')
df_cc = genome.corecruncher_output(
    work_dir,
    '${organism_name}'
)

print(f'CoreCruncher analysis completed')
print('')
print(f'  - Core genes found: {len(df_cc)}')
"""
    
    stub:
    """
    mkdir -p conservation/corecruncher_output/faa
    mkdir -p conservation/corecruncher_output/CC
    mkdir -p conservation/corecruncher_output/core
    touch conservation/core_corecruncher.tsv
    touch conservation/corecruncher_output/summary.txt
    touch conservation/corecruncher_output/families_core.txt
    echo "STUB: CoreCruncher analysis for ${organism_name}"
    """
}
