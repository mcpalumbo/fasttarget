#!/usr/bin/env nextflow

/*
 * Module: Offtarget - Microbiome
 * ===============================
 * Runs BLAST against gut microbiome genomes to identify potential offtargets
 * 
 * Steps:
 * 1. Run BLASTP searches against all microbiome species genomes
 * 2. Parse results with identity and coverage filters
 * 3. Generate normalized scores and counts
 * 
 * Can run in parallel with human and foldseek offtarget modules
 */

process OFFTARGET_MICROBIOME {
    tag "${organism_name}"
    label 'blast_process'
    publishDir "${output_path}", mode: 'copy', pattern: "${organism_name}/offtarget/**"
    
    input:
    path genome_files
    val organism_name
    val output_path
    val databases_path
    val microbiome_catalogues_json
    val cpus
    
    output:
    path "${organism_name}/offtarget/microbiomes/", emit: microbiomes_dir
    path "${organism_name}/offtarget/microbiomes/*/species_blast_results/*_offtarget_norm.tsv", emit: normalized_tables
    path "${organism_name}/offtarget/microbiomes/*/species_blast_results/*_offtarget_counts.tsv", emit: counts_tables
    path "${organism_name}/offtarget/microbiomes/*/species_blast_results/*_genomes_analyzed.tsv", emit: genomes_analyzed_tables
    path "${organism_name}/offtarget/**", emit: all_microbiome_offtarget
    val organism_name, emit: organism_name
    
    script:
    def base_path = workflow.projectDir.parent
    """#!/usr/bin/env python3
    
import sys
import os
import json

# Add parent directory to path to import ftscripts
sys.path.insert(0, '${base_path}')

from ftscripts import offtargets

print('=' * 80)
print('MICROBIOME OFFTARGET ANALYSIS'.center(80))
print('=' * 80)

catalogues = json.loads('''${microbiome_catalogues_json}''')

print('\\nCatalogues:')
for catalogue in catalogues:
    print(
        f"  - {catalogue['name']}: identity={catalogue['identity_filter']}%, "
        f"coverage={catalogue['coverage_filter']}%"
    )
print(f'  - CPUs: ${cpus}')

# Create organism directory structure in work dir
work_dir = os.getcwd()
organism_dir = os.path.join(work_dir, '${organism_name}')
offtarget_dir = os.path.join(organism_dir, 'offtarget')
os.makedirs(offtarget_dir, exist_ok=True)

# Create genome directory and copy genome files
import shutil
genome_dir = os.path.join(organism_dir, 'genome')
os.makedirs(genome_dir, exist_ok=True)

print('Copying genome files...')
for genome_file in os.listdir('.'):
    if genome_file.endswith(('.gbk', '.faa', '.fna', '.fasta', '.gff')):
        src = os.path.join(work_dir, genome_file)
        dst = os.path.join(genome_dir, genome_file)
        if os.path.isfile(src) and not os.path.exists(dst):
            shutil.copy2(src, dst)
            print(f'  Copied: {genome_file}')

for catalogue in catalogues:
    name = catalogue['name']
    identity = float(catalogue['identity_filter'])
    coverage = float(catalogue['coverage_filter'])
    print(f'[1] Running DIAMOND searches against {name}...')
    offtargets.microbiome_offtarget_blast_species(
        '${databases_path}',
        work_dir,
        '${organism_name}',
        name,
        identity,
        coverage,
        ${cpus}
    )
    print(f'[2] Parsing {name} results...')
    result_tables = offtargets.microbiome_species_parse(
        '${databases_path}',
        work_dir,
        '${organism_name}',
        name,
        identity,
        coverage
    )
    print(f'  - Genes analyzed: {len(result_tables[0])}')

print('Microbiome offtarget analysis completed')
"""
    
    stub:
    """
    mkdir -p ${organism_name}/offtarget/microbiomes/human-gut/species_blast_results
    
    # Create dummy species results
    echo -e "gene\thuman_gut_offtarget_norm\ngene1\t0.15" > ${organism_name}/offtarget/microbiomes/human-gut/species_blast_results/human_gut_offtarget_norm.tsv
    echo -e "gene\thuman_gut_offtarget_counts\ngene1\t5" > ${organism_name}/offtarget/microbiomes/human-gut/species_blast_results/human_gut_offtarget_counts.tsv
    echo -e "gene\thuman_gut_genomes_analyzed\ngene1\t4744" > ${organism_name}/offtarget/microbiomes/human-gut/species_blast_results/human_gut_genomes_analyzed.tsv
    
    # Create a dummy individual result
    echo "STUB: Microbiome offtarget for ${organism_name}"
    """
}
