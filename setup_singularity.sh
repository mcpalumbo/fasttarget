#!/bin/bash

# Create directory for Singularity .sif files
mkdir -p singularity_sfi_files

# Pull containers with specific names matching what the code expects
singularity pull --name singularity_sfi_files/sangerpathogens_roary.sif docker://sangerpathogens/roary
singularity pull --name singularity_sfi_files/fpocket_fpocket.sif docker://fpocket/fpocket
singularity pull --name singularity_sfi_files/brinkmanlab_psortb_commandline_1.0.2.sif docker://brinkmanlab/psortb_commandline:1.0.2
singularity pull --name singularity_sfi_files/mcpalumbo_corecruncher_1.sif docker://mcpalumbo/corecruncher:1
singularity pull --name singularity_sfi_files/mcpalumbo_bioperl_1.sif docker://mcpalumbo/bioperl:1
singularity pull --name singularity_sfi_files/mcpalumbo_foldseek_1.sif docker://mcpalumbo/foldseek:1
singularity pull --name singularity_sfi_files/mcpalumbo_p2rank_latest.sif docker://mcpalumbo/p2rank:latest
singularity pull --name singularity_sfi_files/mcpalumbo_metagraphtools_latest.sif docker://mcpalumbo/metagraphtools:latest

echo "Singularity containers have been pulled successfully to singularity_sfi_files/"

