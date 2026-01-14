#!/bin/bash

singularity pull docker://sangerpathogens/roary
singularity pull docker://fpocket/fpocket
singularity pull docker://brinkmanlab/psortb_commandline:1.0.2
singularity pull docker://mcpalumbo/corecruncher:1
singularity pull docker://mcpalumbo/bioperl:1
singularity pull docker://mcpalumbo/foldseek:1
singularity pull docker://mcpalumbo/p2rank:latest
singularity pull docker://mcpalumbo/metagraphtools:latest

echo "Singularity containers have been pulled successfully."