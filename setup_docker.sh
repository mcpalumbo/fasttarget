#!/bin/bash

# Pull Docker images
docker pull sangerpathogens/roary
docker pull ezequieljsosa/fpocket
docker pull brinkmanlab/psortb_commandline:1.0.2
docker pull mcpalumbo/corecruncher:1
docker pull mcpalumbo/bioperl:1

echo "Docker images setup completed."