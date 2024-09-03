#!/bin/bash

# Install Docker 
docker build -f ftscripts/Dockerfile-bioperl -t docker-bioperl .

# Pull Docker images
docker pull sangerpathogens/roary
docker pull ezequieljsosa/fpocket
docker pull brinkmanlab/psortb_commandline:1.0.2

echo "Docker images setup completed."