Developing Docker for

1. exome sequencing

   1. FASTQC
   2. GATK
   3. PICARD

2. Gene Regulatory Network (SCIENICPlus): it is install in the Conda environment
   We have made it for AMD64 and directly pushed to the Docker Hub.
   You can use

### Docker

docker pull abhinavjain1993/wgs
docker pull abhinavjain1993/scenicplus_docker:amd64

# For running docker in the arm64

docker pull abhinavjain1993/scenicplus_docker:latest

# Running the scenicplus, activate conda using

source activate scenicplus

### Singularity

singularity pull docker://abhinavjain1993/wgs
singularity pull docker://abhinavjain1993/scenicplus_docker:amd64
singularity exec –bind <Mounted_directory>:/data --contain-all scenicplus_docker_latest.sif /bin/bash
