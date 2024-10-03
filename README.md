# Developing Docker for

## 1. Whole Genome Sequencing

1. FASTQC
2. GATK
3. PICARD
4. BWA
5. ANNOVAR
6. SAMTOOLS

## 2. Gene Regulatory Network (SCIENICPlus):

it is install in the Conda environment
scenicplus:amd64 is being designed for the linux environment

### Docker on x86_64

docker pull abhinavjain1993/wgs:amd64
docker pull abhinavjain1993/scenicplus:amd64

#### For running docker in the arm64

docker pull abhinavjain1993/scenicplus_docker:latest

#### Running the scenicplus, activate conda using

source activate scenicplus

#### Singularity

singularity pull docker://abhinavjain1993/scenicplus_docker:amd64
singularity exec --bind <Mounted_directory>:/data --contain-all scenicplus_docker_latest.sif /bin/bash

##### Running a WGS files on the server

singularity pull docker://abhinavjain1993/wgs:amd64
singularity exec -bind /diazlab/data3/.abhinav/tools/singularity/data/:/data --contain-all wgs_amd64.sif /bin/bash
fasterq-dump --split-3 DRR001913 -O ./ -t ./tmp/
head -10000 DRR001913_2.fastq > DRR001913_2_10000.fastq
head -10000 DRR001913_1.fastq > DRR001913_1_10000.fastq
java -jar /usr/local/bin/trimmomatic-0.39.jar PE -threads 4 DRR001913_1_10000.fastq DRR001913_2_10000.fastq trim_DRR001913_1_10000.fastq unpaired_DRR001913_1_10000.fastq trim_DRR001913_2_10000.fastq unpaired_DRR001913_2_10000.fastq AVGQUAL:20 SLIDINGWINDOW:5:20 MINLEN:50

###### Download the reference genome
<!-- https://hitchhikersguidetoexomeanalysis.wordpress.com/exome-analysis-exercise/ -->
wget https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

### Added the nextflow files so while running one can directly run using nextflow
nextflow run exome_seq_nextflow.nf --with-apptainer /diazlab/data3/.abhinav/tools/singularity/wgs_amd64.sif

Note
While Running the singularity please keep all files in the same folder where you are working. 