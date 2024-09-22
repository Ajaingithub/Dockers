# Use an official Ubuntu as the base image
# testing new branch
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    wget \
    unzip \
    git \
    ca-certificates \
    bwa \
    samtools \
    python \
    openjdk-17-jdk \
    && rm -rf /var/lib/apt/lists/*

# Install SRA-toolkit
RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xvzf sratoolkit.tar.gz && \
    rm sratoolkit.tar.gz && \
    mv sratoolkit.3.1.1-ubuntu64 /opt/sratoolkit/ && \
    ln -s /opt/sratoolkit/bin/prefetch /usr/local/bin/prefetch && \
    ln -s /opt/sratoolkit/bin/fastq-dump /usr/local/bin/fastq-dump


# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    rm fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    mv FastQC /opt/fastqc && \
    ln -s /opt/fastqc/fastqc /usr/local/bin/fastqc

# Install GATK (version 4.4.0.0)
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && \
    unzip gatk-4.4.0.0.zip && \
    rm gatk-4.4.0.0.zip && \
    mv gatk-4.4.0.0 /opt/gatk && \
    ln -s /opt/gatk/gatk /usr/local/bin/gatk

# Install Picard
RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar && \
    mv picard.jar /opt/picard.jar && \
    ln -s /opt/picard.jar /usr/local/bin/picard.jar

# Install Annovar
# Run wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz && \
#     tar -xvzf annovar.latest.tar.gz && \
#     rm annovar.latest.tar.gz && \
#     mv annovar /opt/ && \
#     ln -s /opt/annovar/annotate_variation.pl /usr/local/bin/annotate_variation.pl \
#     ln -s /opt/annovar/convert2annovar.pl /usr/local/bin/convert2annovar.pl

# Set the working directory
WORKDIR /data

# Default command (optional)
CMD ["bash"]