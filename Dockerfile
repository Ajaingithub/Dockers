# Use an official Ubuntu as the base image
# Make some new Changes very new
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    wget \
    unzip \
    default-jre \
    build-essential \\
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    curl \
    git \
    ca-certificates \
    bwa \
    samtools \
    python \
    && rm -rf /var/lib/apt/lists/*

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

# Set the working directory
WORKDIR /data

# Default command (optional)
CMD ["bash"]