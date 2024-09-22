# DIsha is the best
FROM ubuntu:20.04


ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        bwa \
        samtools \
        wget \
        unzip \
        python \
        openjdk-17-jdk
        
RUN mkdir WGStoolsDir
WORKDIR WGStoolsDir

RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xvzf sratoolkit.tar.gz && rm sratoolkit.tar.gz

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    rm fastqc_v0.12.1.zip
RUN wget -O picard.jar "https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar" 

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && \
    unzip gatk-4.4.0.0.zip && rm gatk-4.4.0.0.zip

RUN wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz && \
    tar -xvzf annovar.latest.tar.gz && \
    rm annovar.latest.tar.gz

# You can specify a default command or entrypoint if needed
CMD ["bash"]
