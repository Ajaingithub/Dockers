#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplefile = '/diazlab/data3/.abhinav/tools/docker/Dockers/exome_seq_docker/sample.csv'
params.reference = '/diazlab/data3/.abhinav/tools/docker/Dockers/exome_seq_docker/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.index_files = '/diazlab/data3/.abhinav/tools/docker/Dockers/exome_seq_docker/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.*'

// Quality Check
process QC {
    publishDir "fastQC", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(fastq1), path(fastq2),path(reference),path(index_files)

    output:
    tuple val(sample_id), path("*.html"),path(reference),path(index_files)

    script:
    """
    fastqc ${fastq1} ${fastq2} -o ./
    """
}

// Removing bad quality Reads
process TRIM {
    publishDir "trimmed", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(fastq1), path(fastq2),path(reference),path(index_files)

    output:
    tuple val(sample_id), path("trim*.fastq"),path(reference),path(index_files)

    script:
    """
    java -jar /usr/local/bin/trimmomatic-0.39.jar PE \
    -threads 4 \
    ${fastq1} ${fastq2} \
    trim_${sample_id}_1_10000.fastq unpaired_${sample_id}_1_10000.fastq \
    trim_${sample_id}_2_10000.fastq unpaired_${sample_id}_2_10000.fastq \
    AVGQUAL:20 \
    SLIDINGWINDOW:5:20 \
    MINLEN:50

    fastqc trim_${sample_id}_1_10000.fastq trim_${sample_id}_2_10000.fastq
    """
}

// Alignment to the reference genome
// before aligning this one should have index the reference genome using bwa index <ref genome>
process ALIGN {
    publishDir "alignment", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(trimmed),path(reference), path(index_files)

    output:
    tuple val(sample_id), path("*bwa_sorted.bam")

    script:
    """
    bwa mem -t 10 ${reference} trim_${sample_id}_2_10000.fastq trim_${sample_id}_2_10000.fastq > ${sample_id}_bwa.sam
    samtools view -S -b ${sample_id}_bwa.sam > ${sample_id}_bwa.bam
    samtools sort ${sample_id}_bwa.bam ${sample_id}_bwa_sorted.bam
    """
}

process DUP_REMOVE {
    publishDir "alignment", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(aligned)

    output:
    tuple val(sample_id), path("*.ba*")

    script:
    """
    
    samtools index ${sample_id}_bwa_sorted.bam
    """
}

workflow {
    samples = Channel
        .fromPath(params.samplefile, type: 'file')
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample_id, file(row.fastq1), file(row.fastq2),file(row.reference),file(row.index_files))}
        
        processed = QC(samples)
        trimmed = TRIM(samples)
        trimmed.view()
        aligned = ALIGN(trimmed)
}
