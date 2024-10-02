#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplefile = '/diazlab/data3/.abhinav/tools/docker/Dockers/exome_seq_docker/sample.csv'

// Performing he fastqc using singularity
process QC {
    publishDir "fastQC", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("*.html")

    script:
    """
    fastqc ${fastq1} ${fastq2} -o ./
    """
}

workflow {
    samples = Channel
        .fromPath(params.samplefile, type: 'file')
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample_id, file(row.fastq1), file(row.fastq2)) }
        
        processed = QC(samples)
}
