#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.singularityImage = "/diazlab/data3/.abhinav/tools/singularity/wgs_amd64.sif"
params.samplefile = '/Users/abjain/Documents/Docker/exome_seq_docker/sample.csv'

// Performing he fastqc using singularity
process QC {
    publishDir "QC", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("*.html")

    // Using the Singularity image
    container params.singularityImage  

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

    processed = PREPROCESSING(samples)

}