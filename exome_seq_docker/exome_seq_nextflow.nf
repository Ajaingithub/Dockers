#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplefile = '/diazlab/data3/.abhinav/tools/docker/Dockers/exome_seq_docker/sample.csv'
params.reference = '/diazlab/data3/.abhinav/tools/singularity/data/resources/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.index = '/diazlab/data3/.abhinav/tools/singularity/data/resources/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.*'
params.gatkindex = '/diazlab/data3/.abhinav/tools/singularity/data/resources/reference/Homo_sapiens.GRCh38.dna.primary_assembly.dict'

// Quality Check
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

// Removing bad quality Reads
process TRIM {
    publishDir "trimmed", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("trim*.fastq")

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
process ALIGN {
    publishDir "alignment", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(trimmed), path(index_files)

    output:
    tuple val(sample_id), path("*_bwa_sorted.bam")

    script:
    """
    bwa mem -t 10 ${params.reference} trim_${sample_id}_1_10000.fastq trim_${sample_id}_2_10000.fastq > ${sample_id}_bwa.sam
    samtools view -S -b ${sample_id}_bwa.sam > ${sample_id}_bwa.bam
    samtools sort ${sample_id}_bwa.bam -o ${sample_id}_bwa_sorted.bam
    """
}

process REM_DUP{
    publishDir "alignment", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(aligned)

    output:
    tuple val(sample_id), path("*_sorted_marked.bam"),path("*_sorted_marked.bam.bai") ,path("*_metric_summary")

    script:
    """
    java -jar /usr/local/bin/picard.jar MarkDuplicates I=${aligned} O=${sample_id}_sorted_marked.bam M=${sample_id}_picard_info.txt REMOVE_DUPLICATES=true AS=true
    samtools index ${sample_id}_sorted_marked.bam
    samtools flagstat ${sample_id}_sorted_marked.bam > ${sample_id}_metric_summary
    """
}

process VAR_CALL{
    publishDir "variant", mode: 'copy', pattern: "*", overwrite: true

    input:
    tuple val(sample_id), path(rem_dup), path(gatkref), path(index),path(gatkindex)

    output:
    tuple val(sample_id), path("*.vcf.gz")

    script:
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${gatkref} -I ${sample_id}_sorted_marked.bam -O ${sample_id}.vcf.gz
    """

}

workflow {
    samples = Channel
        .fromPath(params.samplefile, type: 'file')
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample_id, file(row.fastq1), file(row.fastq2)) }

    processed = QC(samples)
    trimmed = TRIM(samples)

    // Passing index files as input to ALIGN
    aligned = trimmed.map { sample_id, trimmed_fastq -> 
        tuple(sample_id, trimmed_fastq, file(params.index))
    }

    bamfile_ch=ALIGN(aligned)
    duprem_ch=REM_DUP(bamfile_ch)
    
    varcall_ch = duprem_ch.map { sample_id, bam, bai, metrics ->
        tuple(sample_id, bam, file(params.reference), file(params.index), file(params.gatkindex)) 
        }

    varcall_ch.view()

    varcall = VAR_CALL(varcall_ch)
}
