#!/usr/bin/env nextflow

process TRIMMOMATIC_PE {

    container 'ghcr.io/sydneysorbello/trimmomatic:latest'
    label "process_high"
    publishDir "${params.outdir}/trimmomatic"

    input:
    tuple path(read1), path(read2)

    output:
    tuple path("*_1P.fastq.gz"), path("*_2P.fastq.gz"), emit: readp
    tuple path("*_1U.fastq.gz"), path("*_2U.fastq.gz"), emit: readu

    shell:
    def prefix = read1.simpleName
    """
    trimmomatic PE -threads ${task.cpus} \
      $read1 $read2 \
      ${prefix}_1P.fastq.gz ${prefix}_1U.fastq.gz \
      ${prefix}_2P.fastq.gz ${prefix}_2U.fastq.gz \
      SLIDINGWINDOW:4:15 MINLEN:25
    """
}