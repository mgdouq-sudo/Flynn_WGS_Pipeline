#!/usr/bin/env nextflow

process SAMTOOLS_SORTBYNAME {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/samtoolsbam"

    input:
    path(bam)

    output:
    path("*bam"), emit: sorted_bam

    shell:
    def name = bam.simpleName
    """
    samtools sort -n -@ ${task.cpus} -m 2G -o ${name}.namesorted.bam $bam
    """
}