#!/usr/bin/env nextflow

process SAMTOOLS_INDEXSUBSET {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/bam_subset"

    input:
    path(file)

    output:
    path("*bai"), emit: bam_index

    shell:
    """
    samtools index $file
    """
}