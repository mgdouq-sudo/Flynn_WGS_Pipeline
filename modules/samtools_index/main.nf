#!/usr/bin/env nextflow

process SAMTOOLS_INDEX {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/samtoolsbam"

    input:
    path(file)

    output:
    tuple path(file), path("*bai"), emit: bam_index

    shell:
    """
    samtools index $file
    """
}