#!/usr/bin/env nextflow

process SAMTOOLS_FIXMATE {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/samtoolsbam"

    input:
    path(name_sorted_bam)

    output:
    path("*bam"), emit: fixmate_bam

    shell:
    def name = name_sorted_bam.simpleName
    """
    samtools fixmate -m $name_sorted_bam ${name}.fixmate.bam
    """
}