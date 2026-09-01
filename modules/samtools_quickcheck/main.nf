#!/usr/bin/env nextflow

process SAMTOOLS_QUICKCHECK {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_low"
    publishDir params.outdir

    input:
    tuple val(sample), path(cram)

    shell:
    """
    samtools quickcheck -v $cram
    """
}