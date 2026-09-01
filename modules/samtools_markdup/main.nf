#!/usr/bin/env nextflow

process SAMTOOLS_MARKDUP {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/samtoolsbam"

    input:
    path(sorted_bam)

    output:
    path("*markdup.bam"), emit: markdup_bam

    shell:
    def base = sorted_bam.simpleName
    """
    samtools markdup -@ ${task.cpus} $sorted_bam ${base}.markdup.bam
    """
}