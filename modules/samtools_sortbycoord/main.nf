#!/usr/bin/env nextflow

process SAMTOOLS_SORTBYCOORD {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/samtoolsbam"

    input:
    path(fixmate_bam)

    output:
    path("*bam"), emit: sorted_bam

    shell:
    def name = fixmate_bam.simpleName
    """
    samtools sort -@ ${task.cpus} -m 2G -o ${name}.coordsorted.bam $fixmate_bam
    """
}