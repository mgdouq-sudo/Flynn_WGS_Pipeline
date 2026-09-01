#!/usr/bin/env nextflow

process SAMTOOLS_FLAGSTAT {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/flagstat"

    input:
    path(prelim_bam)

    output:
    path("*txt"), emit: flagstat

    shell:
    def name = prelim_bam.simpleName
    """
    samtools flagstat $prelim_bam -@ ${task.cpus} > ${name}.txt
    """
}

process SAMTOOLS_COVERAGE {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/samtools_coverage"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path("*txt"), emit: flagstat

    shell:
    """
    samtools coverage ${bam} > ${sample}.txt
    """
}