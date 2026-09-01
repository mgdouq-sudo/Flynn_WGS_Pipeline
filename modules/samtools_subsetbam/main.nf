#!/usr/bin/env nextflow

process SAMTOOLS_SUBSETBAM {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/bam_subset"

    input:
    tuple path(bam), path(bai)
    
    output:
    path("*bam"), emit: bam

    shell:
    def name = bam.simpleName
    """
    samtools view -b $bam chr17:7660000-7710000 > ${name}_subset.bam
    """
}