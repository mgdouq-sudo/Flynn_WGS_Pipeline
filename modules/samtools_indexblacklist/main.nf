#!/usr/bin/env nextflow

process SAMTOOLS_INDEXBLACKLIST {

    conda '/restricted/projectnb/flynngrp/WGS_25/Sorbello_analysis/Flynn_WGS_Analysis/envs/samtools_env.yml'
    label "process_high"
    publishDir params.outdir

    input:
    path(file)

    output:
    path("*bai"), emit: bam_index

    shell:
    """
    samtools index $file
    """
}