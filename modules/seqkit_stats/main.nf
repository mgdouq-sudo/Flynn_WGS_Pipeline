#!/usr/bin/env nextflow

process SEQKIT_STATS {

    label 'process_low'
    container 'ghcr.io/sydneysorbello/seqkit:latest'
    publishDir "${params.outdir}/seqkit", mode: 'copy'

    input:
    tuple path(read1), path(read2)

    output:
    path("*txt"), emit: seqstat

    shell:
    def name = read1.simpleName
    """
    seqkit stats -j 4 $read1 $read2 > ${name}.txt
    """
    
    }