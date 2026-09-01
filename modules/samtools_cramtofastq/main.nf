#!/usr/bin/env nextflow

process SAMTOOLS_CRAMTOFASTQ {

    container 'ghcr.io/sydneysorbello/samtools:latest'
    label "process_high"
    publishDir "${params.outdir}/fastq"

    input:
    tuple val(sample), path(cram)
    path(ref)

    output:
    tuple path("*out.R1.fastq.gz"), path("*out.R2.fastq.gz"), emit: fastq_pair

    shell:
    def name = cram.simpleName
    """
    samtools fastq -@ ${task.cpus} --reference $ref -1 ${name}.out.R1.fastq.gz -2 ${name}.out.R2.fastq.gz $cram
    """
    
    }