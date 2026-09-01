#!/usr/bin/env nextflow

process SVABA_SVDETECTION {

    container 'ghcr.io/sydneysorbello/svaba:latest'
    label "process_high"
    publishDir "${params.outdir}/svaba"

    input:
    tuple path(bam), path(bai)
    path(ref)
    path(ref_index)
    tuple path(amb), path(ann), path(bwt), path(pac), path(sa)
    path(dict)

    output:
    path("*svaba.sv.vcf"), emit: svaba_sv
    path("*svaba.indel.vcf"), emit: svaba_indel

    script:
    def prefix = bam.simpleName
    """
    svaba run \\
        -p ${task.cpus} \\
        -a ${prefix} \\
        -G ${ref} \\
        -t ${bam} \\
        --override-reference-check
    """
}
