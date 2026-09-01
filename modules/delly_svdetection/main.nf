#!/usr/bin/env nextflow

process DELLY_SVDETECTION {

    container 'ghcr.io/sydneysorbello/delly:latest'
    label "process_high"
    publishDir "${params.outdir}/delly"

    // input:
    // tumor bam file and corresponding index file
    // human reference genome file
    // human reference genome index file
    input:
    tuple path(bam), path(bai)
    path(ref)
    path(ref_index)

    // output unzipped vcf file
    output:
    path("*vcf"), emit: delly_vcf

    // run delly on the tumor sample bam file
    script:
    def name = bam.simpleName
    """
    delly call -g $ref $bam > ${name}_target.delly.vcf
    """

}