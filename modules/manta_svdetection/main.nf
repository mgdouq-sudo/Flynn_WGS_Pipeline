#!/usr/bin/env nextflow

process MANTA_SVDETECTION {

    container 'ghcr.io/sydneysorbello/manta:latest'
    label "process_high"
    publishDir "${params.outdir}/manta", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path(ref)
    path(ref_fai)

    output:
    path("*diploidSV.vcf.gz"), emit: dipvcf_gz
    path("*candidateSV.vcf.gz"), emit: canvcf
    path("*diploidSV.vcf"), emit: dipvcf

    shell:
    def name = bam.simpleName
    """
    configManta.py \
    --bam $bam \
    --referenceFasta $ref \
    --runDir manta_run
    
    manta_run/runWorkflow.py -j ${task.cpus}

    cp manta_run/results/variants/diploidSV.vcf.gz ${name}_target.diploidSV.vcf.gz
    cp manta_run/results/variants/candidateSV.vcf.gz ${name}_target.candidateSV.vcf.gz

    gunzip -c ${name}_target.diploidSV.vcf.gz > ${name}_target.diploidSV.vcf
    """

}