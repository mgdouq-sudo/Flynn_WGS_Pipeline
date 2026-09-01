#!/usr/bin/env nextflow

process SURVIVOR_SVINTEGRATION {

    container 'ghcr.io/sydneysorbello/survivor:latest'
    label "process_high"
    publishDir "${params.outdir}/survivor"

    input:
    tuple val(sample), path(manta_vcf), path(delly_vcf), path(svaba_vcf)
    path(svaba_out)

    output:
    tuple val(sample), path("${sample}.SURVIVOR.vcf"), emit: survivor_vcf

    shell:
    """
    # Create a list of VCFs
    echo "${manta_vcf}" > ${sample}.vcf_list.txt
    echo "${delly_vcf}" >> ${sample}.vcf_list.txt
    echo "${svaba_vcf}" >> ${sample}.vcf_list.txt

    # Merge using SURVIVOR
    SURVIVOR merge ${sample}.vcf_list.txt 50 2 1 1 0 50 ${sample}.SURVIVOR.vcf
    """

}

process SURVIVOR_INSINTEGRATION {
  
    container 'ghcr.io/sydneysorbello/survivor:latest'
    label "process_high"
    publishDir "${params.outdir}/survivor"

    input:
    tuple val(sample), path(manta_vcf), path(delly_vcf), path(svaba_vcf)
    path(svaba_out)

    output:
    tuple val(sample), path("${sample}.INS.SURVIVOR.vcf"), emit: survivor_vcf

    shell:
    """
    # Create a list of VCFs
    echo "${manta_vcf}" > ${sample}.vcf_list.txt
    echo "${delly_vcf}" >> ${sample}.vcf_list.txt
    echo "${svaba_vcf}" >> ${sample}.vcf_list.txt

    # Merge using SURVIVOR
    SURVIVOR merge ${sample}.vcf_list.txt 1000 1 1 0 0 1 ${sample}.INS.SURVIVOR.vcf
    """
  
}