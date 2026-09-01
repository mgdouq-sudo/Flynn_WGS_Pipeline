#!/usr/bin/env nextflow

process TABIX_INDEXVCF {

    conda 'envs/tabix_env.yml'
    label "process_high"
    publishDir "${params.outdir}/survivor"

    input:
    path(vcf)

    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_index

    shell:
    def name = vcf.simpleName
    """
    bcftools sort $vcf -Oz -o ${name}.vcf.gz
    tabix -p vcf ${name}.vcf.gz
    """
}

process TABIX_INDEXDELLYVCF {

    conda 'envs/tabix_env.yml'
    label "process_high"
    publishDir "${params.outdir}/tabix"

    input:
    path(vcf)

    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_index

    shell:
    def name = vcf.simpleName
    """
    bcftools sort $vcf -Oz -o ${name}.vcf.gz
    tabix -p vcf ${name}.vcf.gz
    """
}

process TABIX_INDEXSVABAVCF {

    conda 'envs/tabix_env.yml'
    label "process_high"
    publishDir "${params.outdir}/tabix"

    input:
    path(vcf)

    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_index

    shell:
    def name = vcf.simpleName
    """
    bcftools sort $vcf -Oz -o ${name}.svaba.vcf.gz
    tabix -p vcf ${name}.svaba.vcf.gz
    """
}

