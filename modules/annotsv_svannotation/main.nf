#!/usr/bin/env nextflow

process ANNOTSV_SVANNOTATION {

    conda 'envs/annotsv_env.yml'
    label "process_high"
    publishDir "${params.outdir}/annotsv_1Mb", mode: 'copy'

    // single VCF as input file
    input:
    tuple val(sample), path(vcf)

    // outputs an annotated tsv and unannotated tsv
    output:
    tuple val(sample), path("*.annotsv.tsv"), emit: annotsv

    // run AnnotSV on VCF file
    shell:
    def name = vcf.simpleName
    """
    AnnotSV \
    -SVinputFile $vcf \
    -outputFile ${name}.final.annotsv.tsv \
    -annotationsDir ${params.annotsv_dir} \
    -genomeBuild GRCh38 \
    -outputDir . \
    -overwrite 1
    """

}

process ANNOTSV_FORREPEATS {

    conda 'envs/annotsv_env.yml'
    label "process_high"
    publishDir "${params.outdir}/annotsv_withrepeats", mode: 'copy'

    // single VCF as input file
    input:
    tuple val(sample), path(vcf)

    // outputs an annotated tsv and unannotated tsv
    output:
    tuple val(sample), path("*.annotsv.tsv"), emit: annotsv

    // run AnnotSV on VCF file
    shell:
    def name = vcf.simpleName
    """
    AnnotSV \
    -SVinputFile $vcf \
    -outputFile ${name}.final.annotsv.tsv \
    -annotationsDir ${params.annotsv_dir} \
    -genomeBuild GRCh38 \
    -outputDir . \
    -overwrite 1
    """

}