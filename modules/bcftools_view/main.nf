#!/usr/bin/env nextflow

process BCFTOOLS_VIEW {

    container 'ghcr.io/sydneysorbello/bcftools:latest'
    label "process_high"
    publishDir "${params.outdir}/variant_subset"

    // takes in a zipped vcf and index as input 
    input:
    tuple path(vcf), path(index)

    // outputs a text file with calls within the region of interest
    output:
    path("*txt"), emit: chr17_sv_txt

    // use bcftools to subset the vcf file
    script:
    def name = vcf.simpleName
    """
    bcftools view -r chr17:7668421-7703502 $vcf | grep -v '^#' > ${name}.tp53_tcab1_sv.txt || true
    """

}