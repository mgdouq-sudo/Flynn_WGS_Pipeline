#!/usr/bin/env nextflow

process BEDTOOLS_REMOVE_REPEATS {

    container 'ghcr.io/sydneysorbello/bedtools:latest'
    label "process_high"
    publishDir "${params.outdir}/vcf_filtered_repeatmasker"

    // input: 
    // unzipped sample vcf file
    // unzipped common structural variant file
    input:
    tuple val(sample), path(my_vcf)
    path(repeat_regions)

    // output filtered vcf file unzipped
    output:
    tuple val(sample), path("*.ALLSV.repeatmasker.vcf"), emit: repetitive_filtered_vcf

    // use bedtools to remove the common SVs from the tumor samples
    script:
    """
    bedtools intersect -a ${my_vcf} -b ${repeat_regions} -v -f 0.5 -header > ${sample}.ALLSV.repeatmasker.vcf
    """

}