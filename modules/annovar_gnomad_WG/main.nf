#!/usr/bin/env nextflow

process ANNOVAR_GNOMAD {

    label "process_medium"
    publishDir "${params.outdir}/annovar_output_WG_gnomad", mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path(humandb)
    path(table_annovar)
    path(convert2annovar)
    path(annotate_variant)
    path(coding_change)
    path(retrieve_seq_from_fasta)
    path(variants_reduction)

    output:
    tuple val(sample), path("*.txt"), emit: annovar_txt
    tuple val(sample), path("*vcf"), emit: annovar_vcf

    script:
    """
    perl ${table_annovar} ${vcf} \
    ${humandb}/ \
    -buildver hg38 \
    -out ${sample}_annovar \
    -protocol refGeneWithVer,clinvar_20250721,dbnsfp47a,gnomad41_genome \
    -operation g,f,f,f \
    -remove -polish -vcfinput -nastring .
    """
}