#!/usr/bin/env nextflow

process ANNOVAR {

    label "process_medium"
    publishDir "${params.outdir}/annovar_output", mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path(humandb)
    path(table_annovar)
    path(convert2annovar)
    path(annotate_variant)
    path(coding_change)
    path(retrieve_seq_from_fasta)
    path(variants_reduction)
    path(varannot_out)

    output:
    path("*.txt"), emit: annovar_txt
    path("*vcf"), emit: annovar_vcf

    script:
    def sample = vcf.simpleName
    """
    perl ${table_annovar} ${vcf} \
    ${humandb}/ \
    -buildver hg38 \
    -out ${sample}_annovar \
    -protocol refGeneWithVer,clinvar_20250721,dbnsfp47a \
    -operation g,f,f \
    -remove -polish -vcfinput -nastring .
    """
}