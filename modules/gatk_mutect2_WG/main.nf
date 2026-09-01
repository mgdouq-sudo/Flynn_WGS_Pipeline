#!/usr/bin/env nextflow

process GATK_MUTECT2 {

  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_mutect2_WG", mode: 'copy'

  input:
  tuple val(name), path(bam), path(bai)
  path(ref)
  path(ref_index)
  path(ref_dict)
  path(gnomad_ref)
  path(gnomad_ref_index)
  path(pon_ref)
  path(pon_index)

  output:
  tuple val(name), path("*.mutect2.vcf.gz"), emit: mutect2_vcf
  tuple val(name), path("*.mutect2.vcf.gz.tbi"), emit: mutect2_tbi
  tuple val(name), path("*mutect2.vcf.gz.stats"), emit: mutect2_stats

  script:
  """
  gatk Mutect2 \
   -R ${ref} \
   -I ${bam} \
   --germline-resource ${gnomad_ref} \
   --panel-of-normals ${pon_ref} \
   -O ${name}.mutect2.vcf.gz
  """

}

process GATK_FILTERMUTECTCALLS {

  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_mutect2_filtered_WG", mode: 'copy'

  input:
  tuple val(sample), path(mutect2_vcf)
  tuple val(sample), path(mutect2_tbi)
  tuple val(sample), path(mutect2_stats)
  path(ref)
  path(ref_index)
  path(ref_dict)

  output:
  tuple val(sample), path("*filtered.vcf.gz"), emit: filtered_mutect_vcfgz
  tuple val(sample), path("*filtered.vcf.gz.tbi"), emit: filtered_mutect_tbi
  tuple val(sample), path("*filtered.vcf.gz.filteringStats.tsv"), emit: filtered_mutect_tsv
  tuple val(sample), path("*filtered.vcf"), emit: filtered_mutect_vcf

  script:
  """
  gatk FilterMutectCalls \
   -V ${mutect2_vcf} \
   -R ${ref} \
   -O ${sample}.mutect2.filtered.vcf.gz

  gunzip -k ${sample}.mutect2.filtered.vcf.gz
  """
}

process GATK_VARIANTANNOTATION {

  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_mutect2_variantannotaiton_WG", mode: 'copy'

  input:
  tuple val(sample), path(vcf), path(bam), path(bai)
  path(dbsnp)
  path(dbsnp_idx)
  path(ref)
  path(ref_index)
  path(ref_dict)

  output:
  tuple val(sample), path("*mutect2.variantannotation.vcf"), emit: variantannotated_vcf

  script:
  """
  gatk VariantAnnotator \
   -R ${ref} \
   -I ${bam} \
   -V ${vcf} \
   --output ${sample}.mutect2.variantannotation.vcf \
   --dbsnp ${dbsnp}
  """
}
