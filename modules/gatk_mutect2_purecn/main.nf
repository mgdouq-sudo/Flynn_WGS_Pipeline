#!/usr/bin/env nextflow

process GATK_MUTECT2_PURECN {

  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_mutect2_purecn", mode: 'copy'

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
  tuple val(name), path("*.mutect2.purecn.vcf.gz"), emit: mutect2_vcf
  tuple val(name), path("*.mutect2.purecn.vcf.gz.tbi"), emit: mutect2_tbi
  tuple val(name), path("*.mutect2.purecn.vcf.gz.stats"), emit: mutect2_stats

  script:
  """
  gatk Mutect2 \
   -R ${ref} \
   -I ${bam} \
   --germline-resource ${gnomad_ref} \
   --panel-of-normals ${pon_ref} \
   --genotype-germline-sites true \
   --genotype-pon-sites true \
   -O ${name}.mutect2.purecn.vcf.gz
  """

}