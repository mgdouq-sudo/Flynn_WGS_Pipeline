#!/usr/bin/env nextflow

process SVABA_REWRITE_INSERTIONS {
  
  container 'ghcr.io/sydneysorbello/bcftools:latest'
  label "process_medium"
  publishDir "${params.outdir}/svaba_ins_fixed"

  input:
  path(vcf)
  path(script)

  output:
  path("*ins.fixed.vcf"), emit: svaba_ins_fixed

  script:
  def name = vcf.simpleName
  """
  bash ${script} ${vcf} ${name}.ins.fixed.vcf
  """
}