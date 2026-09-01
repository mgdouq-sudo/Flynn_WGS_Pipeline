process FILTER_COMMON_SVS {

  conda 'envs/r_env.yml'
  label 'process_medium'
  publishDir "${params.outdir}/survivor_concat", mode: 'copy'

  input:
  tuple val(sample), path(annotsv_tsv)
  path(common_bed)
  path(script)

  output:
  tuple val(sample), path("*.overlapping_SV_calls.bed"), emit: overlap
  tuple val(sample), path("*.tumor_only_SV_calls.bed"), emit: tumor_only

  script:
  """
  Rscript ${script} \
    --annotsv=${annotsv_tsv} \
    --common=${common_bed} \
    --outdir=. \
    --tolerance=50
  """
}
