#!/usr/bin/env nextflow

process GATK_PREPROCESS_INTERVALS {
  label 'process_medium'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_10kb/intervals", mode: 'copy'

  input:
  path ref_fasta
  path ref_dict
  path ref_fai
  val contigs

  // Typical WGS binning; tweak sizes as needed
  // For WES, provide a padded targets list via --intervals instead.
  output:
  path "intervals.preprocessed.interval_list", emit: intervals

  script:
  def bin = params.bin_size ?: 10000
  def Ls = contigs.collect { "-L ${it}" }.join(' ')
  """
  gatk PreprocessIntervals \
    -R ${ref_fasta} \
    --bin-length ${bin} \
    --padding 250 \
    --interval-merging-rule OVERLAPPING_ONLY \
    ${Ls} \
    -O intervals.preprocessed.interval_list
  """
}

process GATK_ANNOTATE_INTERVALS {
  label 'process_medium'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_10kb/intervals", mode: 'copy'

  input:
  path ref_fasta
  path intervals

  output:
  path "intervals.annotated.tsv", emit: annotated

  script:
  """
  gatk AnnotateIntervals \
    -R ${ref_fasta} \
    -L ${intervals} \
    -O intervals.annotated.tsv
  """
}

process GATK_COLLECT_READ_COUNTS {
  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_10kb/read_counts", mode: 'copy'

  input:
  tuple path(bam), path(bai)
  path ref_fasta
  path ref_index
  path ref_dict
  path intervals

  output:
  path "*.counts.hdf5", emit: counts_hdf5

  script:
  def sample = bam.baseName.replaceAll(/\.bam$/,'')
  """
  gatk CollectReadCounts \
    -I ${bam} \
    -R ${ref_fasta} \
    -L ${intervals} \
    --format HDF5 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${sample}.counts.hdf5
  """
}

process GATK_CREATE_PON {
  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_10kb/pon", mode: 'copy'

  input:
  path counts_hdf5_files  // a file collection (flattened) from all normals

  output:
  path "*hdf5", emit: pon

  script:
  // Build the repeated --input string
  def inputs = counts_hdf5_files.collect { "--input ${it}" }.join(' ')
  """
  gatk CreateReadCountPanelOfNormals \
    ${inputs} \
    --minimum-interval-median-percentile 10.0 \
    --maximum-zeros-in-sample-percentage 5.0 \
    --maximum-zeros-in-interval-percentage 5.0 \
    --extreme-sample-median-percentile 2.5 \
    --do-impute-zeros true \
    --extreme-outlier-truncation-percentile 0.1 \
    --number-of-eigensamples 20 \
    --maximum-chunk-size 16777216 \
    --output somatic_cnv.pon.hdf5
  """
}
