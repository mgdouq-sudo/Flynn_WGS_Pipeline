#!/usr/bin/env nextflow

process GATK_COLLECT_READ_COUNTS_TUMOR {
  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_10kb/read_counts", mode: 'copy'

  // input:
  // tumor bam files and corresponding index file 
  // human reference file
  // bwa index files
  // genome file dictionary
  // binned windows
  input:
  tuple val(sample), path(bam), path(bai)
  path ref_fasta
  path ref_index
  path ref_dict
  path intervals

  // counts in an hdf5 file format
  output:
  tuple val(sample), path("*.counts.hdf5"), emit: counts_hdf5

  // use gatk to collect tumor read counts
  script:
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

process DENOISE_READ_COUNTS {
  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_10kb/denoise", mode: 'copy'

  // input: 
  // tumor counts in hdf5 format
  // panel of normals
  input:
  tuple val(sample), path(counts_hdf5_file)
  path pon

  output:
  tuple val(sample), path("*standardizedCR.tsv"), emit: standard
  tuple val(sample), path("*denoisedCR.tsv"), emit: denoised

  script:
  """
  set -euo pipefail

  gatk DenoiseReadCounts \
  -I ${counts_hdf5_file} \
  --count-panel-of-normals ${pon} \
  --standardized-copy-ratios ${sample}.standardizedCR.tsv \
  --denoised-copy-ratios ${sample}.denoisedCR.tsv
  """
}

  process COLLECT_ALLELIC_COUNTS_TUMOR {
  label 'process_high'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_10kb/allelic_counts", mode: 'copy'

  input:
  tuple val(sample), path(bam), path(bai)
  path ref
  path ref_index
  path ref_dict
  path snps
  path snps_index

  output:
  tuple val(sample), path("*allelicCounts.tsv"), emit: allelic_counts

  script:
  """
  gatk CollectAllelicCounts \
  -I ${bam} -R ${ref} \
  -L ${snps} -O ${sample}.allelicCounts.tsv
  """
}

process GATK_MODEL_SEGMENTS_TUMOR_ONLY {
    label 'process_high'
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/gatk_10kb/model", mode: 'copy'

    input:
    tuple val(sample), path(denoised), path(allelic_counts)

    output:
    tuple val(sample), path("model_segments_out"), emit: model_dir
    tuple val(sample), path("*cr.seg"), emit: cr_seg
    tuple val(sample), path("*modelFinal.seg"), emit: model_segments_out

    script:
    def ac_arg = allelic_counts ? "--allelic-counts ${allelic_counts}" : ""
    """
    set -euo pipefail
    mkdir -p model_segments_out

    gatk ModelSegments \
      --denoised-copy-ratios ${denoised} \
      ${ac_arg} \
      --output-prefix ${sample} \
      -O model_segments_out

    # standard expected outputs for downstream steps
    ln -s model_segments_out/${sample}.cr.seg ${sample}.cr.seg
    ln -s model_segments_out/${sample}.modelFinal.seg ${sample}.modelFinal.seg
    """
}

process GATK_CALL_COPYRATIO_SEGMENTS {
    label 'process_light'
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/gatk_10kb/called", mode: 'copy'

    input:
    tuple val(sample), path(cr_seg)
    
    output:
    tuple val(sample), path("*called.seg"), emit: called_seg

    script:
    """
    set -euo pipefail

    gatk CallCopyRatioSegments \
     -I ${cr_seg} \
     -O ${sample}.called.seg
    """
}

process GATK_PLOT_DENOISED_COPYRATIOS {
    label 'process_light'
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/gatk_10kb/plots_copyratios", mode: 'copy'

    input:
    tuple val(sample), path(standard)
    tuple val(sample), path(denoised)
    path ref_dict

    output:
    tuple val(sample), path("*"), emit: copyratio_plots

  script:
  """
  set -euo pipefail
  gatk PlotDenoisedCopyRatios \
    --standardized-copy-ratios ${standard} \
    --denoised-copy-ratios ${denoised} \
    --sequence-dictionary ${ref_dict} \
    -O . --output-prefix ${sample}
  """
}

process GATK_PLOT_MODELED_SEGMENTS {
  label 'process_light'
  container 'broadinstitute/gatk:latest'
  publishDir "${params.outdir}/gatk_10kb/plots_segments", mode: 'copy'

  input:
  tuple val(sample), path(denoised), path(model_seg)
  path ref_dict

  output:
  tuple val(sample), path("*"), emit: segment_plots

  script:
  """
  set -euo pipefail
  gatk PlotModeledSegments \
    --denoised-copy-ratios ${denoised} \
    --segments ${model_seg} \
    --sequence-dictionary ${ref_dict} \
    -O . --output-prefix ${sample}
  """
}


