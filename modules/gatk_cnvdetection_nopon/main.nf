#!/usr/bin/env nextflow

// ModelSegments without PoN — uses allelic counts only
// PoN denoising was removed per Dr. Leshchiner's recommendation
// because the TARGET-ALL normals were from a different assay type
// and made the copy number profiles worse

process GATK_MODEL_SEGMENTS_NO_PON {
    label 'process_high'
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/gatk_nopon/model", mode: 'copy'

    input:
    tuple val(sample), path(allelic_counts)

    output:
    tuple val(sample), path("${sample}.cr.seg"), emit: cr_seg
    tuple val(sample), path("${sample}.modelFinal.seg"), emit: model_segments_out
    tuple val(sample), path("${sample}.modelFinal.af.param"), emit: af_param

    script:
    """
    set -euo pipefail
    mkdir -p model_segments_out

    gatk ModelSegments \
      --allelic-counts ${allelic_counts} \
      --output-prefix ${sample} \
      -O model_segments_out

    ln -s model_segments_out/${sample}.cr.seg ${sample}.cr.seg
    ln -s model_segments_out/${sample}.modelFinal.seg ${sample}.modelFinal.seg
    ln -s model_segments_out/${sample}.modelFinal.af.param ${sample}.modelFinal.af.param
    """
}