#!/usr/bin/env nextflow

process PICARD_ADDREADGROUPS {

    container 'broadinstitute/gatk:latest'
    label "process_high"
    publishDir  "${params.outdir}/picard_bam", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("*.rg.bam"), path("*.rg.bam.bai"), emit: bam_rg

    script:
    """
    set -euo pipefail

    gatk AddOrReplaceReadGroups \
      I=${bam} \
      O=${bam.baseName}.rg.bam \
      RGID=${sample} \
      RGSM=${sample} \
      RGPL=ILLUMINA \
      RGLB=lib1 \
      RGPU=${sample}.unit1

    samtools index -@ ${task.cpus} ${bam.baseName}.rg.bam
    """
}
