#!/usr/bin/env nextflow

process PICARD_CREATESEQUENCEDICT {

    conda 'envs/picard_env.yml'
    label "process_high"
    publishDir "${params.refdir}", mode: 'copy'

    input:
    path(ref)

    output:
    path("*.dict"), emit: dict

    script:
    def name = ref.baseName
    """
    picard CreateSequenceDictionary \
      R=$ref \
      O=${name}.dict
    """
}
