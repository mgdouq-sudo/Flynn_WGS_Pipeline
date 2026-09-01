#!/usr/bin/env nextflow

process BWA_INDEX {

    container 'ghcr.io/sydneysorbello/bwa:latest'
    label "process_high"
    publishDir params.refdir

    // input comes from human reference genome
    input:
    path(ref)

    // outputs bwa index files for the reference file
    output:
    tuple path("*amb"), path("*ann"), path("*bwt"), path("*pac"), path("*sa"), emit: bwa_index

    // use bwa to index the human reference genome
    script:
    def name = ref.simpleName
    """
    bwa index -p ${name}.fasta $ref
    """

}