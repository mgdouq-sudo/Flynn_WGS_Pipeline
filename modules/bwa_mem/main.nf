#!/usr/bin/env nextflow

process BWA_MEM {

    // this process requires high computational allocation
    conda 'envs/bwa_env.yml'
    label "process_superhigh"
    publishDir "${params.outdir}/bwa"
    memory '1 TB'

    // inputs:
    // human genome reference file
    // human genome reference file bwa index files
    // first paired end fastq file
    // second paired end fastq file
    input:
    path(ref)
    tuple path(amb), path(ann), path(bwt), path(pac), path(sa)
    tuple path(R1), path(R2)

    // output is a single bam file with the aligned reads
    output:
    path("*bam"), emit: bam

    // align the paired end fastq files to the human reference genome
    script:
    def name = R1.simpleName
    """
    bwa mem -t ${task.cpus} $ref $R1 $R2 | samtools sort -@ ${task.cpus} -o ${name}.bam
    """

}