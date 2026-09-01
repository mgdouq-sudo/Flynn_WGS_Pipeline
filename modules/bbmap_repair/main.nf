#!/usr/bin/env nextflow

process BBMAP_REPAIR {

    // requires high computational allocation
    container 'ghcr.io/sydneysorbello/bbmap:latest'
    label "process_bbmap"
    publishDir "${params.outdir}/bbmap"
    memory '512 GB'

    // takes in paired end fastq files
    input:
    tuple path(read1), path(read2)

    // outputs paired end fastq files containing repaired reads
    // outputs singleton fastq file with reads that are unpaired
    output:
    tuple path("*r1.repaired.fq.gz"), path("*r2.repaired.fq.gz"), emit: bbmap_pair
    path("*.singleton.fq.gz"), emit: singleton

    // call the repair.sh script for data
    script:
    def prefix = read1.simpleName
    """
    repair.sh \\
        -Xmx400G \\
        in1=$read1 \\
        in2=$read2 \\
        out1=${prefix}.r1.repaired.fq.gz \\
        out2=${prefix}.r2.repaired.fq.gz \\
        outs=${prefix}.singleton.fq.gz \\
        overwrite=t \\
        threads=${task.cpus}
    """

}