#!/usr/bin/env nextflow

process XENGSORT_CLASSIFY {

    conda 'envs/xengsort_env.yml'
    label "process_high"
    publishDir "${params.outdir}/xengsort"

    input:
    path(hash)
    path(info)
    tuple path(read1), path(read2)

    output:
    tuple path("*combined.1.fq.gz"), path("*combined.2.fq.gz"), emit: combined_pair
    tuple path("*graft.1.fq.gz"), path("*graft.2.fq.gz"), emit: graft_pair
    tuple path("*both.1.fq.gz"), path("*both.2.fq.gz"), emit: both_pair

    script:
    def base = read1.simpleName
    """
    xengsort classify \
        --index xengsort \
        --fastq $read1 \
        --pairs $read2 \
        --prefix ${base} \
        --threads ${task.cpus}

    # Concatenate graft + both reads for R1 and R2
    cat ${base}-graft.1.fq.gz ${base}-both.1.fq.gz ${base}-ambiguous.1.fq.gz > ${base}-combined.1.fq.gz
    cat ${base}-graft.2.fq.gz ${base}-both.2.fq.gz ${base}-ambiguous.2.fq.gz > ${base}-combined.2.fq.gz
    """
}
