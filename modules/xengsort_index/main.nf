#!/usr/bin/env nextflow

process XENGSORT_INDEX {

    conda 'envs/xengsort_env.yml'
    label "process_high"
    publishDir "${params.outdir}/xengsort", mode: 'copy'

    input:
    path(host_ref)
    path(graft_ref)

    output:
    path("*hash"), emit: hash
    path("*info"), emit: info

    shell:
    """
    xengsort index --index xengsort -H ${host_ref} -G ${graft_ref} -n 4500000000 -k 25 --threads-read ${task.cpus} --threads-split ${task.cpus}
    """

}