#!/usr/bin/env nextflow

process BCFTOOLS_FILTERMANTA {

    container 'ghcr.io/sydneysorbello/bcftools:latest'
    label "process_high"
    publishDir "${params.outdir}/vcf_filtered"

    // takes in a zipped vcf and index as input
    input:
    path(vcf)

    // output as an unzipped vcf file
    output:
    path("*filtered.vcf"), emit: filtered_vcf

    // run filter with bcftools
    shell:
    def sample = vcf.simpleName
    """
    bcftools view -i 'FILTER == "PASS"' -Ov -o ${sample}.manta.filtered.vcf ${vcf}
    """

}

process BCFTOOLS_FILTERDELLY {

    conda 'envs/bcftools_env.yml'
    label "process_high"
    publishDir "${params.outdir}/vcf_filtered"

    // takes in a zipped vcf and index as input
    input:
    path(vcf)

    // output as an unzipped vcf file
    output:
    path("*delly.filtered.vcf"), emit: filtered_vcf

    // run filter with bcftools
    shell:
    def sample = vcf.simpleName
    """
    bcftools view -i '(FMT/DV + FMT/RV) >= 15' -Ov -o ${sample}.delly.filtered.vcf ${vcf}
    """

}

process BCFTOOLS_FILTERSVABA {

    conda 'envs/bcftools_env.yml'
    label "process_high"
    publishDir "${params.outdir}/vcf_filtered"

    // takes in a zipped vcf and index as input
    input:
    path(vcf)

    // output as an unzipped vcf file
    output:
    path("*svaba.filtered.vcf"), emit: filtered_vcf

    // run filter with bcftools
    shell:
    def name = vcf.simpleName
    """
    bcftools view -i 'COUNT(FMT/AD >= 15 || (FMT/DR + FMT/SR) >= 15) > 0' -Ov -o ${name}.svaba.filtered.vcf ${vcf}
    """

}

process BCFTOOLS_FILTERFORINS {

    conda 'envs/bcftools_env.yml'
    label "process_high"
    publishDir "${params.outdir}/vcf_filtered_ins"

    // takes in a zipped vcf and index as input
    input:
    tuple val(sample), path(vcf)

    // output as an unzipped vcf file
    output:
    tuple val(sample), path("*only.ins.vcf"), emit: filtered_vcf

    // run filter with bcftools
    // The if statement creates an empty VCF with just the header when no insertions pass the filter in BCFTOOLS_CONCAT
    script:
    def name = vcf.simpleName
    """
    grep -v "^##contig=<ID=HLA" ${vcf} | grep -v "^HLA" > ${name}.INS.SURVIVOR.noHLA.vcf

    bcftools view -i '(FILTER="PASS") && ((INFO/SUPP="2" || INFO/SUPP="3"))' -Ov -o ${name}.INS.SURVIVOR.filtered_quality.vcf ${name}.INS.SURVIVOR.noHLA.vcf
    bcftools view -i '((INFO/SVTYPE="INS") || (ALT="<INS>") || (ALT~"^<INS"))' -Ov -o ${name}.INS.SURVIVOR.only.ins.vcf ${name}.INS.SURVIVOR.filtered_quality.vcf

    if [ ! -s ${name}.INS.SURVIVOR.only.ins.vcf ]; then
    bcftools view -h ${name}.INS.SURVIVOR.filtered_quality.vcf > ${name}.INS.SURVIVOR.only.ins.vcf
    fi
    """

}

process BCFTOOLS_CONCAT {

    conda 'envs/bcftools_env.yml'
    label "process_high"
    publishDir "${params.outdir}/survivor_concat", mode: 'copy'

    // takes in a zipped vcf and index as input
    input:
    tuple val(sample), path(survivor_1), path(survivor_2)

    // output as an unzipped vcf file
    output:
    tuple path("*SURVIVOR.concat.vcf.gz"), path("*SURVIVOR.concat.vcf.gz.tbi"), emit: vcf_allSV
    // path("*SURVIVOR.concat.vcf"), emit: concat_vcf
    // or this to include sample names for subsequent processes:    
    tuple val(sample), path("*SURVIVOR.concat.vcf"), emit: concat_vcf


    // run filter with bcftools
    script:
    """
    grep -v "^##" ${survivor_2} > no_header.vcf
    grep -v "^#" no_header.vcf > variant_only.vcf

    cat ${survivor_1} variant_only.vcf > unsorted_concat.vcf

    bcftools sort unsorted_concat.vcf -o ${sample}.SURVIVOR.concat.vcf

    bcftools view ${sample}.SURVIVOR.concat.vcf -Oz -o ${sample}.SURVIVOR.concat.vcf.gz
    tabix -p vcf ${sample}.SURVIVOR.concat.vcf.gz
    """
}

process BCFTOOLS_UNZIP {

    conda 'envs/bcftools_env.yml'
    label "process_high"
    publishDir "${params.outdir}/vcf_filtered_commonsv"

    // takes in a zipped vcf and index as input
    input:
    path(common_vcf)
    path(common_vcf_index)

    // output as an unzipped vcf file
    output:
    path("*unzipped.vcf"), emit: unzipped_commonsv

    // run filter with bcftools
    script:
    """
    bcftools view -O v ${common_vcf} -o gnomad.v4.1.sv.sites.unzipped.vcf
    """
}

process BCFTOOLS_REMOVE_COMMONSV {
    conda 'envs/bcftools_env.yml'
    label "process_high"
    publishDir "${params.outdir}/vcf_filtered_commonsv"

    // takes in a zipped vcf and index as input
    input:
    tuple val(sample), path(all_sv), path(common_sv)

    // output as an unzipped vcf file
    output:
    tuple val(sample), path("${sample}.ALLSV.no_common.vcf"), emit: common_filtered_vcf

    // run filter with bcftools
    script:
    """
    set -euo pipefail

    awk 'NR==FNR { key[\$1 FS \$2]; next }
     /^#/ { print; next }
     !(\$1 FS \$2 in key)' ${common_sv} ${all_sv} > ${sample}.ALLSV.no_common.vcf
    """

}

process BCFTOOLS_FILTER_CONTIGS {
    conda 'envs/bcftools_env.yml'
    label "process_high"
    publishDir "${params.outdir}/vcf_final"

    // takes in a zipped vcf and index as input
    input:
    tuple val(sample), path(vcf)

    // output as an unzipped vcf file
    output:
    tuple val(sample), path("*.ALLSV.final.vcf"), emit: final_vcf

    // run filter with bcftools
    script:
    """
    bgzip -f ${vcf}
    tabix -f -p vcf ${vcf}.gz

    bcftools view \
      -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
      -Ov ${vcf}.gz \
      -o ${sample}.ALLSV.final.vcf
    """

}

process BCFTOOLS_FILTERSNV {

    label "process_high"
    conda 'envs/bcftools_env.yml'
    publishDir "${params.outdir}/gatk_mutect2_bcftools_filter_WG"

    input:
    tuple val(sample), path(vcf)
    
    output:
    tuple val(sample), path("*mutect2.filtered.hard.vcf"), emit: filtered_snv

    script:
    """
    bcftools view -i 'FORMAT/DP >= 15 && INFO/MMQ >= 40 && FILTER == "PASS"' -Ov -o ${sample}.mutect2.filtered.hard.vcf ${vcf}
    """

}

process BCFTOOLS_FILTERSIZE {

    conda 'envs/bcftools_env.yml'
    label "process_high"
    publishDir "${params.outdir}/vcf_filtered_size"

    // takes in a zipped vcf and index as input
    input:
    tuple val(sample), path(vcf)

    // output as an unzipped vcf file
    output:
    tuple val(sample), path("*size.filtered.vcf"), emit: filtered_vcf

    // run filter with bcftools
    shell:
    """
    bcftools view -i 'INFO/SVLEN < 1000000 && INFO/SVLEN > -1000000' -Ov -o ${sample}.size.filtered.vcf ${vcf}
    """

}