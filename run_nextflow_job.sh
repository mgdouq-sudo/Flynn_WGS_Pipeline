#!/bin/bash
#$ -P your_project_group
#$ -l h_rt=96:00:00
#$ -l mem_per_core=8G
#$ -pe omp 1
#$ -N nextflow_pipeline
#$ -j y
#$ -o nextflow_pipeline.log

source /share/pkg.8/miniconda/25.3.1/install/etc/profile.d/conda.sh
conda activate nextflow_base

cd /path/to/pipeline/

nextflow run main.nf -resume -profile conda,singularity,cluster
