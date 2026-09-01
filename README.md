# Flynn Lab Pediatric Osteosarcoma Tumor WGS Variant Detection Pipeline

## Overview
This pipeline processes whole-genome sequencing (WGS) data from osteosarcoma xenograft tumor samples to detect structural variants (SVs), copy number variants (CNVs), chromothripsis, and single nucleotide variants (SNVs). The workflow is implemented in Nextflow (v24.04.2) with supporting scripts in R (v4.5.1) and Bash (v5.3).

For more information about the pipeline, please visit the [user guide](https://github.com/BU-BMSIP/Flynn_WGS_Analysis/blob/main/docs/User%20Guide%20Variant%20and%20Chromothripsis%20Detection%20for%20Xenografted%20Osteosarcoma%20Tumor%20Samples.pdf).

The pipeline is designed for samples initially xenografted into mouse models, sequenced with Illumina short-read WGS (paired-end, 60× coverage), and delivered in CRAM format. It automatically handles mouse contamination removal and integrates results from multiple SV callers to reduce false positives. CNV is performed using a Panel of Normals and with further computation to find the absolute copy number. Chromothripsis detection is performed using the SV and CNV analysis. Finally, SNV analysis is performed to investigate genes of interest. 
---
## Pipeline Workflow
1. WGS Data Processing
* QC: xengsort classify sorts reads into graft (human), host (mouse), both, ambiguous, and neither. Only graft, both, and combined datasets proceed.
* Formatting: the final BAM file is sorted and marked for duplicates

2. Structural Variant Detection
* Detection: Three SV detection tools are deployed including Manta, Delly, and SVABA.
* Merge: SV calls are merged using SURVIVOR in order to reduce false positive calls. 

3. Copy Number Variant Detection
* Detection: GATK4 is used to analyze somatic copy number.
* Calculation: DoAbsolute is used to calculate the absolute copy number. 

4. Chromothripsis Analysis
* Formatting: SV and CNV results are formatted into proper data formats. 
* Detection: ShatterSeek is used to detect chromothripsis.

5. Single Nucleotide Variant Analysis
* Detection: GATK Mutect2 detects SNVs in genes of interest.
* Annotation: ANNOVAR is used to annotate SNV calls.
--- 
## Requirements

Software
* Core: Nextflow (v24.04.2), Bash (v5.3), R (v4.5.1)
* Tools: SAMtools (v1.21), BBMap (v39.26), Xengsort (v2.0.8), BWA (v0.7.19), SeqKit (v2.10.0), Trimmomatic (v0.39), BCFtools (v1.22), Picard (v3.4.0), Manta (v1.6.0), Delly (v1.3.3), SvABA (v1.2.0), SURVIVOR (v1.0.7), Tabix (v1.11), bedtools (v2.30.0), AnnotSV (v3.4.6), GATK (v4.6.2.0), GDC client (v2.3.0), ShatterSeek (v1.1), GATK Mutect2 (v4.6.0.0), ANNOVAR (v2025Mar02)

Reference Data
* Human genome: hg38 (UCSC-annotated FASTA)
* Mouse genome: mm39 (UCSC-annotated FASTA)
* GENCODE annotation: hg38 GFF (recommended)
---
## Input Files
* Sample sheet (CSV): Columns for sample_id and cram_path.
* CRAM files: Located in refs/ directory unless otherwise specified.
* Reference genomes: Stored in refs/ directory.
---
## Running the Pipeline
In order to run the pipeline, the sample CRAM files and reference files must be located in the refs/ directory. Additionally, your sample sheet must be labelled 'samplesheet.csv' and must be located in the home directory of the project.

Create or activate the Nextflow Conda environment
If you don’t already have a Nextflow environment:
```
conda create -n nextflow_base nextflow
conda activate nextflow_base
```
Or, if you already have it:
```
conda activate nextflow_base
```
Run the pipeline
```
nextflow run main.nf -profile conda,singularity,cluster
```

