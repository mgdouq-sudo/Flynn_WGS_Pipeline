# Flynn Lab Pediatric Osteosarcoma WGS Pipeline

## Overview
This repository is based on [BU-BMSIP/Flynn_WGS_Analysis](https://github.com/BU-BMSIP/Flynn_WGS_Analysis), originally developed by Joshua Keegan, Sydney Sorbello, Joakin Mori, and Shugo Muratani in the Flynn Lab at Boston University School of Medicine. It includes modifications to support downstream analyses of ALT pathway activation, mutational timing, and mutational profiling in pediatric osteosarcoma.

The pipeline processes whole-genome sequencing (WGS) data from osteosarcoma xenograft tumor samples to detect structural variants (SVs), copy number variants (CNVs), and single nucleotide variants (SNVs). The workflow is implemented in Nextflow (v24.04.2) with supporting scripts in R (v4.5.1) and Bash (v5.3).

The pipeline is designed for samples initially xenografted into mouse models, sequenced with Illumina short-read WGS (paired-end, 60× coverage), and delivered in CRAM format. It automatically handles mouse contamination removal and integrates results from multiple SV callers to reduce false positives.

---

## Pipeline Workflow

### 1. WGS Data Processing
- **QC**: Xengsort classifies reads into graft (human), host (mouse), both, ambiguous, and neither — only graft, both, and combined reads proceed
- **Formatting**: The final BAM file is sorted, duplicate-marked, and subsetted to the TP53/TCAB1 locus

### 2. Structural Variant Detection
- **Detection**: Three SV callers are deployed — Manta, Delly, and SvABA
- **Merging**: SV calls are merged using SURVIVOR to reduce false positives, retaining SVs detected by at least two of three tools within 50 bp
- **Filtering**: Common SVs (gnomAD v4.1), repetitive regions (RepeatMasker), and non-canonical chromosomes are removed; SVs >1Mb are excluded
- **Annotation**: Final SV calls are annotated using AnnotSV

### 3. Copy Number Variant Detection
- **Detection**: GATK4 somatic CNV workflow using a Panel of Normals (PoN) built from 5 TARGET-ALL-P2 normal samples
- **Absolute CN**: DoAbsolute is used to calculate absolute copy number, tumor purity, and ploidy

### 4. Single Nucleotide Variant Analysis
- **Detection**: GATK Mutect2 in tumor-only mode
- **Filtering**: BCFtools retains calls with ≥15 supporting reads, median mapping quality ≥40, and PASS filter
- **Annotation**: GATK VariantAnnotation (dbSNP) and ANNOVAR (pathogenicity)

---

## Downstream Analyses
This pipeline was modified to support the following downstream analyses, which are conducted separately:
- ALT status prediction using telomere feature extraction and Random Forest classification
- Mutational timing analysis using PhylogicNDT
- Mutational profiling including VAF distributions and gene-level recurrence

---

## Requirements

### Software
| Tool | Version |
|------|---------|
| Nextflow | v24.04.2 |
| R | v4.5.1 |
| Bash | v5.3 |
| SAMtools | v1.21 |
| BBMap | v39.26 |
| Xengsort | v2.0.8 |
| BWA | v0.7.19 |
| SeqKit | v2.10.0 |
| Trimmomatic | v0.39 |
| BCFtools | v1.22 |
| Picard | v3.4.0 |
| Manta | v1.6.0 |
| Delly | v1.3.3 |
| SvABA | v1.2.0 |
| SURVIVOR | v1.0.7 |
| Tabix | v1.11 |
| bedtools | v2.30.0 |
| AnnotSV | v3.4.6 |
| GATK | v4.6.2.0 |
| ANNOVAR | v2025Mar02 |
| DoAbsolute | v2.2 |
| GDC client | v2.3.0 |
| PhylogicNDT | - |

### Reference Data
- Human genome: hg38 (UCSC-annotated FASTA)
- Mouse genome: mm39 (UCSC-annotated FASTA)
- GENCODE annotation: hg38 GFF (recommended)
- gnomAD v4.1 (SV sites)
- RepeatMasker hg38 BED
- 1000 Genomes PoN (for Mutect2)
- dbSNP v138

---

## Input Files
- `samplesheet.csv` — sample ID and CRAM file paths
- `matched_normal.csv` — sample ID and matched normal BAM paths
- CRAM files and reference genomes stored in `refs/`

---

## Repository Contents

| Directory/File | Description |
|---|---|
| `main.nf` | Main Nextflow script orchestrating the pipeline |
| `nextflow.config` | Configuration file for resource allocation and execution |
| `samplesheet.csv` | Sample sheet with CRAM paths |
| `matched_normal.csv` | Matched normal BAM paths for PoN construction |
| `run_nextflow_job.sh` | Batch job submission script for HPC cluster |
| `terra_to_hpc.sh` | Script to transfer data from Terra/GCS to HPC |
| `envs/` | Conda environment definitions (.yml) for all tools |
| `modules/` | Nextflow process modules for all pipeline tools |
| `scripts/` | Supporting R and Bash scripts |
| `results/` | Pipeline outputs (auto-generated) |
| `work/` | Nextflow intermediate files (auto-generated, cleared periodically) |

---

## Running the Pipeline

### 1. Set up your environment

Create or activate the Nextflow conda environment:
```bash
conda create -n nextflow_base nextflow
conda activate nextflow_base
```

### 2. Configure paths

Edit `nextflow.config` and set `proj_root` to your project directory and `your_project_group` to your HPC project group.

### 3. Run the pipeline

```bash
nextflow run main.nf -profile conda,singularity,cluster
```

Or submit as a batch job:
```bash
qsub run_nextflow_job.sh
```

---

## Principal Investigator
**Dr. Rachel Flynn** — rlflynn@bu.edu
Departments of Pharmacology, Physiology & Biophysics, and Medicine
Boston University Chobanian & Avedisian School of Medicine

