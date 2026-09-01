#!/bin/bash -l
#$ -P your_project_group
#$ -l download
#$ -l h_rt=24:00:00
#$ -N bucket_to_scc

module load python3/3.10.12
module load google-cloud-sdk

# Set these before running
DEST=/path/to/destination/
GCS_BUCKET=gs://your-bucket-name

mkdir -p "$DEST/crams"
mkdir -p "$DEST/metrics"

echo "Starting CRAM transfer..."
gsutil -m cp -r "$GCS_BUCKET/crams/" "$DEST/crams/"

echo "Starting metrics transfer..."
gsutil -m cp -r "$GCS_BUCKET/metrics/" "$DEST/metrics/"

echo "All done!"
