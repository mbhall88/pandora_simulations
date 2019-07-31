#!/usr/bin/env bash
JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')
LOG_DIR=logs/
MEMORY=4000
PROFILE="lsf"

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY]" \
    -M "$MEMORY" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
        snakemake --profile "$PROFILE" --keep-going 

exit 0
