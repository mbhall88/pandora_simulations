#!/usr/bin/env bash
CLUSTER_CMD=("bsub -n {threads} -R \"select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1]\" -M {resources.mem_mb} -o {cluster.output} -e {cluster.error} -J {cluster.name}")
JOB_NAME=snakemake_master_process
LOG_DIR=logs/
MEMORY=4000

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY]" \
    -M "$MEMORY" \
    -o "$LOG_DIR"/cluster_"$JOB_NAME".o \
    -e "$LOG_DIR"/cluster_"$JOB_NAME".e \
    -J "$JOB_NAME" \
    snakemake --use-singularity \
    --use-conda \
    --cluster-config cluster.yaml \
    --jobs 2000 \
    --keep-going \
    --cluster "${CLUSTER_CMD[@]}" \
    --restart-times 1

exit 0
