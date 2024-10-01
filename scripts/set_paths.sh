#!/bin/bash

export ANALYSIS_DIR=${HOME}/work/analysis

export SEXAQUARK_DIR=${ANALYSIS_DIR}/sexaquark
export OUTPUT_DIR=${ANALYSIS_DIR}/output
mkdir -p ${OUTPUT_DIR}
export SLURM_LOGS_DIR=${ANALYSIS_DIR}/slurm_logs
mkdir -p ${SLURM_LOGS_DIR}

# export DEBUG_DIR=${ANALYSIS_DIR}/debug
# mkdir -p ${DEBUG_DIR}

export SIMS_DIR=${HOME}/some/sims
