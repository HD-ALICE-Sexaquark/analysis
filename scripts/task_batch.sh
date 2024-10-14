#!/bin/bash

# USAGE: source set_paths.sh and execute, options are hardcoded below!

if [[ -z ${ANALYSIS_DIR} || -z ${SEXAQUARK_DIR} || -z ${OUTPUT_DIR} || -z ${SLURM_LOGS_DIR} || -z ${SIMS_DIR} ]]; then
    echo "task_batch.sh :: make sure to \`source set_paths.sh\` or to set:"
    echo "task_batch.sh :: - ANALYSIS_DIR   : where the main analysis repository is"
    echo "task_batch.sh :: - SEXAQUARK_DIR  : where the AliAnalysisTaskSexaquark is"
    echo "task_batch.sh :: - OUTPUT_DIR     : where the output will be stored"
    echo "task_batch.sh :: - SLURM_LOGS_DIR : where the slurm logs will be stored"
    echo "task_batch.sh :: - SIMS_DIR       : where the simulations are"
    exit 1
fi

# OPTIONS
export TASK_MODE="local"
export TASK_LOCAL_PATH="${SIMS_DIR}"
export TASK_IS_MC=1
export TASK_PROD_NAME="LHC23l1a3"
export TASK_SOURCE_V0S="kalman"
export TASK_SIM_SET="A1.8"
export TASK_DO_QA=1
export TASK_READ_LOGS=1

RUN_NUMBERS_FILE=
if [[ ${TASK_PROD_NAME} == "LHC23l1a3" || ${TASK_PROD_NAME} == "LHC20e3a" ]]; then
    RUN_NUMBERS_FILE=${ANALYSIS_DIR}/doc/LHC18qr_pass3_rn.txt
elif [[ ${TASK_PROD_NAME} == "LHC23l1b3" || ${TASK_PROD_NAME} == "LHC20j6a" ]]; then
    RUN_NUMBERS_FILE=${ANALYSIS_DIR}/doc/LHC15o_pass2_rn.txt
fi

DATASET_SHORT=$(basename ${RUN_NUMBERS_FILE})
DATASET_SHORT=${DATASET_SHORT%%_*}
DATASET_SHORT=${DATASET_SHORT#LHC}

SBATCH_SETUP="--partition=main"
SBATCH_SETUP+=" --output=${SLURM_LOGS_DIR}/slurm-%J.out"
SBATCH_SETUP+=" --error=${SLURM_LOGS_DIR}/slurm-%J.err"
SBATCH_SETUP+=" --time=15:00"
SBATCH_SETUP+=" --mem-per-cpu=4000" # in MB

while IFS= read -r RUN_NUMBER; do
    export TASK_OUTPUT_DIR=${OUTPUT_DIR}/${TASK_SIM_SET}_${DATASET_SHORT}_${TASK_MODE}/${RUN_NUMBER}
    mkdir -p ${TASK_OUTPUT_DIR}
    echo ${RUN_NUMBER} > ${TASK_OUTPUT_DIR}/rn.txt

    ln -s ${SEXAQUARK_DIR}/AliAnalysisTaskSexaquark.cxx ${TASK_OUTPUT_DIR}/AliAnalysisTaskSexaquark.cxx
    ln -s ${SEXAQUARK_DIR}/AliAnalysisTaskSexaquark.h ${TASK_OUTPUT_DIR}/AliAnalysisTaskSexaquark.h
    ln -s ${SEXAQUARK_DIR}/runAnalysis.C ${TASK_OUTPUT_DIR}/runAnalysis.C
    ln -s ${SEXAQUARK_DIR}/AddSexaquark.C ${TASK_OUTPUT_DIR}/AddSexaquark.C

    job_name="--job-name=${TASK_SIM_SET}_${DATASET_SHORT}_${TASK_MODE}_${RUN_NUMBER}"
    sbatch ${SBATCH_SETUP} ${job_name} ${ANALYSIS_DIR}/scripts/task_single.sh
done < "${RUN_NUMBERS_FILE}"
