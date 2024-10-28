#!/bin/bash

if [[ $# -ne 3 ]]; then
    echo "ERROR :: task_batch.sh :: wrong number of arguments"
    echo "ERROR :: task_batch.sh :: USAGE: 1. set ANALYSIS_DIR, SEXAQUARK_DIR, OUTPUT_DIR, SLURM_LOGS_DIR and SIMS_DIR"
    echo "ERROR :: task_batch.sh ::        2. \`./task_batch.sh <config_file> <source_of_v0s> <sim_set>\`"
    exit 1
fi

if [[ -z ${ANALYSIS_DIR} || -z ${SEXAQUARK_DIR} || -z ${OUTPUT_DIR} || -z ${SLURM_LOGS_DIR} || -z ${SIMS_DIR} ]]; then
    echo "ERROR :: task_batch.sh :: make sure to set:"
    echo "ERROR :: task_batch.sh :: - ANALYSIS_DIR   : where the main analysis repository is"
    echo "ERROR :: task_batch.sh :: - SEXAQUARK_DIR  : where the AliAnalysisTaskSexaquark is"
    echo "ERROR :: task_batch.sh :: - SIMS_DIR       : where the simulations are"
    echo "ERROR :: task_batch.sh :: - OUTPUT_DIR     : where the output will be stored"
    echo "ERROR :: task_batch.sh :: - SLURM_LOGS_DIR : where the slurm logs will be stored"
    exit 1
fi

echo "INFO :: task_batch.sh :: starting..."

# OPTIONS
CONFIG_FILE=$1
export SOURCE_OF_V0S=$2
export SIMULATION_SET=$3

source ${CONFIG_FILE}
echo "INFO :: task_batch.sh :: CONFIG_FILE      = ${CONFIG_FILE}"
echo "INFO :: task_batch.sh :: ATTEMPT_NAME     = ${ATTEMPT_NAME}"
echo "INFO :: task_batch.sh :: MODE             = ${MODE}"
echo "INFO :: task_batch.sh :: LOCAL_INPUT_PATH = ${LOCAL_INPUT_PATH}"
echo "INFO :: task_batch.sh :: LOCAL_N_DIRS     = ${LOCAL_N_DIRS}"
echo "INFO :: task_batch.sh :: GRID_TEST_MODE   = ${GRID_TEST_MODE}"
echo "INFO :: task_batch.sh :: IS_MC            = ${IS_MC}"
echo "INFO :: task_batch.sh :: PRODUCTION_NAME  = ${PRODUCTION_NAME}"
echo "INFO :: task_batch.sh :: RUN_NUMBERS_LIST = ${RUN_NUMBERS_LIST}"
echo "INFO :: task_batch.sh :: SOURCE_OF_V0S    = ${SOURCE_OF_V0S}"
echo "INFO :: task_batch.sh :: SIMULATION_SET   = ${SIMULATION_SET}"
echo "INFO :: task_batch.sh :: STOP_AFTER       = ${STOP_AFTER}"
echo "INFO :: task_batch.sh :: DO_QA            = ${DO_QA}"
echo "INFO :: task_batch.sh :: READ_SIGNAL_LOGS = ${READ_SIGNAL_LOGS}"
echo "INFO :: task_batch.sh :: CHOOSE_N_EVENTS  = ${CHOOSE_N_EVENTS}"

SBATCH_SETUP="--partition=main"
SBATCH_SETUP+=" --output=${SLURM_LOGS_DIR}/slurm-%J.out"
SBATCH_SETUP+=" --error=${SLURM_LOGS_DIR}/slurm-%J.err"
SBATCH_SETUP+=" --time=30:00"
SBATCH_SETUP+=" --mem-per-cpu=4000" # in MB

while IFS= read -r RUN_NUMBER; do
    export TASK_OUTPUT_DIR=${OUTPUT_DIR}/${ATTEMPT_NAME}/${RUN_NUMBER}
    mkdir -p ${TASK_OUTPUT_DIR}

    echo ${RUN_NUMBER} > ${TASK_OUTPUT_DIR}/rn.txt

    cp ${SEXAQUARK_DIR}/AliAnalysisTaskSexaquark.cxx ${TASK_OUTPUT_DIR}/AliAnalysisTaskSexaquark.cxx
    cp ${SEXAQUARK_DIR}/AliAnalysisTaskSexaquark.h ${TASK_OUTPUT_DIR}/AliAnalysisTaskSexaquark.h
    cp ${SEXAQUARK_DIR}/runAnalysis.C ${TASK_OUTPUT_DIR}/runAnalysis.C
    cp ${SEXAQUARK_DIR}/AddSexaquark.C ${TASK_OUTPUT_DIR}/AddSexaquark.C

    JOB_NAME="--job-name=${ATTEMPT_NAME}_${RUN_NUMBER}"
    echo -n "INFO :: task_batch.sh :: "
    sbatch ${SBATCH_SETUP} ${JOB_NAME} ${ANALYSIS_DIR}/scripts/task_single.sh
done < "${RUN_NUMBERS_LIST}"

echo "INFO :: task_batch.sh :: all jobs submitted."
