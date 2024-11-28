#!/bin/bash

if [[ -z ${ANALYSIS_DIR} || -z ${SOURCE_OF_V0S} || -z ${SIMULATION_SET} || -z ${SIMS_DIR} ]]; then exit 1; fi

export ATTEMPT_NAME="local_data_18qr_test_${SOURCE_OF_V0S}${SIMULATION_SET}"
export MODE="local"
export LOCAL_INPUT_PATH="${SIMS_DIR}"
export LOCAL_N_DIRS=6
export GRID_TEST_MODE=0
export IS_MC=0
export PRODUCTION_NAME="LHC18qr"
export RUN_NUMBERS_LIST="${ANALYSIS_DIR}/doc/LHC18qr_pass3_rn_TEST.txt"
export DO_QA=1
export CHOOSE_N_EVENTS=0
