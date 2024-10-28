#!/bin/bash

if [[ -z ${ANALYSIS_DIR} || -z ${SOURCE_OF_V0S} || -z ${SIMULATION_SET} || -z ${SIMS_DIR} ]]; then exit 1; fi

export ATTEMPT_NAME="local_bkgMC_18qr_test_${SOURCE_OF_V0S}${SIMULATION_SET}"
export MODE="local"
export LOCAL_INPUT_PATH="${SIMS_DIR}"
export LOCAL_N_DIRS=20
export GRID_TEST_MODE=0
export IS_MC=1
export PRODUCTION_NAME="LHC20e3a"
export RUN_NUMBERS_LIST="${ANALYSIS_DIR}/doc/LHC18qr_pass3_rn_TEST.txt"
export STOP_AFTER=""
export DO_QA=1
export READ_SIGNAL_LOGS=0
export CHOOSE_N_EVENTS=0