#!/bin/bash

if [[ -z ${ANALYSIS_DIR} || -z ${SOURCE_OF_V0S} || -z ${SIMULATION_SET} ]]; then exit 1; fi

export ATTEMPT_NAME="grid_signalMC_15o_full_${SOURCE_OF_V0S}${SIMULATION_SET}"
export MODE="grid"
export LOCAL_INPUT_PATH=""
export LOCAL_N_DIRS=0
export GRID_TEST_MODE=0
export IS_MC=1
export PRODUCTION_NAME="LHC23l1b3"
export RUN_NUMBERS_LIST="${ANALYSIS_DIR}/doc/LHC15o_pass2_rn.txt"
export DO_QA=1
export CHOOSE_N_EVENTS=0
