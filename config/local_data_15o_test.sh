#!/bin/bash

if [[ -z ${SOURCE_OF_V0S} ]]; then echo "Missing env. var. SOURCE_OF_V0S" ; return 1; fi
if [[ -z ${DATA_DIR} ]]; then echo "Missing env. var. DATA_DIR" ; return 1; fi
if [[ -z ${ANALYSIS_DIR} ]]; then echo "Missing env. var. ANALYSIS_DIR" ; return 1; fi

export ATTEMPT_NAME="local_data_15o_test_${SOURCE_OF_V0S}"
export MODE="local"
export LOCAL_INPUT_PATH="${DATA_DIR}/2015"
export LOCAL_N_DIRS=6
export GRID_TEST_MODE=0
export IS_MC=0
export PRODUCTION_NAME="LHC15o"
export RUN_NUMBERS_LIST="${ANALYSIS_DIR}/doc/LHC15o_pass2_rn_TEST.txt"
export DO_QA=0
export CHOOSE_N_EVENTS=0
