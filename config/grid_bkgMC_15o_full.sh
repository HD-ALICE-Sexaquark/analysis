#!/bin/bash

if [[ -z ${SOURCE_OF_V0S} ]]; then echo "Missing env. var. SOURCE_OF_V0S" ; return 1; fi
if [[ -z ${ANALYSIS_DIR} ]]; then echo "Missing env. var. ANALYSIS_DIR" ; return 1; fi

export ATTEMPT_NAME="grid_bkgMC_15o_full_${SOURCE_OF_V0S}"
export MODE="grid"
export LOCAL_INPUT_PATH=""
export LOCAL_N_DIRS=0
export GRID_TEST_MODE=0
export IS_MC=1
export PRODUCTION_NAME="LHC20j6a"
export RUN_NUMBERS_LIST="${ANALYSIS_DIR}/doc/LHC15o_pass2_rn.txt"
export DO_QA=1
export CHOOSE_N_EVENTS=0
