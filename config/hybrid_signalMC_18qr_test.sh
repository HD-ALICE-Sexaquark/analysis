#!/bin/bash

if [[ -z ${SOURCE_OF_V0S} ]]; then echo "Missing env. var. SOURCE_OF_V0S" ; return 1; fi
if [[ -z ${SIMULATION_SET} ]]; then echo "Missing env. var. SIMULATION_SET" ; return 1; fi
if [[ -z ${ANALYSIS_DIR} ]]; then echo "Missing env. var. ANALYSIS_DIR" ; return 1; fi

export ATTEMPT_NAME="hybrid_signalMC_18qr_test_${SOURCE_OF_V0S}${SIMULATION_SET}"
export MODE="hybrid"
export LOCAL_INPUT_PATH="/alice/sim/2023"
export LOCAL_N_DIRS=110
export GRID_TEST_MODE=0
export IS_MC=1
export PRODUCTION_NAME="LHC23l1a3"
export RUN_NUMBERS_LIST="${ANALYSIS_DIR}/doc/LHC18qr_pass3_rn_TEST.txt"
export DO_QA=1
export CHOOSE_N_EVENTS=0
