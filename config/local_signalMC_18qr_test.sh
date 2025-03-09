#!/bin/bash

# IMPORTANT: this file should always be sourced by another script

if [[ -z ${LOCAL_SIMS_DIR} ]]; then echo "ERROR :: config_file.sh :: Missing env. var. LOCAL_SIMS_DIR" ; exit 1; fi

export MODE="local"
export PRODUCTION_NAME="LHC23l1a3"
export SIMULATION_SET="A1.8"
export ATTEMPT_NAME="local_signalMC_${SIMULATION_SET}_18qr_test"
export INPUT_PATH="${LOCAL_SIMS_DIR}/${PRODUCTION_NAME}/${SIMULATION_SET}"
export RUN_NUMBER=295585
export LOCAL_N_DIRS=2
export LOCAL_LIMIT_N_EVENTS=10
export GRID_TEST_MODE=0 # not used
export GRID_WORKING_DIR="" # not used
export GRID_CUSTOM_SPLIT=0 # not used
export GRID_CUSTOM_PATTERN="" # not used
