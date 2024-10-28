#!/bin/bash

# !! Don't execute this script directly, it is meant to be used by `task_batch.sh`` !!

cd ${TASK_OUTPUT_DIR}

ANALYSIS_OPTIONS="("
ANALYSIS_OPTIONS+="\"${MODE}\","
ANALYSIS_OPTIONS+="\"${LOCAL_INPUT_PATH}\","
ANALYSIS_OPTIONS+="${LOCAL_N_DIRS},"
ANALYSIS_OPTIONS+="${GRID_TEST_MODE},"
ANALYSIS_OPTIONS+="${IS_MC},"
ANALYSIS_OPTIONS+="\"${PRODUCTION_NAME}\","
ANALYSIS_OPTIONS+="\"rn.txt\","
ANALYSIS_OPTIONS+="\"${SOURCE_OF_V0S}\","
ANALYSIS_OPTIONS+="\"${SIMULATION_SET}\","
ANALYSIS_OPTIONS+="\"${STOP_AFTER}\","
ANALYSIS_OPTIONS+="${DO_QA},"
ANALYSIS_OPTIONS+="${READ_SIGNAL_LOGS}",
ANALYSIS_OPTIONS+="${CHOOSE_N_EVENTS}"
ANALYSIS_OPTIONS+=")"

ALIROOT_COMMAND='aliroot -l -b -q runAnalysis.C'${ANALYSIS_OPTIONS}''
echo "INFO :: task_single.sh :: ${ALIROOT_COMMAND}"
${ALIROOT_COMMAND} 2>&1 | tee analysis.log

rm -v AliAnalysisTaskSexaquark*
rm -v AddSexaquark.C
rm -v runAnalysis.C
