#!/bin/bash

# !! don't execute this script directly, it is meant to be used by task_batch.sh !!

if [[ -z {TASK_MODE}       || -z {TASK_LOCAL_PATH} || -z {TASK_IS_MC} || -z {TASK_PROD_NAME} ||
      -z {TASK_SOURCE_V0S} || -z {TASK_SIM_SET}    || -z {TASK_DO_QA} || -z {TASK_READ_LOGS} ||
      -z {TASK_OUTPUT_DIR} ]]; then
    echo "task_single.sh :: make sure to set all the necessary variables!"
    exit 1
fi

cd ${TASK_OUTPUT_DIR}
aliroot -l -b -q 'runAnalysis.C("'${TASK_MODE}'", "'${TASK_LOCAL_PATH}'", 0, '${TASK_IS_MC}', "'${TASK_PROD_NAME}'", "'${TASK_OUTPUT_DIR}'/rn.txt", "'${TASK_SOURCE_V0S}'", "'${TASK_SIM_SET}'", "", '${TASK_DO_QA}', '${TASK_READ_LOGS}')'
