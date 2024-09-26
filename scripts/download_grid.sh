#!/bin/bash

# download output from the grid
ALIEN_PATH=$1
RUN_NUMBERS_LIST=$2
DEST_DIR=$3

DIRS_ARR=(1 2 3 4 5 6 7 8 9 10 11 12)

# loop over the run numbers
while IFS= read -r RUN_NUMBER
do
    for DIR_N in "${DIRS_ARR[@]}"; do
        DIR_STR=$(printf "%03i" $DIR_N)
        THE_COMMAND="alien_cp alien://${ALIEN_PATH}/${RUN_NUMBER}/${DIR_STR}/AnalysisResults.root file://${DEST_DIR}/AnalysisResults_${RUN_NUMBER}_${DIR_STR}.root"
        echo ${THE_COMMAND}
        ${THE_COMMAND}
    done
done < "${RUN_NUMBERS_LIST}"
