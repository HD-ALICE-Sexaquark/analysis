#!/bin/bash

FILES_PATH=$1
RN_LIST=$2

mkdir -p ${FILES_PATH}/merged
echo "Merging files in ${FILES_PATH}"

while read -r rn; do
    echo "Merging ${rn}"
    hadd -f -k -O ${FILES_PATH}/merged/AnalysisResults_${rn}.root ${FILES_PATH}/AnalysisResults_${rn}_*.root &>> ${FILES_PATH}/merged/merging_rn.log
done < ${RN_LIST}
