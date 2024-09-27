#!/bin/bash

FILES_PATH=$1

echo "Merging files in ${FILES_PATH}"

alihadd -O -s 500000000 -k ${FILES_PATH}/AnalysisResults_merged.root ${FILES_PATH}/AnalysisResults_*.root &> ${FILES_PATH}/merging.log
