#!/bin/bash

# !! do not execute interactively, this scripts depends on `mergers_batch.sh` !!

echo "starting merger for ${THE_PATH} with option ${THE_OPT}"

if [[ ${THE_OPT} == "QA" ]]; then
    root -l -b -q "${ANALYSIS_DIR}/macros/Merger_QA.C(\"${THE_PATH}\", \"${ANALYSIS_DIR}/macros/output\")"
elif [[ ${THE_OPT} == "A" ]]; then
    root -l -b -q "${ANALYSIS_DIR}/macros/Merger_ChannelA.C(\"${THE_PATH}\", \"${ANALYSIS_DIR}/macros/output\")"
elif [[ ${THE_OPT} == "D" ]]; then
    root -l -b -q "${ANALYSIS_DIR}/macros/Merger_ChannelD.C(\"${THE_PATH}\", \"${ANALYSIS_DIR}/macros/output\")"
elif [[ ${THE_OPT} == "E" ]]; then
    root -l -b -q "${ANALYSIS_DIR}/macros/Merger_ChannelE.C(\"${THE_PATH}\", \"${ANALYSIS_DIR}/macros/output\")"
elif [[ ${THE_OPT} == "H" ]]; then
    root -l -b -q "${ANALYSIS_DIR}/macros/Merger_ChannelH.C(\"${THE_PATH}\", \"${ANALYSIS_DIR}/macros/output\")"
fi
