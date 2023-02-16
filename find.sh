#!/bin/bash

#######################################################
#                                                     #
# Script to analyze sim. with SexaquarkFinder macro   #
#                                                     #
#######################################################

# PENDING: different for different channels

# 15.Feb.2023
## A. Bórquez

#############
# Functions #
#############

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--sim" ]]; then
            SIM_SET=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--v0s" ]]; then
            V0S_OPTION=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--rn" ]]; then
            RUN_NUMBER=${arr[$((ic+1))]}
        else
            echo "find.sh :: ERROR: unrecognized argument: ${arr[$((ic))]}."
            print_help
            exit 1
        fi
        ((ic+=2))
    done
}

function print_usage() {
    echo "find.sh :: SCRIPT: find.sh"
    echo "find.sh :: =================="
    echo "find.sh :: "
    echo "find.sh :: USAGE : ./find.sh --sim <sim-set> --v0s <v0s-option> --rn <run-number>"
    echo "find.sh ::         where:"
    echo "find.sh ::         <sim-set> : simulation set"
    echo "find.sh ::                     - only_bkg"
    # echo "find.sh ::                     - only_signal" # PENDING!
    echo "find.sh ::                     - only_V0s"
    echo "find.sh ::                     - signal+bkg"
    echo "find.sh ::         <v0s-option> : source of V0 particles"
    echo "find.sh ::                        - true"
    echo "find.sh ::                        - official"
    echo "find.sh ::                        - custom"
    echo "find.sh ::         <run-number> : choose run number"
    echo "find.sh ::                        (default value: 297595)"
    echo "find.sh ::"
    echo "find.sh :: EXAMPLES :"
    echo "find.sh :: ./find.sh --sim signal+bkg --v0s true --rn 297595"
    echo "find.sh :: ./find.sh --sim only_bkg --v0s official --rn 297595"
    echo "find.sh :: ./find.sh --sim only_V0s --v0s official --rn 297595"
}

##########################################
# Check for possible command-line errors #
##########################################

if [[ -z ${ROOTSYS} ]]; then
    echo "find.sh :: ERROR : please, set ROOT."
    exit 1
fi

if [[ $# -lt 6 ]]; then
    echo "find.sh :: ERROR : insufficient number of arguments, you need to set at least a simulation set, source of V0s and run number."
    echo "find.sh :: " # empty line
    print_usage
    exit 1
fi

########
# Main #
########

# process input
argArray=("$@")
process_args "${argArray[@]}"

# another check for errors
if [[ -z "${SIM_SET}" ]]; then
    echo "find.sh :: ERROR : please, choose a simulation set."
    echo "find.sh :: " # empty line
    print_usage
    exit 1
fi

if [[ -z "${V0S_OPTION}" ]]; then
    echo "find.sh :: ERROR : please, choose a source of V0s."
    echo "find.sh :: " # empty line
    print_usage
    exit 1
fi

# default values
if [[ -z "${RUN_NUMBER}" ]]; then
    RUN_NUMBER="297595"
fi

# after input options, decide further variables
if [[ "${V0S_OPTION}" == "official" ]]; then
    V0S_STR="OfficialV0s"
elif [[ "${V0S_OPTION}" == "custom" ]]; then
    V0S_STR="CustomV0s"
elif [[ "${V0S_OPTION}" == "true" ]]; then
    V0S_STR="TrueV0s"
fi

ANALYSIS_DIR="${HOME}/work/analysis/output/${SIM_SET}/${RUN_NUMBER}"
HAS_SEXAQUARK=1 # valid for "${SIM_SET}" == "only_signal" or "signal+bkg"
if [[ "${SIM_SET}" == "only_bkg" || "${SIM_SET}" == "only_V0s" ]]; then
    HAS_SEXAQUARK=0
fi

echo "find.sh :: initiating..."
echo "find.sh :: - simulation set : ${SIM_SET}"
echo "find.sh :: - source of V0s  : ${V0S_OPTION}"
echo "find.sh :: - run number     : ${RUN_NUMBER}"
echo "find.sh :: " # empty line

# start loop over each file
for ((file_index = 0; file_index < 4; file_index++)); do

    ANALYSIS_FILE=${ANALYSIS_DIR}/AnalysisResults_${V0S_STR}_00${file_index}.root
    SEXAQUARK_FILE=${ANALYSIS_DIR}/SexaquarkResults_${V0S_STR}_00${file_index}.root
    SEXAQUARK_LOG=${ANALYSIS_DIR}/AnalysisResults_${V0S_STR}_00${file_index}.root

    TMP_FILE=${ANALYSIS_DIR}/SexaquarkResults.root
    TMP_LOG=${ANALYSIS_DIR}/sexaquark.log

    # if the input file (the analysis file -- the output of the analysis task) exists
    if [[ -e ${ANALYSIS_FILE} ]]; then
        echo "find.sh :: finding..."
        root -l -b -q 'SexaquarkFinder.cxx("'${ANALYSIS_FILE}'", "'${TMP_FILE}'")' &> ${TMP_LOG}
    fi

    # if the output file exists, rename it
    if [[ -e ${SEXAQUARK_FILE} ]]; then
        echo -n "find.sh :: "; mv -v ${TMP_FILE} ${SEXAQUARK_FILE}
    fi

    # if the log file exists, rename it
    if [[ -e ${SEXAQUARK_LOG} ]]; then
        echo -n "find.sh :: "; mv -v ${TMP_LOG} ${SEXAQUARK_LOG}
    fi

    echo "find.sh :: " # empty line

done

echo "find.sh :: finished"
