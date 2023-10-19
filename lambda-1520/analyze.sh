#!/bin/bash

##############################################################
#                                                            #
# Script to analyze sim. with AliAnalysisTaskLambda1520Lpipi #
# in a single run number                                     #
#                                                            #
##############################################################

# 29.Aug.2023
## A. Bórquez

#############
# Functions #
#############

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--rn" ]]; then
            RUN_NUMBER=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--dir" ]]; then
            DIR_NUMBER=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--n" ]]; then
            N_EVENTS=${arr[$((ic+1))]}
        else
            echo "analyze.sh :: ERROR: unrecognized argument: ${arr[$((ic))]}."
            print_help
            exit 1
        fi
        ((ic+=2))
    done
}

function print_usage() {
    echo "analyze.sh :: SCRIPT: analyze.sh"
    echo "analyze.sh :: =================="
    echo "analyze.sh :: "
    echo "analyze.sh :: USAGE : ./analyze.sh --rn <run-number> --dir <dir-number>"
    echo "analyze.sh ::         where:"
    echo "analyze.sh ::         <run-number> : choose run number"
    echo "analyze.sh ::                        (default value: 297595)"
    echo "analyze.sh ::         <dir-number> : choose specific directory, within a run number dir"
    echo "analyze.sh ::                        (default value: *)"
    echo "analyze.sh ::         <n>          : choose number of events to process"
    echo "analyze.sh ::                        (default value: 0, which means all)"
    echo "analyze.sh ::"
    echo "analyze.sh :: EXAMPLES :"
    echo "analyze.sh :: ./analyze.sh"
    echo "analyze.sh :: ./analyze.sh --rn 297595"
}

##########################################
# Check for possible command-line errors #
##########################################

if [[ -z ${ALIBUILD_WORK_DIR} || -z ${ALIDPG_VERSION} || -z ${ALIPHYSICS_VERSION} ]]; then
    echo "analyze.sh :: ERROR : please, set your ALICE environment and load AliDPG and AliPhysics"
    exit 1
fi

#########################
# Hard-coded parameters #
#########################

CURRENT_DIR=${PWD}
DATA_SET="LHC20e3a"

########
# Main #
########

# process input
argArray=("$@")
process_args "${argArray[@]}"

# default values
if [[ -z "${RUN_NUMBER}" ]]; then
    RUN_NUMBER=297595
fi
if [[ -z "${DIR_NUMBER}" ]]; then
    DIR_NUMBER="*" # wildcard
fi
if [[ -z "${N_EVENTS}" ]]; then
    N_EVENTS=0 # all
fi

# after input options, decide further variables
INPUT_DIR="${HOME}/work/sim/${DATA_SET}/${RUN_NUMBER}/${DIR_NUMBER}"

# create output dir if it doesn't exist
OUT_DIR=${PWD}/../output/lambda-1520/${RUN_NUMBER}
mkdir -p ${OUT_DIR}

echo "analyze.sh :: initiating..."
echo "analyze.sh :: - run number  : ${RUN_NUMBER}"
echo "analyze.sh :: - dir number  : ${DIR_NUMBER}"
if [[ ${N_EVENTS} -eq 0 ]]; then
    echo "analyze.sh :: - n of events : all"
else
    echo "analyze.sh :: - n of events : ${N_EVENTS}"
fi
echo "analyze.sh :: " # empty line

echo "analyze.sh :: copying analysis files"
echo -n "analyze.sh :: "; cp -v ${CURRENT_DIR}/AliAnalysisTaskLambda1520Lpipi.cxx ${OUT_DIR}/
echo -n "analyze.sh :: "; cp -v ${CURRENT_DIR}/AliAnalysisTaskLambda1520Lpipi.h ${OUT_DIR}/
echo -n "analyze.sh :: "; cp -v ${CURRENT_DIR}/runAnalysis.C ${OUT_DIR}/
echo -n "analyze.sh :: "; cp -v ${CURRENT_DIR}/AddTaskLambda1520Lpipi.C ${OUT_DIR}/
echo "analyze.sh :: " # empty line

echo "analyze.sh :: moving into $(readlink -f ${OUT_DIR})"
cd ${OUT_DIR}
echo "analyze.sh :: " # empty line

# start loop, in case of wildcard
for dir in $(readlink -f ${INPUT_DIR}); do

    # extract dir number, remove wildcard
    DIR_NUMBER_NOWC="${dir##*/}"

    echo "analyze.sh :: >> processing ${DIR_NUMBER_NOWC}"
    echo "analyze.sh :: (directory: ${dir})"
    echo "analyze.sh :: " # empty line

    # check existence of required input files
    if [[ ! -e ${dir}/AliESDs.root ]] || [[ ! -e ${dir}/Kinematics.root ]] || [[ ! -e ${dir}/galice.root ]]; then
        echo "analyze.sh :: WARNING: required input files not found, skipping this run..."
        continue
    fi

    # symbolic link to required input files
    echo "analyze.sh :: bringing ${dir}/AliESDs.root"
    ln -s ${dir}/AliESDs.root
    echo "analyze.sh :: bringing ${dir}/Kinematics.root"
    ln -s ${dir}/Kinematics.root
    echo "analyze.sh :: bringing ${dir}/galice.root"
    ln -s ${dir}/galice.root
    echo "analyze.sh :: " # empty line

    echo "analyze.sh :: analyzing..."
    aliroot -l -b -q 'runAnalysis.C(1, '${N_EVENTS}')' 2>&1 | tee analysis.log

    # if the log file exists, move it and rename it
    if [[ -e analysis.log ]]; then
        echo -n "analyze.sh :: "; mv -v analysis.log AnalysisResults_${DIR_NUMBER_NOWC}.log
    fi

    # if the output file exists, move it and rename it
    if [[ -e AnalysisResults.root ]]; then
        echo -n "analyze.sh :: "; mv -v AnalysisResults.root AnalysisResults_${DIR_NUMBER_NOWC}.root
    fi
    echo "analyze.sh :: " # empty line

    # remove sym-linked root files
    echo -n "analyze.sh :: "; rm -v AliESDs.root
    echo -n "analyze.sh :: "; rm -v Kinematics.root
    echo -n "analyze.sh :: "; rm -v galice.root
    echo "analyze.sh :: " # empty line
done

# cleaning
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskLambda1520Lpipi.cxx
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskLambda1520Lpipi_cxx.d
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskLambda1520Lpipi_cxx.so
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskLambda1520Lpipi_cxx_ACLiC_dict_rdict.pcm
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskLambda1520Lpipi.h
echo -n "analyze.sh :: "; rm -v runAnalysis.C
echo -n "analyze.sh :: "; rm -v AddTaskLambda1520Lpipi.C
echo "analyze.sh :: " # empty line

# come back to original dir
echo "analyze.sh :: moving out of ${OUT_DIR}"
cd ${CURRENT_DIR}

echo "analyze.sh :: " # empty line
echo "analyze.sh :: finished"
