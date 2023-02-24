#!/bin/bash

##########################################################
#                                                        #
# Script to analyze sim. with AliAnalysisTaskSexaquark   #
#                                                        #
##########################################################

# 04.Dec.2022
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
        elif [[ "${arr[$ic]}" == "--dir" ]]; then
            DIR_NUMBER=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--channel" ]]; then
            REACTION_CHANNEL=${arr[$((ic+1))]}
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
    echo "analyze.sh :: USAGE : ./analyze.sh --sim <sim-set> --v0s <v0s-option> --rn <run-number> --dir <dir-number> --channel <reaction-channel>"
    echo "analyze.sh ::         where:"
    echo "analyze.sh ::         <sim-set> : simulation set"
    echo "analyze.sh ::                     - only_bkg"
    # echo "analyze.sh ::                     - only_signal" # PENDING!
    echo "analyze.sh ::                     - only_V0s"
    echo "analyze.sh ::                     - signal+bkg"
    echo "analyze.sh ::         <v0s-option> : source of V0 particles"
    echo "analyze.sh ::                        - true"
    echo "analyze.sh ::                        - official"
    echo "analyze.sh ::                        - custom"
    echo "analyze.sh ::         <run-number> : choose run number"
    echo "analyze.sh ::                        (default value: 297595)"
    echo "analyze.sh ::         <dir-number> : choose specific directory, within a run number dir"
    echo "analyze.sh ::                        (default value: *)"
    echo "analyze.sh ::         <reaction-channel> : reaction channel to analyze" # PENDING! will be important later
    echo "analyze.sh ::                              - A (AntiS + N -> AntiL + K0)"
    echo "analyze.sh ::                              - B (AntiS + N -> AntiL + K0 + pim + pip)"
    echo "analyze.sh ::                              - C (AntiS + N -> AntiP + K0 + K0 + pip)"
    echo "analyze.sh ::                              - D (AntiS + P -> AntiL + Kp)"
    echo "analyze.sh ::                              - E (AntiS + P -> AntiL + Kp + pim + pip)"
    echo "analyze.sh ::                              - F (AntiS + P -> AntiP + Kp + K0 + pip)"
    echo "analyze.sh ::                              - G (AntiS + N -> Xip + pim)"
    echo "analyze.sh ::                              (default value: A)"
    echo "analyze.sh ::"
    echo "analyze.sh :: EXAMPLES :"
    echo "analyze.sh :: ./analyze.sh --sim signal+bkg --v0s true --rn 297595 --dir 000"
    echo "analyze.sh :: ./analyze.sh --sim only_bkg --v0s official --rn 297595"
    echo "analyze.sh :: ./analyze.sh --sim only_V0s --v0s official --rn 297595"
}

##########################################
# Check for possible command-line errors #
##########################################

if [[ -z ${ALIBUILD_WORK_DIR} || -z ${ALIDPG_VERSION} || -z ${ALIPHYSICS_VERSION} ]]; then
    echo "analyze.sh :: ERROR : please, set your ALICE environment and load AliDPG and AliPhysics"
    exit 1
fi

if [[ $# -lt 4 ]]; then
    echo "analyze.sh :: ERROR : insufficient number of arguments, you need to set at least a simulation set and a source of V0s."
    echo "analyze.sh :: " # empty line
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
    echo "analyze.sh :: ERROR : please, choose a simulation set."
    echo "analyze.sh :: " # empty line
    print_usage
    exit 1
fi

if [[ -z "${V0S_OPTION}" ]]; then
    echo "analyze.sh :: ERROR : please, choose a source of V0s."
    echo "analyze.sh :: " # empty line
    print_usage
    exit 1
fi

# default values
if [[ -z "${REACTION_CHANNEL}" ]]; then
    REACTION_CHANNEL="A"
fi
if [[ -z "${RUN_NUMBER}" ]]; then
    RUN_NUMBER="297595"
fi
if [[ -z "${DIR_NUMBER}"  ]]; then
    DIR_NUMBER="*" # wildcard
fi

# after input options, decide further variables
if [[ "${V0S_OPTION}" == "official" ]]; then
    V0S_STR="OfficialV0s"
elif [[ "${V0S_OPTION}" == "custom" ]]; then
    V0S_STR="CustomV0s"
elif [[ "${V0S_OPTION}" == "true" ]]; then
    V0S_STR="TrueV0s"
fi

if [[ "${REACTION_CHANNEL}" == "A" ]]; then
    CHANNEL_STR="AntiS-N_AntiL-K0"
elif [[ "${REACTION_CHANNEL}" == "B" ]]; then
    CHANNEL_STR="AntiS-N_AntiL-K0-pim-pip"
elif [[ "${REACTION_CHANNEL}" == "C" ]]; then
    CHANNEL_STR="AntiS-N_AntiP-K0-K0-pip"
elif [[ "${REACTION_CHANNEL}" == "D" ]]; then
    CHANNEL_STR="AntiS-P_AntiL-Kp"
elif [[ "${REACTION_CHANNEL}" == "E" ]]; then
    CHANNEL_STR="AntiS-P_AntiL-Kp-pim-pip"
elif [[ "${REACTION_CHANNEL}" == "F" ]]; then
    CHANNEL_STR="AntiS-P_AntiP-Kp-K0-pip"
elif [[ "${REACTION_CHANNEL}" == "G" ]]; then
    CHANNEL_STR="AntiS-N_Xip-pim"
fi

if [[ "${SIM_SET}" == "only_bkg" ]]; then
    INPUT_DIR="${HOME}/work/sim/LHC22b9a1/${RUN_NUMBER}/${DIR_NUMBER}"
    HAS_SEXAQUARK=0
elif [[ "${SIM_SET}" == "only_signal" ]]; then
    INPUT_DIR="PENDING"
    HAS_SEXAQUARK=1
elif [[ "${SIM_SET}" == "only_V0s" ]]; then
    INPUT_DIR="${HOME}/work/sim/LHC22b7d/${RUN_NUMBER}/${DIR_NUMBER}"
    HAS_SEXAQUARK=0
elif [[ "${SIM_SET}" == "signal+bkg" ]]; then
    INPUT_DIR="${HOME}/work/sim/signal-v1+bkg/${CHANNEL_STR}/${RUN_NUMBER}/${DIR_NUMBER}"
    HAS_SEXAQUARK=1
fi

echo "analyze.sh :: initiating..."
echo "analyze.sh :: - simulation set   : ${SIM_SET}"
echo "analyze.sh :: - source of V0s    : ${V0S_OPTION}"
echo "analyze.sh :: - run number       : ${RUN_NUMBER}"
echo "analyze.sh :: - dir number       : ${DIR_NUMBER}"
echo "analyze.sh :: - reaction channel : ${REACTION_CHANNEL} (${CHANNEL_STR})"
echo "analyze.sh :: " # empty line

# create output dir if it doesn't exist
CURRENT_DIR=${PWD}

OUT_DIR=${PWD}/output/${SIM_SET}/${RUN_NUMBER}
mkdir -p ${OUT_DIR}

echo "analyze.sh :: copying analysis files"
echo -n "analyze.sh :: "; cp -v AliAnalysisTaskSexaquark.cxx ${OUT_DIR}/
echo -n "analyze.sh :: "; cp -v AliAnalysisTaskSexaquark.h ${OUT_DIR}/
echo -n "analyze.sh :: "; cp -v runAnalysis.C ${OUT_DIR}/
echo -n "analyze.sh :: "; cp -v AddSexaquark.C ${OUT_DIR}/
echo "analyze.sh :: " # empty line

echo "analyze.sh :: moving into ${OUT_DIR}"
cd ${OUT_DIR}
echo "analyze.sh :: " # empty line

# start loop, in case of wildcard
for dir in $(readlink -f ${INPUT_DIR}); do

    # extract dir number, remove wildcard
    DIR_NUMBER_NOWC="${dir##*/}"

    echo "analyze.sh :: >> processing ${DIR_NUMBER_NOWC}"
    echo "analyze.sh :: (directory: ${dir})"
    echo "analyze.sh :: " # empty line

    # symbolic link to required input files
    echo "analyze.sh :: bringing ${dir}/AliESDs.root"
    ln -s ${dir}/AliESDs.root
    echo "analyze.sh :: bringing ${dir}/Kinematics.root"
    ln -s ${dir}/Kinematics.root
    echo "analyze.sh :: bringing ${dir}/galice.root"
    ln -s ${dir}/galice.root
    echo "analyze.sh :: " # empty line

    echo "analyze.sh :: analyzing..."
    aliroot -l -b -q 'runAnalysis.C(1, '${HAS_SEXAQUARK}', "'${V0S_OPTION}'", "'${REACTION_CHANNEL}'")' &> analysis.log

    # if the log file exists, move it and rename it
    if [[ -e analysis.log ]]; then
        echo -n "analyze.sh :: "; mv -v analysis.log AnalysisResults_${V0S_STR}_${DIR_NUMBER_NOWC}.log
    fi

    # if the output file exists, move it and rename it
    if [[ -e AnalysisResults.root ]]; then
        echo -n "analyze.sh :: "; mv -v AnalysisResults.root AnalysisResults_${V0S_STR}_${DIR_NUMBER_NOWC}.root
    fi
    echo "analyze.sh :: " # empty line

    # remove sym-linked root files
    echo -n "analyze.sh :: "; rm -v AliESDs.root
    echo -n "analyze.sh :: "; rm -v Kinematics.root
    echo -n "analyze.sh :: "; rm -v galice.root
    echo "analyze.sh :: " # empty line
done

# cleaning
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskSexaquark.cxx
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskSexaquark_cxx.d
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskSexaquark_cxx.so
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskSexaquark_cxx_ACLiC_dict_rdict.pcm
echo -n "analyze.sh :: "; rm -v AliAnalysisTaskSexaquark.h
echo -n "analyze.sh :: "; rm -v runAnalysis.C
echo -n "analyze.sh :: "; rm -v AddSexaquark.C
echo "analyze.sh :: " # empty line

# come back to original dir
echo "analyze.sh :: moving out of ${OUT_DIR}"
cd ${CURRENT_DIR}

echo "analyze.sh :: " # empty line
echo "analyze.sh :: finished"
