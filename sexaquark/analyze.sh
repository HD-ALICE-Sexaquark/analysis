#!/bin/bash

##########################################################
#                                                        #
# Script to analyze sim. with AliAnalysisTaskSexaquark   #
#                                                        #
##########################################################

# 25.Jan.2024
## A. Bórquez

#############
# Hardcoded #
#############

INPUT_DIR= # find and modify below...

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
        elif [[ "${arr[$ic]}" == "--period" ]]; then
            RUN_PERIOD=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--rn" ]]; then
            RUN_NUMBER=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--dir" ]]; then
            DIR_NUMBER=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--n" ]]; then
            N_EVENTS=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--help" ]]; then
            print_usage
            exit 0
        else
            echo "analyze.sh :: ERROR: unrecognized argument: ${arr[$((ic))]}."
            print_usage
            exit 1
        fi
        ((ic+=2))
    done
}

function print_usage() {
    echo "analyze.sh :: SCRIPT: analyze.sh"
    echo "analyze.sh :: =================="
    echo "analyze.sh :: "
    echo "analyze.sh :: USAGE : ./analyze.sh --sim <sim-set> --v0s <v0s-option> --period <run-period> --rn <run-number> --dir <dir-number> --n <n>"
    echo "analyze.sh ::         where:"
    echo "analyze.sh ::         <sim-set> : simulation set, should be formatted as <reaction_channel><sexaquark_mass>"
    echo "analyze.sh ::                     where the five possible sexaquark masses are 1.73, 1.8, 1.87, 1.94, 2.01, and the reaction channels:"
    echo "analyze.sh ::                     - A : AntiSexaquark + Neutron -> AntiLambda, K0"
    echo "analyze.sh ::                     - D : AntiSexaquark + Proton -> AntiLambda, K+"
    echo "analyze.sh ::                     - E : AntiSexaquark + Proton -> AntiLambda, K+, pi-, pi+"
    echo "analyze.sh ::                     - H : AntiSexaquark + Proton -> AntiProton, K+, K+, pi0"
    echo "analyze.sh ::                     (default value: A1.8)"
    echo "analyze.sh ::         <v0s-option> : source of V0 particles"
    echo "analyze.sh ::                        - offline"
    echo "analyze.sh ::                        - online"
    echo "analyze.sh ::                        - custom"
    echo "analyze.sh ::                        - kalman"
    echo "analyze.sh ::         <run-period> : choose run period"
    echo "analyze.sh ::                        - LHC15o"
    echo "analyze.sh ::                        - LHC18q"
    echo "analyze.sh ::                        - LHC18r (default)"
    echo "analyze.sh ::         <run-number> : choose run number"
    echo "analyze.sh ::                        (default value: 297595)"
    echo "analyze.sh ::         <dir-number> : choose specific directory, within a run number dir"
    echo "analyze.sh ::                        (default value: *)"
    echo "analyze.sh ::         <n>          : choose number of events"
    echo "analyze.sh ::                        (default value: 0, which means all)"
    echo "analyze.sh ::"
    echo "analyze.sh :: EXAMPLES :"
    echo "analyze.sh :: ./analyze.sh --sim A1.8 --v0s online --period LHC18r --rn 297595 --dir 000"
}

##########################################
# Check for possible command-line errors #
##########################################

if [[ -z ${ALIBUILD_WORK_DIR} || -z ${ALIDPG_VERSION} || -z ${ALIPHYSICS_VERSION} ]]; then
    echo "analyze.sh :: ERROR : please, set your ALICE environment and load AliDPG and AliPhysics"
    exit 1
fi

if [[ $# -lt 1 ]]; then
    echo "analyze.sh :: ERROR : insufficient number of arguments, you need to at least set the source of V0s."
    exit 1
fi

########
# Main #
########

# process input
argArray=("$@")
process_args "${argArray[@]}"

# check for input errors
if [[ -z "${V0S_OPTION}" ]]; then
    echo "analyze.sh :: ERROR : please, choose a source of V0s."
    echo "analyze.sh :: " # empty line
    print_usage
    exit 1
fi

# set default values
if [[ -z "${SIM_SET}" ]]; then
    SIM_SET="A1.8"
fi
if [[ -z "${RUN_PERIOD}" ]]; then
    RUN_PERIOD="LHC18r"
fi
if [[ -z "${RUN_NUMBER}" ]]; then
    RUN_NUMBER="297595"
fi
if [[ -z "${DIR_NUMBER}"  ]]; then
    DIR_NUMBER="*" # wildcard
fi
if [[ -z "${N_EVENTS}"  ]]; then
    N_EVENTS=0 # all
fi

# derivate more options
REACTION_ID=${SIM_SET:0:1}
SEXAQUARK_MASS=${SIM_SET:1}
if [[ "${REACTION_ID}" == "A" ]]; then
    CHANNEL_STR="AntiS-N_AntiL-K0"
elif [[ "${REACTION_ID}" == "D" ]]; then
    CHANNEL_STR="AntiS-P_AntiL-Kp"
elif [[ "${REACTION_ID}" == "E" ]]; then
    CHANNEL_STR="AntiS-P_AntiP-Kp-K0-pip"
elif [[ "${REACTION_ID}" == "H" ]]; then
    CHANNEL_STR="AntiS-N_Xip-pim"
fi

echo "analyze.sh :: initiating..."
echo "analyze.sh :: - run period             : ${RUN_PERIOD}"
echo "analyze.sh :: - run number             : ${RUN_NUMBER}"
echo "analyze.sh :: - simulation set         : ${SIM_SET}"
echo "analyze.sh ::   -> reaction id         : ${REACTION_ID}"
echo "analyze.sh ::   -> reaction channel    : ${CHANNEL_STR}"
echo "analyze.sh ::   -> inj. sexaquark mass : ${SEXAQUARK_MASS}"
echo "analyze.sh :: - source of V0s          : ${V0S_OPTION}"
echo "analyze.sh :: - dir number             : ${DIR_NUMBER}"
if [[ ${N_EVENTS} -eq 0 ]]; then
    echo "analyze.sh :: - n of events            : all"
else
    echo "analyze.sh :: - n of events            : ${N_EVENTS}"
fi
echo "analyze.sh :: " # empty line

# set input dir
INPUT_DIR=${HOME}/work/simulations/samples/${SIM_SET}_S14/${RUN_NUMBER}

# create output dir if it doesn't exist
CURRENT_DIR=${PWD}
OUT_DIR=${PWD}/output/${RUN_PERIOD}_${RUN_NUMBER}_${SIM_SET}
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
    aliroot -l -b -q 'runAnalysis.C(1, "'${RUN_PERIOD}'", "'${V0S_OPTION}'", "'${SIM_SET}'", '${N_EVENTS}')' 2>&1 | tee analysis.log

    # if the log file exists, move it and rename it
    if [[ -e analysis.log ]]; then
        echo -n "analyze.sh :: "; mv -v analysis.log AnalysisResults_${V0S_OPTION}V0s_${DIR_NUMBER_NOWC}.log
    fi

    # if the output file exists, move it and rename it
    if [[ -e AnalysisResults.root ]]; then
        echo -n "analyze.sh :: "; mv -v AnalysisResults.root AnalysisResults_${V0S_OPTION}V0s_${DIR_NUMBER_NOWC}.root
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
