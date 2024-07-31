#!/bin/bash

if [[ -z ${ANALYSIS_DIR} || -z ${SEXAQUARK_DIR} || -z ${OUTPUT_DIR} || -z ${SIMS_DIR} ]]; then
    echo "submit_many.sh :: make sure to first set:"
    echo "submiy_many.sh :: - ANALYSIS_DIR  : where the main repository is"
    echo "submit_many.sh :: - SEXAQUARK_DIR : where the AliAnalysisTaskSexaquark is"
    echo "submit_many.sh :: - OUTPUT_DIR    : where the output will be stored"
    echo "submit_many.sh :: - SIMS_DIR      : where the simulations are"
    exit 1
fi

mkdir -p ${OUTPUT_DIR}/slurm_logs
export SBATCH_SETUP="--partition=main"
SBATCH_SETUP+=" --output=${OUTPUT_DIR}/slurm_logs/slurm-%J.out"
SBATCH_SETUP+=" --error=${OUTPUT_DIR}/slurm_logs/slurm-%J.err"
SBATCH_SETUP+=" --time=05:00"
SBATCH_SETUP+=" --mem-per-cpu=2500"

export IS_MC=1
export MC_TYPE="signal"
export V0S_OPTION="kalman"
export RC="A"
export SM=1.8
export READ_LOGS=1
export DO_QA=1
export REWEIGHT_PT=1
export REWEIGHT_RADIUS=1

CURRENT_DIR=$(pwd)

DATASETS=( # "LHC15o_pass2"
           "LHC18q_pass3"
           "LHC18r_pass3"
           )

for ((i=0; i<${#DATASETS[@]}; i++)); do
    export DS=${DATASETS[$i]}

    RUN_NUMBERS=( $(cat "${ANALYSIS_DIR}/doc/${DS}_rn.txt") )

    for ((j=0; j<${#RUN_NUMBERS[@]}; j++)); do
        export RN=${RUN_NUMBERS[$j]}
        # export RN=${RUN_NUMBERS[0]} # CONTROL

        # determine SIMSET and RUN_DIR
        if [[ ${MC_TYPE} == "signal" ]]; then
            if [[ ${DS:0:6} == "LHC15o" ]]; then SS="LHC23l1b3"
            elif [[ ${DS:0:6} == "LHC18q" || ${DS:0:6} == "LHC18r" ]]; then SS="LHC23l1a3"; fi
            export RUN_DIR=${SIMS_DIR}/${SS}/${RC}${SM}/${RN}
        elif [[ ${MC_TYPE} == "bkg" ]]; then
            if [[ ${DS:0:6} == "LHC15o" ]]; then SS="LHC20j6a"
            elif [[ ${DS:0:6} == "LHC18q" || ${DS:0:6} == "LHC18r" ]]; then SS="LHC20e3a"; fi
            export RUN_DIR=${SIMS_DIR}/${SS}/${RN}
        fi

        sbatch ${SBATCH_SETUP} single_rn.sh
    done
done
