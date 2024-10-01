#!/bin/bash

if [[ -z ${ANALYSIS_DIR} || -z ${OUTPUT_DIR} || -z ${SLURM_LOGS_DIR} ]]; then
    echo "mergers_batch.sh :: make sure to \`source set_paths.sh\` or to set:"
    echo "mergers_batch.sh :: - ANALYSIS_DIR   : where the main repository is"
    echo "mergers_batch.sh :: - OUTPUT_DIR     : where the analysis task output is stored"
    echo "mergers_batch.sh :: - SLURM_LOGS_DIR : where the slurm logs will be stored"
    exit 1
fi

# reaction_channels=("A" "D" "E" "H")
reaction_channels=("H")
# data_periods=("15o" "18qr")
data_periods=("18qr")
# sexaquark_masses=(1.73 1.8 1.87 1.94 2.01)
sexaquark_masses=(1.73 1.8 1.87 1.94)

export SBATCH_SETUP="--partition=main"
SBATCH_SETUP+=" --output=${SLURM_LOGS_DIR}/slurm-%J.out"
SBATCH_SETUP+=" --error=${SLURM_LOGS_DIR}/slurm-%J.err"
SBATCH_SETUP+=" --time=30:00"
SBATCH_SETUP+=" --mem-per-cpu=120000"

# QA

# for reaction_id in ${reaction_channels[@]}; do
    # for sexaquark_mass in ${sexaquark_masses[@]}; do
        # for data_period in ${data_periods[@]}; do
            # ### temporary ###
            # if [[ ${reaction_id} == "D" && ${sexaquark_mass} == 2.01 && ${data_period} == "18qr" ]]; then continue; fi
            # if [[ ${reaction_id} == "E" && ${sexaquark_mass} == 2.01 && ${data_period} == "18qr" ]]; then continue; fi
            # if [[ ${reaction_id} == "H" && ${sexaquark_mass} == 2.01 && ${data_period} == "18qr" ]]; then continue; fi
            # ### temporary ###
            # export THE_PATH="${OUTPUT_DIR}/${reaction_id}${sexaquark_mass}_${data_period}_latest"
            # export THE_OPT="QA"
            # job_name="--job-name=QA_${reaction_id}${sexaquark_mass}_${data_period}"
            # sbatch ${job_name} ${SBATCH_SETUP} ${ANALYSIS_DIR}/scripts/mergers_single.sh
        # done
    # done
# done

# QA -- Bkg

# for data_period in ${data_periods[@]}; do
#     export THE_PATH="${OUTPUT_DIR}/A_bkg_${data_period}_latest"
#     export THE_OPT="QA"
#     job_name="--job-name=QA_bkg_${data_period}"
#     sbatch ${job_name} ${SBATCH_SETUP} ${ANALYSIS_DIR}/scripts/mergers_single.sh
# done

# V0s + Sexaquarks/KaonPairs

for reaction_id in ${reaction_channels[@]}; do
    for sexaquark_mass in ${sexaquark_masses[@]}; do
        for data_period in ${data_periods[@]}; do
            ### temporary ###
            if [[ ${reaction_id} == "D" && ${sexaquark_mass} == 2.01 && ${data_period} == "18qr" ]]; then continue; fi
            if [[ ${reaction_id} == "E" && ${sexaquark_mass} == 2.01 && ${data_period} == "18qr" ]]; then continue; fi
            if [[ ${reaction_id} == "H" && ${sexaquark_mass} == 2.01 && ${data_period} == "18qr" ]]; then continue; fi
            ### temporary ###
            export THE_PATH="${OUTPUT_DIR}/${reaction_id}${sexaquark_mass}_${data_period}_latest"
            export THE_OPT=${reaction_id}
            job_name="--job-name=${reaction_id}${sexaquark_mass}_${data_period}"
            sbatch ${job_name} ${SBATCH_SETUP} ${ANALYSIS_DIR}/scripts/mergers_single.sh
        done
    done
done
