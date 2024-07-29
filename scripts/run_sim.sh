#!/bin/bash

#SBATCH --partition=main
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=3000

if [[ -z ${MC_TYPE} || -z ${RC} || -z ${SM} || -z ${RN} || -z ${CURRENT_DIR} || -z ${RUN_DIR} || -z ${DS} ]]; then
    echo "run_sim.sh :: one or more env. vars. are undefined!!"
    exit 1
fi

V0S_OPTION="kalman"

MC_SUFFIX="sig"
ANALYSIS_DIR=${CURRENT_DIR}/../sexaquark
OUTPUT_DIR=${ANALYSIS_DIR}/output_${MC_SUFFIX}_${RC}${SM}
READ_LOGS=1
mkdir -p ${OUTPUT_DIR}

for DN in {001..006}; do # CONTROL N DIRS

    if [[ ! -e ${RUN_DIR}/${DN}/AliESDs.root ||
        ! -e ${RUN_DIR}/${DN}/Kinematics.root ||
        ! -e ${RUN_DIR}/${DN}/galice.root ]]; then
        continue
    fi

    if [[ ! -e ${RUN_DIR}/${DN}/sim.log && ${MC_TYPE} == "signal" ]]; then
        continue
    fi

    OUTPUT_SUBDIR=${OUTPUT_DIR}/${RN}_${DN}
    mkdir -p ${OUTPUT_SUBDIR}

    cp ${ANALYSIS_DIR}/AliAnalysisTaskSexaquark.cxx ${OUTPUT_SUBDIR}/
    cp ${ANALYSIS_DIR}/AliAnalysisTaskSexaquark.h ${OUTPUT_SUBDIR}/
    cp ${ANALYSIS_DIR}/AddSexaquark.C ${OUTPUT_SUBDIR}/
    cp ${ANALYSIS_DIR}/runAnalysis.C ${OUTPUT_SUBDIR}/
    cp ${ANALYSIS_DIR}/blastwave.root ${OUTPUT_SUBDIR}/
    cp ${ANALYSIS_DIR}/radius-weights.root ${OUTPUT_SUBDIR}/

    cd ${OUTPUT_SUBDIR}

    ln -s ${RUN_DIR}/${DN}/AliESDs.root
    ln -s ${RUN_DIR}/${DN}/Kinematics.root
    ln -s ${RUN_DIR}/${DN}/galice.root
    if [[ ${MC_TYPE} == "signal" ]]; then  ln -s ${RUN_DIR}/${DN}/sim.log; fi

    aliroot -l -b -q 'runAnalysis.C(1, "'${DS}'", "'${V0S_OPTION}'", "'${RC}${SM}'", '${READ_LOGS}', 1, 1, 1, 0)'
done

cd ${CURRENT_DIR}
