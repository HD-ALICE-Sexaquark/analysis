#!/bin/bash

####### ###### ##### #### ### ## #
#  Execute all DNs of an specific RN
 # However, when given an argument, it will execute only that DN
### ## #

if [[ -z ${ANALYSIS_DIR}    || -z ${SEXAQUARK_DIR} || -z ${OUTPUT_DIR} || -z ${SIMS_DIR} || -z ${RN}          ||
      -z ${RUN_DIR}         || -z ${MC_TYPE}       || -z ${IS_MC}      || -z ${DS}       || -z ${V0S_OPTION}  ||
      -z ${RC}              || -z ${SM}            || -z ${READ_LOGS}  || -z ${DO_QA}    || -z ${REWEIGHT_PT} ||
      -z ${REWEIGHT_RADIUS} ]]; then
    echo "run_single.sh :: make sure to first set all the options!"
    exit 1
fi

if [[ $# -eq 1 ]]; then
    DN_INI=$(printf "%d" $1)
    DN_END=$(printf "%d" $1)
else
    DN_INI=1
    DN_END=6
fi

for ((DN_i=${DN_INI}; DN_i<=${DN_END}; DN_i++)); do

    DN=$(printf "%03d" ${DN_i})

    # check if files exist
    if [[ ! -e ${RUN_DIR}/${DN}/AliESDs.root ||
          ! -e ${RUN_DIR}/${DN}/Kinematics.root ||
          ! -e ${RUN_DIR}/${DN}/galice.root ]]; then
        continue
    fi

    if [[ ! -e ${RUN_DIR}/${DN}/sim.log && ${IS_MC} == 1 && ${MC_TYPE} == "signal" ]]; then
        continue
    fi

    OUTPUT_SUBDIR=${OUTPUT_DIR}/sexaquark_${MC_TYPE}_${RC}${SM}/${RN}_${DN}
    if [[ $# -eq 1 && -e ${OUTPUT_SUBDIR} ]]; then
        rm -rfv ${OUTPUT_SUBDIR}
    fi
    mkdir -p ${OUTPUT_SUBDIR}

    cp ${SEXAQUARK_DIR}/AliAnalysisTaskSexaquark.cxx ${OUTPUT_SUBDIR}/
    cp ${SEXAQUARK_DIR}/AliAnalysisTaskSexaquark.h ${OUTPUT_SUBDIR}/
    cp ${SEXAQUARK_DIR}/AddSexaquark.C ${OUTPUT_SUBDIR}/
    cp ${SEXAQUARK_DIR}/runAnalysis.C ${OUTPUT_SUBDIR}/
    if [[ ${REWEIGHT_PT} == 1 ]]; then cp ${SEXAQUARK_DIR}/blastwave.root ${OUTPUT_SUBDIR}/; fi
    if [[ ${REWEIGHT_RADIUS} == 1 ]]; then cp ${SEXAQUARK_DIR}/radius-weights.root ${OUTPUT_SUBDIR}/; fi

    cd ${OUTPUT_SUBDIR}

    ln -s ${RUN_DIR}/${DN}/AliESDs.root
    ln -s ${RUN_DIR}/${DN}/Kinematics.root
    ln -s ${RUN_DIR}/${DN}/galice.root
    if [[ ${MC_TYPE} == "signal" ]]; then ln -s ${RUN_DIR}/${DN}/sim.log; fi

    aliroot -l -b -q 'runAnalysis.C('${IS_MC}', "'${DS}'", "'${V0S_OPTION}'", "'${RC}${SM}'", '${READ_LOGS}', '${DO_QA}', '${REWEIGHT_PT}', '${REWEIGHT_RADIUS}', 0)' &> analysis.slurm-${SLURM_JOB_ID}.log
done
