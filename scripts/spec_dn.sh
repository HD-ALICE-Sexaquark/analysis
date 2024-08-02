#!/bin/bash

#### ### ## #
# Execute an specific RN and DN
## #

# hardcoded options
IS_MC=1
MC_TYPE="signal"
V0S_OPTION="kalman"
RC="A"
SM=1.8
READ_LOGS=1
DO_QA=1
REWEIGHT_PT=1
REWEIGHT_RADIUS=1

if [[ $# -ne 1 ]]; then
    echo "spec_dn.sh :: usage: ./spec_dn.sh <RN>_<DN>"
    exit 1
fi

RN=${1%_*}
DN=${1#*_}

if [[ -z ${ANALYSIS_DIR} || -z ${SEXAQUARK_DIR} || -z ${DEBUG_DIR} || -z ${SIMS_DIR} ]]; then
    echo "spec_dn.sh :: make sure to first set:"
    echo "spec_dn.sh :: - ANALYSIS_DIR  : where the main repository is"
    echo "spec_dn.sh :: - SEXAQUARK_DIR : where the AliAnalysisTaskSexaquark is"
    echo "spec_dn.sh :: - DEBUG_DIR     : where the debug output will be stored"
    echo "spec_dn.sh :: - SIMS_DIR      : where the simulations are"
    exit 1
fi

OUTPUT_SUBDIR=${DEBUG_DIR}/sexaquark_${MC_TYPE}_${RC}${SM}/${RN}_${DN}
if [[ -e ${OUTPUT_SUBDIR} ]]; then
    rm -rfv ${OUTPUT_SUBDIR}
fi
mkdir ${OUTPUT_SUBDIR}

# determine DATASET based on RN, via ripgrep
DS_FILE=$(rg -l ${RN} ${ANALYSIS_DIR}/doc/)
DS_FILE=$(basename ${DS_FILE})
DS=${DS_FILE/%_*}

# determine SIMSET and RUN_DIR
if [[ ${MC_TYPE} == "signal" ]]; then
    if [[ ${DS:0:6} == "LHC15o" ]]; then SS="LHC23l1b3"
    elif [[ ${DS:0:6} == "LHC18q" || ${DS:0:6} == "LHC18r" ]]; then SS="LHC23l1a3"; fi
    RUN_DIR=${SIMS_DIR}/${SS}/${RC}${SM}/${RN}
elif [[ ${MC_TYPE} == "bkg" ]]; then
    if [[ ${DS:0:6} == "LHC15o" ]]; then SS="LHC20j6a"
    elif [[ ${DS:0:6} == "LHC18q" || ${DS:0:6} == "LHC18r" ]]; then SS="LHC20e3a"; fi
    RUN_DIR=${SIMS_DIR}/${SS}/${RN}
fi

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

aliroot -l -b -q 'runAnalysis.C('${IS_MC}', "'${DS}'", "'${V0S_OPTION}'", "'${RC}${SM}'", '${READ_LOGS}', '${DO_QA}', '${REWEIGHT_PT}', '${REWEIGHT_RADIUS}', 0)' 2>&1 | tee analysis.log
