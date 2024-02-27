#!/bin/bash

# Analyze gen. purpose sim. with AliAnalysisTaskLambda1520Lpipi

TOPDIR=/home/ceres/borquez/some/sims

SIM_SET=("LHC20e3a")
# RUN_NUMBERS=(295585 295586 295588 295589 295610)
RUN_NUMBERS=(295585)
DATA_SET=("LHC18q")

mkdir -p test_gp
cp AliAnalysisTaskLambda1520Lpipi.cxx test_gp/
cp AliAnalysisTaskLambda1520Lpipi.h test_gp/
cp AddTask_Lambda1520Lpipi.C test_gp/
cp runAnalysis.C test_gp/

cd test_gp

for SS in ${SIM_SET[@]}; do
    for ((i=0; i<${#RUN_NUMBERS[@]}; i++)); do

        RN=${RUN_NUMBERS[i]}
        DS=${DATA_SET[i]}

        RUNDIR=${TOPDIR}/${SS}/${RN}

        for SIMDIR in ${RUNDIR}/*/; do
        # SIMDIR=${RUNDIR}/001

            DN=$(basename ${SIMDIR})

            if [[ ! -e ${SIMDIR}/AliESDs.root ||
                ! -e ${SIMDIR}/Kinematics.root ||
                ! -e ${SIMDIR}/galice.root ]]; then
                continue
            fi

            ln -s ${SIMDIR}/AliESDs.root
            ln -s ${SIMDIR}/Kinematics.root
            ln -s ${SIMDIR}/galice.root

            aliroot -l -b -q 'runAnalysis.C(1, "'${DS}'", 1)' 2>&1 | tee analysis.log

            mv -v AnalysisResults.root AnalysisResults_${RN}_${DN}.root
            mv -v analysis.log analysis_${RN}_${DN}.log

            rm -v AliESDs.root Kinematics.root galice.root

        done
    done
done

rm -v AliAnalysisTaskLambda1520Lpipi*
rm -v AddTask_Lambda1520Lpipi.C runAnalysis.C

cd ..
