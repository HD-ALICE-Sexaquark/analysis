#!/bin/bash

#########################################
#                                       #
#  Master Script to send batch jobs     #
#    to analyze Sexaquark simulations   #
#                                       #
#########################################

# 80 run numbers
RN_ARRAY_2015=(246994 246991 246989 246984 246982 246980 246948 246945 246928 246851
               246847 246846 246845 246844 246810 246809 246808 246807 246805 246804
               246766 246765 246763 246760 246759 246758 246757 246751 246750 246495
               246493 246488 246487 246434 246431 246428 246424 246276 246275 246272
               246271 246225 246222 246217 246185 246182 246181 246180 246178 246153
               246152 246151 246148 246115 246113 246089 246087 246053 246052 246049
               246048 246042 246037 246036 246012 246003 246001 245963 245954 245952
               245949 245923 245833 245831 245829 245705 245702 245700 245692 245683
               )

# 91 run numbers
RN_ARRAY_2018=(297595 297590 297588 297558 297544 297542 297541 297540 297537 297512
               297483 297479 297452 297451 297450 297446 297442 297441 297415 297414
               297413 297406 297405 297380 297379 297372 297367 297366 297363 297336
               297335 297333 297332 297317 297315 297312 297311 297310 297278 297222
               297221 297218 297196 297195 297193 297133 297132 297129 297128 297124
               297123 297119 297118 297117 297085 297035 297031 296966 296941 296938
               296935 296934 296932 296931 296930 296903 296900 296899 296894 296852
               296851 296850 296848 296839 296838 296836 296835 296799 296794 296793
               296790 296787 296786 296785 296784 296781 296752 296694 296693 296691
               296690
               )

V0_OPT="true" # "true" "custom" "official"
V0_STR="TrueV0s" # "TrueV0s" "CustomV0s" "OfficialV0s"

NICE_OPT=0

for CURRENT_RN in ${RN_ARRAY_2018[@]}; do
    # create tmux session
    THIS_SESSION="${CURRENT_RN}_${V0_OPT}"
    tmux new -d -s ${THIS_SESSION}

    # prepare environment
    tmux send -t ${THIS_SESSION} "own-alice" ENTER
    tmux send -t ${THIS_SESSION} "source set_env.sh" ENTER

    # stage 1: run AliAnalysisTask
    # tmux send -t ${THIS_SESSION} 'nice -n '${NICE_OPT}' ./analyze.sh --sim signal+bkg --v0s '${V0_OPT}' --rn '${CURRENT_RN}'' ENTER

    # stage 2: run SexaquarkFinder
    ANALYSIS_DIR="output/signal+bkg/${CURRENT_RN}"
    for analysis_file in $(readlink -f ${ANALYSIS_DIR}/AnalysisResults_${V0_STR}*.root); do
        sexaquark_file="${analysis_file/AnalysisResults/SexaquarkResults}"
        tmux send -t ${THIS_SESSION} 'nice -n '${NICE_OPT}' root -l -b -q '\''SexaquarkFinder.cxx("'${analysis_file}'", "'${sexaquark_file}'")'\''' ENTER
    done

    # finish
    tmux send -t ${THIS_SESSION} "exit" ENTER
done
