#!/bin/bash

# Check out the output

if [[ -z ${OUTPUT_DIR} ]]; then
    echo "cout_all.sh :: make sure to first set OUTPUT_DIR"
    exit 1
fi

summarize_log() {
    local log_file=$1
    ln_seg_fault=$(grep -n "segmentation violation" ${1} | grep -o '^[0-9]*')
    ln_hint=$(grep -n "The lines below might hint at the cause of the crash" ${1} | grep -o '^[0-9]*')

    # read every file line by line, from ln_hint to EOF
    output=$(tail -n +${ln_hint} ${1} | grep "AnalysisTaskSexaquark" | awk -F '[ :()]+' '
             {
                 if (NR == 1) {
                     printf "%s:%s", $5, $NF
                 } else {
                     printf ", %s:%s\n", $5, $NF
                 }
             }
             ')
    echo $output
}

N_TOTAL_DN=$(ls -1 ${OUTPUT_DIR}/sexaquark_signal_A1.8 | wc -l)
N_TOTAL_AR=$(ls -lrth ${OUTPUT_DIR}/sexaquark_signal_A1.8/*/AnalysisResults.root | wc -l)

N_FAILED=0
for dir in ${OUTPUT_DIR}/sexaquark_signal_A1.8/*/; do
    if [[ ! -e ${dir}/AnalysisResults.root ]]; then
        log_file=$(ls -1 ${dir}/analysis.slurm-*.log)
        if [[ N_FAILED -eq 0 ]]; then
            echo "$(basename $dir) $(summarize_log $log_file)" > failed_dirs.txt
        else
            echo "$(basename $dir) $(summarize_log $log_file)" >> failed_dirs.txt
        fi
        N_FAILED=$((N_FAILED + 1))
    fi
done

echo ""
echo "Summary of Output"
echo "================="
echo ">> N Successful Dir Numbers = ${N_TOTAL_AR}"
echo ">> N Failed Dir Numbers     = ${N_FAILED}"
echo ">> N Total Dir Numbers      = ${N_TOTAL_DN}"
echo ""
