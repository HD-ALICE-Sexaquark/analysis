#!/bin/bash

if [[ -z ${OUTPUT_DIR} ]]; then
    echo "seff_all.sh :: make sure to first set OUTPUT_DIR"
    exit 1
fi

convert_to_mb() {
    local size=$1
    echo $size | numfmt --from=iec --to-unit=1048576
}

# get all job ids
JOBIDS=( $(find ${OUTPUT_DIR} -name "slurm-*.out" | grep -oP 'slurm-\K[0-9]+') )

# get the max and avg memory and time
MAX_MEM=0
AVG_MEM=0

for JOBID in ${JOBIDS[@]}; do

    SEFF_OUT=$(seff ${JOBID})
    MEM=$(echo ${SEFF_OUT} | grep -oP 'Memory Utilized: \K[0-9.]+ [A-Z]+' | sed 's/ //' | sed 's/B//')
    MEM=$(convert_to_mb $MEM)

    if [[ ${MEM} -gt ${MAX_MEM} ]]; then
        MAX_MEM=${MEM}
    fi

    AVG_MEM=$((AVG_MEM + MEM))
done

echo "Max memory: $MAX_MEM MB"
echo "Avg memory: $(($AVG_MEM / ${#JOBIDS[@]})) MB"
