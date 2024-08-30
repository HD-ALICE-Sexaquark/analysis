#!/bin/bash

# usage:
#   ./sanitize_local_files.sh /home/ceres/borquez/some/sims/LHC23l1a3
#   ./sanitize_local_files.sh /home/ceres/borquez/some/sims/LHC23l1b3

# get input path from first argument
input_path=$1

# get the list of directories in the input path
directories=$(ls -d ${input_path}/*/*/*/)

# loop through each directory
for directory in ${directories}; do
    # get the SIM_SET, RUN_NUMBER and DIR_NUMBER from the directory name
    sim_set=$(basename $(dirname $(dirname ${directory})))
    run_number=$(basename $(dirname ${directory}))
    dir_number=$(basename ${directory})

    # check if the directory has the files AliESDs.root, Kinematics.root, galice.root and sim.log
    necessary_files=("AliESDs.root" "Kinematics.root" "galice.root" "sim.log")
    echo "${sim_set} ${run_number} ${dir_number}"
    for file in ${necessary_files[@]}; do
        if [[ ! -e $directory/$file ]]; then
            echo ">> ${file} not found"
            production_name=$(basename ${input_path})
            alien_path=/alice/sim/2023/${production_name}/${sim_set}/${run_number}/${dir_number}
            alien.py cp alien://${alien_path}/${file} file://${directory}/${file}
        fi
    done

    if [[ ! -e ${directory}/AliESDs.root || ! -e ${directory}/Kinematics.root || ! -e ${directory}/galice.root || ! -e ${directory}/sim.log ]]; then
        echo "${sim_set} ${run_number} ${dir_number} : tried to copy all files but some are still missing, deleting dir"
        rm -rf ${directory}
    else
        echo "${sim_set} ${run_number} ${dir_number} : it's fine now!"
    fi
done
