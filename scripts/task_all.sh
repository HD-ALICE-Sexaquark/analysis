#!/bin/bash

# USAGE: ./task_all.sh <CONFIG_FILE>

reaction_channels=("A")
sexaquark_masses=("1.73 1.8 1.87 1.94 2.01")

for reaction_channel in ${reaction_channels[@]}; do
    for sexaquark_mass in ${sexaquark_masses[@]}; do
        ./task_batch.sh $1 "kalman" "${reaction_channel}${sexaquark_mass}"
    done
done
