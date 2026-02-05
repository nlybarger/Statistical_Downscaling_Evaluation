#!/bin/bash

# Check if the template script exists
if [ ! -f "template_WT_function_calls.py" ]; then
    echo "Error: Template script 'template_WT_function_calls.py' not found!"
    exit 1
fi

# Loop through domain, exp, and aper options
for i in {0..8}; do
    for j in {0..6}; do
    # for j in 7; do
        # Create a filename based on parameters
        file_name="obs_cwtpr_${i}_${j}.py"

        # Copy content from the template script to a new file
        cp "template_WT_function_calls.py" "$file_name"

        sed -i "s/regind = 0/regind = $i/" $file_name
        sed -i "s/obsind = 0/obsind = $j/" $file_name
        sed -i "s/^    # call_wtpr_obs()/    call_wtpr_obs()/" $file_name
                
        echo "Generated script: $file_name"
    done
done


# Check if the template script exists
if [ ! -f "template_WT_function_calls_batch.pbs" ]; then
    echo "Error: Template script 'template_WT_function_calls_batch.pbs' not found!"
    exit 1
fi

for i in {0..8}; do
    for j in {0..6}; do
    # for j in 7; do
    # for j in 0 1 2 3 4 5 7; do
        # Create a filename based on parameters
        batch_name="obs_cwtpr_batch_${i}_${j}.pbs"

        # Copy content from the template script to a new file
        cp "template_WT_function_calls_batch.pbs" "$batch_name"

        sed -i "s/0_func_call/${i}_${j}_obs_cwtpr/" $batch_name
        sed -i "s/func_call_0.py/obs_cwtpr_${i}_${j}.py/" $batch_name
        sed -i "s/walltime=24:00:00/walltime=06:00:00/" $batch_name
        sed -i "s/mem=125GB/mem=150GB/" $batch_name

        echo "Generated script: $batch_name"

        qsub $batch_name
    done
done
