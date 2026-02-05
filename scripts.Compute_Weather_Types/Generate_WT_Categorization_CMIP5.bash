#!/bin/bash

# Check if the template script exists
if [ ! -f "template_WT_function_calls.py" ]; then
    echo "Error: Template script 'template_WT_function_calls.py' not found!"
    exit 1
fi

# Loop through domain, exp, and aper options
for i in {0..8}; do
    for j in {0..5}; do
    # for j in 0 1 3 4 5; do
        # Create a filename based on parameters
        file_name="c5_cwtpr_${i}_${j}.py"
                
        # Copy content from the template script to a new file
        cp "template_WT_function_calls.py" "$file_name"

        sed -i "s/regind = 0/regind = $i/" $file_name
        sed -i "s/imod = 0/imod = $j/" $file_name
        sed -i "s/^    # call_wtpr_cmip5()/    call_wtpr_cmip5()/" $file_name
                
        echo "Generated script: $file_name"
    done
done


# Check if the template script exists
if [ ! -f "template_WT_function_calls_batch.pbs" ]; then
    echo "Error: Template script 'template_WT_function_calls_batch.pbs' not found!"
    exit 1
fi

for i in {0..8}; do
    for j in {0..5}; do
    # for j in 0 1 3 4 5; do
        # Create a filename based on parameters
        batch_name="c5_cwtpr_batch_${i}_${j}.pbs"

        # Copy content from the template script to a new file
        cp "template_WT_function_calls_batch.pbs" "$batch_name"

        sed -i "s/0_func_call/${i}_${j}_cmip5_cwtpr/" $batch_name
        sed -i "s/func_call_0.py/c5_cwtpr_${i}_${j}.py/" $batch_name
        sed -i "s/walltime=24:00:00/walltime=08:00:00/" $batch_name
        sed -i "s/mem=125GB/mem=200GB/" $batch_name

        echo "Generated script: $batch_name"

        qsub $batch_name
    done
done
