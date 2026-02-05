#!/bin/bash

# Check if the template script exists
if [ ! -f "template_preproc_function_calls.py" ]; then
    echo "Error: Template script 'template_preproc_function_calls.py' not found!"
    exit 1
fi

# Loop through domain, exp, and aper options
# for i in {1..6}; do
for i in 6; do
    # for j in {1981..2016}; do
    for j in 2000; do
        
    # Create a filename based on parameters
        file_name="obs_icargrid_${i}_${j}.py"
            
        # Copy content from the template script to a new file
        cp "template_preproc_function_calls.py" "$file_name"
            
        sed -i "s/obsind = 0/obsind = $i/" $file_name
        sed -i "s/yr0 = 1980/yr0 = $j/" $file_name
        sed -i "s/^    # call_preproc_obs()/    call_preproc_obs()/" $file_name
    
        echo "Generated script: $file_name"
    done
done

# Check if the template script exists
if [ ! -f "template_preproc_function_calls_batch.pbs" ]; then
    echo "Error: Template script 'template_preproc_function_calls_batch.pbs' not found!"
    exit 1
fi

# for i in {1..6}; do
for i in 6; do
    rm ./${i}*_obs_icg.o*
    # for j in {1981..2016}; do
    for j in 2000; do
        # Create a filename based on parameters
        batch_name="obs_icg_batch_${i}_${j}.pbs"

        # Copy content from the template script to a new file
        cp "template_preproc_function_calls_batch.pbs" "$batch_name"

        sed -i "s/0_func_call/${i}_${j}_obs_icg/" $batch_name
        sed -i "s/func_call_0.py/obs_icargrid_${i}_${j}.py/" $batch_name
        sed -i "s/walltime=24:00:00/walltime=08:00:00/" $batch_name
        sed -i "s/mem=125GB/mem=54GB/" $batch_name

        qsub $batch_name
        echo "Generated script: $batch_name"
    done
done

