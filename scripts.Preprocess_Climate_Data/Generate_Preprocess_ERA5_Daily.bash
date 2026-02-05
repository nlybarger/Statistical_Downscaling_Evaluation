#!/bin/bash

# Check if the template script exists
if [ ! -f "template_preproc_function_calls.py" ]; then
    echo "Error: Template script 'template_preproc_function_calls.py' not found!"
    exit 1
fi

# Loop through domain, exp, and aper options
for iy in {1981..2016}; do
    for im in {1..12}; do
        im_padded=$(printf "%02d" $im)
        # Create a filename based on parameters
        file_name="e5d_${iy}_${im_padded}.py"

        # Copy content from the template script to a new file
        cp "template_preproc_function_calls.py" "$file_name"
        sed -i "s/year = 1980/year = ${iy}/" $file_name
        sed -i "s/month = 1/month = ${im}/" $file_name
        sed -i "s/^    # call_preproc_era5()/    call_preproc_era5()/" $file_name

        echo "Generated script: $file_name"
    done
done


# Check if the template script exists
if [ ! -f "template_preproc_function_calls_batch.pbs" ]; then
    echo "Error: Template script 'template_preproc_function_calls_batch.pbs' not found!"
    exit 1
fi

for iy in {1981..2016}; do
    for im in {1..12}; do
        # Create a filename based on parameters
        im_padded=$(printf "%02d" $im)
        batch_name="e5d_batch_${iy}_${im_padded}.pbs"

        # Copy content from the template script to a new file
        cp "template_preproc_function_calls_batch.pbs" "$batch_name"

        sed -i "s/0_func_call/${iy}_${im_padded}_era5_concat_daily/" $batch_name
        sed -i "s/func_call_0.py/e5d_${iy}_${im_padded}.py/" $batch_name
        sed -i "s/walltime=24:00:00/walltime=04:00:00/" $batch_name
        sed -i "s/mem=125GB/mem=25GB/" $batch_name

        echo "Generated script: $batch_name"
        qsub $batch_name
    done
done
