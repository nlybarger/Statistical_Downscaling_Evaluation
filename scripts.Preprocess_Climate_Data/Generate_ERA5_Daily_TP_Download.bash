#!/bin/bash

# Check if the template script exists
if [ ! -f "template_preproc_function_calls.py" ]; then
    echo "Error: Template script 'template_preproc_function_calls.py' not found!"
    exit 1
fi

# Loop through domain, exp, and aper options
for i in {1981..2016}; do
    # Create a filename based on parameters
    file_name="e5d_${i}.py"
            
            # Copy content from the template script to a new file
    cp "template_preproc_function_calls.py" "$file_name"

    sed -i "s/year = 1980/year = ${i}/" $file_name
    sed -i "s/^    # call_era5_download_daily_tp()/    call_era5_download_daily_tp()/" $file_name

    echo "Generated script: $file_name"
done


# Check if the template script exists
if [ ! -f "template_preproc_function_calls_batch.pbs" ]; then
    echo "Error: Template script 'template_preproc_function_calls_batch.pbs' not found!"
    exit 1
fi

for i in {1981..2016}; do
    # Create a filename based on parameters
    batch_name="e5d_batch_${i}.pbs"

    # Copy content from the template script to a new file
    cp "template_preproc_function_calls_batch.pbs" "$batch_name"

    sed -i "s/0_func_call/${i}_era5_download_tp/" $batch_name
    sed -i "s/func_call_0.py/e5d_${i}.py/" $batch_name
    sed -i "s/walltime=24:00:00/walltime=24:00:00/" $batch_name
    sed -i "s/mem=125GB/mem=20GB/" $batch_name
    # sed -i "s/conda activate geoanalysis/conda activate cdsapi/" $batch_name

    echo "Generated script: $batch_name"

    qsub $batch_name
done


