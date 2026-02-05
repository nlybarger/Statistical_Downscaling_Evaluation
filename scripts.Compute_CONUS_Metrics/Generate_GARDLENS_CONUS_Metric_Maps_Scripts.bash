#!/bin/bash

# Check if the template script exists
if [ ! -f "template_call_compute_conus_metrics.py" ]; then
    echo "Error: Template script 'template_call_compute_conus_metrics.py' not found!"
    exit 1
fi

# Loop through all GARDLENS ensemble members
for i in {0..99}; do
    file_name="GARDLENS_metcomp_${i}.py"
    # Copy content from the template script to a new file
    cp "template_call_compute_conus_metrics.py" "$file_name"
    
    sed -i "s/ensind = 0/ensind = ${i}/" $file_name
    sed -i "s/^    # compute_conus_metric_maps_gardlens()/    compute_conus_metric_maps_gardlens()/" $file_name

    echo "Generated script: $file_name"
done

# Check if the template script exists
if [ ! -f "template_call_compute_conus_metrics_batch.pbs" ]; then
    echo "Error: Template script 'template_call_compute_conus_metrics_batch.pbs' not found!"
    exit 1
fi

for i in {0..99}; do
    batch_name="GARDLENS_metcomp_batch_${i}.pbs"

    # Copy content from the template script to a new file
    cp "template_call_compute_conus_metrics_batch.pbs" "$batch_name"

    ###== Replace default values with custom parameters
    # Batch job name
    sed -i "s/0_func_call/${i}_GARDLENS_metcomp/" $batch_name
    # Python script to call
    sed -i "s/func_call_0.py/GARDLENS_metcomp_${i}.py/" $batch_name
    # Resources
    sed -i "s/walltime=24:00:00/walltime=06:00:00/" $batch_name
    sed -i "s/mem=125GB/mem=125GB/" $batch_name

    echo "Generated script: $batch_name"
    # qsub $batch_name
done
