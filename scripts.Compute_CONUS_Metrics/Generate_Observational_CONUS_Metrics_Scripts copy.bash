#!/bin/bash

# Check if the template script exists
if [ ! -f "template_call_compute_conus_metrics.py" ]; then
    echo "Error: Template script 'template_call_compute_conus_metrics.py' not found!"
    exit 1
fi

for i in {0..6}; do
# for i in 0; do
    # Create a filename based on parameters
    file_name="obs_metcomp_${i}.py"
            
    # Copy content from the template script to a new file
    cp "template_call_compute_conus_metrics.py" "$file_name"
            
    sed -i "s/obsind = 0/obsind = $i/" $file_name
    sed -i "s/^    # compute_conus_metric_maps_obs()/    compute_conus_metric_maps_obs()/" $file_name

    echo "Generated script: $file_name"
done

# Check if the template script exists
if [ ! -f "template_call_compute_conus_metrics_batch.pbs" ]; then
    echo "Error: Template script 'template_call_compute_conus_metrics_batch.pbs' not found!"
    exit 1
fi

for i in {0..6}; do
# for i in 0; do
    # rm ${i}_obs_metcomp.o*
    # Create a filename based on parameters
    batch_name="obs_metcomp_batch_${i}.pbs"

    # Copy content from the template script to a new file
    cp "template_call_compute_conus_metrics_batch.pbs" "$batch_name"

    sed -i "s/0_func_call/${i}_obs_metcomp/" $batch_name
    sed -i "s/func_call_0.py/obs_metcomp_${i}.py/" $batch_name
    sed -i "s/walltime=24:00:00/walltime=06:00:00/" $batch_name
    sed -i "s/mem=125GB/mem=125GB/" $batch_name

    qsub $batch_name
    echo "Generated script: $batch_name"
done

# qsub obs_metcomp_batch_0.pbs

