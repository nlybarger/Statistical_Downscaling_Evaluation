#!/bin/bash

# Check if the template script exists
if [ ! -f "template_call_compute_conus_metrics.py" ]; then
    echo "Error: Template script 'template_call_compute_conus_metrics.py' not found!"
    exit 1
fi

# Loop through downscaling methods, models, scenarios, and time periods
# for i in {0..6}; do
for i in 6; do
    for j in {0..5}; do
    # for j in 3; do
        for k in {0..2}; do
        # for k in 0; do
            # historical period is only 1 time period
            if [ $k -eq 0 ]; then
                # Create a filename based on parameters
                file_name="c5_metcomp_${i}_${j}_${k}_0.py"
                
                # Copy content from the template script to a new file
                cp "template_call_compute_conus_metrics.py" "$file_name"
                
                sed -i "s/imeth = 0/imeth = $i/" $file_name
                sed -i "s/imod = 0/imod = $j/" $file_name
                sed -i "s/iscen = 0/iscen = $k/" $file_name
                sed -i "s/^    # compute_conus_metric_maps_cmip5_ds()/    compute_conus_metric_maps_cmip5_ds()/" $file_name

                echo "Generated script: $file_name"
            else
                # rcp45 and rcp85 have 3 time periods each
                for l in {0..2}; do
                    # Create a filename based on parameters
                    file_name="c5_metcomp_${i}_${j}_${k}_${l}.py"
                    
                    # Copy content from the template script to a new file
                    cp "template_call_compute_conus_metrics.py" "$file_name"
                    
                    sed -i "s/imeth = 0/imeth = $i/" $file_name
                    sed -i "s/imod = 0/imod = $j/" $file_name
                    sed -i "s/iscen = 0/iscen = $k/" $file_name
                    sed -i "s/itim = 0/itim = $l/" $file_name
                    sed -i "s/^    # compute_conus_metric_maps_cmip5_ds()/    compute_conus_metric_maps_cmip5_ds()/" $file_name

                    echo "Generated script: $file_name"
                done
            fi
        done
    done
done

# Check if the template script exists
if [ ! -f "template_call_compute_conus_metrics_batch.pbs" ]; then
    echo "Error: Template script 'template_call_compute_conus_metrics_batch.pbs' not found!"
    exit 1
fi

# for i in {0..6}; do
for i in 6; do
    for j in {0..5}; do
    # for j in 3; do
        for k in {0..2}; do
        # for k in 0; do
            if [ $k -eq 0 ]; then
                # Create a filename based on parameters
                batch_name="c5_metcomp_batch_${i}_${j}_${k}_0.pbs"

                # Copy content from the template script to a new file
                cp "template_call_compute_conus_metrics_batch.pbs" "$batch_name"

                sed -i "s/0_func_call/${i}_${j}_${k}_0_c5_metcomp/" $batch_name
                sed -i "s/func_call_0.py/c5_metcomp_${i}_${j}_${k}_0.py/" $batch_name
                sed -i "s/walltime=24:00:00/walltime=06:00:00/" $batch_name
                # MACA in particular needs a ton of memory
                sed -i "s/mem=125GB/mem=275GB/" $batch_name

                qsub $batch_name
                echo "Generated script: $batch_name"
            else
                for l in {0..2}; do
                    # Create a filename based on parameters
                    batch_name="c5_metcomp_batch_${i}_${j}_${k}_${l}.pbs"

                    # Copy content from the template script to a new file
                    cp "template_call_compute_conus_metrics_batch.pbs" "$batch_name"

                    sed -i "s/0_func_call/${i}_${j}_${k}_${l}_c5_metcomp/" $batch_name
                    sed -i "s/func_call_0.py/c5_metcomp_${i}_${j}_${k}_${l}.py/" $batch_name
                    sed -i "s/walltime=24:00:00/walltime=06:00:00/" $batch_name
                    # MACA in particular needs a ton of memory
                    sed -i "s/mem=125GB/mem=275GB/" $batch_name

                    qsub $batch_name
                    echo "Generated script: $batch_name"
                done
            fi
        done
    done
done