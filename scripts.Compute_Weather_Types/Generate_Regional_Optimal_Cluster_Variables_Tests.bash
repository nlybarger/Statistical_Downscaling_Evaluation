#!/bin/bash

# Check if the template script exists
if [ ! -f "template_WT_function_calls.py" ]; then
    echo "Error: Template script 'template_WT_function_calls.py' not found!"
    exit 1
fi

# Loop through domain, exp, and aper options
for i in {0..8}; do
    # Create a filename based on parameters
    file_name="ocv_$i.py"
            
            # Copy content from the template script to a new file
    cp "template_WT_function_calls.py" "$file_name"
            
            # Uncomment and replace DOMAIN, EXP, and APER lines with appropriate values
    sed -i "s/regind = 0/regind = $i/" $file_name
    sed -i "s/^    # call_Optimal_Cluster_Test()/    call_Optimal_Cluster_Test()/" $file_name
            
    echo "Generated script: $file_name"
done


# Check if the template script exists
if [ ! -f "template_WT_function_calls_batch.pbs" ]; then
    echo "Error: Template script 'template_WT_function_calls_batch.pbs' not found!"
    exit 1
fi

for i in {0..8}; do
    # Create a filename based on parameters
    batch_name="ocv_batch_$i.pbs"

            # Copy content from the template script to a new file
    cp "template_WT_function_calls_batch.pbs" "$batch_name"

            # Uncomment and replace DOMAIN, EXP, and APER lines with appropriate values
    sed -i "s/0_func_call/${i}_ocv_test/" $batch_name
    sed -i "s/func_call_0.py/ocv_$i.py/" $batch_name
    sed -i "s/ocv0/ocv$i/" $batch_name
    sed -i "s/mem=125GB/mem=10GB/" $batch_name
    sed -i "s/walltime=24:00:00/walltime=02:00:00/" $batch_name

    echo "Generated script: $batch_name"
    qsub $batch_name
done


