#!/bin/bash

# Check if the template script exists
if [ ! -f "template_preproc_function_calls.py" ]; then
    echo "Error: Template script 'template_preproc_function_calls.py' not found!"
    exit 1
fi

# Loop through domain, exp, and aper options
for i in {0..5}; do
# for i in 0 2 3 4; do
    # Create a filename based on parameters
    file_name="ppc5_$i.py"
            
            # Copy content from the template script to a new file
    cp "template_preproc_function_calls.py" "$file_name"
            
            # Uncomment and replace DOMAIN, EXP, and APER lines with appropriate values
    sed -i "s/imod = 0/imod = $i/" $file_name
    sed -i "s/^    # call_preproc_cmip5()/    call_preproc_cmip5()/" $file_name

            
    echo "Generated script: $file_name"
done


# Check if the template script exists
if [ ! -f "template_preproc_function_calls_batch.pbs" ]; then
    echo "Error: Template script 'template_preproc_function_calls_batch.pbs' not found!"
    exit 1
fi

for i in {0..5}; do
# for i in 2; do
    # Create a filename based on parameters
    batch_name="ppc5_batch_$i.pbs"

            # Copy content from the template script to a new file
    cp "template_preproc_function_calls_batch.pbs" "$batch_name"

            # Uncomment and replace DOMAIN, EXP, and APER lines with appropriate values
    sed -i "s/preproc_CMIP5_0/preproc_CMIP5_$i/" $batch_name
    sed -i "s/ppc5_0.py/ppc5_$i.py/" $batch_name


    sed -i "s/0_func_call/${i}_preproc_CMIP5/" $batch_name
    sed -i "s/func_call_0.py/ppc5_${i}.py/" $batch_name
    sed -i "s/walltime=24:00:00/walltime=03:00:00/" $batch_name
    sed -i "s/mem=125GB/mem=200GB/" $batch_name

    qsub $batch_name
    echo "Generated script: $batch_name"

done


