#!/bin/bash

# Check if the template script exists
if [ ! -f "template_preproc_function_calls.py" ]; then
    echo "Error: Template script 'template_preproc_function_calls.py' not found!"
    exit 1
fi

# Loop through domain, exp, and aper options
for i in {0..99}; do
# for i in 0; do
   # for j in 0 1; do
   i_padded=$(printf "%02d" $i)
   # Create a filename based on parameters
   file_name="preproc_cesm2_${i_padded}.py"
   
   # Copy content from the template script to a new file
   cp "template_preproc_function_calls.py" "$file_name"

   sed -i "s/ensind = 0/ensind = ${i}/" $file_name
   sed -i "s/^    # call_preproc_CESM2_LE()/    call_preproc_CESM2_LE()/" $file_name
   # sed -i "s/^    # call_combine_cesm2le_forWT()/    call_combine_cesm2le_forWT()/" $file_name

   echo "Generated script: $file_name"
   # done
done

# Check if the template script exists
if [ ! -f "template_preproc_function_calls_batch.pbs" ]; then
    echo "Error: Template script 'template_preproc_function_calls_batch.pbs' not found!"
    exit 1
fi

for i in {0..99}; do
# for i in 0; do
   # for j in 0 1; do
   i_padded=$(printf "%02d" $i)
   # Create a filename based on parameters
   batch_name="preproc_cesm2_batch_${i_padded}.pbs"

   # Copy content from the template script to a new file
   cp "template_preproc_function_calls_batch.pbs" "$batch_name"

   sed -i "s/0_func_call/${i_padded}_preproc_CESM2/" $batch_name
   sed -i "s/func_call_0.py/preproc_cesm2_${i_padded}.py/" $batch_name
   sed -i "s/walltime=24:00:00/walltime=00:30:00/" $batch_name
   sed -i "s/mem=125GB/mem=20GB/" $batch_name

   echo "Generated script: $batch_name"
   qsub $batch_name

done

