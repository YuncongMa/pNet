#!/bin/sh

####################################################################################
# Yuncong Ma, 1/9/2023
# It is to compile fNMF_QC_SC_FH.m
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/compile_fNMF_QC_SC_FH.sh
####################################################################################

dir_main=/cbica/home/mayun/Projects/NMF
dir_code=${dir_main}/scripts/MATLAB
dir_compile_function=${dir_main}/scripts/MATLAB
dir_MATLAB_package=/cbica/home/mayun/Projects/MATLAB'~'${dir_main}/scripts/MATLAB
function_name=fNMF_QC_SC_FH
dir_output=${dir_code}

file_log=${dir_main}/sge_output/collect/compile_fNMF_QC_SC_FH.log
if test -f "${file_log}"
then
    rm $file_log
fi

# Run matlab 
jid=$(qsub -terse -j y -pe threaded 4 -l h_vmem=10g -o ${file_log} ${dir_code}/run_compile_matlab.sh -cf ${dir_compile_function} -dp ${dir_MATLAB_package} -fn ${function_name} -do ${dir_output});