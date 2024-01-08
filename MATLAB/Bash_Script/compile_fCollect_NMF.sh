#!/bin/sh

####################################################################################
# Yuncong Ma, 10/31/2022
# It is to compile fCollect_NMF.m
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/compile_fCollect_NMF.sh
####################################################################################

dir_main=/cbica/home/mayun/Projects/NMF
dir_code=${dir_main}/scripts/MATLAB
dir_compile_function=${dir_main}/scripts/MATLAB
dir_MATLAB_package=/cbica/home/mayun/Projects/MATLAB'~'${dir_main}/scripts/MATLAB
function_name=fCollect_NMF
dir_output=${dir_code}

file_log=${dir_main}/sge_output/collect/compile_fCollect_NMF.log
if test -f "${file_log}"
then
    rm $file_log
fi

# Run matlab 
jid=$(qsub -terse -j y -pe threaded 4 -l h_vmem=4g -o ${file_log} ${dir_code}/run_compile_matlab.sh -cf ${dir_compile_function} -dp ${dir_MATLAB_package} -fn ${function_name} -do ${dir_output});
