#!/bin/sh

####################################################################################
# Yuncong Ma, 10/10/2022
# It is to compile fVisualize_NMF_Step3_HCP_FS_single.m
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/compile_step3_visualize.sh
####################################################################################

dir_main=/cbica/home/mayun/Projects/NMF
dir_code=${dir_main}/scripts/MATLAB
dir_compile_function=${dir_code}
dir_MATLAB_package=/cbica/home/mayun/Projects/MATLAB'~'${dir_main}/Toolbox/Code_mvNMF_l21_ard_v3_debug
function_name=fVisualize_NMF_Step3_HCP_FS_single
dir_output=${dir_code}

file_log=${dir_main}/sge_output/step3/compile_step3_visualize.log
if test -f "${file_log}"
then
    rm $file_log
fi

# Run matlab 
jid=$(qsub -terse -j y -pe threaded 4 -l h_vmem=10g -o ${file_log} ${dir_code}/run_compile_matlab_figure.sh -cf ${dir_compile_function} -dp ${dir_MATLAB_package} -fn ${function_name} -do ${dir_output});