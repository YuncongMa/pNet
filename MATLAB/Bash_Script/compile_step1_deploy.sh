#!/bin/sh

####################################################################################
# Yuncong Ma, 12/12/2022
# It is to compile deployFuncInit_surf_hcp_tp.m
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/compile_step1_deploy.sh
####################################################################################

dir_main=/cbica/home/mayun/Projects/NMF
dir_code=${dir_main}/scripts/MATLAB
dir_compile_function=${dir_main}/scripts/MATLAB
dir_MATLAB_package=/cbica/home/mayun/Projects/MATLAB~${dir_main}/Toolbox/Code_mvNMF_l21_ard_v3_debug~${dir_main}/Toolbox/Two_Factor_Clustering
function_name=deployFuncInit_surf_hcp_tp
dir_output=${dir_code}

file_log=${dir_main}/sge_output/step1/compile_step1_deploy.log
if test -f "${file_log}"
then
    rm $file_log
fi

# Run matlab 
jid=$(qsub -terse -j y -pe threaded 4 -l h_vmem=10g -o ${file_log} ${dir_code}/run_compile_matlab.sh -cf ${dir_compile_function} -dp ${dir_MATLAB_package} -fn ${function_name} -do ${dir_output});
