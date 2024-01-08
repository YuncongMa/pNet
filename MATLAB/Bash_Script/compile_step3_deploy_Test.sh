#!/bin/sh

####################################################################################
# Yuncong Ma, 1/2/2022
# It is to compile deployFuncMvnmfL21p1_func_surf_hcp_single_gV.m or similar functions
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/compile_step3_deploy_Test.sh
####################################################################################

# For use as a function
parse()
{
    # Default setting
    function_name=deployFuncMvnmfL21p1_func_surf_hcp_single_gV_12
    file_log=/cbica/home/mayun/Projects/NMF/sge_output/step3/compile_step3_deploy.log

    while [ -n "$1" ];
    do
        case $1 in
            -name)
                function_name=$2;
            shift 2;;
            -log)
                file_log=$2;
            shift 2;;
            -*)
                echo "ERROR:no such option $1"
            exit;;
            *)
            break;;
        esac
    done
}

parse $*

dir_main=/cbica/home/mayun/Projects/NMF
dir_code=${dir_main}/scripts/MATLAB
dir_compile_function=${dir_main}/scripts/MATLAB
dir_MATLAB_package=/cbica/home/mayun/Projects/MATLAB~${dir_main}/Toolbox/Code_mvNMF_l21_ard_v3_debug
dir_output=${dir_code}

if test -f "${file_log}"
then
    rm $file_log
fi

# Run matlab 
jid=$(qsub -terse -j y -pe threaded 4 -l h_vmem=10g -o ${file_log} ${dir_code}/run_compile_matlab.sh -cf ${dir_compile_function} -dp ${dir_MATLAB_package} -fn ${function_name} -do ${dir_output});
