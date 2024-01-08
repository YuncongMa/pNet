#!/bin/sh

####################################################################################
# Yuncong Ma, 10/15/2022
# It is to generate figures for robust initilization results stored in init.mat
# It uses Yuncong's MATLAB package to do visualization
# bash /cbica/home/mayun/Projects/NMF/scripts/PNC_GS/visualize_step2_All_HCP.sh
####################################################################################

# For use as a function
parse()
{
    # Default setting
    K=3
    dataset=HCP/HCP_1200
    subset=/All
    dir_log=/cbica/home/mayun/Projects/NMF/sge_output/step1
    N_Fold=2
    fold=1
    file_log_main=/cbica/home/mayun/Projects/NMF/sge_output/step1/submit_step1_All_HCP_FS.log
    
    while [ -n "$1" ];
    do
        case $1 in
            -k)
                K=$2;
            shift 2;;
            -ds)
                dataset=$2;
            shift 2;;
            -ss)
                subset=$2;
            shift 2;;
            -nf)
                N_Fold=$2;
            shift 2;;
            -f)
                fold=$2;
            shift 2;;
            -log)
                file_log_main=$2;
            shift 2;;
            -dlog)
                dir_log=$2;
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

# Start bash processing

echo -e "\nRunning : Visualization 2" >> $file_log_main
echo -e "Start time : `date +%F-%H:%M:%S`\n" >> $file_log_main

# Path
dir_main=/cbica/home/mayun/Projects/NMF
dir_MATLAB_code=/cbica/home/mayun/Projects/MATLAB
dir_result=${dir_main}/Result/${dataset}${subset}/FN_${K}/${N_Fold}F${fold}
file_fs_template=${dir_main}/Template/HCP/My_Template/HCP_FS.mat
file_init_robust=${dir_result}/init_robust/init.mat
dir_figure=${dir_result}/Figure
mkdir -p ${dir_figure}

file_log=${dir_log}/step2/Visualize.log
mkdir -p ${dir_log}/step2
if test -f "$file_log"
then
    rm ${file_log} 
fi

# Run compiled MATLAB function
file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
MATLAB_function=fVisualize_NMF_Step2_HCP_FS
MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh

jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log} \
${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${file_fs_template}'~'${file_init_robust}'~'${dir_figure} -log ${file_log})

echo "$jid" >> $file_log_main