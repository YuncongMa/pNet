#!/bin/sh

####################################################################################
# Yuncong Ma, 12/13/2022
# It is to generate figures for individual NMF results stored in init.mat
# It uses Yuncong's MATLAB package to do visualization
# Each subject has multiple scans
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/visualize_step3_All_HCP_2F.sh
####################################################################################

# For use as a function
parse()
{
    # Default setting
    K=3
    dataset=HCP/HCP_1200
    subset=/All
    N_Fold=2
    fold=1
    dir_log=/cbica/home/mayun/Projects/NMF/sge_output/step1
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

echo -e "\nRunning : Visualization 3" >> $file_log_main
echo -e "Start time : `date +%F-%H:%M:%S`\n" >> $file_log_main

# Path
dir_main=/cbica/home/mayun/Projects/NMF
dir_MATLAB_code=/cbica/home/mayun/Projects/MATLAB
file_fs_template=${dir_main}/Template/HCP/My_Template/HCP_FS.mat
dir_sbjList=${dir_main}/Subject_List/${dataset}
file_sbjIDList=${dir_sbjList}/Subject_Folder.txt

sbjID=$(cat $file_sbjIDList |tr "\n" " ")
sbjFolderarray=($sbjID)
let N=${#sbjFolderarray[@]}-1
# let N=19 # for fast testing

echo -e "Submit job:\n" >> $file_log_main
for i in $(seq 0 $N)
do
    dir_result=${dir_main}/Result/${dataset}${subset}/FN_${K}/${N_Fold}F${fold}/Individual/${sbjFolderarray[i]}
    file_init_robust=$(find $dir_result -type f -name 'final_UV.mat')
    dir_figure=${dir_result}
    
    # log file
    file_log=${dir_log}/step3/visualize_step3_All_HCP_${i}.log
    if test -f "$file_log"
    then
        rm ${file_log}
    fi
    
    # Run compiled MATLAB function
    file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
    MATLAB_function=fVisualize_NMF_Step3_HCP_FS_single
    MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh
    
    jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log} \
        ${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${file_fs_template}'~'${file_init_robust}'~'${dir_figure} -log ${file_log})
    
    echo "$jid" >> $file_log_main
done