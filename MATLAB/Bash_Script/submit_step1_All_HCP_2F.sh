#!/bin/sh

####################################################################################
# Yuncong Ma, 10/11/2022
# This script is to do step 1 of NMF for timepoints at high global signal
# It is to use FreeSurfer data with bootstrap for initialization
# Only use timepoints noted in TP_subjID.mat (Time_Point)
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/submit_step1_All_HCP_2F.sh
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
    file_log_main=/cbica/home/mayun/Projects/NMF/sge_output/step1/submit_step1_All_HCP_2F.log

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

echo -e "\nRunning : Process 1" >> $file_log_main
echo -e "Start time : `date +%F-%H:%M:%S`\n" >> $file_log_main

# Path
dir_main=/cbica/home/mayun/Projects/NMF
dir_code=${dir_main}/scripts/HCP_SGC
dir_sbjList=${dir_main}/Subject_List/${dataset}
dir_NMF_code=${dir_main}/Toolbox

# File
dir_raw=/cbica/projects/HCP_Data_Releases/HCP_1200
file_surfL=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.L.inflated.32k_fs_LR.surf.gii
file_surfR=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.R.inflated.32k_fs_LR.surf.gii
file_surfML=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.L.atlasroi.32k_fs_LR.shape.gii
file_surfMR=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.R.atlasroi.32k_fs_LR.shape.gii

file_gNb=${dir_main}/Data/HCP/gNb.mat
file_wb=${dir_main}/Toolbox/workbench/bin_rh_linux64/wb_command

# MATLAB
file_run_matlab=${dir_main}/scripts/MATLAB/run_matlab_function.sh
dir_MATLAB_package=${dir_main}/Toolbox/Code_mvNMF_l21_ard_v3_debug
MATLAB_function=deployFuncInit_surf_hcp_tp
file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh

# Prepare bootstrap file list for timepoints
N_BS=50

# Initialization of NMF
# Parameter for running run_step1*.sh
spaR=1
vxI=0
ard=1
iterNum=1000
timeNum=1200
alpha=2
beta=10
resId="HCP"

dir_result=${dir_main}/Result/${dataset}${subset}/FN_${K}/${N_Fold}F${fold}
mkdir -p $dir_result
dir_out=${dir_result}/init_bs50_100
mkdir -p ${dir_out}
mkdir -p ${dir_log}/step1
echo -e "Submit jobs:\n" >> $file_log_main
for i in $(seq 1 $N_BS)
do
    file_sbjList=${dir_sbjList}/Surface_FIX_${N_Fold}F${fold}_2BS100_${i}.txt
    file_timepoint=0
    
    list_id=${i}
    file_log=${dir_log}/step1/All_${i}.log
    mkdir -p ${dir_log}/step1
    if test -f "$file_log"
    then
        rm ${file_log}
    fi
    
    #jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=30G -o ${file_log} \
    #    ${file_run_matlab} -m ${dir_MATLAB_package} -f ${MATLAB_function} \
    #    -p ${file_sbjList}'~'${file_wb}'~'${file_surfL}'~'${file_surfR}'~'${file_surfML}'~'${file_surfMR}'~'${file_gNb}'~'${dir_out}'~'${list_id}'~'${spaR}'~'${vxI}'~'${ard}'~'${iterNum}'~'${K}'~'${timeNum}'~'${alpha}'~'${beta}'~'${resId}'~'${file_timepoint} \
    #   -log ${file_log})
    jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=100G -o ${file_log} \
        ${file_run_matlab_compiled} -f ${MATLAB_bash} \
        -p ${file_sbjList}'~'${file_wb}'~'${file_surfL}'~'${file_surfR}'~'${file_surfML}'~'${file_surfMR}'~'${file_gNb}'~'${dir_out}'~'${list_id}'~'${spaR}'~'${vxI}'~'${ard}'~'${iterNum}'~'${K}'~'${timeNum}'~'${alpha}'~'${beta}'~'${resId}'~'${file_timepoint} \
    -log ${file_log})

    echo "$jid" >> $file_log_main
done
