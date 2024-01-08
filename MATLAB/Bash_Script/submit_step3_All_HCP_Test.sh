#!/bin/sh

####################################################################################
# Yuncong Ma, 1/23/2023
# This code is to run individual NMF using Hongming Li's NMF
# Each subject has multiple scans
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/submit_step3_All_HCP_Test.sh
####################################################################################

# For use as a function
parse()
{
    # Default setting
    K=50
    dataset=HCP/HCP_1200
    subset=/All
    dir_log=/cbica/home/mayun/Projects/NMF/sge_output/step3
    file_log_main=/cbica/home/mayun/Projects/NMF/sge_output/step3/submit_step3_All_HCP.log
    
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

echo -e "\nRunning : Process 3 " >> $file_log_main
echo -e "Start time : `date +%F-%H:%M:%S`\n" >> $file_log_main

# Path
dir_main=/cbica/home/mayun/Projects/NMF
dir_NMF_code=${dir_main}/Toolbox
dir_result=${dir_main}/Result/${dataset}${subset}/FN_${K}_Test
dir_init_robust=${dir_result}/init_robust
file_init=${dir_init_robust}/init.mat
dir_sbjList=${dir_main}/Subject_List/${dataset}
file_sbjList=${dir_sbjList}/Surface_FIX.txt
file_sbjIDList=${dir_sbjList}/Subject_Folder.txt
file_gNb=${dir_main}/Data/HCP/gNb.mat

# Set parameter
spaR=1
vxI=0
ard=0 # 1 for strong sparsity
eta=0 # 1 for strong sparsity
iterNum=100
alphaS21=20
alphaL=10
calcGrp=0
parforOn=0
resId="HCP"

# file
sbjfile=$(cat $file_sbjList |tr "\n" " ")
sbjfilearray=($sbjfile)
sbjID=$(cat $file_sbjIDList |tr "\n" " ")
sbjFolderarray=($sbjID)
let N=${#sbjFolderarray[@]}-1
#let N=9 # for fast testing

dir_raw=/cbica/projects/HCP_Data_Releases/HCP_1200
file_surfL=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.L.inflated.32k_fs_LR.surf.gii
file_surfR=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.R.inflated.32k_fs_LR.surf.gii
file_surfML=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.L.atlasroi.32k_fs_LR.shape.gii
file_surfMR=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.R.atlasroi.32k_fs_LR.shape.gii

file_gNb=${dir_main}/Data/HCP/gNb.mat
file_wb=${dir_main}/Toolbox/workbench/bin_rh_linux64/wb_command

echo -e "Submit job:\n" >> $file_log_main
mkdir -p ${dir_log}/step3
for i in $(seq 0 $N)
do
    dir_out=${dir_result}/Individual/${sbjFolderarray[i]}
    mkdir -p ${dir_out}
    sbjFileSingle=${dir_out}/Subject_File.txt
    if test -f "$sbjFileSingle"
    then 
        rm $sbjFileSingle
    fi
    let p1=i # start from 0
    echo ${sbjfilearray[$p1]} >> $sbjFileSingle

    sbjTCFolder=$dir_out
    timepointFile=0

    # log file
    file_log=${dir_log}/step3/submit_step3_All_HCP_FS_${i}.log
    if test -f "$file_log"
    then
        rm ${file_log} 
    fi

    # Run compiled MATLAB function
    file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
    MATLAB_function=deployFuncMvnmfL21p1_func_surf_hcp_single_gV
    MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh

    jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log} \
        ${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${sbjFileSingle}'~'${sbjTCFolder}'~'${file_wb}'~'${file_gNb}'~'${dir_out}'~'${resId}'~'${file_init}'~'${K}'~'${alphaS21}'~'${alphaL}'~'${vxI}'~'${spaR}'~'${ard}'~'${eta}'~'${iterNum}'~'${calcGrp}'~'${parforOn} \
        -log ${file_log})
    echo "$jid" >> $file_log_main

    sleep 2 # to avoid crush on matlab runtime
done