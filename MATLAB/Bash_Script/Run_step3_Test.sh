#!/bin/sh

####################################################################################
# Yuncong Ma, 1/18/2023
# This code is to run individual NMF using Hongming Li's NMF
# Use scans failed QC control to test different NMF methods
# Submit job
# dir_main=/cbica/home/mayun/Projects/NMF
# qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${dir_main}/sge_output/streamline/Run_step3_Test.log ${dir_main}/scripts/HCP_SGC/Run_step3_Test.sh
####################################################################################

# Setup
K=50
dataset=HCP/HCP_1200
subset=/All
# Test folder
Test_Number=26
dir_result=/cbica/home/mayun/Projects/NMF/Result/${dataset}${subset}/FN_${K}_Test_${Test_Number}
# Function for step3
MATLAB_function=deployFuncMvnmfL21p1_func_surf_hcp_single_gV_17


# Start bash processing

# Path
dir_main=/cbica/home/mayun/Projects/NMF
dir_NMF_code=${dir_main}/Toolbox
# clean up everthing in the folder
rm -rf $dir_result
mkdir -p ${dir_result}
mkdir -p ${dir_result}/init_robust
mkdir -p ${dir_result}/Individual
mkdir -p ${dir_result}/sge_output/step3
mkdir -p ${dir_result}/sge_output/collect
# Copy init.mat
file_init=${dir_result}/init_robust/init.mat
scp ${dir_main}/Result/${dataset}${subset}/FN_${K}/init_robust/init.mat $file_init
# File list
dir_sbjList=${dir_main}/Subject_List/${dataset}
#file_sbjList=${dir_sbjList}/Surface_FIX_Failed_QC_K50.txt
#file_sbjIDList=${dir_sbjList}/Subject_Folder_Failed_QC_K50.txt
file_sbjList=${dir_sbjList}/Surface_FIX.txt
file_sbjIDList=${dir_sbjList}/Subject_Folder.txt
# gNb
file_gNb=${dir_main}/Data/HCP/gNb.mat

# main log file
file_log_main=${dir_result}/sge_output/step3/Run_step3_Test.log
if test -f "$file_log_main"
then
    rm $file_log_main
fi

echo -e "\nRunning : Process 3 " >> $file_log_main
echo -e "Start time : `date +%F-%H:%M:%S`\n" >> $file_log_main

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


dir_raw=/cbica/projects/HCP_Data_Releases/HCP_1200
file_surfL=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.L.inflated.32k_fs_LR.surf.gii
file_surfR=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.R.inflated.32k_fs_LR.surf.gii
file_surfML=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.L.atlasroi.32k_fs_LR.shape.gii
file_surfMR=${dir_raw}/100307/MNINonLinear/fsaverage_LR32k/100307.R.atlasroi.32k_fs_LR.shape.gii

file_gNb=${dir_main}/Data/HCP/gNb.mat
file_wb=${dir_main}/Toolbox/workbench/bin_rh_linux64/wb_command

# Compile MATLAB function
echo -e "Compile MATLAB function ${MATLAB_function}" >> $file_log_main
MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh

file_log=${dir_result}/sge_output/step3/compile_step3_deploy.log
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
jid=$(qsub -terse -j y -pe threaded 4 -l h_vmem=10g -o ${file_log} ${dir_code}/run_compile_matlab.sh -cf ${dir_compile_function} -dp ${dir_MATLAB_package} -fn ${MATLAB_function} -do ${dir_output});
let Flag=0
let N_Flag=0
Flag=$(qstat | grep $jid | wc -l)
while [ "$Flag" -gt "$N_Flag" ]
do
    Flag=$(qstat | grep $jid | wc -l)
    sleep 20
done
echo -e "File Compiled\n" >> $file_log_main

# Run step 3
echo -e "Submit job:\n" >> $file_log_main
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
    dir_log=${dir_result}/sge_output
    file_log=${dir_log}/step3/Run_step3_Test_${i}.log
    if test -f "$file_log"
    then
        rm ${file_log}
    fi
    
    # Run compiled MATLAB function
    file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
    
    jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log} \
        ${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${sbjFileSingle}'~'${sbjTCFolder}'~'${file_wb}'~'${file_gNb}'~'${dir_out}'~'${resId}'~'${file_init}'~'${K}'~'${alphaS21}'~'${alphaL}'~'${vxI}'~'${spaR}'~'${ard}'~'${eta}'~'${iterNum}'~'${calcGrp}'~'${parforOn} \
    -log ${file_log})
    echo "$jid" >> $file_log_main
    
    sleep 2 # to avoid crush on matlab runtime
done

# check the completion of step 3
let Flag=0
let N_Flag=$N
while [ "$Flag" -lt "$N_Flag" ]
do
    dir=$(find ${dir_result}/Individual -name "final_UV.mat" -type f)
    dirarray=($dir)
    let Flag=${#dirarray[@]}-1
    if  [ "$Flag" -lt "$N_Flag" ]
    then
        sleep 60
    fi
done
echo -e "Process 3 is finished\n" >> $file_log_main


# QC
file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
MATLAB_function=fNMF_QC_SC_FH
MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh

file_wb=${dir_main}/Toolbox/workbench/bin_rh_linux64/wb_command

file_groupV=${dir_result}/init_robust/init.mat
dir_individual=${dir_result}/Individual
file_result=${dir_result}/QC_SC_FH.mat

# log file
file_log=${dir_result}/sge_output/collect/NMF_QC_Job.log
mkdir -p ${dir_result}/sge_output/collect
if test -f "$file_log"
then
    rm ${file_log}
fi

dir_sbjList=${dir_main}/Subject_List/HCP/HCP_1200
file_scanList=${dir_sbjList}/Surface_FIX.txt
file_folderList=${dir_sbjList}/Subject_Folder.txt

jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log} \
        ${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${file_scanList}'~'${file_folderList}'~'${file_groupV}'~'${dir_individual}'~'${file_result}'~HCP_Surface~'${file_wb} \
        -log ${file_log})

#jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log} \
#    ${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${file_groupV}'~'${dir_individual}'~'${file_result} \
#-log ${file_log})
echo -e "Submit QC job: ${jid}\n" >> $file_log_main

# check the completion of QC
let Flag=0
let N_Flag=0
Flag=$(qstat | grep $jid | wc -l)
while [ "$Flag" -gt "$N_Flag" ]
do
    Flag=$(qstat | grep $jid | wc -l)
    sleep 20
done
echo -e "QC completed\n" >> $file_log_main

exit