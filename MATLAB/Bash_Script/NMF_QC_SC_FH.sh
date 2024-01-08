#!bin/sh

####################################################################################
# Yuncong Ma, 1/15/2022
# Calculate QC for individual FC in different NMF and collect results into multiple MAT files
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/NMF_QC_SC_FH.sh
####################################################################################

# Path List for different NMF processes
#list_dataset=('/HCP/HCP_1200/All/FN_17' '/HCP/HCP_1200/All/FN_17/2F1' '/HCP/HCP_1200/All/FN_17/2F2')
#list_dataset=('/HCP/HCP_1200/All/FN_50' '/HCP/HCP_1200/All/FN_50/2F1' '/HCP/HCP_1200/All/FN_50/2F2')
list_dataset=('/HCP/HCP_1200/All/FN_50_Test_13' '/HCP/HCP_1200/All/FN_50_Test_19' '/HCP/HCP_1200/All/FN_50_Test_20' '/HCP/HCP_1200/All/FN_50_Test_21')
#list_dataset=('/HCP/HCP_1200/All/FN_50_Test_11' '/HCP/HCP_1200/All/FN_50_Test_13' '/HCP/HCP_1200/All/FN_50_Test_15' '/HCP/HCP_1200/All/FN_50_Test_17')

# Path
dir_main=/cbica/home/mayun/Projects/NMF
dir_script=${dir_main}/scripts/HCP_SGC
dir_result=${dir_main}/Result

# Subject file
dir_sbjList=${dir_main}/Subject_List/HCP/HCP_1200
file_scanList=${dir_sbjList}/Surface_FIX.txt
file_folderList=${dir_sbjList}/Subject_Folder.txt

# Run compiled MATLAB function
dir_MATLAB_package=/cbica/home/mayun/Projects/MATLAB'~'${dir_main}/scripts/MATLAB'~'${dir_main}/Toolbox/Code_mvNMF_l21_ard_v3_debug
file_run_matlab=${dir_main}/scripts/MATLAB/run_matlab_function.sh
file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
MATLAB_function=fNMF_QC_SC_FH
MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh



file_wb=${dir_main}/Toolbox/workbench/bin_rh_linux64/wb_command

let N=${#list_dataset[@]}
let N=N-1
for i in $(seq 0 $N)
do
    dataset=${list_dataset[i]}
    echo "QC results from ${dataset}"
    file_groupV=${dir_result}${dataset}/init_robust/init.mat
    dir_individual=${dir_result}${dataset}/Individual
    file_result=${dir_result}${dataset}/QC_SC_FH.mat
    
    # log file
    file_log=${dir_main}/sge_output/collect/NMF_QC_SC_FH_Job${i}.log
    mkdir -p ${dir_main}/sge_output/collect
    if test -f "$file_log"
    then
        rm ${file_log}
    fi
    
    #jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=20G -o ${file_log} \
    #    ${file_run_matlab} -m ${dir_MATLAB_package} -f ${MATLAB_function} -p ${dir_individual}'~'${file_result} \
    #   -log ${file_log})
    
    jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log} \
        ${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${file_scanList}'~'${file_folderList}'~'${file_groupV}'~'${dir_individual}'~'${file_result}'~HCP_Surface~'${file_wb} \
        -log ${file_log})
done

