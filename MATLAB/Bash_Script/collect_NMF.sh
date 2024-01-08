#!bin/sh

####################################################################################
# Yuncong Ma, 1/17/2023
# Calculate individual FC in different NMF and collect results into multiple MAT files
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/collect_NMF.sh
####################################################################################

# Path List for different NMF processes
#list_dataset=('/HCP/HCP_1200/All/FN_17' '/HCP/HCP_1200/All/FN_17/2F1' '/HCP/HCP_1200/All/FN_17/2F2')
#list_dataset=('/HCP/HCP_1200/All/FN_50' '/HCP/HCP_1200/All/FN_50/2F1' '/HCP/HCP_1200/All/FN_50/2F2')
#list_dataset=('/HCP/HCP_1200/All/FN_50' '/HCP/HCP_1200/All/FN_50_Test')
#list_dataset=('/HCP/HCP_1200/All/FN_50_Test_11' '/HCP/HCP_1200/All/FN_50_Test_12' '/HCP/HCP_1200/All/FN_50_Test_13' '/HCP/HCP_1200/All/FN_50_Test_14' '/HCP/HCP_1200/All/FN_50_Test_15')
list_dataset=('/HCP/HCP_1200/All/FN_50_Test_22' '/HCP/HCP_1200/All/FN_50_Test_23')

# Path
dir_main=/cbica/home/mayun/Projects/NMF
dir_script=${dir_main}/scripts/HCP_SGC
dir_result=${dir_main}/Result

# Run compiled MATLAB function
dir_MATLAB_package=/cbica/home/mayun/Projects/MATLAB'~'${dir_main}/scripts/MATLAB
file_run_matlab=${dir_main}/scripts/MATLAB/run_matlab_function.sh
file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
MATLAB_function=fCollect_NMF
MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh

let N=${#list_dataset[@]}
let N=N-1
for i in $(seq 0 $N)
do
    dataset=${list_dataset[i]}
    echo "Collecting results from ${dataset}"
    dir_individual=${dir_result}${dataset}/Individual
    file_result=${dir_result}${dataset}/Individual_Result.mat
    
    # log file
    file_log=${dir_main}/sge_output/collect/collect_NMF_Job${i}.log
    mkdir -p ${dir_main}/sge_output/collect
    if test -f "$file_log"
    then
        rm ${file_log}
    fi
    
    #jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=20G -o ${file_log} \
    #    ${file_run_matlab} -m ${dir_MATLAB_package} -f ${MATLAB_function} -p ${dir_individual}'~'${file_result} \
    #   -log ${file_log})
    
    jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=100G -o ${file_log} \
        ${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${dir_individual}'~'${file_result} \
        -log ${file_log})
done

