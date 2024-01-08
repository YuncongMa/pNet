#!/bin/sh

####################################################################################
# Yuncong Ma, 10/18/2022
# This code is to choose the best initialization for Hongming Li's NMF
# It only extract init.mat that is finished in pre-selected folders
# It could be used to run on both volume and surface data
# It will overwrite existing init_file_list.txt and log file
# Select useful init results based on the average std of each FN
# bash /cbica/home/mayun/Projects/NMF/scripts/PNC_SGC/submit_step2_All_HCP.sh
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

echo -e "\nRunning : Process 2 " >> $file_log_main
echo -e "Start time : `date +%F-%H:%M:%S`\n" >> $file_log_main

# Path
dir_main=/cbica/home/mayun/Projects/NMF
dir_NMF_code=${dir_main}/Toolbox
dir_result=${dir_main}/Result/${dataset}${subset}/FN_${K}/${N_Fold}F${fold}
dir_init_robust=${dir_result}/init_robust

mkdir -p ${dir_init_robust}

dir_init=/cbica/home/mayun/Projects/NMF/Result/${dataset}${subset}/FN_${K}/${N_Fold}F${fold}
dir_init_bs=${dir_init}/init_bs50_100

# init_file_list.txt
file_init=${dir_init}/init_file_list.txt
if test -f "$file_init"
then
    echo -e "Delete existing init_file_list.txt\n" >> ${file_log_main}
    rm ${file_init}
fi

# Find all exisiting init.mat
# Report
dir=$(ls -l ${dir_init_bs} | awk '/^d/ {print $NF}') # $NF output last field
let count=0
let count2=0
for i in $dir
do
    file_init_bs=${dir_init_bs}/${i}/init.mat
    if test -f "$file_init_bs"
    then
        echo ${file_init_bs} >> ${file_init}
        echo "Add init.mat in folder "${i} >> ${file_log_main}
        let count++
    else
        echo "Skip init.mat in folder "${i} >> ${file_log_main}
        let count2++
    fi
done

echo -e "\nNumber of init.mat is $count" >> ${file_log_main}
echo -e "Number of unexisted init.mat is $count2 \n" >> ${file_log_main}
echo -e "\n" >> ${file_log_main}

# log
file_log=$dir_log/step2/submit_step2_All_HCP_2F.log
mkdir -p $dir_log/step2
if test -f "$file_log"
then
    rm ${file_log}
fi

# fSelect_NMF_Init
# Run compiled MATLAB function
dir_MATLAB_package=${dir_main}/scripts/MATLAB
file_run_matlab=${dir_main}/scripts/MATLAB/run_matlab_function.sh
file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
MATLAB_function=fSelect_NMF_Init
MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh
file_init2=${dir_init}/init_file_list_selected.txt
file_init_robust=${dir_result}/init_robust/init.mat

if test -f "$file_init2"
then
    rm $file_init2
fi
if test -f "$file_init_robust"
then
    rm $file_init_robust
fi
file_log_select=$dir_log/step2/step2_select.log
#jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log_select} \
#${file_run_matlab} -m ${dir_MATLAB_package} -f ${MATLAB_function} -p ${file_init}'~'${file_init2} -log ${file_log_select})
jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log_select} \
${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${file_init}'~'${file_init2} -log ${file_log_select})

# selRobustInit
# Run compiled MATLAB function
file_run_matlab_compiled=${dir_main}/scripts/MATLAB/run_matlab_compiled.sh
MATLAB_function=selRobustInit
MATLAB_bash=${dir_main}/scripts/MATLAB/run_${MATLAB_function}.sh

while ! test -f "$file_init2"
do
    sleep 10
done
jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=20G -o ${file_log} \
        ${file_run_matlab_compiled} -f ${MATLAB_bash} -p ${file_init2}'~'${K}'~'${dir_init_robust} -log ${file_log})

echo -e "Submit job:\n$jid" >> $file_log_main