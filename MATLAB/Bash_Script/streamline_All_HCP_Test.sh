#!/bin/sh

####################################################################################
# Yuncong Ma, 1/23/2023
# A streamline control file to ensure each Process is finished to run the next one
# Allow global settings of K, dataset and subset
# Works for FreeSurfer format
# Add groupV for individualized FN computation
# Do not run this code file using bash!
# This file requires qsub to run:
# dir_main=/cbica/home/mayun/Projects/NMF
# qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${dir_main}/sge_output/streamline/NMF_All_HCP.log ${dir_main}/scripts/HCP_SGC/streamline_All_HCP_Test.sh
####################################################################################

# Process to run
Process=
Process[1]=0
Process[2]=0
Process[3]=1
# Visualization to run
Visualization=
Visualization[2]=0
Visualization[3]=0

# Parameter
K=50
let N_BS=50
let N_Visualization=20

# Path
dataset=HCP/HCP_1200
subset=/All
dir_main=/cbica/home/mayun/Projects/NMF
dir_script=${dir_main}/scripts/HCP_SGC
dir_result=${dir_main}/Result/${dataset}${subset}/FN_${K}_Test
dir_init_robust=${dir_result}/init_robust

# log
dir_log=${dir_result}/sge_output
mkdir -p $dir_log
file_log_main=${dir_log}/NMF_All_HCP.log
if test -f "$file_log_main"
then
    rm $file_log_main
fi

# Subject File
file_sbjID=${dir_main}/Subject_List/${dataset}/SubjectID.txt
sbjID=$(cat $file_sbjID |tr "\n" " ")
sbjIDarray=($sbjID)
let N_File=${#sbjIDarray[@]}

# Wait time for streamline
Wait_Process=
Wait_Process[2]=3600
Wait_Process[3]=60
Wait_Visualization=
Wait_Visualization[2]=60
Wait_Visualization[3]=600


echo -e "\nRunning streamline_All_HCP.sh" >> $file_log_main
echo -e "Start time : `date +%F-%H:%M:%S`\n" >> $file_log_main

echo -e "Setting: \n" >> $file_log_main
echo "dataset = $dataset" >> $file_log_main
echo "subset = $subset" >> $file_log_main
echo -e "K = $K\n" >> $file_log_main

####################################################################################
############################### Streamline Process #################################
####################################################################################

### Process 1 ###
if [ "${Process[1]}" == "1" ]
then
    echo -e "\nProcess 1\n" >> $file_log_main

    echo -e "Submit Process 1 at `date +%F-%H:%M:%S`\n" >> $file_log_main
    bash $dir_script/submit_step1_All_HCP_Test.sh -k $K -ds $dataset -ss $subset -log $file_log_main -dlog $dir_log
fi


### Process 2 ###
if [ "${Process[2]}" == "1" ]
then
    echo -e "\nProcess 2\n" >> $file_log_main
    let Flag=0
    let N_Flag=N_BS
    dir_Process2=${dir_result}/init_bs50_100
    while [ "$Flag" -lt "$N_Flag" ]
    do
        echo -e "Waiting for Process 2 to start at `date +%F-%H:%M:%S`\n" >> $file_log_main
        
        dir=$(find $dir_Process2 -name "init.mat" -type f)
        dirarray=($dir)
        let Flag=${#dirarray[@]}
        echo -e "Found $Flag jobs finished in a total number of $N_Flag\n" >> $file_log_main
        if [ "$Flag" -lt "$N_Flag" ]
        then     
            sleep ${Wait_Process[2]}
        fi
    done

    echo -e "Submit Process 2 at `date +%F-%H:%M:%S`\n"
    bash $dir_script/submit_step2_All_HCP_Test.sh -k $K -ds $dataset -ss $subset -log $file_log_main -dlog $dir_log
fi

### visualize 2 ###
if [ "${Visualization[2]}" == "1" ]
then
    echo -e "\nVisualization 2\n" >> $file_log_main
    let Flag=0
    let N_Flag=1
    file_Process3=${dir_init_robust}/init.mat
    while [ "$Flag" -lt "$N_Flag" ]
    do
        echo "Waiting for Visualization 2 to start at `date +%F-%H:%M:%S`" >> $file_log_main
        
        if test -f "$file_Process3"
        then
            let Flag=1
        fi
        echo -e "Found $Flag jobs finished in a total number of $N_Flag\n" >> $file_log_main
        if [ "$Flag" -lt "$N_Flag" ]
        then     
            sleep ${Wait_Visualization[2]}
        fi
    done
    
    echo -e "Submit Visualization 2 at `date +%F-%H:%M:%S`\n" >> $file_log_main
    bash $dir_script/visualize_step2_All_HCP_Test.sh -k $K -ds $dataset -ss $subset -log $file_log_main -dlog $dir_log
fi

### Process 3 ###
if [ "${Process[3]}" == "1" ]
then
    echo -e "\nProcess 3\n" >> $file_log_main
    let Flag=0
    let N_Flag=1
    file_Process3=${dir_init_robust}/init.mat
    while [ "$Flag" -lt "$N_Flag" ]
    do
        echo "Waiting for Process 3 to start at `date +%F-%H:%M:%S`" >> $file_log_main
        
        if test -f "$file_Process3"
        then
            let Flag=1
        fi
        echo -e "Found $Flag jobs finished in a total number of $N_Flag\n" >> $file_log_main
        if [ "$Flag" -lt "$N_Flag" ]
        then     
            sleep ${Wait_Process[3]}
        fi
    done
    
    echo -e "Submit Process 3 at `date +%F-%H:%M:%S`\n" >> $file_log_main
    bash $dir_script/submit_step3_All_HCP_Test.sh -k $K -ds $dataset -ss $subset -log $file_log_main -dlog $dir_log
fi

### Visualization 3 ###
if [ "${Visualization[3]}" == "1" ]
then
    echo -e "\nVisualization 3\n" >> $file_log_main
    let Flag=0
    let N_Flag=N_File
    dir_Visualize3=${dir_result}/Individual
    while [ "$Flag" -lt "$N_Flag" ]
    do
        echo "Waiting for Visualization 3 to start at `date +%F-%H:%M:%S`" >> $file_log_main
        
        dir=$(find $dir_Visualize3 -name "final_UV.mat" -type f)
        dirarray=($dir)
        let Flag=${#dirarray[@]}
        echo -e "Found $Flag jobs finished in a total number of $N_Flag\n" >> $file_log_main
        if [ "$Flag" -lt "$N_Flag" ]
        then     
            sleep ${Wait_Visualization[3]}
        fi
    done

    echo -e "Submit Visualization 3 at `date +%F-%H:%M:%S`\n" >> $file_log_main
    bash $dir_script/visualize_step3_All_HCP_Test.sh -k $K -ds $dataset -ss $subset -log $file_log_main -dlog $dir_log
fi




####################################################################################
############################### Check job status ###################################
####################################################################################

echo -e "\n--  Checking job status --\n" >> $file_log_main

if [ "${Process[1]}" == "1" ]
then
    let Flag=0
    let N_Flag=N_BS
    dir_Process2=${dir_result}/init_bs50_100
    while [ "$Flag" -lt "$N_Flag" ]
    do
        dir=$(find $dir_Process2 -name "init.mat" -type f)
        dirarray=($dir)
        let Flag=${#dirarray[@]}
        if  [ "$Flag" -lt "$N_Flag" ]
        then
            sleep ${Wait_Process[2]}
        fi
    done
    echo -e "Process 1 is finished\n" >> $file_log_main
fi

if [ "${Process[2]}" == "1" ]
then
    let Flag=0
    let N_Flag=1
    file_Process3=${dir_init_robust}/init.mat
    while [ "$Flag" -lt "$N_Flag" ]
    do
        if test -f "$file_Process3"
        then
            let Flag=1
        fi
        if  [ "$Flag" -lt "$N_Flag" ]
        then
            sleep ${Wait_Process[2]}
        fi
    done 
    echo -e "Process 2 is finished\n" >> $file_log_main
fi

if [ "${Visualization[2]}" == "1" ]
then
    let Flag=0
    let N_Flag=1
    dir_Process2=${dir_result}/Figure
    while [ "$Flag" -lt "$N_Flag" ]
    do
        dir=$(find $dir_Process2 -name "All.jpg" -type f)
        dirarray=($dir)
        let Flag=${#dirarray[@]}
        if  [ "$Flag" -lt "$N_Flag" ]
        then
            sleep ${Wait_Visualization[2]}
        fi
    done
    echo -e "Visualization 2 is finished\n" >> $file_log_main
fi

if [ "${Process[3]}" == "1" ]
then
    let Flag=0
    let N_Flag=N_File
    dir_Process3=${dir_result}/Individual
    while [ "$Flag" -lt "$N_Flag" ]
    do   
        dir=$(find $dir_Process3 -name "final_UV.mat" -type f)
        dirarray=($dir)
        let Flag=${#dirarray[@]}
        if  [ "$Flag" -lt "$N_Flag" ]
        then
            sleep ${Wait_Process[3]}
        fi
    done
    echo -e "Process 3 is finished\n" >> $file_log_main
fi

if [ "${Visualization[3]}" == "1" ]
then
    let Flag=0
    let N_Flag=N_Visualization
    dir_Visualize3=${dir_result}/Individual
    while [ "$Flag" -lt "$N_Flag" ]
    do  
        dir=$(find $dir_Visualize3 -name "All.jpg" -type f)
        dirarray=($dir)
        let Flag=${#dirarray[@]}
        if  [ "$Flag" -lt "$N_Flag" ]
        then
            sleep ${Wait_Visualization[3]}
        fi
    done
    echo -e "Visualization 3 is finished\n" >> $file_log_main
fi

echo -e "\nAll done at `date +%F-%H:%M:%S`" >> $file_log_main