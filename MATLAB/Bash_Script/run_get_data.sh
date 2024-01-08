#$ /bin/bash
#$ -a

####################################################################################
# Yuncong Ma, 10/10/2022
# Get subject file from PNC dataset in bbl_pnc_zaixu
# Prepare sub files for NMF initialization
# bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/run_get_data.sh
####################################################################################

echo -e "\nRunning run_get_data.sh"
echo -e "Start time : `date +%F-%H:%M:%S`\n"

# Path
dataset=HCP/HCP_1200
dir_main=/cbica/home/mayun/Projects/NMF
dir_code=${dir_main}/scripts/HCP_SGC
dir_sbjList=${dir_main}/Subject_List/${dataset}
mkdir -p $dir_sbjList

# Parameter
N_BS=50
N_Sub_BS=100

# Raw data and masks of combined resting state and task
dir_raw=/cbica/projects/HCP_Data_Releases/HCP_1200

# Matlab function to do the job
file_run_matlab=${dir_main}/scripts/MATLAB/run_matlab_function.sh
dir_MATLAB_package=${dir_main}/scripts/MATLAB
MATLAB_function=fGet_Data_HCP_1200

# log file
file_log=${dir_main}/sge_output/step0/get_data_HCP_1200.log
mkdir -p ${dir_main}/sge_output/step0
if [ -f "$file_log" ]; then
    rm $file_log
fi

# Run matlab
jid=$(qsub -terse -j y -pe threaded 1 -l h_vmem=10G -o ${file_log} \
${file_run_matlab} -m ${dir_MATLAB_package} -f ${MATLAB_function} -p ${dir_raw}'~'${dir_sbjList}'~'${N_BS}'~'${N_Sub_BS} \
    -log ${file_log})


