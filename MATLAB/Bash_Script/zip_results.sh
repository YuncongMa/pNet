#!/bin/sh

####################################################################################
# Yuncong Ma, 12/6/2022
# backup code for NMF including MATLAB, scripts and Toolbox
# Terminal code: bash /cbica/home/mayun/Projects/NMF/scripts/HCP_SGC/zip_results.sh
# Submit job
# dir_main=/cbica/home/mayun/Projects/NMF
# qsub -terse -j y -pe threaded 1 -l h_vmem=20G -o ${dir_main}/sge_output/collect/zip_results.log ${dir_main}/scripts/HCP_SGC/zip_results.sh
####################################################################################

echo -e "\nRunning backup.sh on                       : `hostname`"
echo -e "Start time                                   : `date +%F-%H:%M:%S`\n"

# Destination
dir_backup=/cbica/home/mayun/Projects/NMF/Result/HCP/HCP_1200/All/FN_50
mkdir -p $dir_backup

file_backup=${dir_backup}"/Result.zip"
if test -f "$file_backup"
then
    rm -r $file_backup
fi

# Folder to be backuped
dir_data=/cbica/home/mayun/Projects/NMF/Result/HCP/HCP_1200/All/FN_50/Individual

# Backup
zip -r $file_backup $dir_data
