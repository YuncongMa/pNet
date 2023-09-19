# Yuncong Ma, 9/19/2023
# pNet
# Provide examples of running the whole workflow of pNet

#########################################
# Packages

# Example
from Example import Example

# Module
# This script builds the five modules of pNet
# Functions for modules of pNet
from Data_Input import *
from FN_Computation import *
from FN_Computation_torch import *
from Computation_Environment import *
from Quality_Control import *


def run_workflow(dir_pnet_result: str, dataType: str, dataFormat: str,
                 file_scan: str, file_subject_ID=None, file_subject_folder=None, file_group=None, scan_info='Manual',
                 K=17, Combine_Scan=False,
                 Compute_gFN=True, samplingMethod='Subject', sampleSize=10, nBS=50,
                 maxIter=1000, minIter=30, meanFitRatio=0.1, error=1e-6, normW=1,
                 Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
                 Parallel=False, Computation_Mode='CPU', N_Thread=1, dataPrecision='double'):
    """
    run_workflow(dir_pnet_result: str, dataType: str, dataFormat: str, file_scan: str, file_subjectID: str,
    file_scan: str, file_subject_ID: str, file_subject_folder: str, file_group: str, scan_info='Manual',
    K=17, Combine_Scan=False,
    Compute_gFN=True, samplingMethod='Subject', sampleSize=10, nBS=50, maxIter=1000, minIter=30,
    meanFitRatio=0.1, error=1e-6, normW=1, Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
    Parallel=False, Computation_Mode='CPU', N_Thread=1, dataPrecision='double')
    Run the workflow of pFN, including Data Input, FN Computation, and Quality Control

    :param dir_pnet_result: directory of the pNet result folder
    :param dataType: 'Surface', 'Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'Volume (*.nii, *.nii.gz, *.mat)'
    :param file_scan: a txt file that stores directories of all fMRI scans
    :param file_subject_ID: a txt file that store subject ID information corresponding to fMRI scan in file_scan
    :param file_subject_folder: a txt file that store subject folder names corresponding to fMRI scans in file_scan
    :param file_group: a txt file that store group information corresponding to fMRI scan in file_scan
    :param scan_info: 'Automatic' or 'Manual', 'Manual' requires manual input of file_subject_ID, file_subject_folder and file_group
    :param K: number of FNs
    :param Combine_Scan: False or True, whether to combine multiple scans for the same subject
    :param Compute_gFN: True or False, whether to compute gFNs from the provided data or load a precomputed gFN set
    :param samplingMethod: 'Subject' or 'Group_Subject'. Uniform sampling based subject ID, or group and then subject ID
    :param sampleSize: number of subjects selected for each bootstrapping run
    :param nBS: number of runs for bootstrap
    :param maxIter: maximum iteration number for multiplicative update
    :param minIter: minimum iteration in case fast convergence
    :param meanFitRatio: a 0-1 scaler, exponential moving average coefficient, used for the initialization of U when using group initialized V
    :param error: difference of cost function for convergence
    :param normW: 1 or 2, normalization method for W used in Laplacian regularization
    :param Alpha: hyper parameter for spatial sparsity
    :param Beta: hyper parameter for Laplacian sparsity
    :param alphaS: internally determined, the coefficient for spatial sparsity based Alpha, data size, K, and gNb
    :param alphaL: internally determined, the coefficient for Laplacian sparsity based Beta, data size, K, and gNb
    :param vxI: flag for using the temporal correlation between nodes (vertex, voxel)
    :param ard: 0 or 1, flat for combining similar clusters
    :param eta: a hyper parameter for the ard regularization term
    :param nRepeat: Any positive integer, the number of repetition to avoid poor initialization
    :param Parallel: False or True, whether to enable parallel computation
    :param Computation_Mode: 'CPU'
    :param N_Thread: positive integers, used for parallel computation
    :param dataPrecision: 'double' or 'single'

    Yuncong Ma, 9/19/2023
    """

    # setup all sub-folders in the pNet result folder
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT = setup_result_folder(dir_pnet_result)

    # ============== Data Input ============== #
    # setup_dataInput
    setup_dataInput(dir_pnet_dataInput, dataType=dataType, dataFormat=dataFormat)
    setup_scan_info(dir_pnet_dataInput=dir_pnet_dataInput, file_scan=file_scan, file_subject_ID=file_subject_ID, file_subject_folder=file_subject_folder, file_group=file_group,
                    scan_info=scan_info, Combine_Scan=Combine_Scan)
    # ============================================= #

    # ============== FN Computation ============== #
    # setup parameters for FN computation
    setup_NMF_setting(dir_pnet_result, K=K, Combine_Scan=Combine_Scan, Compute_gFN=Compute_gFN, samplingMethod=samplingMethod, sampleSize=sampleSize, nBS=nBS, maxIter=maxIter,
                      minIter=minIter, meanFitRatio=meanFitRatio, error=error, normW=normW, Alpha=Alpha, Beta=Beta, alphaS=alphaS, alphaL=alphaL, vxI=vxI, ard=ard, eta=eta, nRepeat=nRepeat,
                      Parallel=Parallel, Computation_Mode=Computation_Mode, N_Thread=N_Thread, dataPrecision=dataPrecision)
    # perform FN computation
    run_FN_Computation(dir_pnet_result)
    # ============================================= #

    # ============== Quality Control ============== #
    # perform quality control
    run_quality_control(dir_pnet_result)
    # ============================================= #




