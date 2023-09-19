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


def run_workflow(dir_pnet_result: str, dataType: str, dataFormat: str, file_scan: str, file_subjectID: str):

    # setup all sub-folders in the pNet result folder
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT = setup_result_folder(dir_pnet_result)

    # ============== Data Input ============== #
    # setup_dataInput

    # ============================================= #

    # ============== FN Computation ============== #
    # setup parameters for FN computation

    # ============================================= #

    # ============== Quality Control ============== #
    # perform quality control

    # ============================================= #




