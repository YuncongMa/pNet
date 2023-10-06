# Yuncong Ma, 9/25/2023
# pNet
# This script provides the highest level organization of pNet
# It provides workflows of pNet, and examples
# It includes five modules of pNet, underlying functions

#########################################
# Packages
import os

# path of pNet
dir_python = os.path.dirname(os.path.abspath(__file__))
dir_pNet = os.path.dirname(dir_python)
dir_brain_template = os.path.join(dir_pNet, 'Brain_Template')
dir_example = os.path.join(dir_pNet, 'Example')


# Example
from Example import Example
from Brain_Template import Brain_Template

# Module
# This script builds the five modules of pNet
# Functions for modules of pNet
from Data_Input import *
from FN_Computation import *
from FN_Computation_torch import *
from Computation_Environment import *
from Quality_Control import *
from Quality_Control_torch import *
from Workflow import *



