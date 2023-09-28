# This is a brief description of the Python version of pNet
The Python version is designed to carry out the pNet setup and computation, consistent to the MATLAB version with GUI.
It has Numpy and PyTorch versions for easy code development and fast computation.
pNet.py is the main function to import all required packages and functions

# Deployment
This Python package requires additional toolboxes which can be found in the Install_Toolbox.py
Specifically, it requires:
1. numpy
2. pytorch
3. scipy
4. nibabel
5. h5py

It is recommended to use 'import pNet' to avoid conflicting function names

# Workflow
The workflow code is in Workflow.py
This streamlines the setup and computation of pNet, including modules below
1. Data Input
2. FN Computation
3. Quality Control.


# Module
1. Data Input: organize fMRI scans based on subject information, prepare brain template for subsequent FN computation and visualization
2. FN Computation: setup and carry out the FN computation, generating group-level and personalized FNs
3. Quality Control: ensures that pFNs have the highest spatial similarity to their group-level counterparts

# Additionals
Example: It includes example results of pNet on HCP, UKBB and PNC datasets
Brain Template: It includes prepared brain template files for HCP, MNI space

# Example script
Several simple examples are provided in Example_Workflow.py
