# This is a brief description of the Python version of pNet


* It is designed to carry out the pNet setup and computation, consistent to the MATLAB version with GUI. 
* It features a step-by-step guidance to configure a customized script to run a workflow.
* It provides both Numpy and PyTorch versions for easy code development and fast computation.

## Deployment

It is recommended to use Anaconda (https://www.anaconda.com) to create an environment. <br />
**Create a new conda environment for pNet**
```
conda env create --name pnet -f environment.yml
```
**Or install required tools in conda for an existing environment**
```
pip install -r requirements.txt
```

**Download pNet**

Change <User's directory> to a desired directory
```
git clone https://github.com/YuncongMa/pNet <User's directory>
```

## Step-by-step guidance
This Python version of pNet provides a step-by-step guidance to customize a Python script for running a desired workflow <br />
**Terminal code**
```
# run a step-by-step terminal guide to set up a Python script for a customized pNet workflow
python Workflow_guide.py
```

**Video tutorial**

<img src="https://github.com/YuncongMa/pNet/assets/20191790/c083cb5b-764c-4c3f-a6ab-00c072c43f73" width="800">

**Generated script**

<img src="https://github.com/YuncongMa/pNet/assets/20191790/53139c7d-d5c0-4aaa-bda5-65a0f72c7a08" width="800">



## Workflow
The workflow code is in Workflow.py
This streamlines the setup and computation of pNet, including modules below
1. Data Input (Data_Input.py)
2. FN Computation (FN_Computation.py)
3. Quality Control (Quality_Control.py)

```
# This is a comprehensive setting for pNet workflow 
pNet.workflow(dir_pnet_result: str,
                 file_scan: str,
                 dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)',
                 file_subject_ID=None, file_subject_folder=None, file_group=None,
                 file_Brain_Template=None,
                 file_surfL=None, file_surfR=None, file_maskL=None, file_maskR=None,
                 file_mask_vol=None, file_overlayImage=None,
                 maskValue=0,
                 file_surfL_inflated=None, file_surfR_inflated=None,
                 K=17, Combine_Scan=False,
                 file_gFN=None,
                 samplingMethod='Subject', sampleSize=10, nBS=50,
                 maxIter=1000, minIter=30, meanFitRatio=0.1, error=1e-6, normW=1,
                 Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
                 Parallel=False, Computation_Mode='CPU_Torch', N_Thread=1,
                 dataPrecision='double'):
```

```
# This is a minimal version of pNet workflow for fast deployment
pNet.workflow_simple(dir_pnet_result: str,
                        dataType: str, dataFormat: str,
                        file_scan: str,
                        file_Brain_Template: str,
                        K=17,
                        Combine_Scan=False,
                        file_gFN=None)
```

## Module
1. **Data Input** <br />
Organize fMRI scans based on subject information, prepare brain template for subsequent FN computation and visualization
2. **FN Computation** <br />
Setup and carry out the FN computation, generating group-level and personalized FNs
3. **Quality Control** <br />
Ensures that pFNs have the highest spatial similarity to their group-level counterparts <br />

Visualization and statistics modules will come soon

## Additionals
### Example: It includes example results of pNet on HCP, UKBB and PNC datasets
```
# HCP formatted data, surface only
pNet.Example.HCP_Surf.dir_pnet
```
### Brain Template: It includes prepared brain template files for HCP, MNI space
```
# load the prepared brain template file for HCP formatted surface data
pNet.Brain_Template.HCP_surf
```
```
# Get the file directory of the prepared brain template for HCP formatted surface data
pNet.Brain_Template.file_HCP_surf
```

## Example script
Several simple examples are provided in Example_Workflow.py <br />

**Example for HCP format**
```
# This example is to perform pNet using surface-based fMRI in HCP format
# 1. Specify the result folder directory in dir_pnet_result
# 2. Provide a txt formatted scan list file, such as Scan_List.txt
# 3. Use a prepared brain template file provided in pNet
# 4. Choose the desired number of FNs

# Setup
# data type is surface, fixed for this example
dataType = 'Surface'
# data format is HCP surface, usually in CIFTI format, but can also be store as a 2D matrix in MAT file, fixed
dataFormat = 'HCP Surface (*.cifti, *.mat)'
# directory of the pNet result folder. Change <User> to a desired directory
dir_pnet_result = '<User>/Test_FN17_HCP_Workflow'
# a txt file storing directory of each fMRI scan file, required to provide
file_scan = '/Volumes/Scratch_0/pNet/Example/HCP_Surface/Data/Scan_List.txt'
# a built-in brain template file, made for the HCP surface data (59412 vertices for cortical gray matter), optional to change
file_Brain_Template = pNet.Brain_Template.file_HCP_surf
# number of FNs, can be changed to any positive integer number
K = 17

# Run pNet workflow
pNet.workflow_simple(
        dir_pnet_result=dir_pnet_result,
        dataType=dataType,
        dataFormat=dataFormat,
        file_scan=file_scan,
        file_Brain_Template=file_Brain_Template,
        K=K
)
```
**Example for NIFTI format**
```
# This example is to perform pNet using volume-based fMRI in NIFTI format, with co-registration to MNI space
# 1. Specify the result folder directory in dir_pnet_result
# 2. Provide a txt formatted scan list file, such as Scan_List.txt
# 3. Use a prepared brain template file provided in pNet
# 4. Choose the desired number of FNs

# Setup
# data type is volume, fixed for this example
dataType = 'Volume'
# data format is NIFTI, which stores a 4D matrix, fixed for this example
dataFormat = 'Volume (*.nii, *.nii.gz, *.mat)'
# directory of the pNet result folder. Change <User> to a desired directory
dir_pnet_result = '<User>/Test_FN17_UKBB_Workflow'
# a txt file storing directory of each fMRI scan file, required to provide
file_scan = '/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB/Data_Input/Scan_List.txt'
# a built-in brain template file, MNI standard space (2mm isotropic), made for the HCP surface data, optional to change
file_Brain_Template = pNet.Brain_Template.file_MNI_vol
# number of FNs, can be changed to any positive integer number
K = 17

# Run pNet workflow
pNet.workflow_simple(
        dir_pnet_result=dir_pnet_result,
        dataType=dataType,
        dataFormat=dataFormat,
        file_scan=file_scan,
        file_Brain_Template=file_Brain_Template,
        K=K
)
```

## Help
All the functions in pNet come with detailed description of each input <br />
Example help comment <br />

<img src="https://github.com/YuncongMa/pNet/assets/20191790/2b040a86-2212-4a48-9081-8d89983babf8" width="500">

