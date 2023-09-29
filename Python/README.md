# This is a brief description of the Python version of pNet
The Python version is designed to carry out the pNet setup and computation, consistent to the MATLAB version with GUI. <br />
It has Numpy and PyTorch versions for easy code development and fast computation. <br />
pNet.py is the main function to import all required packages and functions. <br />

## Deployment
This Python package requires additional toolboxes
Specifically, it requires:
1. numpy
2. pytorch
3. scipy
4. nibabel
5. h5py
```
# Install required toolboxes
python Install_Toolbox.py
```
```
# Add pNet
import pNet
```

## Workflow
The workflow code is in Workflow.py
This streamlines the setup and computation of pNet, including modules below
1. Data Input
2. FN Computation
3. Quality Control.

```
# This is a comprehensive setting for 
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
                 Compute_gFN=True, file_gFN=None,
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
                        Compute_gFN=True,
                        file_gFN=None)
```

## Module
1. **Data Input** <br />
Organize fMRI scans based on subject information, prepare brain template for subsequent FN computation and visualization
2. **FN Computation** <br />
Setup and carry out the FN computation, generating group-level and personalized FNs
3. **Quality Control** <br />
Ensures that pFNs have the highest spatial similarity to their group-level counterparts

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
# Get the file directory the prepared brain template for HCP formatted surface data
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
dir_pnet_result = '<User>/Test_FN17_HCP_Workflow'  # Change <User> to a desired directory
file_scan = '/Volumes/Scratch_0/pNet/Example/HCP_Surface/Data/Scan_List.txt'
file_Brain_Template = pNet.Brain_Template.file_HCP_surf
K = 17

# Run pNet workflow
pNet.workflow_simple(
        dir_pnet_result=dir_pnet_result,
        dataType='Surface',
        dataFormat='HCP Surface (*.cifti, *.mat)',
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
dir_pnet_result = '/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB_Workflow'
file_scan = '/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB/Data_Input/Scan_List.txt'
file_Brain_Template = pNet.Brain_Template.file_MNI_vol
K = 17

# Run pNet workflow
pNet.workflow_simple(
        dir_pnet_result=dir_pnet_result,
        dataType='Volume',
        dataFormat='Volume (*.nii, *.nii.gz, *.mat)',
        file_scan=file_scan,
        file_Brain_Template=file_Brain_Template,
        K=K
)
```

## Step-by-step guidance
This Python version of pNet offers a step-by-step guided setup to generate a Python script for running a desired workflow <br />
```
pNet.workflow_guide()
```
Example terminal interaction <br />
<img src="https://github.com/YuncongMa/pNet/assets/20191790/cc0fcfef-67c3-4375-a586-1eb7d3a51b10" width="400">


## Help
All the functions in pNet come with detailed description of each input <br />
Example help comment <br />
<img src="https://github.com/YuncongMa/pNet/assets/20191790/fd0a95b7-9eb6-46e3-9503-eb43b13f63b7" width="500">

