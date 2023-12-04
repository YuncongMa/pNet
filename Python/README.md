# This is a brief description for the Python version of pNet


* It is designed to carry out the pNet setup and computation, consistent to the MATLAB version with GUI. 
* It features a step-by-step guidance to configure a customized script to run a workflow.
* It provides both Numpy and PyTorch versions for easy code development and fast computation.

## Deployment

This Python version is tested in Anaconda (https://www.anaconda.com) to create an environment. <br />
**Create a new conda environment for pNet**
```
conda env create --name pnet -f environment.yml
```
**Or install required tools in conda using requirements.txt**
```
pip install -r requirements.txt
```
**Or install required tools directly**
```
pip install numpy scipy scikit-learn pandas h5py
conda install -c conda-forge nibabel
conda install pytorch::pytorch torchvision torchaudio -c pytorch
conda install -c conda-forge vtk
pip install matplotlib surfplot
python -m pip install -U scikit-image
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

<p align="center">
<img src="https://github.com/YuncongMa/pNet/assets/20191790/38efac5d-631a-47d8-a34f-fd686f84602d" width="800">
</p>

**Generated script**

<p align="center">
<img src="https://github.com/YuncongMa/pNet/assets/20191790/96df4f6d-a24c-4d7c-afb8-482cba5ad682" width="800">
</p>

## Workflow
The workflow code is in Workflow.py
This streamlines the setup and computation of pNet, including modules below
1. Data Input (Data_Input.py)
2. FN Computation (FN_Computation.py)
3. Visualization (Visualization.py)
4. Quality Control (Quality_Control.py)

```
# This is a comprehensive setting for pNet workflow 
pNet.workflow(dir_pnet_result: str,
             file_scan: str,
             dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)',
             file_subject_ID=None, file_subject_folder=None, file_group_ID=None,
             file_Brain_Template=None,
             templateFormat='HCP',
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
             dataPrecision='double',
             outputFormat='Both')
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

**Server mode**

pNet offers a server mode for its Python version to deploy customized workflows with minimal knowledge about bash jobs
This server mode requires only additional setups for a few server commands to submit bash jobs, as shown below
1. submit_command = 'qsub -terse -j y'
2. thread_command = '-pe threaded '
3. memory_command = '-l h_vmem='
4. log_command = '-o '

Example Python script
https://github.com/YuncongMa/pNet/blob/main/Python/Example_Workflow_Server.py

## Module
1. **Data Input** <br />
Organize fMRI scans based on subject information, prepare brain template for subsequent FN computation and visualization
2. **FN Computation** <br />
Setup and carry out the FN computation, generating group-level and personalized FNs
3. **Visualization** <br />
Provide preconfigured visualization for surface and volume data types
4. **Quality Control** <br />
Ensures that pFNs have the highest spatial similarity to their group-level counterparts <br />

Statistics module will come soon



## Additionals

### Data format
The Python version also supports fMRI data in HCP, MGH, MGZ, NIFTI, MAT. <br />
It supports surface, volume, surface-volume data types. Surface-volume data type is based on the grayordinate of HCP data format which stores cortical surface and subcortical volume data in one file (https://www.jch.com/jch/notes/RestingState2014/grayordinates.html).

### Example: It includes example results of pNet on HCP, UKBB and PNC datasets
```
# HCP formatted data, surface only
pNet.Example.HCP_Surf.dir_pnet
```

### Brain Template: It includes prepared brain template files for HCP, MNI space
Brain templates are stored at https://github.com/YuncongMa/pNet/tree/main/Brain_Template
Users can access brain template files in .mat and .json.zip formats, or from the python version of pNet
```
# Get the file directory of the prepared brain template for HCP formatted surface data
pNet.Brain_Template.file_HCP_surf
```

## Script examples
Several simple examples are provided in Example_Workflow.py <br />

**HCP format**
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
**NIFTI format**
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

### Visualization examples
The Python version of pNet offers preconfigured visualizations consistent to the MATLAB version. However, the 3D surface display is slightly different due to different lighting and material settings between MATLAB and matplotlib.

**Surface type**

<img src= "https://github.com/YuncongMa/pNet/assets/20191790/0b69d0c9-d7f5-4a00-9ed3-bad6873afdf0"  width="800">

**Volume type**

<img src= "https://github.com/YuncongMa/pNet/assets/20191790/ca55d5c4-8757-4d2e-9818-57c6e9b0a7e4"  width="800">

## Help
All the functions in pNet come with detailed description of each input <br />
Example help comment <br />

<img src="https://github.com/YuncongMa/pNet/assets/20191790/2b040a86-2212-4a48-9081-8d89983babf8" width="500">

