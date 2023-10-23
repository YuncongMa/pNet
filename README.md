
<p align="center">
<img src="https://github.com/YuncongMa/pNet/assets/20191790/dccff00b-424c-47fb-a3fe-351e2a75bd11"  width="800">
</p>

# pNet <br /> A toolbox for personalized functional network modeling <br />
<p align="center">
<img width="600" src="https://github.com/YuncongMa/pNet/assets/20191790/ff466c7f-42a1-4020-9691-311b6aaf4537">
</p>

* This toolbox is designed to extract personalized functional networks from fMRI data using spatial-regularized non-negative matrix factorization (**NMF**) method (https://doi.org/10.1016/j.neuroimage.2017.05.004) <br />
* It works with both MATLAB and Python (https://github.com/YuncongMa/pNet/tree/main/Python), and comes with a user-friendly GUI interface in MATLAB, as well as a step-by-step guide to setup a customized workflow in Python.<br />
* It supports both surface and volume based fMRI data, as well as the grayordinate which combines the two types. <br />
* It supports multiple fMRI formats, including CIFTI, MGH, MGZ, NIFTI, and MAT.<br />
* It supports multi-cohort datasets and precomputed group-level FNs.<br />


# GUI
The MATLAB version features an intuitive user interface to configure the workflow, carry out the computation and check the results.

<p align="center">
<img src="https://github.com/YuncongMa/pNet/assets/20191790/f5168794-63d2-4076-92c9-8ca3d78e478f" width="800">
</p>

(A) The data input module for loading fMRI scans and brain template files. <br />
(B) A module to setup computation parameters for both the group-level and individualized FNs. <br />
(C) Surface-based visualization of both group and personalized FNs (k=17) using HCP S1200 dataset, with left panel showing a binarized atlas generated from the group FNs and the right panel showing five views of three personalized FNs. All the color bars of intensity maps were set from the maximum value of the map to its half value. <br />
(D) Volume-based visualization of both group and personalized FNs (k=17) from Zhen’s multi-cohort iSTAGING study, with left panel showing a binarized functional atlas generated from the group FNs and the right panel showing a three-slice view of three personalized FNs. <br />
(E) Surface-based visualization of a personalized FN. <br />
(F) Surface-based visualization of the maximum t value (two sample t-test, FDR correction, p-value=0.001) of sex differences of individualized FNs (k=50) of the HCP S1200 dataset. <br />
(G) A module for quality control, showing one scan with two FNs mismatched to their group-level counterparts. <br />


# Step-by-step guidance
The Python version features a customized workflow using command line
This will generate a Python script to run a workflow

<p align="center">
<img src="https://github.com/YuncongMa/pNet/assets/20191790/38efac5d-631a-47d8-a34f-fd686f84602d" width="800">
</p>

# Manual
A PDF manual can be accessed below and from the MATLAB GUI <br />
https://github.com/YuncongMa/pNet/blob/main/Help_Document.pdf


# Modules and workflow
pNet consists of five main modules and a customizable workflow.

<p align="center">
<img src="https://github.com/YuncongMa/pNet/assets/20191790/f83fa661-3a2f-4aae-9e70-bdbea3b2fdbb"  width="600">
</p>


* Data input module loads preprocessed fMRI data, with corresponding brain template files to extract coordination information and facilitate visualization. 
* FN computation module is to carry out computation for group-level and personalized FNs. 
* Visualization module can provide both preconfigured and interative displays
* Quality control module will ensure the reliability of the whole computation. 
* Additionally, statistical analyses can be employed to investigate differences in personalized FNs between groups or their relation to behavioral traits.



# Data Structure
pNet organizes data structure with respect to the six panels of its GUI version: Data Input, FN Computation, Group FN, Personalized FN, Statistics, Quality Control.

<p align="center">
<img src="https://github.com/YuncongMa/pNet/assets/20191790/a72263c5-0a0d-4507-89f9-f22ba80c9853" width="800">
</p>

Folder are noted in boxes with black borders above the black horizontal line, and files are in boxes without black border, separated into three categories (setting, data, and figure files) by the dashed lines. Boxes are color coded for different modules. In the first three modules, the setting files are present. In personalized FN folder, subfolders are named by subject information. In statistics folder, subfolders are named by their method information. 


# Installation
1. The MATLAB version requires MATLAB no older than 2021A. Its GUI version can be accessed by openning pNet.mlapp, or install the pNet.mlappinstall into MATLAB APPS. Compiled pNet is located in folder Compiled_Runtime. It requires the free MATLAB runtime software at https://www.mathworks.com/products/compiler/matlab-runtime.html. Details can be found in the README.txt in each subfolder containing the compiled pNet.
2. The Python version is available at https://github.com/YuncongMa/pNet/tree/main/Python.
3. Additional installation information can be found in the help document (Help_Document.pdf).
4. Built-in brain template files are stored in subfolder "Brain_Template".

# Hardware requirement
It is recommended to run pNet with CPU at least 4 cores, memory (RAM) at least 16GB, and free disk storage at least 10GB. Large dataset and parallel computation will require more CPU cores, memory space, and disk space. CPU and memory inf, as well as estimated memory usage can be found in the parallel section.The pNet UI requires screen displaying resolution at least 1280x960. It is recommended to use an external montior with reosolution at least 1920x1080. 2k and 4k displays are desirable.


# Brain Template
Five built-in brain templates are available in the subfolder "Brain_Template"
1. HCP surface: subfolder "HCP_Surface". It contains 3D mesh shapes (vertices and faces) and brain masks for two hemishperes.
2. FreeSurfer fsaverage5: "FreeSurfer_fsaverage5". It is similar to HCP surface .
3. MNI volume space: subfolder "MNI_Volume". It contains two MATLAB files "Brain_Mask.mat" and "Overlay_Image.mat".
4. HCP surface-volume: It contains both cortical surface information, and subcortical volume.
5. HCP volume: It is similar to MNI volume space.

# Example
We provided three examples for users to learn how to run the pFN computation, and navigate through the precomputed results. For each example, we provides examples of simulated fMRI data, and precomputed results of real fMRI data. Examples can be accessed on Google Drive: https://drive.google.com/drive/folders/1xkCy-0WqYvPA9ooq8txTdc0GsGC-YXMq?usp=share_link.   <br />

1.	HCP surface data. (2 subject, 2 scans per subject, 400 volumes per scan) in .mat format. Precomputed results of 10 subject data are stored in subfolder "Test_FN17".
2.	PNC surface data (2 subject, 1 scan per subject, 300 volumes per scan) in .mgh format. Precomputed results of 10 subject data are stored in subfolder "Test_FN17".
3.	UKBB volume data (2 subject, 1 scan per subject, 200 volumes per scan) in .nii.gz format. Precomputed results of 10 subject data are stored in subfolder "Test_FN17".



The data organization is as below.

<p align="center">
<img src="https://user-images.githubusercontent.com/20191790/223320985-12aed4d4-6bef-4b23-a9a2-ff67f79eced3.jpg" width="600">
</p>

# MATLAB and system Compatibility
Note: it is recommended to run pNet as a MATLAB APP to allow for maximum compatibility. Although we are striving for maximizing the compatibility of pNet on different MATLAB versions and operation systems, it is difficult to conduct comprehensive tests on all cases. Our primary test environment is MATLAB 2022B and macOS Ventura 13.

<p align="center">
<img src="https://user-images.githubusercontent.com/20191790/226790097-18b601fa-84ca-4ee6-aab6-19c8322ffc7d.jpg" width="600">
</p>





# Reference
[1] Cui, Z. (2020). Individual variation in functional topography of association networks in youth. Neuron, 106(2), 340-353. <br />
[2] Cui, Z. (2022). Linking Individual Differences in Personalized Functional Network Topography to Psychopathology in Youth. Biological Psychiatry. <br />
[3] Li, H. and Fan, Y., 2016, April. Individualized brain parcellation with integrated funcitonal and morphological information. In 2016 IEEE 13th International Symposium on Biomedical Imaging (ISBI) (pp. 992-995). IEEE. <br />
[4] Li H. (2017) Large-scale sparse functional networks from resting state fMRI. Neuroimage 156:1-13. <br />
[5] Shanmugan, S. (2022). Sex differences in the functional topography of association networks in youth. Proceedings of the National Academy of Sciences, 119(33), e2110416119. <br />
[6] Zhou, Z. (2023). Multiscale functional connectivity patterns of the aging brain learned from harmonized rsfMRI data of the multi-cohort iSTAGING study. NeuroImage, 269, 119911. <br />

