![Screenshot 2023-03-06 at 2 23 52 PM](https://user-images.githubusercontent.com/20191790/223320724-67f18e53-cda3-41a4-9215-620c991877f4.jpg)


# pNet <br /> A toolbox for personalized functional network modeling <br />

This toolbox is designed to provide a user-friendly interface to perform personalized functional network computation <br />
Major functions:
1. search and organize fMRI files
2. compute or load pre-computed group-level functional network (FN) results
3. compute personalized FNs, with quality control
4. perform statistical analysis to investigate the correlation between FN loading and behavior measurement
5. provide both pre-computed and interactive visualization for group and personalized FNs and statistical results
6. provide MATLAB functions and scripts to carry out computation and visualization of group and personalized FNs
7. provide Python functions and scripts to carry out the computation
8. come with a help document that can be accessed directly on the app UI

# Module
![Screenshot 2023-06-29 at 10 57 59 AM](https://github.com/YuncongMa/pNet/assets/20191790/9546d31e-aa9a-4a9c-8558-a875c92da27e)






# GUI workflow

![Screenshot 2023-06-29 at 10 58 25 AM](https://github.com/YuncongMa/pNet/assets/20191790/8bd6a580-dd3f-4003-b92e-487dd5c3bf11)



It consists of six major modules, including data loading, FN computation, group-level FN display, personalized FN display, statistics, and quality control. The six modules correspond to six major UI panels in pNet, and six subfolders in the pNet result folder. The pNet UI is designed to follow the left-to-right workflow direction to setup each step. Curved arrows show alternative workflows. The FN computation module (colored in light red) can use precomputed group-level FNs to derive personalized FNs as one option. Precomputed visualization of group-level and personalized FNs will be used in the two display modules (group FN and personalized FN). pNet features synchronized display to facilitate the comparisons between group-level and personalized FNs, as well as statistical results. The dashed double-sided arrows denote this display synchronization between three modules.

# Snapshot
![Screenshot 2023-03-06 at 4 46 53 PM](https://user-images.githubusercontent.com/20191790/223241753-b5a0a300-480a-4397-8585-5874f91c6590.jpg)

(A) Welcome page for the toolbox. <br />
(B) A module for loading fMRI scans and brain template files. <br />
(C) A module to setup computation parameters for both the group-level and individualized FNs. <br />
(D) Surface-based visualization of both group and personalized FNs (k=17) using HCP S1200 dataset, with left panel showing a binarized atlas generated from the group FNs and the right panel showing five views of three personalized FNs. All the color bars of intensity maps were set from the maximum value of the map to its half value. <br />
(E) Volume-based visualization of both group and personalized FNs (k=17) from Zhen’s multi-cohort iSTAGING study, with left panel showing a binarized functional atlas generated from the group FNs and the right panel showing a three-slice view of three personalized FNs. <br />
(F) Surface-based visualization of the maximum t value (two sample t-test, FDR correction, p-value=0.001) of sex differences of individualized FNs (k=50) of the HCP S1200 dataset. <br />
(G) A module for quality control, showing one scan with two FNs mismatched to their group-level counterparts. <br />

# Data Structure
![Screenshot 2023-06-29 at 10 58 52 AM](https://github.com/YuncongMa/pNet/assets/20191790/a72263c5-0a0d-4507-89f9-f22ba80c9853)



Folder are noted in boxes with black borders above the black horizontal line, and files are in boxes without black border, separated into three categories (setting, data, and figure files) by the dashed lines. Boxes are color coded for different modules. In the first three modules, the setting files are present. In personalized FN folder, subfolders are named by subject information. In statistics folder, subfolders are named by their method information. 

# Guideline
1. A help document ([Help_Document.pdf](https://github.com/YuncongMa/pNet/blob/main/Help_Document.pdf)) is provided in the root directory and can be accessed in the app UI. 
2. Help videos are provided on a Google drive folder: https://drive.google.com/drive/folders/1aVU99aJlFJZ5kKw2KW3wkwZKv6TBkgu_?usp=share_link

# Installation
1. Download the whole package and run the pNet.mlapp in MATLAB with version not older than 2021A. Or install the pNet.mlappinstall into MATLAB APPS. Also we provide several versions of compiled pNet (in folder Compiled_Runtime) for different operation systems. It requires the free MATLAB runtime software at https://www.mathworks.com/products/compiler/matlab-runtime.html.  Details can be found in the README.txt in each subfolder containing the compiled pNet.
2. The Python version requires version=3.8, with additional tools installed via the script 'Install_Toolbox.py'. pNet offers computation code built with either NumPy or PyTorch. The PyTorch version has significantly faster speed and capability on GPU accleration.
3. Additional installation information can be found in the help document (Help_Document.pdf).
4. HCP data (.cifti) requires additional pre-compiled packages to read, please follow the section below "Additional Package"
5. Brain template files are stored in subfolder "Brain_Template". It includeds templates for three different brain formats. Please check the section below "Brain Template".
6. Several examples can be downloaded from Google Drive, please follow the section below "Example".

# Hardware requirement
It is recommended to run pNet with CPU at least 4 cores, memory (RAM) at least 16GB, and free disk storage at least 10GB. Large dataset and parallel computation will require more CPU cores, memory space, and disk space. CPU and memory inf, as well as estimated memory usage can be found in the parallel section.The pNet UI requires screen displaying resolution at least 1280x960. It is recommended to use an external montior with reosolution at least 1920x1080. 2k and 4k displays are desirable.


# Brain Template
Note: three widely used brain templates are available in the toolbox subfolder "Brain_Template"
1. HCP surface space: subfolder "HCP_Surface". It contains a single MATLAB file "Brain_Surface.mat" which includes the brain shapes, medial wall index.
2. FreeSurfer fsaverage5: "FreeSurfer_ fsaverage5". It contains a single MATLAB file "Brain_Surface.mat" which includes the brain shapes, medial wall index.
3. MNI volume space: subfolder "MNI_Volume". It contains two MATLAB files "Brain_Mask.mat" and "Overlay_Image.mat".

# Example
We provided three examples for users to learn how to run the pFN computation, and navigate through the precomputed results. For each example, we provides examples of simulated fMRI data, and precomputed results of real fMRI data. Examples can be accessed on Google Drive: https://drive.google.com/drive/folders/1xkCy-0WqYvPA9ooq8txTdc0GsGC-YXMq?usp=share_link.   <br />

1.	HCP surface data. (2 subject, 2 scans per subject, 400 volumes per scan) in .mat format. Precomputed results of 10 subject data are stored in subfolder "Test_FN17".
2.	PNC surface data (2 subject, 1 scan per subject, 300 volumes per scan) in .mgh format. Precomputed results of 10 subject data are stored in subfolder "Test_FN17".
3.	UKBB volume data (2 subject, 1 scan per subject, 200 volumes per scan) in .nii.gz format. Precomputed results of 10 subject data are stored in subfolder "Test_FN17".



The data organization is as below.

![Screenshot 2023-03-06 at 11 29 38 PM](https://user-images.githubusercontent.com/20191790/223320985-12aed4d4-6bef-4b23-a9a2-ff67f79eced3.jpg)

# MATLAB and system Compatibility
Note: it is recommended to run pNet as a MATLAB APP to allow for maximum compatibility. Although we are striving for maximizing the compatibility of pNet on different MATLAB versions and operation systems, it is difficult to conduct comprehensive tests on all cases. Our primary test environment is MATLAB 2022B and macOS Ventura 13.

![Screenshot 2023-03-21 at 10 53 56 PM](https://user-images.githubusercontent.com/20191790/226790097-18b601fa-84ca-4ee6-aab6-19c8322ffc7d.jpg)





# Reference
[1] Cui, Z. (2020). Individual variation in functional topography of association networks in youth. Neuron, 106(2), 340-353. <br />
[2] Cui, Z. (2022). Linking Individual Differences in Personalized Functional Network Topography to Psychopathology in Youth. Biological Psychiatry. <br />
[3] Li, H. and Fan, Y., 2016, April. Individualized brain parcellation with integrated funcitonal and morphological information. In 2016 IEEE 13th International Symposium on Biomedical Imaging (ISBI) (pp. 992-995). IEEE. <br />
[4] Li H. (2017) Large-scale sparse functional networks from resting state fMRI. Neuroimage 156:1-13. <br />
[5] Shanmugan, S. (2022). Sex differences in the functional topography of association networks in youth. Proceedings of the National Academy of Sciences, 119(33), e2110416119. <br />
[6] Zhou, Z. (2023). Multiscale functional connectivity patterns of the aging brain learned from harmonized rsfMRI data of the multi-cohort iSTAGING study. NeuroImage, 269, 119911. <br />

