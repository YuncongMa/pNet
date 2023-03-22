![Screenshot 2023-03-06 at 2 23 52 PM](https://user-images.githubusercontent.com/20191790/223320724-67f18e53-cda3-41a4-9215-620c991877f4.jpg)


# pNet <br /> A toolbox for personalized functional network modeling <br />

This toolbox is designed to provide a user-friendly interface to perform personalized functional network computation <br />
Major functions:
1. search and organize fMRI files
2. compute or load pre-computed group-level functional network (FN) results
3. compute personalized FNs, with quality control
4. perform statistical analysis to investigate the correlation between FN loading and behavior measurement
5. provide both pre-computed and interactive visualization for group and personalized FNs and statistical resutls
6. provide MATLAB functions and scripts to carry out computation and visualization of group and personalized FNs
7. come with a help document that can be accessed directly on the app UI

# Snapshot
![Screenshot 2023-03-06 at 4 46 53 PM](https://user-images.githubusercontent.com/20191790/223241753-b5a0a300-480a-4397-8585-5874f91c6590.jpg)

(A) Welcome page for the toolbox. <br />
(B) A module for loading fMRI scans and brain template files. <br />
(C) A module to setup computation parameters for both the group-level and individualized FNs. <br />
(D) Surface-based visualization of both group and personalized FNs (k=17) using HCP S1200 dataset, with left panel showing a binarized atlas generated from the group FNs and the right panel showing five views of three personalized FNs. All the color bars of intensity maps were set from the maximum value of the map to its half value. <br />
(E) Volume-based visualization of both group and personalized FNs (k=17) from Zhen’s multi-cohort iSTAGING study, with left panel showing a binarized functional atlas generated from the group FNs and the right panel showing a three-slice view of three personalized FNs. <br />
(F) Surface-based visualization of the maximum t value (two sample t-test, FDR correction, p-value=0.001) of sex differences of individualized FNs (k=50) of the HCP S1200 dataset. <br />
(G) A module for quality control, showing one scan with two FNs mismatched to their group-level counterparts. <br />

# Guideline
1. A help document (Help_Document.pdf) is provided in the root directory and can be accessed in the app UI. It can also be accessed from Google drive: https://drive.google.com/file/d/1ArWRLM72tTgQbAo0P84v-qkoIYYwdSrF/view?usp=share_link
2. Help videos are provided on a Google drive folder: https://drive.google.com/drive/folders/1aVU99aJlFJZ5kKw2KW3wkwZKv6TBkgu_?usp=share_link

# Installation
1. Download the whole package and open the pNet.mlapp in MATLAB with version not older than 2021A. Or add the whole folder to the path, and run pNet in the command window.
2. Additional installation information can be found in the help document (Help_Document.pdf).
3. HCP data (.cifti) requires additional pre-compiled packages to read, please follow the section below "Additional Package"
4. Brain template files are stored in subfolder "Brain_Template". It includeds templates for three different brain formats. Please check the section below "Brain Template".
5. Several examples can be downloaded from Google Drive, please follow the section below "Example".

# Hardware requirement
It is recommended to run pNet with CPU at least 4 cores, memory (RAM) at least 16GB, and free disk storage at least 10GB. Large dataset and parallel computation will require more CPU cores, memory space, and disk space. CPU and memory inf, as well as estimated memory usage can be found in the parallel section.

# Additional Package
1. To allow the input of HCP data in cifti format, download workbench from HCP website: https://www.humanconnectome.org/software/get-connectome-workbench

# Brain Template
Note: three widely used brain templates are available in the toolbox subfolder "Brain_Template"
1. HCP surface space: subfolder "HCP_Surface". It contains a single MATLAB file "Brain_Surface.mat" which includes the brain shapes, medial wall index.
2. FreeSurfer fsaverage5: "FreeSurfer_ fsaverage5". It contains a single MATLAB file "Brain_Surface.mat" which includes the brain shapes, medial wall index.
3. MNI volume space: subfolder "MNI_Volume". It contains two MATLAB files "Brain_Mask.mat" and "Overlay_Image.mat".

# Example
Note: fMRI data can easily take more than 1GB space, so we provide examples on Google Drive: https://drive.google.com/drive/folders/1xkCy-0WqYvPA9ooq8txTdc0GsGC-YXMq?usp=share_link  <br />

1. HCP surface data (10 subjects, 2 scans per subject, 400 volumes per scan) in .mat format. pNet results are stored in subfolder "Test_FN17".
2. PNC surface data (10 subjects, 1 scan per subject, 555 volumes per scan) in .mgh format. pNet results are stored in subfolder "Test_FN17".
3. UKBB volume data (10 subjects, 1 scan per subject, 490 volumes per scan) in .nii.gz format. pNet results are stored in subfolder "Test_FN17".

The data organization is as below.

![Screenshot 2023-03-06 at 11 29 38 PM](https://user-images.githubusercontent.com/20191790/223320985-12aed4d4-6bef-4b23-a9a2-ff67f79eced3.jpg)

# MATLAB and system Compatibility
Note: it is recommended to run pNet as a MATLAB function to allow for maximum compatibility. Although we are striving for maximizing the compatibility of pNet on different MATLAB versions and operation systems, it is difficult to conduct comprehensive tests on all cases. Our primary test environment is MATLAB 2022B and macOS Ventura 13.

![Screenshot 2023-03-21 at 10 53 56 PM](https://user-images.githubusercontent.com/20191790/226790097-18b601fa-84ca-4ee6-aab6-19c8322ffc7d.jpg)





# Reference
[1] Cui, Z. (2020). Individual variation in functional topography of association networks in youth. Neuron, 106(2), 340-353. <br />
[2] Cui, Z. (2022). Linking Individual Differences in Personalized Functional Network Topography to Psychopathology in Youth. Biological Psychiatry. <br />
[3] Li H. (2017) Large-scale sparse functional networks from resting state fMRI. Neuroimage 156:1-13. <br />
[4] Pines, A. R. (2022). Dissociable multi-scale patterns of development in personalized brain networks. Nature communications, 13(1), 1-15. <br />
[5] Shanmugan, S. (2022). Sex differences in the functional topography of association networks in youth. Proceedings of the National Academy of Sciences, 119(33), e2110416119. <br />
[6] Zhou, Z. (2023). Multiscale functional connectivity patterns of the aging brain learned from harmonized rsfMRI data of the multi-cohort iSTAGING study. NeuroImage, 269, 119911. <br />

