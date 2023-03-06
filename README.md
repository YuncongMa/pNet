# pNet
# A toolbox for personalized functional network modeling

![Welcome](https://user-images.githubusercontent.com/20191790/222938004-af056d30-1ddd-4e35-86ed-bcb0c31b7094.jpg)

This toolbox is designed to provide a user-friendly interface to perform personalized functional network computation <br />
Major functions:
1. search and organize fMRI files
2. compute or load pre-computed group-level functional network (FN) results
3. compute personalized FNs, with quality control
4. perform statistical analysis to investigate the correlation between FN loading and behavior measurement
5. provide both pre-computed and interactive visualization for group and personalized FNs and statistical resutls
6. provide MATLAB functions and scripts to carry out computation and visualization of group and personalized FNs
7. come with built-in help documents and video

# Snapshots
![Screenshot 2023-03-06 at 10 31 39 AM](https://user-images.githubusercontent.com/20191790/223155797-ab8f8b7c-e528-4e2a-a6f0-946e7f9c2694.jpg)

(A) Welcome page for the toolbox. <br />
(B) A module for loading fMRI scans and brain template files. <br />
(C) A module to setup computation parameters for both the group-level and individualized FNs. <br />
(D) Surface-based visualization of both group and personalized FNs (k=17) using HCP S1200 dataset, with left panel showing a binarized atlas generated from the group FNs and the right panel showing five views of three personalized FNs. All the color bars of intensity maps were set from the maximum value of the map to its half value. <br />
(E) Volume-based visualization of both group and personalized FNs (k=17) from Zhen’s multi-cohort iSTAGING study, with left panel showing a binarized atlas generated from the group FNs and the right panel showing a three-slice view of three personalized FNs. <br />
(F) Surface-based visualization of the maximum t value (two sample t-test, FDR correction, p-value=0.001) of sex differences of individualized FNs (k=50) of the HCP S1200 dataset. <br />
(G) A module for quality control, showing one scan with two FNs mismatched to their group-level counterparts. <br />


# Additional package needed
1. To allow the input of HCP data in cifti format, downlaod workbench from HCP website: https://www.humanconnectome.org/software/get-connectome-workbench

# Reference
[1] Cui, Z. (2020). Individual variation in functional topography of association networks in youth. Neuron, 106(2), 340-353. <br />
[2] Cui, Z. (2022). Linking Individual Differences in Personalized Functional Network Topography to Psychopathology in Youth. Biological Psychiatry. <br />
[3] Li H. (2017) Large-scale sparse functional networks from resting state fMRI. Neuroimage 156:1-13. <br />
[4] Pines, A. R. (2022). Dissociable multi-scale patterns of development in personalized brain networks. Nature communications, 13(1), 1-15. <br />
[5] Shanmugan, S. (2022). Sex differences in the functional topography of association networks in youth. Proceedings of the National Academy of Sciences, 119(33), e2110416119. <br />
[6] Zhou, Z. (2023). Multiscale functional connectivity patterns of the aging brain learned from harmonized rsfMRI data of the multi-cohort iSTAGING study. NeuroImage, 269, 119911. <br />

