# This is a brief introduction of the MATLAB version
The MATLAB version offers both GUI and scripts.

## GUI
The MATLAB version features an intuitive user interface (/MATLAB/GUI/pNet.mlapp) to configure the workflow, carry out the computation and check the results. <br />
Video tutorial is availabel at https://github.com/YuncongMa/pNet/blob/main/MATLAB_Old/Video_Tutorial/Demo_Video.gif

<p align="center">
<img src="https://github.com/YuncongMa/pNet/assets/20191790/f5168794-63d2-4076-92c9-8ca3d78e478f" width="800">
</p>

(A) The data input module for loading fMRI scans and brain template files. <br />
(B) A module to setup computation parameters of our spatial-reguarlized NMF method. <br />
(C) Surface-based visualization of both group and personalized FNs (k=17) using HCP S1200 dataset, with left panel showing a binarized atlas generated from the group FNs and the right panel showing five views of three personalized FNs. All the color bars of intensity maps were set from the maximum value of the map to its half value. <br />
(D) Volume-based visualization of both group and personalized FNs (k=17) from Zhen’s multi-cohort iSTAGING study, with left panel showing a binarized functional atlas generated from the group FNs and the right panel showing a three-slice view of three personalized FNs. <br />
(E) Surface-based visualization of a personalized FN. <br />
(F) Surface-based visualization of the maximum t value (two sample t-test, FDR correction, p-value=0.001) of sex differences of individualized FNs (k=50) of the HCP S1200 dataset. <br />
(G) A module for quality control, showing one scan with two FNs mismatched to their group-level counterparts. <br />

## Manual
A comprehensive manual about the GUI can be found at https://github.com/YuncongMa/pNet/blob/main/MATLAB/MATLAB_Manual.pdf

# Installation
1. The MATLAB version requires MATLAB no older than 2021A. Its GUI version can be accessed by openning pNet.mlapp, or install the pNet.mlappinstall into MATLAB APPS. Compiled pNet is located in folder Compiled_Runtime. It requires the free MATLAB runtime software at https://www.mathworks.com/products/compiler/matlab-runtime.html. Details can be found in the README.txt in each subfolder containing the compiled pNet.
2. Additional installation information can be found in the help document (Help_Document.pdf).

# MATLAB and system Compatibility
Note: it is recommended to run pNet as a MATLAB APP to allow for maximum compatibility. Although we are striving for maximizing the compatibility of pNet on different MATLAB versions and operation systems, it is difficult to conduct comprehensive tests on all cases. Our primary test environment is MATLAB 2022B and macOS Ventura 13.

<p align="center">
<img src="https://user-images.githubusercontent.com/20191790/226790097-18b601fa-84ca-4ee6-aab6-19c8322ffc7d.jpg" width="600">
</p>

# Hardware requirement
It is recommended to run pNet with CPU at least 4 cores, memory (RAM) at least 16GB, and free disk storage at least 10GB. Large dataset and parallel computation will require more CPU cores, memory space, and disk space. CPU and memory inf, as well as estimated memory usage can be found in the parallel section.The GUI requires screen displaying resolution at least 1280x960. It is recommended to use an external montior with reosolution at least 1920x1080. 2k and 4k displays are desirable.
