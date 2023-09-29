# Install required packages
# Yuncong Ma, 9/29/2023

import os

os.system('pip install numpy')
os.system('pip install scipy')
os.system('pip install h5py')
os.system('conda install -c conda-forge nibabel')
os.system('conda install pytorch::pytorch torchvision torchaudio -c pytorch')
