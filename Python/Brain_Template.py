# Yuncong Ma, 8/25/2023
# Template information

#########################################
# Packages
import os
from Data_Input import load_json_setting
#########################################

# Get the directory of pNet based the location of this file
current_file_path = os.path.abspath(__file__)
dir_python_package = os.path.dirname(current_file_path)
dir_pNet = os.path.dirname(dir_python_package)

#########################################
dir_Template = os.path.join(dir_pNet, 'Brain_Template')
# Organize example into a class variable

class Brain_Template:

    HCP_surf = load_json_setting(os.path.join(dir_Template, 'HCP_Surface', 'Brain_Template.json'))
    FS_surf = load_json_setting(os.path.join(dir_Template, 'FreeSurfer_fsaverage5', 'Brain_Template.json'))
    MNI_vol = load_json_setting(os.path.join(dir_Template, 'MNI', 'Brain_Template.json'))

