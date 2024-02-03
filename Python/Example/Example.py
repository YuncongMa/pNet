# Yuncong Ma, 9/27/2023
# Example data information

#########################################
# Packages
import os

#########################################

# Get the directory of pNet based the location of this file
current_file_path = os.path.abspath(__file__)
dir_python_package = os.path.dirname(current_file_path)
dir_pNet = os.path.dirname(dir_python_package)

#########################################
# Examples for data and pnet results
# Preprocessed fMRI data in surface and volume formats for testing
dir_example = os.path.join(dir_pNet, '')
# HCP surface data
dir_example_hcp_surf = os.path.join(dir_example, 'HCP_Surface')
dir_example_data_hcp_surf = os.path.join(dir_example_hcp_surf, 'Data')
dir_example_pnet_hcp_surf = os.path.join(dir_example_hcp_surf, 'Test_FN17')
# PNC surface data
dir_example_pnc_surf = os.path.join(dir_example, 'PNC_Surface')
dir_example_data_pnc_surf = os.path.join(dir_example_pnc_surf, 'Data')
dir_example_pnet_pnc_surf=os.path.join(dir_example_pnc_surf, 'Test_FN17')
# UKBB volume data
dir_example_ukbb_vol = os.path.join(dir_example, 'UKBB_Volume')
dir_example_data_ukbb_vol = os.path.join(dir_example_ukbb_vol, 'Data')
dir_example_pnet_ukbb_vol = os.path.join(dir_example_ukbb_vol, 'Test_FN17')

# Organize example into a class variable


class Example:

    # HCP formatted data, surface only
    class HCP_surf:
        dir = dir_example_hcp_surf
        dir_data = dir_example_data_hcp_surf
        dir_pnet = dir_example_pnet_hcp_surf

        # load_gFN(K, dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)')

    # MGH formatted data, surface only
    class PNC_surf:
        dir = dir_example_pnc_surf
        dir_data = dir_example_data_pnc_surf
        dir_pnet = dir_example_pnet_pnc_surf

    # NIFTI formatted data, volume type, co-registered in MNI space
    class UKBB_vol:
        dir = dir_example_ukbb_vol
        dir_data = dir_example_data_ukbb_vol
        dir_pnet = dir_example_pnet_ukbb_vol
