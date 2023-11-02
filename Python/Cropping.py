# Yuncong Ma, 11/2/2023
# Python version of cropping functions written in MATLAB

import numpy as np
from scipy.ndimage import label


def fTruncate_Image_3D_4D(Image_3D_4D, Voxel_Size, Extend):
    """
    Truncate 3D or 4D image based on mask and extend parameters.

    Parameters:
    Image_3D_4D : ndarray
        The input 3D or 4D image.
    Voxel_Size : list or tuple of length 3
        The size of each voxel in the image.
    Extend : list or tuple of length 3
        The extend parameters in each dimension, could be positive, negative, or inf.

    Returns:
    Image_3D_4D_Truncated : ndarray
        The truncated image.
    Center : list of length 3
        The center of the FOV after truncation.
    Crop_Parameter : dict
        Dictionary containing 'FOV_Old' and 'FOV' as keys.
    """
    Size = Image_3D_4D.shape
    Crop_Parameter = {'FOV_Old': [[1, Size[0]], [1, Size[1]], [1, Size[2]]]}

    if Image_3D_4D.ndim == 4:
        Mask = Image_3D_4D[:, :, :, 0] > 0
    else:
        Mask = Image_3D_4D > 0

    FOV = np.zeros((3, 2), dtype=int)
    for dim in range(3):
        temp = np.any(np.any(Mask, axis=(dim+1)%3), axis=(dim+2)%3)
        FOV[dim, :] = [np.argmax(temp) + 1, Size[dim] - np.argmax(temp[::-1])]

    for dim in range(3):
        if not np.isfinite(Extend[dim]):
            FOV[dim, :] = [1, Size[dim]]
        else:
            FOV[dim, :] = [max(1, FOV[dim, 0] - Extend[dim]), min(Size[dim], FOV[dim, 1] + Extend[dim])]

    if Extend[2] < 0:
        Mask[:,:,set(range(1, Size[2]+1)) - set(range(FOV[2, 0], FOV[2, 1]+1))] = False
        for dim in range(3):
            temp = np.any(np.any(Mask, axis=(dim+1)%3), axis=(dim+2)%3)
            FOV[dim, :] = [np.argmax(temp) + 1, Size[dim] - np.argmax(temp[::-1])]
            if not np.isfinite(Extend[dim]):
                FOV[dim, :] = [1, Size[dim]]
            else:
                FOV[dim, :] = [max(1, FOV[dim, 0] - Extend[dim]), min(Size[dim], FOV[dim, 1] + Extend[dim])]

    Crop_Parameter['FOV'] = FOV

    Center = Voxel_Size * (np.mean(FOV, axis=1) - (np.array(Size[:3]) + 1) / 2)

    if Image_3D_4D.ndim == 4:
        Image_3D_4D_Truncated = Image_3D_4D[FOV[0, 0]-1:FOV[0, 1], FOV[1, 0]-1:FOV[1, 1], FOV[2, 0]-1:FOV[2, 1], :]
    else:
        Image_3D_4D_Truncated = Image_3D_4D[FOV[0, 0]-1:FOV[0, 1], FOV[1, 0]-1:FOV[1, 1], FOV[2, 0]-1:FOV[2, 1]]

    return Image_3D_4D_Truncated, Center, Crop_Parameter


def fApply_Cropped_FOV(Image_3D_4D, Crop_Parameter):
    """
    Apply predefined Crop_Parameter to crop a 3D or 4D image.

    Parameters:
    Image_3D_4D : ndarray
        The input 3D or 4D image.
    Crop_Parameter : dict
        Dictionary containing 'FOV_Old' and 'FOV' as keys.

    Returns:
    Image_3D_4D : ndarray
        The cropped image.
    """
    Size = Image_3D_4D.shape
    FOV_Old = Crop_Parameter['FOV_Old']
    FOV = Crop_Parameter['FOV']

    if not (FOV_Old[0][1] - FOV_Old[0][0] == Size[0] - 1 and
            FOV_Old[1][1] - FOV_Old[1][0] == Size[1] - 1 and
            FOV_Old[2][1] - FOV_Old[2][0] == Size[2] - 1):
        print('Error in fApply_Cropped_FOV: the old FOV does not match to the data')
        return

    # Adjusting indices to be 0-based instead of 1-based
    FOV = np.array(FOV) - 1

    if Image_3D_4D.ndim == 3:
        Image_3D_4D = Image_3D_4D[FOV[0, 0]:FOV[0, 1] + 1,
                                  FOV[1, 0]:FOV[1, 1] + 1,
                                  FOV[2, 0]:FOV[2, 1] + 1]
    else:
        Image_3D_4D = Image_3D_4D[FOV[0, 0]:FOV[0, 1] + 1,
                                  FOV[1, 0]:FOV[1, 1] + 1,
                                  FOV[2, 0]:FOV[2, 1] + 1,
                                  :]
    return Image_3D_4D


def fInverse_Crop_EPI_Image_3D_4D(EPI_Image_3D_4D, Crop_Parameter):
    """
    An inverse function for cropping a 3D or 4D EPI image.

    Parameters:
    EPI_Image_3D_4D : ndarray
        The input 3D or 4D EPI image.
    Crop_Parameter : dict
        Dictionary containing 'FOV_Old' and 'FOV' as keys.

    Returns:
    EPI_Image_3D_4D : ndarray
        The inverse cropped image.
    """
    if 'FOV_Old' not in Crop_Parameter or 'FOV' not in Crop_Parameter:
        print('Error in Crop_Parameter: invalid settings of Crop_Parameter')
        return

    FOV_Old = np.array(Crop_Parameter['FOV_Old']) - 1  # Adjust indices for 0-based Python indexing
    FOV = np.array(Crop_Parameter['FOV']) - 1

    Data_shape = list(FOV_Old[:, 1] - FOV_Old[:, 0] + 1)
    if EPI_Image_3D_4D.ndim == 4:
        Data_shape.append(EPI_Image_3D_4D.shape[3])

    Data = np.zeros(Data_shape, dtype=EPI_Image_3D_4D.dtype)
    Data[FOV[0, 0]:FOV[0, 1] + 1,
         FOV[1, 0]:FOV[1, 1] + 1,
         FOV[2, 0]:FOV[2, 1] + 1,
         ...] = EPI_Image_3D_4D

    EPI_Image_3D_4D = Data
    return EPI_Image_3D_4D


def fMass_Center(Image_2D_3D):
    """
    To find the center of the largest connected region in a 2D or 3D image.

    Parameters:
    Image_2D_3D : ndarray
        The input 2D or 3D image.

    Returns:
    Center : ndarray or None
        The center of the largest connected region.
        If an error occurs, None is returned.
    """
    # Binarize the image
    Image_2D_3D = Image_2D_3D > 0

    if Image_2D_3D.ndim < 2:
        print('Error in fMass_Center: wrong size of Image_2D_3D')
        return None

    # Find the largest connected component
    labeled_array, num_features = label(Image_2D_3D)
    component_sizes = np.bincount(labeled_array.ravel())
    largest_component = component_sizes.argmax()
    Image_2D_3D = labeled_array == largest_component

    # Calculate the center of mass
    Center = np.array([np.sum(coords * Image_2D_3D) / np.sum(Image_2D_3D) for coords in np.indices(Image_2D_3D.shape)])

    return Center
