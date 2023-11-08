# Yuncong Ma, 11/6/2023
# Visualization module of pNet

#########################################
# Packages
import numpy as np
import os
import re
import time
import gc

import pandas as pd
import matplotlib
import surfplot
from brainspace.mesh.mesh_creation import build_polydata
from PIL import Image

# other functions of pNet
import pNet
from Data_Input import *
from FN_Computation import setup_pFN_folder
from Cropping import fApply_Cropped_FOV, fInverse_Crop_EPI_Image_3D_4D, fTruncate_Image_3D_4D, fMass_Center
import scipy
from collections import defaultdict

# reduce the memory leakage issue in macOS
matplotlib.use('TkAgg')  # Use the Tkinter backend

# =============== basic functions =============== #


def prepare_BSPolyData(vertices: np.ndarray, faces: np.ndarray):
    """
    prepare a BSPolyData class

    :param vertices: 2D matrix [N, 3], N is the number of vertices
    :param faces: 2D matrix, [N, 3], N is the number of triangle faces, vertex index starts from 0
    :return: bspolydata, a data class for using brain space surface plot

    Yuncong Ma, 11/1/2023
    """

    bspolydata = build_polydata(points=vertices, cells=faces)

    return bspolydata


def color_theme(theme: str,
                parameter: tuple,
                darker=False):
    """
    Get predefined color functions based color theme and additional parameters.
    This function matches to Yuncong Ma's MATLAB function fColor_Theme

    :param theme: 'Seed_Map', 'Seed_Map_3_Positive'
    :param parameter: parameters used for that specific theme. It is a tuple or scalar containing of 1-4 real numbers
    :param darker: False, or a real number in (0, 1)
    :return: color_function, a 2D ndarray matrix [N, 4], with first column as values, and the rest three as RGB colors

    Examples:
    color_function = color_theme('Seed_Map',(min_CC,Threshold,max_CC))
    color_function = color_theme('Seed_Map_2',(min_CC,max_CC))
    color_function = color_theme('Seed_Map_3',(min_CC,max_CC))
    color_function = color_theme('Seed_Map_Positive',(min_CC,Threshold,max_CC))
    color_function = color_theme('Seed_Map_Positive_Only',(min_CC,Threshold,max_CC))
    color_function = color_theme('Seed_Map_3_Positive',(min_CC,max_CC))
    color_function = color_theme('FC',max_CC)
    color_function = color_theme('Atlas',Max_Index)
    color_function = color_theme('Gray_Scale',(Min_Value,Min_Gray,Max_Value,Max_Gray))
    color_function = color_theme('Rainbow',(Min_Value,Max_Value))

    Yuncong Ma, 10/27/2023
    """

    if theme == 'Seed_Map':
        min_CC = parameter[0]
        Threshold = parameter[1]
        max_CC = parameter[2]
        color_function = np.array(((-max_CC,0,1,1),
                     (-Threshold,0,(Threshold-min_CC)/(max_CC-min_CC),1),
                     (-Threshold,0,0,0),
                     (0,0,0,0),
                     (Threshold,0,0,0),
                     (Threshold,1,(Threshold-min_CC)/(max_CC-min_CC),0),
                     (max_CC,1,1,0)), dtype=np.float32)

    elif theme == 'Seed_Map_2':
        min_CC = parameter[0]
        max_CC = parameter[1]
        color_function = np.array((
            (-max_CC,0,1,1),
            (-min_CC,0,0,1),
            (-min_CC,0,0,0),
            (0,0,0,0),
            (min_CC,0,0,0),
            (min_CC,1,0,0),
            (max_CC,1,1,0)), dtype=np.float32)

    elif theme == 'Seed_Map_3':
        min_CC = parameter[0]
        max_CC = parameter[1]
        color_function = np.array((
            (-max_CC,0,1,1),
            (-min_CC*0.7 - max_CC*0.3,0,0,1),
            (-min_CC,0,0,0.5),
            (-min_CC,0,0,0),
            (0,0,0,0),
            (min_CC,0,0,0),
            (min_CC,0.5,0,0),
            (min_CC*0.7 + max_CC*0.3,1,0,0),
            (max_CC,1,1,0)), dtype=np.float32)

    elif theme == 'Seed_Map_Positive':
        min_CC = parameter[0]
        Threshold = parameter[1]
        max_CC = parameter[2]
        color_function = np.array((
            (0,0,0,0),
            (Threshold,0,0,0),
            (Threshold,1,(Threshold-min_CC)/(max_CC-min_CC),0),
            (max_CC,1,1,0)), dtype=np.float32)

    elif theme == 'Seed_Map_Positive_Only':
        min_CC = parameter[0]
        Threshold = parameter[1]
        max_CC = parameter[2]
        color_function = np.array((
            (Threshold,0,0,0),
            (Threshold,1,(Threshold-min_CC)/(max_CC-min_CC),0),
            (max_CC,1,1,0)), dtype=np.float32)

    elif theme == 'Seed_Map_3_Positive':
        min_CC = parameter[0]
        max_CC = parameter[1]
        color_function = np.array((
            (min_CC,0,0,0),
            (min_CC,0.5,0,0),
            (min_CC*0.7 + max_CC*0.3,1,0,0),
            (max_CC,1,1,0)), dtype=np.float32)

    elif theme == 'FC':
        max_CC = parameter[0]
        color_function = np.array((
            (-max_CC,0,1,1),
            (0,0,0,1),
            (0,1,0,0),
            (max_CC,1,1,0)), dtype=np.float32)

    elif theme == 'Atlas':
        Color_Table = load_matlab_single_array(os.path.join(pNet.dir_python, 'Color_Table.mat'))
        Max_Index = parameter[0]
        color_function = np.zeros((2*Max_Index+1, 4), dtype=np.float32)
        color_function[0, 0] = 0.5
        for i in range(1, 2, 2*Max_Index):
            color_function[i, :] = (i/2-0.5, Color_Table[i, :])
            color_function[i+1, :] = (i/2+0.5, Color_Table[i, :])

    elif theme == 'Gray_Scale':
        Min_Value = parameter[0]
        Min_Gray = parameter[1]
        Max_Value = parameter[2]
        Max_Gray = parameter[3]
        color_function = np.array((
            (Min_Value, Min_Gray, Min_Gray, Min_Gray),
            (Max_Value, Max_Gray, Max_Gray, Max_Gray)), dtype=np.float32)

    elif theme == 'Rainbow':
        Min_Value = parameter[0]
        Max_Value = parameter[1]
        color_function = np.array((
                      (Min_Value, 75.0/255.0, 0, 130.0/255.0),
                      (0.2*(Max_Value-Min_Value)+Min_Value,0,0,1),
                      (0.4*(Max_Value-Min_Value)+Min_Value,0,1,0),
                      (0.6*(Max_Value-Min_Value)+Min_Value,1,1,0),
                      (0.8*(Max_Value-Min_Value)+Min_Value,1,127.0/255.0,0),
                      (Max_Value,1,0,0)), dtype=np.float32)

    else:
        raise ValueError('Unknown color theme: ' + theme)

    if isinstance(darker, float) or isinstance(darker, np.ndarray):
        color_function[:, 1:4] = darker * color_function[:, 1:4]

    return color_function


def prepare_color_map(map_name=None or str,
                      color_function=None or np.ndarray,
                      N=256,
                      alpha=1,
                      black_transparent=True):
    """
    Prepare a color map using color function

    :param map_name: a color map name for matplotlib.colors.ListedColormap
    :param color_function: np.ndarray [N, 4], N is the number of predefined color changing points
    :param N: digitized colors
    :param alpha: transparency
    :param black_transparent: True or False, True is to set black color to transparent
    :return: cmap: a structure created by matplotlib.colors.ListedColormap

    Yuncong Ma, 10/30/2023
    """

    # use built-in color maps in matplotlib
    if isinstance(map_name, str) and len(map_name) > 0:
        cmap = matplotlib.pyplot.get_cmap(map_name)
        cmapV = cmap(np.arange(cmap.N))
        cmapV[:, -1] = alpha
        cmap = matplotlib.colors.ListedColormap(cmapV)
        return cmap

    # number of digitized colors
    N_cf = color_function.shape[0]

    # resample color_map to a color function with fixed step size and length
    color_map = np.zeros((N, 4), dtype=np.float32)
    for i in range(N):
        value = float(i) / float(N-1) * (color_function[-1, 0] - color_function[0, 0]) + color_function[0, 0]
        if value <= color_function[0, 0]:
            color_map[i, 0:3] = color_function[0, 1:4]
            continue
        elif value >= color_function[-1, 0]:
            color_map[i, 0:3] = color_function[-1, 1:4]
            continue

        for j in range(N_cf-1):
            if color_function[j, 0] <= value <= color_function[j+1, 0]:
                if color_function[j, 0] == value:
                    color_map[i, 0:3] = color_function[j, 1:4]
                    break
                elif value == color_function[j+1, 0]:
                    color_map[i, 0:3] = color_function[j+1, 1:4]
                    break
                weight = (value - color_function[j, 0]) / (color_function[j+1, 0] - color_function[j, 0])
                color_map[i, 0:3] = (1-weight) * color_function[j, 1:4] + weight * color_function[j+1, 1:4]
                break

    # add alpha
    color_map[:, 3] = alpha

    # Set black color to transparent
    if black_transparent is True:
        color_map[:, 3] = color_map[:, 3] * (np.sum(color_map[:, 0:3], axis=1) > 0)

    # convert to matplotlib
    cmap = matplotlib.colors.ListedColormap(color_map)

    return cmap


def colorize(value_map: np.ndarray, color_function: np.ndarray):
    """
    colorize a 2D value map into a 2D image [X, Y, 3]

    :param value_map:
    :param color_function:
    :return: image_rgb: np.ndarray, np.float32, 0-1

    Yuncong Ma, 11/3/2023
    """

    # Ensure Color_Function is valid
    if color_function.shape[1] != 4 or color_function.shape[0] < 2 or not np.all(color_function[:, 0] == np.sort(color_function[:, 0])):
        raise ValueError('check color_function')

    # Flatten the matrix and initialize the RGB image
    flat_matrix = value_map.flatten()
    image_rgb = np.zeros((flat_matrix.size, 3), dtype=np.float32)

    # Handle values below the first threshold
    mask = flat_matrix < color_function[0, 0]
    image_rgb[mask] = color_function[0, 1:4]

    # Handle values above the last threshold
    mask = flat_matrix >= color_function[-1, 0]
    image_rgb[mask] = color_function[-1, 1:4]

    # Handle values within the color thresholds
    for p in range(color_function.shape[0]-1):
        mask = (color_function[p, 0] <= flat_matrix) & (flat_matrix < color_function[p+1, 0])
        if color_function[p, 0] == color_function[p+1, 0]:
            weight = np.ones(flat_matrix.shape)
        else:
            weight = (flat_matrix - color_function[p, 0]) / (color_function[p+1, 0] - color_function[p, 0])
        for c in range(3):  # Iterate over RGB channels
            image_rgb[mask, c] = (weight[mask] * color_function[p+1, c + 1] +
                                  (1 - weight[mask]) * color_function[p, c + 1])

    # Reshape to original dimensions with RGB channels, unit8
    image_rgb = image_rgb.reshape((*value_map.shape, 3)).astype(np.float32)
    return image_rgb


def assemble_image(file_list_image: tuple, file_output_assembled=None or str, organization=(0, 10), interval=(50, 5), background=(0, 0, 0)):
    """

    :param file_list_image: a tuple of image directories or image matrices, images must have the same size
    :param file_output_assembled: output file directory, can be None to get image matrix as output
    :param organization: number of rows and columns, default is (0, 10) means to set 10 columns with automatic row number
    :param interval: (50, 5) in default, meaning the interval is 50 by 5 pixels
    :param background: (0, 0, 0) in default, meaning the background color is black
    :return: image_assembled (M, N, 3) matrix, if file_output_assembled is None

    Yuncong Ma, 11/3/2023
    """

    N_image = len(file_list_image)
    if isinstance(organization, tuple):
        organization = np.array(organization)
    if organization[0] == 0 and organization[1] > 0:
        organization[0] = np.ceil(float(N_image)/float(organization[1]))
    elif organization[1] == 0 and organization[0] > 0:
        organization[1] = np.ceil(float(N_image)/float(organization[0]))
    elif organization[0] == 0 and organization[1] == 0:
        organization[0] = np.ceil(np.sqrt(float(N_image)))
        organization[1] = np.ceil(float(N_image)/float(organization[0]))

    count = 0
    ps = np.array((0, 0))
    for x in range(organization[0]):
        for y in range(organization[1]):
            if count < N_image:
                if isinstance(file_list_image[count], str):
                    image_sub = np.array(Image.open(file_list_image[count]))
                elif isinstance(file_list_image[count], tuple):
                    image_sub = np.array(file_list_image[count])
                elif isinstance(file_list_image[count], np.ndarray):
                    image_sub = file_list_image[count]
                else:
                    raise ValueError('file_list_image must be either a tuple of image directories or image matrices')

            if x == 0 and y == 0:
                image_assembled = np.zeros((image_sub.shape[0] * organization[0] + (organization[0]-1)*interval[0], image_sub.shape[1] * organization[1] + (organization[1]-1)*interval[1], 3), dtype=np.uint8)
                image_assembled[0:image_sub.shape[0], 0:image_sub.shape[1], :] = image_sub
            else:
                if count < N_image:
                    image_assembled[ps[0]:ps[0]+image_sub.shape[0], ps[1]:ps[1]+image_sub.shape[1], :] = image_sub

            if x < organization[0] - 1:
                image_assembled[ps[0]+image_sub.shape[0]:ps[0]+image_sub.shape[0]+interval[0], ps[1]:ps[1]+image_sub.shape[1], :] = \
                    np.reshape(background, (1, 1, 3))
            if y < organization[1] - 1:
                image_assembled[ps[0]:ps[0]+image_sub.shape[0], ps[1]+image_sub.shape[1]:ps[1]+image_sub.shape[1]+interval[1], :] = \
                    np.reshape(background, (1, 1, 3))
            if x < organization[0] - 1 and y < organization[1] - 1:
                image_assembled[ps[0]:ps[0]+image_sub.shape[0]+interval[0], ps[1]+image_sub.shape[1]:ps[1]+image_sub.shape[1]+interval[1], :] = \
                    np.reshape(background, (1, 1, 3))

            ps[1] = ps[1]+interval[1] + image_sub.shape[1]
            count += 1
        ps[0] = ps[0] + interval[0] + image_sub.shape[0]
        ps[1] = 0

    if file_output_assembled is None:
        return image_assembled
    else:
        image = Image.fromarray(image_assembled, 'RGB')
        image.save(file_output_assembled)
        image.close()


# =============== Surface data type =============== #


def plot_brain_surface(brain_map: np.ndarray,
                       mesh: dict,
                       mask: np.ndarray,
                       color_function: np.ndarray,
                       file_output=None or str,
                       orientation='medial',
                       view_angle=1.5,
                       mask_color=(0.2, 0.2, 0.2),
                       brain_color=(0.5, 0.5, 0.5),
                       background=(0, 0, 0),
                       figure_size=(500, 400),
                       dpi=25):

    # Prepare BSPolyData for using its plot
    polyData = prepare_BSPolyData(mesh['vertices'], mesh['faces'] - 1)
    p = surfplot.Plot(surf_lh=polyData, zoom=view_angle, views=orientation, background=background, brightness=1, size=figure_size)

    # brain surface and mask layer
    map_mask = (mask == 0).astype(np.float32)
    color_function_brain = np.array(((0, brain_color[0], brain_color[1], brain_color[2]), (1, mask_color[0], mask_color[1], mask_color[2])), dtype=np.float32)
    p.add_layer(map_mask, cmap=prepare_color_map(color_function=color_function_brain),
                color_range=(0, 1), cbar=None, zero_transparent=False)

    # map layer
    # convert map to the mesh surface space based on the brain mask
    Nv = sum(mask > 0)
    map_2 = np.zeros(mask.shape, dtype=np.float32)
    ps = np.where(mask > 0)[0].astype(int)
    for i in range(int(Nv)):
        map_2[ps[i]] = brain_map[i]
    color_range = (color_function[0, 0], color_function[-1, 0])

    p.add_layer(map_2, cmap=prepare_color_map(color_function=color_function), color_range=color_range, cbar=None, zero_transparent=True)

    # save or return
    if file_output is not None:
        # build the figure
        fig = p.build()
        fig.savefig(file_output, dpi=dpi, bbox_inches="tight", facecolor=background)
    else:
        p = p.render()
        p._check_offscreen()
        image_rgb = p.to_numpy(transparent_bg=True, scale=(2, 2))
        del p
        gc.collect()
        return image_rgb


def merge_mesh_LR(mesh_LR: dict, offset=np.array((90, 0, 0))):
    """
    Merge two meshes into one

    :param mesh_LR: a dict containing 'L' and 'R', in which vertices and faces are stored
    :param offset: np.ndarray, 3D vector for the distance between the centers of left and right hemispheres
    :return: mesh, a dict of vertices and faces

    Yuncong Ma, 11/6/2023
    """

    mesh = {'vertices': np.concatenate((mesh_LR['L']['vertices'], mesh_LR['R']['vertices'] + offset), axis=0),
            'faces': np.concatenate((mesh_LR['L']['faces'], mesh_LR['R']['faces']+mesh_LR['L']['vertices'].shape[0]), axis=0)}
    return mesh


def merge_mask_LR(mask_LR: dict):
    """
    Merge masks of two hemispheres for surface type

    :param mask_LR: a dict containing 'L' and 'R', each contains a 1D vector, np.ndarray
    :return: a 1D vector, np.ndarray

    Yuncong Ma, 11/6/2023
    """
    mask = np.concatenate((mask_LR['L'], mask_LR['R']), axis=0)
    return mask


def plot_FN_brain_surface_5view(brain_map: np.ndarray,
                                brain_template,
                                file_output=None or str,
                                threshold=99,
                                color_function=None,
                                background=(0, 0, 0),
                                figure_organization=(0.6, 1.2, 1, 0.6),
                                view_angle=(1.35, 1.4),
                                hemisphere_offset=90,
                                figure_title=None,
                                title_font_dic=dict(fontsize=20, fontweight='bold'),
                                figure_size=(10, 50),
                                dpi=50):

    # check NaN in brain_map
    brain_map[np.isnan(brain_map)] = 0

    # settings for subplot
    fig, axs = matplotlib.pyplot.subplots(nrows=8, ncols=1, figsize=figure_size)

    # set color function
    if color_function is None:
        threshold_value = np.percentile(np.abs(brain_map), threshold)
        color_range = np.array((threshold_value/2, threshold_value))
        color_function = color_theme('Seed_Map_3_Positive', color_range)
    else:
        color_range = (color_function[0, 0], color_function[-1, 0])

    # sub figure organization
    H = 4*figure_organization[2]+figure_organization[1]+figure_organization[0] + figure_organization[3]
    H_T = figure_organization[0]/H  # height of title
    H_D = figure_organization[1]/H  # height of dorsal view
    H_S = figure_organization[2]/H  # height of saggital view
    H_C = figure_organization[3]/H  # height of color bar

    # number of used vertices in left hemisphere
    Nv_L = int(sum(brain_template['Brain_Mask']['L'] > 0))

    # title
    axs[0].set_position((0, 4*H_S+H_D+H_C, 1, H_T))
    axs[0].figure_size = (int(dpi*figure_size[0]), int(dpi*H_T*figure_size[1]))
    axs[0].facecolor = background
    axs[0].axis('off')

    # dorsal view
    image_rgb = plot_brain_surface(brain_map, mesh=merge_mesh_LR(brain_template['Shape'], offset=np.array((hemisphere_offset, 0, 0))), mask=merge_mask_LR(brain_template['Brain_Mask']),
                           color_function=color_function,
                           orientation='dorsal', view_angle=view_angle[0], file_output=None, background=(1, 1, 1),
                           figure_size=(int(dpi*figure_size[0]), int(dpi*H_D*figure_size[1])), dpi=dpi)

    axs[1].figure_size = (int(dpi*figure_size[0]), int(dpi*H_D*figure_size[1]))
    axs[1].set_position((0, 4*H_S+H_C, 1, H_D))
    axs[1].imshow(image_rgb)
    axs[1].axis('off')
    axs[1].set_title(label=figure_title, loc='center', pad=140, fontsize=150, fontweight='bold', fontname='Arial', color=(1, 1, 1))

    # saggital views
    # 1st
    image_rgb = plot_brain_surface(brain_map[0:Nv_L], mesh=brain_template['Shape']['L'], mask=brain_template['Brain_Mask']['L'],
                           color_function=color_function,
                           orientation='lateral', view_angle=view_angle[1], file_output=None,
                           figure_size=(int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1])), dpi=dpi)

    axs[2].figure_size = (int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1]))
    axs[2].set_position((0, 3*H_S+H_C, 1, H_S))
    axs[2].imshow(image_rgb)
    axs[2].axis('off')

    # 2nd
    image_rgb = plot_brain_surface(brain_map[0:Nv_L], mesh=brain_template['Shape']['L'], mask=brain_template['Brain_Mask']['L'],
                           color_function=color_function,
                           orientation='medial', view_angle=view_angle[1], file_output=None, dpi=dpi)

    axs[3].figure_size = (int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1]))
    axs[3].set_position((0, 2*H_S+H_C, 1, H_S))
    axs[3].imshow(image_rgb)
    axs[3].axis('off')
    # 3rd
    image_rgb = plot_brain_surface(brain_map[Nv_L:], mesh=brain_template['Shape']['R'], mask=brain_template['Brain_Mask']['R'],
                           color_function=color_function,
                           orientation='medial', view_angle=view_angle[1], file_output=None, dpi=dpi)

    axs[4].figure_size = (int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1]))
    axs[4].set_position((0, 1*H_S+H_C, 1, H_S))
    axs[4].imshow(image_rgb)
    axs[4].axis('off')
    # 4th
    image_rgb = plot_brain_surface(brain_map[Nv_L:], mesh=brain_template['Shape']['R'], mask=brain_template['Brain_Mask']['R'],
                           color_function=color_function,
                           orientation='lateral', view_angle=view_angle[1], file_output=None, dpi=dpi)
    axs[5].figure_size = (int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1]))
    axs[5].set_position((0, H_C, 1, H_S))
    axs[5].imshow(image_rgb)
    axs[5].axis('off')

    # color bar
    axs[6].figure_size = (int(dpi*figure_size[0]), int(dpi*H_C*figure_size[1]))
    colorbar_width = 0.8
    colorbar_height = 0.2
    colorbar_pad = 0.7
    colorbar_step = 100
    colorbar_ratio = 10
    colorbar_scale = 100
    axs[6].set_position(((1-colorbar_width)/2, H_C*colorbar_pad, colorbar_width, H_C*colorbar_height))
    X = np.tile(np.arange(color_range[0], color_range[1], (color_range[1] - color_range[0])/colorbar_step), (colorbar_ratio, 1))
    axs[6].imshow(X, cmap=prepare_color_map(color_function=color_function))
    axs[6].axis('off')
    cb_tick_pad = 25
    cb_value_round = True
    cb_name = 'Loading (%)'
    cb_name_pad = 45
    if cb_value_round is True:
        color_range = np.round(color_range * colorbar_scale)
        cb_tick = (str(int(color_range[0])), str(int(color_range[1])))
    else:
        color_range = color_range * colorbar_scale
        cb_tick = (str(color_range[0]), str(color_range[1]))

    # add label
    axs[6].text(0, cb_tick_pad, cb_tick[0],
                ha='left', va='bottom',
                fontdict=dict(fontsize=80, fontweight='bold', color=(1, 1, 1), fontname='Arial'))
    axs[6].text(colorbar_step, cb_tick_pad, cb_tick[1],
                ha='right', va='bottom',
                fontdict=dict(fontsize=80, fontweight='bold', color=(1, 1, 1), fontname='Arial'))
    axs[6].text(colorbar_step/2, cb_name_pad, cb_name,
                ha='center', va='bottom',
                fontdict=dict(fontsize=80, fontweight='bold', color=(1, 1, 1), fontname='Arial'))

    # add a block to maintain the preconfigured figure size
    axs[7].figure_size = (int(dpi*figure_size[0]), int(dpi*H_C*figure_size[1]))
    axs[7].set_position(((1-colorbar_width)/2, 0, colorbar_width, H_C))
    axs[7].axis('off')

    # save fig
    if file_output is None:
        return fig, axs
    else:
        fig.savefig(file_output, dpi=dpi, bbox_inches="tight", facecolor=background)
        # Clear the axes
        for i in range(8):
            axs[i].cla()
        matplotlib.pyplot.close(fig)
        matplotlib.pyplot.cla()
        matplotlib.pyplot.clf()
        del fig, axs
        gc.collect()


# =============== Volume data type =============== #


def plot_voxel_map_3view(Anatomy: np.ndarray, Voxel_Map: np.ndarray, center: np.ndarray, color_function: np.ndarray,
                         rotation=np.array((0, 0, 0)), organization=np.array((0, 0, 0)),
                         background=np.array((0, 0, 0))):

    # Check if Anatomy and Voxel_Map have the same size
    if Anatomy.shape != Voxel_Map.shape:
        raise ValueError('Anatomy and Voxel_Map should have the same size')

    Dimension = Anatomy.shape

    # Create the 2D anatomical views by rotation and slicing
    Anatomy2D = [
        np.rot90(Anatomy[:, :, center[2]], rotation[0]),
        np.rot90(np.squeeze(Anatomy[:, center[1], :]), rotation[1]),
        np.rot90(np.squeeze(Anatomy[center[0], :, :]), rotation[2])
    ]

    # Normalize and replicate the anatomical views for RGB channels
    for i in range(3):
        range_values = np.percentile(Anatomy2D[i][Anatomy2D[i] > 0], [1, 99])
        Anatomy2D[i] = (Anatomy2D[i] - range_values[0]) / np.diff(range_values) * 0.8
        Anatomy2D[i] = np.repeat(Anatomy2D[i][:, :, np.newaxis], 3, axis=2)

    # Create the 2D voxel map views
    Voxel_Map2D = [
        np.rot90(Voxel_Map[:, :, center[2]], rotation[0]),
        np.rot90(np.squeeze(Voxel_Map[:, center[1], :]), rotation[1]),
        np.rot90(np.squeeze(Voxel_Map[center[0], :, :]), rotation[2])
    ]

    # Create masks and colorize
    Mask = [None] * 3
    for i in range(3):
        Voxel_Map2D[i] = colorize(Voxel_Map2D[i], color_function)
        Mask[i] = np.repeat(np.sum(Voxel_Map2D[i], axis=2) == 0, 3).reshape(Voxel_Map2D[i].shape)

    # Assemble the figure based on the organization
    if organization in ([[1, 3], [2, 0]], [[2, 3], [1, 0]]):
        image_rgb = np.zeros((*Dimension[:2], 3)) + background
        for i, val in np.ndenumerate(np.array(organization)):
            if val:
                image_rgb[i[0]*Dimension[0]:(i[0]+1)*Dimension[0], i[1]*Dimension[1]:(i[1]+1)*Dimension[1], :] = Anatomy2D[val-1] * Mask[val-1] + Voxel_Map2D[val-1] * (1 - Mask[val-1])
    elif np.array(organization).shape == (3,):
        image_rgb = np.zeros((Dimension[0]*3, Dimension[1], 3), dtype=np.float32) + background
        for i, val in enumerate(organization.flatten()):
            image_rgb[i*Dimension[0]:(i+1)*Dimension[0], :, :] = Anatomy2D[val-1] * Mask[val-1] + Voxel_Map2D[val-1] * (1 - Mask[val-1])
    else:
        raise ValueError('unsupported Organization settings')

    image_rgb[image_rgb < 0] = 0
    image_rgb[image_rgb > 1] = 1

    return image_rgb


def large_3view_center(weight_map: np.ndarray):
    """
    Get the view center for a 3D map, maximizing the content in each view

    :param weight_map: a 3D matrix, either 0-1 or weighted
    :return: Center ndarray: [3, 1]

    Yuncong Ma, 11/6/2023
    """

    Center = np.zeros(3, dtype=np.int32)

    X_size = np.sum(np.sum(weight_map, axis=2), axis=1)
    Center[0] = X_size.argmax()
    Y_size = np.sum(np.sum(weight_map, axis=2), axis=0)
    Center[1] = Y_size.argmax()
    Z_size = np.sum(np.sum(weight_map, axis=1), axis=0)
    Center[2] = Z_size.argmax()

    return Center


def plot_FN_brain_volume_3view(brain_map: np.ndarray,
                               brain_template,
                               file_output=None or str,
                               threshold=99,
                               color_function=None,
                               view_center='cluster_center',
                               figure_organization=(0.4, 3, 0.6),
                               background=(0, 0, 0),
                               figure_title=None,
                               title_font_dic=dict(fontsize=20, fontweight='bold'),
                               interpolation='nearest',
                               figure_size=(10, 40),
                               dpi=250
                               ):

    # check NaN in brain_map
    brain_map[np.isnan(brain_map)] = 0

    # set color function
    if color_function is None:
        threshold_value = np.percentile(np.abs(reshape_FN(brain_map, dataType='Volume', Brain_Mask=brain_template['Brain_Mask'])), threshold)
        color_range = np.array((threshold_value/2, threshold_value))
        color_function = color_theme('Seed_Map_3_Positive', color_range)
    else:
        color_range = (color_function[0, 0], color_function[-1, 0])

    if not brain_map.shape == brain_template['Brain_Mask'].shape:
        raise ValueError('the brain_map must have the same image size as the Brain_Mask in brain_template')

    # upsampling
    upsampling = np.round(brain_template['Overlay_Image'].shape[0] / brain_map.shape[0])
    if np.sum(np.abs(np.array(brain_map.shape) * upsampling - np.array(brain_template['Overlay_Image'].shape))) > 0:
        raise ValueError('the Overlay_Image does NOT have an integer upsampling scale to the brain map')

    if interpolation == 'nearest':
        Map = scipy.ndimage.zoom(brain_map, upsampling, order=0)  # 'nearest' interpolation is order=0
    elif interpolation == 'spline-3':
        Map = scipy.ndimage.zoom(brain_map, upsampling, order=3)
    else:
        raise ValueError('Unknown ')

    Brain_Mask = scipy.ndimage.zoom(brain_template['Brain_Mask'], upsampling, order=0)
    Brain_Mask_2, _, Crop_Parameter = fTruncate_Image_3D_4D(Brain_Mask, Voxel_Size=np.array((1, 1, 1)), Extend=np.array((2, 2, 2)))

    Overlay_Image = brain_template['Overlay_Image']
    Overlay_Image_2 = fApply_Cropped_FOV(Overlay_Image, Crop_Parameter)
    Map_2 = fApply_Cropped_FOV(Map, Crop_Parameter)

    Max_Dim = np.max(Map_2.shape)
    Crop_Parameter['FOV_Old'] = [[1, Max_Dim]] * 3
    Crop_Parameter['FOV'] = np.array([[1, s] for s in Map_2.shape]) + np.tile(np.round((np.array([Max_Dim] * 3) - np.array(Map_2.shape)) / 2), (2, 1)).T
    Crop_Parameter['FOV'] = np.array(Crop_Parameter['FOV'], dtype=np.int32)

    Map_2 = fInverse_Crop_EPI_Image_3D_4D(Map_2, Crop_Parameter)
    Brain_Mask_2 = fInverse_Crop_EPI_Image_3D_4D(Brain_Mask_2, Crop_Parameter)
    Overlay_Image_2 = fInverse_Crop_EPI_Image_3D_4D(Overlay_Image_2, Crop_Parameter)

    if isinstance(view_center, np.ndarray):
        Center = view_center
    elif view_center == 'cluster_center':
        Center = large_3view_center(Map_2)

    # Get three images
    rotation = np.array((1, 1, 1))
    organization = np.array((2, 1, 0))
    image_rgb = plot_voxel_map_3view(Overlay_Image_2, Map_2, center=Center, rotation=rotation, organization=organization, color_function=color_function)

    # settings for subplot
    fig, axs = matplotlib.pyplot.subplots(nrows=4, ncols=1, figsize=figure_size)

    # sub figure organization
    H = np.sum(np.array(figure_organization))
    H_T = figure_organization[0]/H  # height of title
    H_V = figure_organization[1]/H  # height of 3 views
    H_C = figure_organization[2]/H  # height of color bar

    # title
    axs[0].set_position((0, H_V+H_C, 1, H_T))
    axs[0].figure_size = (int(dpi*figure_size[0]), int(dpi*H_T*figure_size[1]))
    axs[0].facecolor = background
    axs[0].axis('off')

    # 3 views
    axs[1].figure_size = (int(dpi*figure_size[0]), int(dpi*H_V*figure_size[1]))
    axs[1].set_position((0, H_C, 1, H_V))
    axs[1].imshow(image_rgb)
    axs[1].axis('off')
    axs[1].set_title(label=figure_title, loc='center', pad=30, fontsize=150, fontweight='bold', fontname='Arial', color=(1, 1, 1))

    # color bar
    axs[2].figure_size = (int(dpi*figure_size[0]), int(dpi*H_C*figure_size[1]))
    colorbar_width = 0.8
    colorbar_height = 0.2
    colorbar_pad = 0.7
    colorbar_step = 100
    colorbar_ratio = 10
    colorbar_scale = 100
    axs[2].set_position(((1-colorbar_width)/2, H_C*colorbar_pad, colorbar_width, H_C*colorbar_height))
    X = np.tile(np.arange(color_range[0], color_range[1], (color_range[1] - color_range[0])/colorbar_step), (colorbar_ratio, 1))
    axs[2].imshow(X, cmap=prepare_color_map(color_function=color_function))
    axs[2].axis('off')
    cb_tick_pad = 25
    cb_value_round = True
    cb_name = 'Loading (%)'
    cb_name_pad = 45
    if cb_value_round is True:
        color_range = np.round(color_range * colorbar_scale)
        cb_tick = (str(int(color_range[0])), str(int(color_range[1])))
    else:
        color_range = color_range * colorbar_scale
        cb_tick = (str(color_range[0]), str(color_range[1]))

    axs[2].text(0, cb_tick_pad, cb_tick[0],
                ha='left', va='bottom',
                fontdict=dict(fontsize=80, fontweight='bold', color=(1, 1, 1), fontname='Arial'))
    axs[2].text(colorbar_step, cb_tick_pad, cb_tick[1],
                ha='right', va='bottom',
                fontdict=dict(fontsize=80, fontweight='bold', color=(1, 1, 1), fontname='Arial'))
    axs[2].text(colorbar_step/2, cb_name_pad, cb_name,
                ha='center', va='bottom',
                fontdict=dict(fontsize=80, fontweight='bold', color=(1, 1, 1), fontname='Arial'))

    # add a block to maintain the preconfigured figure size
    axs[3].figure_size = (int(dpi*figure_size[0]), int(dpi*H_C*figure_size[1]))
    axs[3].set_position(((1-colorbar_width)/2, 0, colorbar_width, H_C))
    axs[3].axis('off')

    # save fig
    if file_output is None:
        return fig, axs
    else:
        fig.savefig(file_output, dpi=dpi, bbox_inches="tight", facecolor=background)
        # Clear the axes
        for i in range(4):
            axs[i].cla()
        fig.clf()
        matplotlib.pyplot.close(fig)

    return


# =========== Surface-volume data type =========== #


# =========== Module =========== #


def setup_Visualization(file_figure, ):
    return


def run_gFN_Visualization(dir_pnet_result: str):
    """
    Run preconfigured visualizations for gFNs

    :param dir_pnet_result: directory of the pnet result folder
    :return:

    Yuncong Ma, 11/1/2023
    """

    # get directories of sub-folders
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, _, _ = setup_result_folder(dir_pnet_result)

    # load settings for data input and FN computation
    if not os.path.isfile(os.path.join(dir_pnet_dataInput, 'Setting.json')):
        raise ValueError('Cannot find the setting json file in folder Data_Input')
    settingDataInput = load_json_setting(os.path.join(dir_pnet_dataInput, 'Setting.json'))

    setting = {'Data_Input': settingDataInput}

    # load basic settings
    dataType = setting['Data_Input']['Data_Type']
    dataFormat = setting['Data_Input']['Data_Format']

    gFN = load_matlab_single_array(os.path.join(dir_pnet_gFN, 'FN.mat'))
    brain_template = load_brain_template(os.path.join(dir_pnet_dataInput, 'Brain_Template.json.zip'))

    if dataType == 'Surface' and dataFormat == 'HCP Surface (*.cifti, *.mat)':
        K = gFN.shape[1]
        file_output = [os.path.join(dir_pnet_gFN, str(int(i+1))+'.jpg') for i in range(K)]
        for i in range(K):
            figure_title = 'FN '+str(int(i+1))
            plot_FN_brain_surface_5view(gFN[:, i], brain_template, color_function=None, file_output=file_output[i], figure_title=figure_title)

    elif dataType == 'Volume':
        K = gFN.shape[3]
        file_output = [os.path.join(dir_pnet_gFN, str(int(i+1))+'.jpg') for i in range(K)]
        for i in range(K):
            figure_title = 'FN '+str(int(i+1))
            plot_FN_brain_volume_3view(gFN[:, :, :, i], brain_template, color_function=None, file_output=file_output[i], figure_title=figure_title)

    elif dataType == 'Surface-Volume' and dataFormat == 'HCP Surface-Volume (*.cifti)':
        K = gFN.shape[1]
        file_output = [os.path.join(dir_pnet_gFN, str(int(i+1))+'.jpg') for i in range(K)]
        for i in range(K):
            figure_title = 'FN '+str(int(i+1))
            plot_FN_brain_surface_volume_5view(gFN[:, i], brain_template, color_function=None, file_output=file_output[i], figure_title=figure_title)

    # output an assembled image
    file_output_assembled = os.path.join(dir_pnet_gFN, 'All.jpg')
    assemble_image(file_output, file_output_assembled, interval=(50, 5), background=(0, 0, 0))

    return


def run_pFN_Visualization(dir_pnet_result: str):
    """
    Run preconfigured visualizations for pFNs

    :param dir_pnet_result: directory of the pnet result folder
    :return:

    Yuncong Ma, 11/6/2023
    """

    # get directories of sub-folders
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, _, _ = setup_result_folder(dir_pnet_result)

    # load settings for data input and FN computation
    if not os.path.isfile(os.path.join(dir_pnet_dataInput, 'Setting.json')):
        raise ValueError('Cannot find the setting json file in folder Data_Input')
    settingDataInput = load_json_setting(os.path.join(dir_pnet_dataInput, 'Setting.json'))

    setting = {'Data_Input': settingDataInput}

    # load basic settings
    dataType = setting['Data_Input']['Data_Type']
    dataFormat = setting['Data_Input']['Data_Format']

    brain_template = load_brain_template(os.path.join(dir_pnet_dataInput, 'Brain_Template.json'))

    # setup folders in Personalized_FN
    list_subject_folder = setup_pFN_folder(dir_pnet_result)
    N_Scan = len(list_subject_folder)
    for scan in range(1, N_Scan+1):
        # print(f'Start to visualize pFNs for {i}-th folder: {list_subject_folder[i-1]}', file=logFile_FNC, flush=True)
        dir_pnet_pFN_indv = os.path.join(dir_pnet_pFN, list_subject_folder[scan-1])
        pFN = load_matlab_single_array(os.path.join(dir_pnet_pFN_indv, 'FN.mat'))

        if dataType == 'Surface' and dataFormat == 'HCP Surface (*.cifti, *.mat)':
            K = pFN.shape[1]
            file_output = [os.path.join(dir_pnet_pFN_indv, str(int(i+1))+'.jpg') for i in range(K)]
            for i in range(K):
                figure_title = 'FN '+str(int(i+1))
                brain_map = pFN[:, i]
                plot_FN_brain_surface_5view(brain_map, brain_template, color_function=None, file_output=file_output[i], figure_title=figure_title)

        elif dataType == 'Volume':
            K = pFN.shape[3]
            file_output = [os.path.join(dir_pnet_pFN_indv, str(int(i+1))+'.jpg') for i in range(K)]
            for i in range(K):
                figure_title = 'FN '+str(int(i+1))
                brain_map = pFN[:, :, :, i]
                plot_FN_brain_volume_3view(brain_map, brain_template, color_function=None, file_output=file_output[i], figure_title=figure_title)

        # output an assembled image
        file_output_assembled = os.path.join(dir_pnet_pFN_indv, 'All.jpg')
        assemble_image(file_output, file_output_assembled, interval=(50, 5), background=(0, 0, 0))

    return


def run_Visualization(dir_pnet_result: str):
    """
    Run preconfigured visualizations for gFNs and pFNs

    :param dir_pnet_result: directory of the pnet result folder
    :return:

    Yuncong Ma, 11/6/2023
    """

    run_gFN_Visualization(dir_pnet_result)
    run_pFN_Visualization(dir_pnet_result)

    return
