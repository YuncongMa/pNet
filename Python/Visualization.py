# Yuncong Ma, 10/30/2023
# Visualization module of pNet

#########################################
# Packages
import numpy as np
import os
import re
import time

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import surfplot
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import LightSource
import vtk
from brainspace.vtk_interface.wrappers import BSPolyData
from brainspace.vtk_interface.wrappers import BSPolyDataMapper
import nilearn as nl

from brainspace.mesh.mesh_creation import build_polydata
from matplotlib.font_manager import FontProperties

# other functions of pNet
import pNet
from Data_Input import *


# =============== Surface data type =============== #


def prepare_BSPolyData(vertices: np.ndarray, faces: np.ndarray):
    """
    prepare a BSPolyData class

    :param vertices: 2D matrix [N, 3], N is the number of vertices
    :param faces: 2D matrix, [N, 3], N is the number of triangle faces, vertex index starts from 0
    :return: bspolydata, a data class for using brain space surface plot

    Yuncong Ma, 10/26/2023
    """

    bspolydata = build_polydata(points=vertices, cells=faces)

    if not isinstance(bspolydata, BSPolyData):
        raise ValueError('Cannot generate BSPolyData for a 3D surface mesh')

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
        color_function = np.arra((
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
        color_function = np.arra((
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
        color_function = np.arra((
            (0,0,0,0),
            (Threshold,0,0,0),
            (Threshold,1,(Threshold-min_CC)/(max_CC-min_CC),0),
            (max_CC,1,1,0)), dtype=np.float32)

    elif theme == 'Seed_Map_Positive_Only':
        min_CC = parameter[0]
        Threshold = parameter[1]
        max_CC = parameter[2]
        color_function = np.arra((
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
        cmap = plt.get_cmap(map_name)
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


def plot_brain_surface(brain_map: np.ndarray,
                       mesh: dict,
                       mask: np.ndarray,
                       color_function: np.ndarray,
                       file_output=None or str,
                       orientation='medial',
                       view_angle=1.5,
                       mask_color=(0.2, 0.2, 0.2),
                       brain_color=(0.5, 0.5, 0.5),
                       background_color=(0,0,0),
                       figure_size=(500,400),
                       dpi=25):

    # Prepare BSPolyData for using its plot
    polyData = prepare_BSPolyData(mesh['vertices'], mesh['faces'] - 1)
    p = surfplot.Plot(surf_lh=polyData, zoom=view_angle, views=orientation, background=background_color, brightness=1, size=figure_size)

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
    #max_value = np.percentile(brain_map, map_threshold)
    #map_2[np.abs(map_2) < max_value/2] = 0
    color_range = (color_function[0, 0], color_function[-1, 0])
    p.add_layer(map_2, cmap=prepare_color_map(color_function=color_function), color_range=color_range, cbar=None, zero_transparent=True)

    # save or return
    if file_output is not None:
        # build the figure
        fig = p.build()
        fig.savefig(file_output, dpi=dpi, bbox_inches="tight", facecolor=background_color)
    else:
        return p


def merge_mesh_LR(mesh_LR: dict, offset=np.array((90, 0, 0))):
    mesh = {'vertices': np.concatenate((mesh_LR['L']['vertices'], mesh_LR['R']['vertices'] + offset), axis=0),
            'faces': np.concatenate((mesh_LR['L']['faces'], mesh_LR['R']['faces']+mesh_LR['L']['vertices'].shape[0]), axis=0)}
    return mesh


def merge_mask_LR(mask_LR: dict):
    mask = np.concatenate((mask_LR['L'], mask_LR['R']), axis=0)
    return mask


def colorize(value_map: np.ndarray, color_function: np.ndarray):
    color_map = ()
    return color_map


def plot_FN_brain_surface_5view(brain_map: np.ndarray,
                             brain_template,
                             file_output: str,
                             threshold=99,
                             color_function=None,
                             background_color=(0,0,0),
                             figure_organization=(0.6, 1.2, 1, 0.6),
                             view_angle=(1.35, 1.4),
                             hemisphere_offset=90,
                             figure_title=None,
                             title_font_dic=dict(fontsize=20, fontweight='bold'),
                             figure_size=(10, 50),
                             dpi=50):

    # settings for subplot
    fig, axs = matplotlib.pyplot.subplots(nrows=7+1, ncols=1, figsize=figure_size)

    # set color function
    if color_function is None:
        threshold_value = np.percentile(np.abs(brain_map), threshold)
        color_range = np.array((threshold_value/2, threshold_value))
        color_function = color_theme('Seed_Map_3_Positive', color_range)
    else:
        color_range = (color_function[0,0], color_function[-1,0])

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
    axs[0].facecolor = background_color
    axs[0].axis('off')

    # dorsal view
    p = plot_brain_surface(brain_map, mesh=merge_mesh_LR(brain_template['Shape'], offset=np.array((hemisphere_offset, 0, 0))), mask=merge_mask_LR(brain_template['Brain_Mask']),
                           color_function=color_function,
                           orientation='dorsal', view_angle=view_angle[0], file_output=None, background_color=(1,1,1),
                           figure_size=(int(dpi*figure_size[0]), int(dpi*H_D*figure_size[1])), dpi=dpi)
    p = p.render()
    p._check_offscreen()
    x = p.to_numpy(transparent_bg=True, scale=(2, 2))
    axs[1].figure_size = (int(dpi*figure_size[0]), int(dpi*H_D*figure_size[1]))
    axs[1].set_position((0, 4*H_S+H_C, 1, H_D))
    axs[1].imshow(x)
    axs[1].axis('off')
    axs[1].set_title(label=figure_title, loc='center', pad=140, fontsize=150, fontweight='bold', color=(1, 1, 1))

    # saggital views
    # 1st
    p = plot_brain_surface(brain_map[0:Nv_L], mesh=brain_template['Shape']['L'], mask=brain_template['Brain_Mask']['L'],
                           color_function=color_function,
                           orientation='lateral', view_angle=view_angle[1], file_output=None,
                           figure_size=(int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1])), dpi=dpi)
    p = p.render()
    p._check_offscreen()
    x = p.to_numpy(transparent_bg=True, scale=(2, 2))
    axs[2].figure_size = (int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1]))
    axs[2].set_position((0, 3*H_S+H_C, 1, H_S))
    axs[2].imshow(x)
    axs[2].axis('off')

    # 2nd
    p = plot_brain_surface(brain_map[0:Nv_L], mesh=brain_template['Shape']['L'], mask=brain_template['Brain_Mask']['L'],
                           color_function=color_function,
                           orientation='medial', view_angle=view_angle[1], file_output=None, dpi=dpi)
    p = p.render()
    p._check_offscreen()
    x = p.to_numpy(transparent_bg=True, scale=(2, 2))
    axs[3].figure_size = (int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1]))
    axs[3].set_position((0, 2*H_S+H_C, 1, H_S))
    axs[3].imshow(x)
    axs[3].axis('off')
    # 3rd
    p = plot_brain_surface(brain_map[Nv_L:], mesh=brain_template['Shape']['R'], mask=brain_template['Brain_Mask']['R'],
                           color_function=color_function,
                           orientation='medial', view_angle=view_angle[1], file_output=None, dpi=dpi)
    p = p.render()
    p._check_offscreen()
    x = p.to_numpy(transparent_bg=True, scale=(2, 2))
    axs[4].figure_size = (int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1]))
    axs[4].set_position((0, 1*H_S+H_C, 1, H_S))
    axs[4].imshow(x)
    axs[4].axis('off')
    # 4th
    p = plot_brain_surface(brain_map[Nv_L:], mesh=brain_template['Shape']['R'], mask=brain_template['Brain_Mask']['R'],
                           color_function=color_function,
                           orientation='lateral', view_angle=view_angle[1], file_output=None, dpi=dpi)
    p = p.render()
    p._check_offscreen()
    x = p.to_numpy(transparent_bg=True, scale=(2, 2))
    axs[5].figure_size = (int(dpi*figure_size[0]), int(dpi*H_S*figure_size[1]))
    axs[5].set_position((0, H_C, 1, H_S))
    axs[5].imshow(x)
    axs[5].axis('off')

    # color bar
    axs[6].figure_size = (int(dpi*figure_size[0]), int(dpi*H_C*figure_size[1]))
    colorbar_width = 0.8
    colorbar_height = 0.2
    colorbar_pad = 0.6
    colorbar_step = 100
    colorbar_ratio = 10
    colorbar_scale = 100
    axs[6].set_position(((1-colorbar_width)/2, H_C*colorbar_pad, colorbar_width, H_C*colorbar_height))
    X = np.tile(np.arange(color_range[0], color_range[1], (color_range[1] - color_range[0])/colorbar_step), (colorbar_ratio, 1))
    axs[6].imshow(X, cmap=prepare_color_map(color_function=color_function))
    axs[6].axis('off')
    cb_font_pad = 30
    cb_value_round = True
    if cb_value_round is True:
        color_range = np.round(color_range * colorbar_scale)
        cb_tick = (str(int(color_range[0])), str(int(color_range[1])))
    else:
        color_range = color_range * colorbar_scale
        cb_tick = (str(color_range[0]), str(color_range[1]))

    axs[6].text(0, cb_font_pad, cb_tick[0],
                ha='left', va='bottom',
                fontdict=dict(fontsize=100, fontweight='bold', color=(1, 1, 1)))
    axs[6].text(colorbar_step, cb_font_pad, cb_tick[1],
                ha='right', va='bottom',
                fontdict=dict(fontsize=100, fontweight='bold', color=(1, 1, 1)))

    # add a block to maintain the preconfigured figure size
    axs[7].figure_size = (int(dpi*figure_size[0]), int(dpi*H_C*figure_size[1]))
    axs[7].set_position(((1-colorbar_width)/2, 0, colorbar_width, H_C))
    axs[7].axis('off')

    # save fig
    fig.savefig(file_output, dpi=dpi, bbox_inches="tight", facecolor=background_color)


def assemble_image(file_output: tuple, file_output_assembled: str, organization=(0, 10), interval=(50, 5), background_color=(0, 0, 0)):

    return


def plot_pFN():
    return


def setup_Visualization(file_figure, ):
    return


def run_gFN_Visualization(dir_pnet_result: str):

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

    if dataType == 'Surface' and dataFormat == 'HCP Surface (*.cifti, *.mat)':
        gFN = load_matlab_single_array(os.path.join(dir_pnet_gFN, 'FN.mat'))
        brain_template = load_brain_template(os.path.join(dir_pnet_dataInput, 'Brain_Template.json'))
        K = gFN.shape[1]
        file_output = [os.path.join(dir_pnet_gFN, str(int(i+1))+'.jpg') for i in range(K)]
        for i in range(K):
            figure_title = 'FN '+str(int(i+1))
            plot_FN_brain_surface_5view(gFN[:, i], brain_template, color_function=None, file_output=file_output[i], figure_title=figure_title)

        file_output_assembled = os.path.join(dir_pnet_gFN, 'All.jpg')
        assemble_image(file_output, file_output_assembled, interval=(50, 5), background_color=(0, 0, 0))

    return


def run_pFN_Visualization():
    return


def run_Visualization():
    return

# =============== Volume data type =============== #




# =========== Surface-volume data type =========== #
