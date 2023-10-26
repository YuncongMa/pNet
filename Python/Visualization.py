# Yuncong Ma, 10/26/2023
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

# other functions of pNet
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


def plot_surface_5view(FN: np.ndarray,
                       brain_template,
                       file_output: str,
                       threshold=98,
                       color_fun='Automatic',
                       figure_title=None):

    # settings for subplot
    fig, axs = matplotlib.pyplot.subplots(nrows=5+1, ncols=1, figsize=(5, 30))
    #fig = plt.figure()
    axs[1] = fig.add_subplot(projection='3d')

    # use if individual color scale for each brain
    threshold_value = np.percentile(np.abs(FN), threshold)

    # plot each sub view
    # Create Poly3DCollection and add to the plot
    vertices = brain_template['Shape']['L']['vertices']
    faces = brain_template['Shape']['L']['faces'] - 1  # change to Python index
    poly3d = [[vertices[vert_id] for vert_id in face] for face in faces]
    ls = LightSource(azdeg=315, altdeg=45)
    poly3d_collection = Poly3DCollection(poly3d, facecolors='cyan', linewidths=1, edgecolors=None, alpha=1, lightsource=ls)
    axs[1].add_collection3d(poly3d_collection)
    axs[1].view_init(elev=0, azim=180, roll=0)

    # Define axis properties
    axs_lim = np.max(np.abs(vertices), axis=0) * 1.1
    axs[1].set_xlim([-axs_lim[0], axs_lim[0]])
    axs[1].set_ylim([-axs_lim[1], axs_lim[1]])
    axs[1].set_zlim([-axs_lim[2], axs_lim[2]])

    axs[1].axis('off')
    axs[1].grid(False)
    axs[1].set_xticks([])
    axs[1].set_yticks([])
    axs[1].set_zticks([])  # For 3D plots

    # save
    fig.savefig(file_output, dpi=250, bbox_inches="tight")


def color_theme(theme: str, parameter: tuple, alpha=1, darker=False):

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

    elif theme == 'Seed_Map_3_Positive':
        min_CC = parameter[0]
        max_CC = parameter[1]
        color_function = np.array((
            (min_CC,0,0,0),
            (min_CC,0.5,0,0),
            (min_CC*0.7 + max_CC*0.3,1,0,0),
            (max_CC,1,1,0)), dtype=np.float32)

    else:
        raise ValueError('Unknown color theme: ' + theme)

    return color_function


def convert_color_function(color_function: np.ndarray, N=256, alpha=1):
    # resample color_map to a color function with fixed step size and length
    N_cf = color_function.shape[0]
    N = 256
    color_map = np.zeros((N, 4), dtype=np.float32)

    for i in range(N):
        value = float(i) / float(N-1) * (color_function[-1, 0] - color_function[0, 0]) + color_map[0, 0]
        if value < color_function[0, 0]:
            color_map[i, 0:3] = color_function[0, 0]
            continue
        elif value > color_function[-1, 0]:
            color_map[i, 0:3] = color_function[-1, 0]
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

    color_map[:, 3] = alpha

    cmap = matplotlib.colors.ListedColormap(color_map)

    return cmap


def plot_brain_surface(brain_map: np.ndarray,
                       brain_template,
                       file_output: str,
                       threshold=98,
                       color_fun='Automatic',
                       figure_title=None,
                       background_color='black'):

    # use if individual color scale for each brain
    threshold_value = np.percentile(np.abs(brain_map), threshold)

    # Prepare BSPolyData for using its plot
    polyData = prepare_BSPolyData(brain_template['Shape']['L']['vertices'], brain_template['Shape']['L']['faces'] - 1)

    # brain surface
    p = surfplot.Plot(surf_lh=polyData, zoom=1, views='lateral')

    # convert map to the mesh surface space based on the brain mask
    Nv_L = sum(brain_template['Brain_Mask']['L'] > 0)
    map_2 = np.zeros(brain_template['Brain_Mask']['L'].shape, dtype=np.float32)
    ps_L = np.where(brain_template['Brain_Mask']['L']>0)[0].astype(int)
    for i in range(int(Nv_L)):
        map_2[ps_L[i]] = brain_map[i]

    # map layer
    max_value = np.percentile(brain_map, 99)
    color_range = (max_value/2, max_value)
    color_function = color_theme('Seed_Map_3_Positive', color_range)
    print(color_function)
    p.add_layer(map_2, cmap=convert_color_function(color_function), color_range=color_range)

    # given that both colorbars have the same range.
    fig = p.build()

    # save
    fig.savefig(file_output, dpi=500, bbox_inches="tight", facecolor=background_color)


def plot_pFN():
    return


def setup_Visualization():
    return


def run_gFN_Visualization():
    return


def run_pFN_Visualization():
    return


def run_Visualization():
    return

# =============== Volume data type =============== #




# =========== Surface-volume data type =========== #
