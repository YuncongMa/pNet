# Yuncong Ma, 10/19/2023
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

# other functions of pNet
from Data_Input import *


# =============== Surface data type =============== #


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
