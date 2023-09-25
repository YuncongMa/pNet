# Yuncong Ma, 9/8/2023

import numpy as np
import scipy
# import pNet
import json

import pNet

np_float, np_eps = pNet.set_data_precision('double')
a = np.inf
b = 100
print(np.abs(a - b) / np.maximum(a, 100))