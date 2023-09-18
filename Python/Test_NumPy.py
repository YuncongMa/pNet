# Yuncong Ma, 9/8/2023

import numpy as np
import scipy
import pNet
import json

temp=np.random.rand(109,91,109,17)
mask = np.random.randint(0,2,(109,91,109))
temp=pNet.reshape_FN(temp, 'Volume', mask)
print(temp.shape)
temp=pNet.reshape_FN(temp, 'Volume', mask)
print(temp.shape)
