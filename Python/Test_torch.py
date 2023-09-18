# Test some PyTorch functions

import torch
import numpy as np

X = torch.rand((5,3))
Y = torch.rand((5,3))

print(torch.min(X, dim=0).values)
print(torch.max(X, dim=0)[0])
