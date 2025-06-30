#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 11:19:51 2025

@author: julia
"""

import numpy as np
from scipy.io import loadmat

mat      = loadmat(f"{'/p/julia/ROMS/eccofs/OBS/SST_STATS/sst_statistics.mat'}")
# %%

std_harm = mat['std_harm']
np.savez_compressed('/p/julia/ROMS/eccofs/OBS_python/CombineObs/sst_statistics.npz', std_harm=std_harm)