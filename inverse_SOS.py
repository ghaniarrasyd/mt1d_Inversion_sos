"""Inversion process using SOS Algorithm"""

import numpy as np
import random as r
from MT_1D import FWD_MT1D, misfit_MT
import matplotlib.pyplot as plt

# Create synthetic data
r_layer = [100,10,1000] #resistivity of each layer
t_layer = [600,1200] #thickness of each layer
period = np.logspace(-3,3,55) #period from 10^-3 to 10^3 dividen into 50 points
ndata = len(period)
rho_data, phase_data = FWD_MT1D(r_layer, t_layer, period)