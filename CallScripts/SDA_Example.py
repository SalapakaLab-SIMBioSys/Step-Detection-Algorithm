# -*- coding: utf-8 -*-

# =============================================================================
# Description
# =============================================================================

 # Code for running and validating the SDA

# =============================================================================


# =============================================================================
# Imports
# =============================================================================
import random
import numpy as np
import matplotlib.pyplot as plt
import control as ctrl


import os

import sys
sys.path.append('../')

# my toolkits
from APIs.SDACore.CostFunctions import CostFunctions
from APIs.SDACore.SDAHistogram import SDAHistogram
from APIs.SDACore.SDA import SDA
from APIs.Validation.ErrorMetrics import ErrorMetrics
from APIs.Validation.SampleSignals import SampleSignals
from APIs.SDACore.SDADynamics import SDADynamics
# =============================================================================

#%%

# =============================================================================
# Generating noisy stepping data 
# =============================================================================
d = 10 # No. of steps
dt = 0.1 # Sampling interval
N = 1000 # total number of samples 
L = 1 # Step size
iterations = 50 # No. of SDA iterations
resolution = 0.05 # Resolution of the SDA

# Noise standard deviation
sigma = 0.3162 #0.3162 = SNR of 5 for signal of size 1 and a duty ratio of 0.5 

# Sensor dynamics
SDynObj = SDADynamics()
num = [1]
den = [0.1, 0.8, 1, 1]
sysd = SDynObj.get_DT_sys(num, den, dt)

include_dynamics = 0

sigObj = SampleSignals()
seed = 2
x = sigObj.get_step_sample(step_size = L, no_steps = d, no_samples = N, bidirectional = 0, seed = seed)
noise = sigObj.get_gaussian_noise(std = sigma, no_samples = N, seed = seed)


if(include_dynamics == 1):
    [T, x_sensor] = sigObj.sensor_dynamics(num, den, dt, x)
    x_sensor = np.expand_dims(x_sensor, axis = 0)
else:
    x_sensor = x
    sysd = None

y = x_sensor + noise  # Noisy measurements

plt.figure()
plt.plot(y.T)
plt.plot(x.T)
plt.plot(x_sensor.T)

# =============================================================================


#%%

SDAObj = SDA()

#[step_estimate, est_array, hist_array]= SDAObj.SDA_standard(sigma, y, resolution = resolution,  iterations = iterations, make_plots = 0, diagnostics = 1)

[step_estimate, est_array, hist_array]= SDAObj.SDA_dynamics(sigma, y, resolution = resolution,  iterations = iterations, hist_tol = 0.001, make_plots = 0, diagnostics = 1, sysd = sysd)

SDAObj.plot_SDA_stages(x.T, est_array, hist_array, no_zeros = 1)

#%%










    
    