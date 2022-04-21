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

# Noise standard deviation
sigma = 0.3162 #0.3162 Corresponds to an SNR of 5 for signal of size 1 and a duty ratio of 0.5 


# =============================================================================
# Generating noisy stepping data 
# =============================================================================
d = 10 # No. of steps
dt = 0.1 # Sampling interval in s
N = 1000 # total number of samples 
L = 1 # Step size
iterations = 50 # No. of SDA iterations
resolution = 0.05 # Resolution of the SDA
seed = 3 # For reproducibility (change seed for a new signal)
include_dynamics = True # Choose to simulate sensor dynamics 

# Simulating sensor dynamics
SDynObj = SDADynamics()
num = [1] # Specify numerator of sensor transfer function
den = [0.1, 0.8, 1, 1] # Specify denominator of sensor transfer function
sysd = SDynObj.get_DT_sys(num, den, dt)


# Generating a sample stepping signal
SIGObj = SampleSignals()
# Stepping signal
x = SIGObj.get_step_sample(step_size = L, no_steps = d, no_samples = N, bidirectional = 0, seed = seed)
# Noise
noise = SIGObj.get_gaussian_noise(std = sigma, no_samples = N, seed = seed)

if(include_dynamics == True):
    [T, x_sensor] = SIGObj.sensor_dynamics(num, den, dt, x)
    x_sensor = np.expand_dims(x_sensor, axis = 0)
else:
    x_sensor = x
    sysd = None

y = x_sensor + noise  # Noisy measurements

# =============================================================================

# To use the SDA to estimate the stepping signal from data, 
# you would need the standard deviation of noise in the signal.
# Here the standard deviation of noise is denoted by 'sigma'.  
# If you believe your data is curropted by sensor dynamics, 
# specify the sensor dynamics. Here this is denoted by 'sysd'.  
SDAObj = SDA()

[step_estimate, est_array, hist_array]= SDAObj.SDA_dynamics(sigma, y, diagnostics = 0, sysd = sysd)

SDAObj.plot_SDA_stages(x.T, y, est_array, hist_array)











    
    