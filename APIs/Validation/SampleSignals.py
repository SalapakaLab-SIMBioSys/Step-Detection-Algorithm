# -*- coding: utf-8 -*-


# =============================================================================
# Description
# =============================================================================

# Code for sample stepping signals

# =============================================================================
 
 
# =============================================================================
# Imports
# ============================================================================= 
import numpy as np
import matplotlib.pyplot as plt

import control as ctrl

import sys
sys.path.append('../../')
# =============================================================================

# =============================================================================
# Class containing cost functions 
# =============================================================================
class SampleSignals():
    
    
    #--------------------------------------------------------------------------        
    # Init
    #--------------------------------------------------------------------------        
    def __init__(self):
        a = 1
   
    #--------------------------------------------------------------------------       
    
    
    #--------------------------------------------------------------------------        
    # Get step signal with randomly placed steps (no noise added)
    #--------------------------------------------------------------------------      
    def get_step_sample(self, step_size = 1, no_steps = 10, no_samples = 1000, bidirectional = 1, bidirec_prob = 0.5, seed = None):
    
        # Initialize both randomizations with seed if a valid value is given 
        if(seed is not None):
            np.random.seed(seed)
        
        seed_arr = np.random.randint(2**31, size = 2)
        
        L = step_size
        d = no_steps
        N = no_samples
        
        np.random.seed(seed_arr[0])
        s = np.floor(np.random.rand(no_steps)*N)
        s = s.astype(int)
        
        x = np.matrix(np.zeros((1,N)))
        
        np.random.seed(seed_arr[1])
        if(bidirectional == 1 ):
            x[0,s] = L*np.sign(bidirec_prob-np.random.rand(1,d))
        else:
            x[0,s] = L
            
        x = np.cumsum(x)   # Underlying signal
    
        step_signal = x
        
        return(step_signal)
    #--------------------------------------------------------------------------  
    

    #--------------------------------------------------------------------------        
    # Get step signal with randomly placed steps (no noise added) with multiple step sizes
    #--------------------------------------------------------------------------      
    def get_m_step_sample(self, step_sizes = [1], no_steps = 10, no_samples = 1000, bidirectional = 1, bidirec_prob = 0.5, step_size_prob = None, seed = None):
        
        # Initialize both randomizations with seed if a valid value is given 
        if(seed is not None):
            np.random.seed(seed)
        
        seed_arr = np.random.randint(2**31, size = 2)
        
        L = step_sizes
        d = no_steps
        N = no_samples
        
        np.random.seed(seed_arr[0])
        s = np.floor(np.random.rand(no_steps)*N)
        s = s.astype(int)
        
        x = np.matrix(np.zeros((1,N)))
        
        
        np.random.seed(seed_arr[1])
        for ii in range(d):
            
            L_local = np.random.choice(L, p = step_size_prob)
            if(bidirectional == 1):        
                x[0,s[ii]] = L_local*np.sign(bidirec_prob-np.random.rand(1))
            else:
                x[0,s[ii]] = L_local
        
        x = np.cumsum(x)   # Underlying signal
        
        step_signal = x            
    
        return(step_signal)
    #--------------------------------------------------------------------------  
    
    #--------------------------------------------------------------------------        
    # Get Gaussian noise (for additive noise)
    #--------------------------------------------------------------------------      
    def get_gaussian_noise(self, std = 1, mean = 0, no_samples = 1000, step_signal = None, seed = None):
        
        if(seed is not None):
            np.random.seed(seed)
        
        if(step_signal is None):
            noise = std*np.random.randn(no_samples) + mean
            
        else:
            shape_signal = np.shape(step_signal)
            noise = std*np.random.randn(shape_signal[0],shape_signal[1]) + mean
        
        return(noise)
    #--------------------------------------------------------------------------  
    
    
    #--------------------------------------------------------------------------        
    # Get gamma noise (for additive noise)
    #--------------------------------------------------------------------------      
    def get_gamma_noise(self, shape = 1, scale = 1, no_samples = 1000, step_signal = None, seed = None):
        
        if(seed is not None):
            np.random.seed(seed)
        
        if(step_signal is None):
            noise = np.random.gamma(shape, scale, size = no_samples)
            
        else:
            shape_signal = np.shape(step_signal)
            noise = np.random.randn(shape, scale, size = shape_signal)
        
        return(noise)
    #--------------------------------------------------------------------------   
    
    
    #--------------------------------------------------------------------------        
    # simulate sensor dynamics (linear transfer function output)
    #--------------------------------------------------------------------------    
    def sensor_dynamics(self, num, den, dt, U):

        
        G = ctrl.tf(num,den)
        
        inp_len = np.size(U)
        
        T = np.arange(0,inp_len)*dt
        
        [T, yout, xout] = ctrl.forced_response(G, T, U)

        return(T, yout)
    #--------------------------------------------------------------------------        
    
    
# =============================================================================    
    
'''
#%%

sigObj = SampleSignals()

x = sigObj.get_step_sample()

noise = sigObj.get_gaussian_noise(std = 0.2, step_signal=x, seed = 10)
    
noise2 = sigObj.get_gaussian_noise(std = 0.2, step_signal=x, seed = 10)
y = x + noise

plt.figure()

plt.plot(y.T)
plt.plot(x.T)  
plt.plot(noise.T)
plt.plot(noise2.T)

#%%
x1 = sigObj.get_step_sample(seed = 1)
x2 = sigObj.get_step_sample(seed = 1)

plt.figure()
plt.plot(x1.T)
plt.plot(x2.T)


#%%
g1 = sigObj.get_gamma_noise(shape = 1,scale = 1)

y = x + g1

plt.figure()
plt.plot(y.T)
plt.plot(x.T)  
plt.plot(g1.T)

'''