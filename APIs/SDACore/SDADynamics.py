# -*- coding: utf-8 -*-
"""
@author: siva_
"""

# =============================================================================
# Description
# =============================================================================

 # Code for simulating dynamics 

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
class SDADynamics():
    
    
    #--------------------------------------------------------------------------        
    # Init
    #--------------------------------------------------------------------------        
    def __init__(self):
        a = 1
   
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

    
    
    #--------------------------------------------------------------------------
    # Convert continuous time system to discrete time system
    #--------------------------------------------------------------------------
    def get_DT_sys(self, num, den, dt):
                
        sysc = ctrl.tf2ss(num,den)
        
        sysd = ctrl.c2d(sysc, dt, method = 'zoh')
        
        return(sysd)
    #--------------------------------------------------------------------------
    
    
    #--------------------------------------------------------------------------
    # Get one step update for a given discrete time system, states, inputs
    #--------------------------------------------------------------------------    
    def get_1step_update(self, sysd, x_k_1, u_k):
    
        x_k = (sysd.A*x_k_1) + (sysd.B*u_k)
        y_k = (sysd.C*x_k) + (sysd.D*u_k) 
        
        return(x_k, y_k)
    #--------------------------------------------------------------------------

    
    #--------------------------------------------------------------------------
    # Given an input, compute the output sequence
    #--------------------------------------------------------------------------    
    def sensor_dynamics_DT(self, sysd, x_k_0, u_k):
    
        u_k = np.squeeze(np.array(u_k))
        
        seq_len = np.size(u_k)
        
        [x_k, y_k] = self.get_1step_update(sysd, x_k_0, u_k[0])
        
        y_k_arr = []
        y_k_arr.append(y_k)
        for ii in range(0,seq_len-1,1):
            [x_k, y_k] = self.get_1step_update(sysd, x_k, u_k[ii+1])
            y_k_arr.append(y_k)
        
        y_k_arr = np.squeeze(np.array(y_k_arr))
        
        return(y_k_arr)
    #--------------------------------------------------------------------------

    
    
        
    
    
    
# =============================================================================








