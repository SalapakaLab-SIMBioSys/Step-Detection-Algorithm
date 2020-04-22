# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:20:49 2019

@author: siva_
"""

# =============================================================================
# Description
# =============================================================================

# Code for Alternate step detection methods

# =============================================================================
 
 
# =============================================================================
# Imports
# ============================================================================= 
import numpy as np
import matlab.engine


import sys
sys.path.append('../../')

from APIs.AlternateTools.tvr import tv
# =============================================================================
#%%
# =============================================================================
# Class containing cost functions 
# =============================================================================
class AlternateMethods():
    
    #--------------------------------------------------------------------------        
    # Init
    #--------------------------------------------------------------------------        
    def __init__(self):        
        
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(r'../../APIs/AlternateTools/HMMsForMotors_2/HMMsForMotors 2/')
        self.eng.addpath(r'../../APIs/AlternateTools/StepFitOriginal/')
        self.eng.addpath(r'../../APIs/AlternateTools/Kerssemakers/')
    #--------------------------------------------------------------------------    

    #--------------------------------------------------------------------------
    # Total Variation Regularization
    #--------------------------------------------------------------------------
    def total_variation_regularization(self,  noise_std, y):                            
        
        # TV as MAP estimate (Algorithm 3.1 )
        #diff_y = np.diff(y)
        #std_diff_y = np.std(diff_y)
        #beta = std_diff_y/np.sqrt(2)        
        #delta = (noise_std**2)/beta

        
        delta = (noise_std**2)*10 #19.3069
        
        x_estimate = tv(delta, y.T)
                
        return(x_estimate)
    #--------------------------------------------------------------------------    
    
    #--------------------------------------------------------------------------
    # HMM step fit (Siggworth)
    #--------------------------------------------------------------------------
    def HMMsForMotors(self, y):                            
        
        y = matlab.double(y[0].tolist())
        
        x_estimate = self.eng.GetStepsHMM1(y)       
        
        x_estimate = np.asarray(x_estimate)
                
        return(x_estimate)
    #--------------------------------------------------------------------------     
    
    
    #--------------------------------------------------------------------------
    # SDA Original (Tanuj Aggarwal, Matlab implementation)
    #--------------------------------------------------------------------------
    def SDA_original(self, y, Ts = None, tau = 0):                            
        
        y = matlab.double(y[0].tolist())

        if(Ts == None):
            # No sensor dynamics
            x_estimate = self.eng.stepfit(y)       
        else:
            # Estimate with sensor dynamics
            x_estimate = self.eng.stepfit(y,'Fs',1/Ts,'tau',tau,'verbose',0)
        
        x_estimate = np.asarray(x_estimate)
                
        return(x_estimate)
    #--------------------------------------------------------------------------       
    
    
    #--------------------------------------------------------------------------
    # Kerssemakers step fit (Kerssemakers' method, Matlab implementation)
    #--------------------------------------------------------------------------    
    def Kerssemakers(self, y):
        # Make sure y is a column vector.
        
        
        y = matlab.double(y[0].tolist())
        
        x_estimate = self.eng.call_Kerssemakers(y)
        
        x_estimate = np.asarray(x_estimate)
        
        return(x_estimate)
    #--------------------------------------------------------------------------       

    
# =============================================================================    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    