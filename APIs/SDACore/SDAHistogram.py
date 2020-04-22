# -*- coding: utf-8 -*-

# =============================================================================
# Description
# =============================================================================

# Code for SDA's custom histogram and binning steps
 
# =============================================================================


# =============================================================================
# Imports
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../../')

# =============================================================================



# =============================================================================
# Class containing error metric functions
# =============================================================================
class SDAHistogram():

    #--------------------------------------------------------------------------        
    # Init
    #--------------------------------------------------------------------------        
    def __init__(self, bin_boundaries = None):
        self.bin_boundaries = np.array(bin_boundaries)
   
    #--------------------------------------------------------------------------    
    
    
    #--------------------------------------------------------------------------        
    # Reset bin_boundaries
    #--------------------------------------------------------------------------        
    def set_boundaries(self, bin_boundaries):
        self.bin_boundaries = np.array(bin_boundaries)
        return()
    #--------------------------------------------------------------------------
    
    
    #--------------------------------------------------------------------------        
    # Set bin boundaries
    #--------------------------------------------------------------------------        
    def get_boundaries(self, min_val, max_val, no_bins):
        
        bin_size = np.abs(max_val-min_val)/no_bins
        bin_boundaries = np.arange(min_val, max_val+bin_size, bin_size)
        
        return(bin_boundaries)
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------        
    # Get boundary_decisions for upper semi-continuity
    #--------------------------------------------------------------------------        
    def get_USC_decisions(self, bin_heights):
        
        # Decision variables for an upper semi-continuous histogram
        bin_heights_right = np.append(bin_heights[1::],0)            
        # if boundary_decisions = 1, then assign boundary to the bin to the left of the boundary
        boundary_decisions = (bin_heights >= bin_heights_right)
        boundary_decisions = np.append(False,boundary_decisions)
        
        return(boundary_decisions)
    #--------------------------------------------------------------------------
    
    
    #--------------------------------------------------------------------------        
    # Check if given histogram is USC
    #--------------------------------------------------------------------------        
    def is_hist_USC(self, bin_heights, boundary_decisions):
                
        # inputs, convert to numpy arrays
        bin_heights = np.array(bin_heights)         
        
        # Boundary decisions for USC
        USC_decisions = self.get_USC_decisions(bin_heights)
        
        is_USC = np.array_equal(USC_decisions, boundary_decisions)
        
        return(is_USC)
    #--------------------------------------------------------------------------    
    
    
    #--------------------------------------------------------------------------        
    # Evaluate histogram function at points 'xvals'
    # Assumes upper semi-continuous histograms by default
    #--------------------------------------------------------------------------  
    def SDA_hist(self, xvals, bin_heights, bin_boundaries = None, boundary_decisions = None):
        
        # Update optional variables
        if(np.any(bin_boundaries) == None):
            bin_boundaries = self.bin_boundaries
        boundary_decisions = boundary_decisions 
        
        # inputs, convert to numpy arrays
        bin_boundaries = np.array(bin_boundaries)
        bin_heights = np.array(bin_heights) 
        
        xvals = np.array(xvals)
        
        # Assume support for the histogram is between max and min of bin_boundaries 
        if(np.min(bin_boundaries) > -np.inf):
            bin_boundaries = np.append(-np.inf, bin_boundaries)        
            bin_heights = np.append(0,bin_heights)
        
        if(np.max(bin_boundaries) < np.inf):
            bin_boundaries = np.append(bin_boundaries, np.inf)
            bin_heights = np.append(bin_heights,0)
                
        
        if(boundary_decisions == None):
            # Decision variables for an upper semi-continuous histogram
            boundary_decisions = self.get_USC_decisions(bin_heights)
        else:            
            boundary_decisions = np.array(boundary_decisions)
            boundary_decisions = np.append(False,boundary_decisions)
            boundary_decisions = np.append(boundary_decisions,True)
        
        
        # Determine histogram values
        SDA_hist_vals = np.zeros(np.shape(xvals))
        
        right_locs = np.searchsorted(bin_boundaries, xvals, side = 'right')
        left_locs = np.searchsorted(bin_boundaries, xvals, side = 'left')
        
        # If check_rl == false, then the point is on the boundary, else in between
        check_rl = np.equal(right_locs, left_locs)
        
        # histogram values for points in between 
        in_betweens = right_locs[check_rl]
        SDA_hist_vals[check_rl] = bin_heights[in_betweens - 1]
        
        # histogram values for points on the boundary
        # If boundary_decision for boundary == 1, set value to the bin on the left of the boundary
        
        # Boundary that the value lies on
        val_boundary = left_locs[np.logical_not(check_rl)]
        val_boundary_decisions = boundary_decisions[val_boundary]
        
        SDA_hist_vals[np.logical_not(check_rl)] = bin_heights[val_boundary - val_boundary_decisions]
        
        
        return(SDA_hist_vals)    
    #--------------------------------------------------------------------------  
    
    
    
    #--------------------------------------------------------------------------        
    # Build histogram by binning, given xvals, boundaries and boundary_decisions
    #--------------------------------------------------------------------------        
    def get_bin_heights(self, hist_data, bin_boundaries, boundary_decisions):
                
        bin_boundaries = np.array(bin_boundaries)
        boundary_decisions = np.array(boundary_decisions)
        
        # If current boundaries are finite, add inf as endpoints
        if(np.min(bin_boundaries) > -np.inf):
            bin_boundaries = np.append(-np.inf, bin_boundaries)        
            boundary_decisions = np.append(False, boundary_decisions)
            
        if(np.max(bin_boundaries) < np.inf):
            bin_boundaries = np.append(bin_boundaries, np.inf)
            boundary_decisions = np.append(boundary_decisions, True)
            
        boundary_decisions[0] = False
        boundary_decisions[-1] = True
        
        boundary_decisions = boundary_decisions.astype('bool')
        
        # Binning step
        bin_heights = np.zeros(np.size(bin_boundaries)-1)
        
        # Could improve efficiency of the implementation
        for ii in range(0,np.size(bin_boundaries)-1,1):
            
            bin_check = (hist_data > bin_boundaries[ii]) & (hist_data < bin_boundaries[ii+1])
            bin_check = bin_check | ((hist_data == bin_boundaries[ii]) & np.logical_not(boundary_decisions[ii])) \
            | ((hist_data == bin_boundaries[ii+1]) & boundary_decisions[ii+1])
            
            bin_count = np.sum(bin_check)
            bin_heights[ii] = bin_count

        
        return(bin_heights, bin_boundaries)
    #--------------------------------------------------------------------------    
    
    
    #--------------------------------------------------------------------------        
    # Build upper semi-continuous histogram via binning & checking for USC
    #--------------------------------------------------------------------------        
    def build_USC_hist(self, hist_data, bin_boundaries, density = None, boundary_decisions = None):    
        
        # Update optional variables
        density = density or False
        
        boundary_decisions = np.zeros(np.shape(bin_boundaries))
        [bin_heights, bin_boundaries] = self.get_bin_heights(hist_data, bin_boundaries, boundary_decisions)
        
        is_USC = False
        while is_USC == False:
            
            boundary_decisions = self.get_USC_decisions(bin_heights)
            [bin_heights, bin_boundaries] = self.get_bin_heights(hist_data, bin_boundaries, boundary_decisions)
            is_USC = self.is_hist_USC(bin_heights, boundary_decisions)
            
        # Normalization assumes equal bin heights
        if(density == True):
            bin_heights = bin_heights/(np.sum(bin_heights))
            
        return(bin_heights, bin_boundaries)
    #--------------------------------------------------------------------------
    
    
    #--------------------------------------------------------------------------
    # Plot histogram using barchart
    #--------------------------------------------------------------------------
    def plot_hist(self, bin_heights, bin_boundaries, no_zeros = 0 , bin_width = None, suppress_plots = 0):
        
        
        # Update optional variables              
        bin_width = bin_width or False        
        if(bin_width == False):
            bin_width = np.min(np.abs(np.diff(bin_boundaries))) # bin widths for plotting
        

        # inputs, convert to numpy arrays
        bin_boundaries = np.array(bin_boundaries)
        bin_heights = np.array(bin_heights) 


        # If plotting without the zero step size bin
        if(no_zeros == 1):
            index_left = np.searchsorted(bin_boundaries, 0, side = 'left')
            index_right = np.searchsorted(bin_boundaries, 0, side = 'right')
            
            if(index_left == index_right):
                bin_heights[index_left-1] = 0
            else:
                if(bin_heights[index_left-1]>bin_heights[index_left]):
                    bin_heights[index_left-1] = 0
                elif(bin_heights[index_left-1]<bin_heights[index_left]):
                    bin_heights[index_left] = 0
                else:
                    bin_heights[index_left-1] = 0
                    bin_heights[index_left] = 0
    
        # Check for and remove infs in the bin_boundaries just for plotting
        bin_range = np.min(np.abs(np.diff(bin_boundaries)))
        
        if(bin_boundaries[0]==-np.inf):
            bin_boundaries[0] = bin_boundaries[1] - bin_range
            
        if(bin_boundaries[-1]== np.inf):
            bin_boundaries[-1] = bin_boundaries[-2] + bin_range
    

        # Assuming bin_boundaries have one extra entry than the bin_heights, and have no infs
        bin_centers = (bin_boundaries[0:-1] + bin_boundaries[1::])/2
                    
        if(suppress_plots == 0):
            plt.figure()        
            plt.bar(bin_centers, bin_heights, width = bin_width, color=(0.1, 0.1, 0.1, 0.2),  edgecolor='blue', align = 'center')
            plt.xticks(bin_boundaries)
            plt.show()        
        
        return(bin_centers, bin_heights, bin_width, bin_boundaries)
    #--------------------------------------------------------------------------
    
    
    #--------------------------------------------------------------------------
    # Plot histogram using barchart
    #--------------------------------------------------------------------------
    def get_hist_entropy(self, bin_heights, bin_boundaries):
        finite_support_flag = 1
        
        if( (np.min(bin_boundaries)==-np.inf) and (bin_heights[0]!=0) ):
            finite_support_flag = 0
            
        
        if( (np.max(bin_boundaries)== np.inf) and (bin_heights[-1]!=0) ):
            finite_support_flag = 0
            
        if(finite_support_flag == 1):
            
            
            
            if( (np.min(bin_boundaries)==-np.inf) and (bin_heights[0]==0) ):    
                bin_boundaries = bin_boundaries[1::]
                bin_heights = bin_heights[1::]
                
                
            if( (np.max(bin_boundaries)== np.inf) and (bin_heights[-1]==0) ):
                bin_boundaries = bin_boundaries[0:-1]
                bin_heights = bin_heights[0:-1]
                
            bin_widths = np.diff(bin_boundaries)            
            prob = np.multiply(bin_widths, bin_heights)
            hist_mass = np.sum(prob)
            prob = prob/hist_mass
            
            zero_check = (prob != 0)
            
            entropy = - np.sum(np.multiply(prob[zero_check], np.log(prob[zero_check])))
            
        else:
            entropy = -1
            print('INVALID DENSITY: Histogram has infinite mass, not normalizable. ')
            return(entropy)    
        
        
        return(entropy)
    #--------------------------------------------------------------------------
    
    
# =============================================================================    
    

#%%
'''
bin_boundaries = [-np.inf, 0,1,2,3,4, np.inf]
bin_heights = [0,1,2,1,0,0] 
SDAHistObj = SDAHistogram(bin_boundaries) 

SDAHistObj.plot_hist(bin_heights, bin_boundaries)


bin_width = 1


# inputs, convert to numpy arrays
bin_boundaries = np.array(bin_boundaries)
bin_heights = np.array(bin_heights) 


# Check for and remove infs in the bin_boundaries just for plotting
bin_range = np.min(np.abs(np.diff(bin_boundaries)))

if(bin_boundaries[0]==-np.inf):
    bin_boundaries[0] = bin_boundaries[1] - bin_range
    
if(bin_boundaries[-1]== np.inf):
    bin_boundaries[-1] = bin_boundaries[-2] + bin_range

bin_centers = (bin_boundaries[0:-1] + bin_boundaries[1::])/2

                
plt.figure()        
plt.bar(bin_centers, bin_heights, width = bin_width, color=(0.1, 0.1, 0.1, 0.2),  edgecolor='blue', align = 'center')
plt.xticks(bin_boundaries)
plt.show()      
    

entropy = SDAHistObj.get_hist_entropy(bin_heights, bin_boundaries)
print(entropy)
'''
#%%
'''
bin_boundaries = [1,2,3,4]
bin_heights = [111,-111,222]


SDAHistObj = SDAHistogram(bin_boundaries) 
xvals = [3, 3.5, 4.5, 2.3, 2, -10, 120, 1, np.inf]

SDA_hist_vals = SDAHistObj.SDA_hist(xvals, bin_heights)

boundary_decisions = [1, 0, 1, 0]
SDA_hist_vals_2 = SDAHistObj.SDA_hist(xvals, bin_heights, boundary_decisions = boundary_decisions)

#bin_boundaries = SDAHistObj.get_boundaries(-10, 10, 10)


hist_data = [2, 2.3, 2.3, 1, 1, 1.3, 3.5, 500, 500, 500]
boundary_decisions = np.zeros(np.shape(bin_boundaries))

[bin_heights, bin_boundaries] = SDAHistObj.get_bin_heights(hist_data, bin_boundaries, boundary_decisions)
is_USC = SDAHistObj.is_hist_USC(bin_heights, boundary_decisions)


# Keeping binning till USC histogram is formed
hist_data = [2, 2.3, 2.3, 1, 1, 1.3, 3.5, 500, 500, 500]

[bin_heights, boundary_decisions] = SDAHistObj.build_USC_hist(hist_data, bin_boundaries) 
'''

    

    
    