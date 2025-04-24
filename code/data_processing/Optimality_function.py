# -*- coding: utf-8 -*-
"""
This script is an adapted version of the Max van Gerrevink's python script. Changed all functions to vectorized execution for Optimality calculation.

Changes made by Nils Rietze on 7 April 2025.

---- original docstring below:

Created on 14 March 2020
This script is an adapted version of Sander Veraverbeke's matlab sript.
It calculates the optimality based on differences in pre and post-fire 
spectral values. Optimality is defined by the direction and magnitude of 
pixel displacements in the bi-spectral space (Verstraete & Pinty, 1996)

This script uses 1-D array data as an inputfile, Optimality output is dependent on input size.
All input files must have the same shape.

Required packages: numpy version 1.20.2

To do this in your own script use the following code:
    root_folder = 'C:\VU_Amsterdam\Research_project\scripts'  # Setting rootfolder where you save this file
    os.chdir(root_folder) # Working directory
    from OPTIMALITY_function import Spectral_optimality as OPT
    
By Max van Gerrevink
"""
#%% Packages
import numpy as np
#%% Optimality function

def Spectral_optimality(xpre, ypre, xpost, ypost):
       
    NoData = -9999
    xpre[np.isnan(xpre)] = NoData # Replacing nans with the no data value
    ypre[np.isnan(ypre)] = NoData # Replacing nans with the no data value
    xpost[np.isnan(xpost)] = NoData # Replacing nans with the no data value
    ypost[np.isnan(ypost)] = NoData # Replacing nans with the no data value
    
    # Optimality equation
    bperp = xpre + ypre
    bperp = np.where(bperp==0.02,np.nan,bperp)

    with np.errstate(divide='ignore'):
        apost = ypost / xpost
        apost = np.where(apost==1,np.nan,apost)
        
    with np.errstate(divide='ignore'):
        xideal = bperp / (apost + 1)

    with np.errstate(divide='ignore'):
        yideal = xideal * apost

    with np.errstate(divide='ignore'):
        dinsens = np.sqrt(((xpost-xideal)**2) + ((ypost-yideal)**2))

    with np.errstate(divide='ignore'):
        dsens = np.sqrt(((xpost-xpre)**2) + ((ypost-ypre)**2))

    with np.errstate(divide='ignore'):
        Opt = 1 - (dinsens/dsens)
        Opt = np.where(Opt<0,0,Opt)   
    
    # Get optimality
    return Opt