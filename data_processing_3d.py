# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 16:38:59 2023

@author: naina
"""
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import json
from datetime import date, datetime, timedelta, time
import pickle
#from mayavi import mlab
import wrf

# In[]
def velocity_magnitude_compute(UTS,VTS,WTS):
    #WTS_destag = wrf.destagger(WTS,0,meta=True)
    U_square = (UTS)**2
    V_square = (VTS)**2
    W_square = (WTS)**2
    
    UMAG = np.sqrt(U_square+V_square+W_square)
    print(W_square.shape)
    print(UMAG.shape)
    
    return(UMAG)
# In[]
def vorticty_compute(UTS,VTS,WTS):  
    #WTS_destag = wrf.destagger(WTS,0,meta=True)
       
    dU_dx, dU_dy, dU_dz = np.gradient(UTS)
    dV_dx, dV_dy, dV_dz = np.gradient(VTS)
    dW_dx, dW_dy, dW_dz = np.gradient(WTS)
    vorticity_x = dW_dy-dV_dz
    vorticity_y = dU_dz-dW_dx
    vorticity_z = dV_dx-dU_dy
    vorticity_MAG = np.sqrt(vorticity_x**2+vorticity_y**2+vorticity_z**2)
    
    return(vorticity_MAG)
    
# In[]
def QCriterion_compute(UTS,VTS,WTS):
    #WTS_destag = wrf.destagger(WTS,0,meta=True)
    
    dU_dx, dU_dy, dU_dz = np.gradient(UTS)
    dV_dx, dV_dy, dV_dz = np.gradient(VTS)
    dW_dx, dW_dy, dW_dz = np.gradient(WTS)
    
    del_V = np.array([[dU_dx, dU_dy, dU_dz], [dV_dx, dV_dy, dV_dz], [dW_dx, dW_dy, dW_dz]])
    del_V_T = np.array([[dU_dx, dV_dx, dW_dx], [dU_dy, dV_dy, dW_dy], [dU_dz, dV_dz, dW_dz]])

    S = 0.5*(del_V + del_V_T)
    W = 0.5*(del_V - del_V_T)

    Q_criterion = np.array(0.5*(np.square(np.linalg.norm(W, axis=(0,1)))-np.square(np.linalg.norm(S, axis=(0,1)))))

    return(Q_criterion)

# In[]
def create_data_3d (nc_data, tstart, tend, tstamp, qoi):
    pickled_data = {}
    
    # In[] Start time stamp to identify the beginning of a data set 
    pickled_data['start_time'] = tstart
    
    # In[]
    # Extract data size
    dX = np.size(nc_data['west_east'])
    dY = np.size(nc_data['south_north'])
    pickled_data.update({'dX': dX, 'dY': dY})
    
    # In[]
    # Find number of time stamps, z-locations, and axial locations for sampling
    n_time_stamps = nc_data.dims['Time']
    pickled_data.update({'n_time_stamps': n_time_stamps})
    WTS_destag = wrf.destagger(nc_data.WTS,1,meta=True)
    Q_crit = {}
    Vort = {}
    
    # In[]
    for time_stamp in range (n_time_stamps):
        U_0 = nc_data.UTS.isel(Time=time_stamp)
        V_0 = nc_data.VTS.isel(Time=time_stamp)
        W_0 = WTS_destag.isel(Time=time_stamp)
        Vort[time_stamp] = vorticty_compute(U_0,V_0,W_0)
        Q_crit[time_stamp] = QCriterion_compute(U_0,V_0,W_0)
    
    # In[]
    pickled_data[qoi[0]] = WTS_destag
    pickled_data[qoi[1]] = Vort
    pickled_data[qoi[2]] = Q_crit
    
    return pickled_data
    
    