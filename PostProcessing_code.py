#!/usr/bin/env python
# coding: utf-8

# In[]:
# Import Modules
from data_processing import *
from data_processing_3d import *
from plotting import *
from plottingIsosurface import *

import os
import sys
#import wrf
import numpy as np
import math
import xarray as xr
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib import gridspec
import json
import pickle
from datetime import date, datetime, timedelta, time

plt.style.use('seaborn-white')

# In[]
# Load NetCDF file containing tsout data and define plot location
netCDF_file_loc = 'C:/Users/naina/OneDrive/Documents/LLNL_Internship_2022/AlphaVentus-main/test_output_small'
wrf_domain6 = xr.open_dataset(netCDF_file_loc)

proc_data_loc = 'C:/Users/naina/OneDrive/Documents/LLNL_Internship_2022/AlphaVentus-main/'

'''# In[]:
# reference time (based on NetCDF file name) and sampling time interval
ref_time = netCDF_file_loc.split('_')
dt = 10 # sec

# In[] pickle file name
#start_time = datetime.fromisoformat(ref_time[-2]+ '_' + ref_time[-1]) - timedelta(seconds = dt)
start_time = datetime.fromisoformat(ref_time[-2]+ '_' + '00:00:00')
#start_time_stamp = start_time.isoformat('_') #.split('_')[1]
#pickle_filename = f"pickled_%s.pkl"%start_time_stamp
pickle_filename = "plot_testing.pkl"
pickle_file = os.path.join(proc_data_loc, pickle_filename)

# In[]:
# Variables of interest
plane = ['xy','yz','xz']

qoi_units_map = {'UTS': 'm/s',
                 'VTS': 'm/s',
                 'WTS': 'm/s',
                 'TKETS': 'm2/s2'
                }
qoi_list = list(qoi_units_map.keys())
#qoi_list = ['UTS']
qoi_units = list(qoi_units_map.values())

qoi_plot_map = {'UMAG': 'm/s',
                'TKE_RES': 'm2/s2',
                'TKE_SGS': 'm2/s2',
                'TKE_TOT': 'm2/s2',
                'W'   : 'm/s'
               }
qoi_plot_avg_map = {'UMAG_AVG': 'm/s',
                    'TKE_RES_AVG': 'm2/s2',
                    'TKE_SGS_AVG': 'm2/s2',
                    'TKE_TOT_AVG': 'm2/s2',
                    'W_AVG'   : 'm/s'
                   }

qoi_range_map = {'UMAG'    : [8.0, 20.0],
                 'UMAG_AVG': [8.0, 20.0],
                 'TKE_RES' : [0.0, 10.0],
                 'TKE_RES_AVG' : [0.0, 10.0],
                 'TKE_SGS' : [0.0, 0.6],
                 'TKE_SGS_AVG' : [0.0, 0.6],
                 'W'       : [-1.35, 1.35],
                 'W_AVG'   : [-1.35, 1.35]
    }
'''
z_plane_locs = [25]

vert_line_locs = [1] #D

#frac_time = 0.04 # Fraction of time series to use
frac_time = 1.00 # Fraction of time series to use

'''# In[]
pickled_data = create_data (wrf_domain6, qoi_units_map, z_plane_locs, vert_line_locs, ref_time, dt, start_time, frac_time)
with open(pickle_file, 'wb') as pickle_file_handle:
    pickle.dump(pickled_data, pickle_file_handle)

# In[]
plot_contours_instantaneous(pickle_file, proc_data_loc, qoi_plot_map, qoi_range_map, [250, 350], [250, 350])
#plot_contours_instantaneous(pickle_file, proc_data_loc, qoi_plot_map, qoi_range_map)

# In[]
plot_contours_time_avg(pickle_file, proc_data_loc, qoi_plot_avg_map, qoi_range_map, [250, 350], [250, 350])
plot_contours_time_avg(pickle_file, proc_data_loc, qoi_plot_avg_map, qoi_range_map)
'''
# In[]
# 3D iso-surfaces Computation

threedim_qoi = ['W_destag','Vorticity','Qcriterion']
t_start = wrf_domain6['Time'][1]
t_end = wrf_domain6['Time'][-1]
t_stamp = int(t_end)
n_time = np.size(np.array(wrf_domain6['Time']))
switch = 0

pickle_filename_3d = "Data_3d.pkl"
pickle_file_3d = os.path.join(proc_data_loc, pickle_filename_3d)

# In[]
pickled_data_3d = create_data_3d(wrf_domain6, t_start, t_end, t_stamp, threedim_qoi)
with open(pickle_file_3d, 'wb') as pickle_file_handle:
    pickle.dump(pickled_data_3d, pickle_file_handle)

# In[]
#Plot 3d iso-surfaces
iso_values = list(np.linspace(0.005,0.009,4))
plot_isosurfaces_instant(pickle_file_3d, threedim_qoi[2], t_stamp, iso_values, z_plane_locs)
#plot_isosurfaces_animate(pickle_file_3d, threedim_qoi[2], t_start, t_end, iso_values, z_plane_locs)


# In[]
'''
# Loop over QoI of interest to create plots
for (qoi, qoi_unit) in zip(qoi_list, qoi_units):
    qoi_data, qoi_dim_names, qoi_dim = getQoI(wrf_domain6, qoi)

    # In[]    
    
    # Plot contours at the desired times and spatial locations
    slt = sampling_loc_time(qoi_dim_names, qoi_dim)
    for plane_ind in range(0,1):
        pl = plane[plane_ind]
        if pl=='xy':
            plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[2]])
            contourPlotSpaceTime(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[1],slt,plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
        elif pl=='yz':
            plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[2]], wrf_domain6[qoi_dim_names[1]])
            contourPlotSpaceTime(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[3],slt,plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
        elif pl=='xz':
            plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[1]])
            contourPlotSpaceTime(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[2],slt,plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
    '''
    
    # In[]  
    '''
    pl = 'xy'
    vert_loc = [20, 25]
    plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[2]])
    contourPlotInstantaneous(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[1], vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
    '''
    
    # In[]  
    '''
    pl = 'xy'
    vert_loc = [20]
    plane_cols, plane_rows = np.meshgrid(wrf_domain6[qoi_dim_names[3]], wrf_domain6[qoi_dim_names[2]])
    contourPlotTimeAvg(wrf_domain6, qoi_data, qoi, qoi_unit, qoi_dim_names[1], vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt)
    '''
    
    # In[]
    '''
    [we_ind, sn_ind] = [300, 300]
    linePlotTimeAvg(wrf_domain6, qoi_data, qoi, qoi_unit, we_ind, sn_ind, plot_loc, ref_time, dt)
    '''