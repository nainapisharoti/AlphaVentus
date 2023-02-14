# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 19:18:45 2023

@author: naina
"""

# In[]
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import json
from datetime import date, datetime, timedelta, time
import pickle
from mayavi import mlab
import wrf

# In[]
def plot_isosurfaces_instant(pickle_file_name, qoi, tstamp, a, zloc):

    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)

    # In[]    
    QOI = pickled_data_read[qoi]
    DX = pickled_data_read['dX']
    DY = pickled_data_read['dY']
    scalars = QOI[tstamp] #QOI at desired time step
    z = int(zloc[0])
    
    #mlab.options.offscreen = True
    
    src = mlab.pipeline.scalar_field(scalars)
    obj = mlab.pipeline.iso_surface(src, colormap = 'rainbow', contours=a,extent=[z-20,z+20,0,DX-1,0,DY-1], opacity=0.5, vmin=0.005,vmax=0.01)
    clb = mlab.colorbar(obj)
    clb.number_of_colors = 20
    mlab.view(0,270, distance=200) #(180,270: view from top)
    mlab.orientation_axes()
    mlab.savefig(f"%s.png"%(qoi))

# In[]
def plot_isosurfaces_animate(pickle_file_name, qoi, tstart, tend, a, zloc):
   
    with open(pickle_file_name, 'rb') as pickle_file_handle:
        pickled_data_read = pickle.load(pickle_file_handle)
    
    # In[]    
    QOI = pickled_data_read[qoi]
    DX = pickled_data_read['dX'] #number of grid points in x-direction
    DY = pickled_data_read['dY'] #number of grid points in y-direction
    scalars = QOI[1] #QOI at first time step
    z = int(zloc[0])
    t = np.linspace(int(tstart), int(tend), int(tend))
    
    #Iso-surface at first time step
    src = mlab.pipeline.scalar_field(scalars)
    obj = mlab.pipeline.iso_surface(src, colormap = 'rainbow', contours=a,extent=[z-20,z+20,0,DX-1,0,DY-1], opacity=0.5, vmin=0.005,vmax=0.01)
    clb = mlab.colorbar(obj)
    clb.number_of_colors = 20
    mlab.view(0,270, distance=200) #(180,270: view from top)
    mlab.orientation_axes()
    mlab.savefig(f"%s.png"%(qoi))

    for i in range(np.size(t)):
        T= int(t[i])
        print(T)
        obj.mlab_source.scalars = QOI[T]
        mlab.title(f"%s: time = %d"%(qoi,T), height=0.9,line_width =2, size = 0.5)
        filename = f"%s_%s.png"%(qoi,T)
        mlab.savefig(filename)
        mlab.clf

   