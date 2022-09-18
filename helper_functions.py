#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 03:13:43 2022

@author: jha3
"""
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import json

from datetime import date, datetime, timedelta, time
# In[7]:


def getQoI (netCDF_data, qoi):
    qoi_data = netCDF_data[qoi]
    qoi_dim_names = qoi_data.dims
    qoi_dim = qoi_data.shape
    
    #shape format in : [time, z, y, x]
    print('QoI: {}, \nDimension names: {}, \nDimensions: {} '.format(qoi, qoi_dim_names, qoi_dim))
    return qoi_data, qoi_dim_names, qoi_dim
    


# ### Provide desired time and location for contour plots

# In[8]:


def sampling_loc_time(qoi_dim_names, qoi_dim):
    [time_dim, bottom_top_dim, south_north_dim, west_east_dim] = qoi_dim
    
    desired_time = [0, time_dim-1]
    #desired_time = list(np.arange(0, time_dim, 1))
    #hub_height = int(wrf_domain6['HUB_HEIGHT'].isel(Time=0)[0])
    hub_height = 100
    #bottom_top = [0.05*hub_height, 0.1*hub_height, 0.5*hub_height,hub_height] #bottom_top_locs
    bottom_top = [15,20,25]
    south_north = [int(west_east_dim*0.5), int(west_east_dim*0.75)] #south_north_locs
    west_east = [int(bottom_top_dim*0.25), int(south_north_dim*0.5), int(south_north_dim*0.75)] #west_east_locs

    #sampling location dictionary
    sampling_loc_time = '{"%s" : %s, "%s" : %s, "%s" : %s, "%s" : %s}'%(qoi_dim_names[0],desired_time, qoi_dim_names[1],bottom_top, qoi_dim_names[2], south_north, qoi_dim_names[3],west_east)
    slt = json.loads(sampling_loc_time)
    
    return slt


# ### Plot Contours

# In[9]:


def contourPlotSpaceTime(nc_data, qoi_data, qoi, qoi_unit, plane_header, slt, plane_rows, plane_cols, pl, plot_loc, ref_time, dt):
    desired_time = slt['Time']
    print('Desired Time : {}'.format(desired_time))
    desired_space = slt[plane_header]
    print('Desired Space ({}) : {}'.format(plane_header, desired_space))
    
    DX = nc_data.DX
    DY = nc_data.DY
    
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
   
    nSpace = np.size(desired_space)    
    nTime = np.size(desired_time) 
   
    qoi_time = []
    for time_count in range(nTime):
        qoi_space = []
        for space_count in range(nSpace):
            slice_loc = '{"%s" : %d}'%(plane_header, desired_space[space_count])
            slice_loc_dict = json.loads(slice_loc)
            arr = qoi_data.isel(Time=desired_time[time_count]).isel(slice_loc_dict)
            axes_name = arr.dims
            qoi_space.append(np.array(arr))
                    
        qoi_time.append(qoi_space)
       
    #print('qoi_time ', qoi_time)
    
        
    cmap_name = 'rainbow'
     
    if nTime==1:
        fig = plt.figure(figsize = (20, 4))
        fig.subplots_adjust(hspace=0.2, wspace=0.25)
        #fig.figure(figsize = (20, 4))
        time_inst_space = qoi_time[0]
        current_time = start_time + timedelta(seconds = desired_time[0]*dt)
        current_time_stamp = current_time.isoformat('_').split('_')[1]
        
        for space_count in range(nSpace):
            space = time_inst_space[space_count]
            nx = int(round(np.shape(space)[0],-2))
            ny = int(round(np.shape(space)[1],-2))
            ratio = ny/nx
            z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = desired_time[0]).isel(bottom_top_stag = desired_space[space_count])
            z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = desired_time[0]).isel(bottom_top_stag = desired_space[space_count] + 1)
            z = 0.5*(z1 + z2) 
            #print(ratio)
            ax = fig.add_subplot(nTime, nSpace, space_count+1)
            fig.set_figheight(4*nTime)
            fig.set_figwidth(5*nSpace*ratio)
            cont = ax.contourf(plane_cols*DX,plane_rows*DY, space, 20, cmap=cmap_name)
            clb = fig.colorbar(cont, ax=ax)
            clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
            ax.set_xlabel(f"%s [m]"%axes_name[1], fontsize=14)
            ax.set_ylabel(f"%s [m]"%axes_name[0], fontsize=14)
            ax.tick_params(axis='x', labelsize=14)
            ax.tick_params(axis='y', labelsize=14)
            ax.set_title(f"%s, z = %.2f m" %(current_time_stamp, z), fontsize=14)
            ax.set_xlim([250*DX,350*DX])
            ax.set_ylim([250*DY,350*DY])
        fig.tight_layout()
    else:
        fig, ax = plt.subplots(nTime, nSpace, figsize = (20,4.5), sharex='col', sharey='row')#figsize = (20,4.5)
        fig.subplots_adjust(hspace=0.5, wspace=0.1)
    
        for time_count in range(nTime):
            time_inst_space = qoi_time[time_count]
            current_time = start_time + timedelta(seconds = desired_time[time_count]*dt)
            current_time_stamp = current_time.isoformat('_').split('_')[1]

            for space_count in range(nSpace):
                space = time_inst_space[space_count]
                nx = int(round(np.shape(space)[0],-2))
                ny = int(round(np.shape(space)[1],-2))
                ratio = ny/nx
                z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = desired_time[time_count]).isel(bottom_top_stag = desired_space[space_count])
                z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = desired_time[time_count]).isel(bottom_top_stag = desired_space[space_count] + 1)
                z = 0.5*(z1 + z2) 
                #print(ratio)
                cont = ax[time_count,space_count].contourf(plane_cols*DX,plane_rows*DY, space, 20, cmap=cmap_name)
                clb = fig.colorbar(cont, ax=ax[time_count,space_count])
                clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
                fig.set_figheight(4.5*nTime)#figsize = (4.5*nTime)
                fig.set_figwidth(5*nSpace*ratio)
                ax[time_count,space_count].set_xlabel(f"%s [m]"%axes_name[1], fontsize=14)
                ax[time_count,space_count].set_ylabel(f"%s [m]"%axes_name[0], fontsize=14)
                ax[time_count,space_count].tick_params(axis='x', labelsize=14)
                ax[time_count,space_count].tick_params(axis='y', labelsize=14)
                ax[time_count,space_count].set_title(f"%s, z = %.2f m" %(current_time_stamp, z), fontsize=14)
                ax[time_count,space_count].set_xlim([250*DX,350*DX])
                ax[time_count,space_count].set_ylim([250*DY,350*DY])
        fig.tight_layout()
    
    filename = "%s_plane_%s"%(pl, qoi)
    plt.savefig(os.path.join(plot_loc, filename), bbox_inches='tight')       
            
    
# In[]
def contourPlotInstantaneous(nc_data, qoi_data, qoi, qoi_unit, plane_header, vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt):   
    DX = nc_data.DX
    DY = nc_data.DY
    
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
   
    nSpace = np.size(vert_loc)    
    nTime = nc_data.dims['Time']
   
    qoi_time = []
    for time_count in range(int(nTime)):
        qoi_space = []
        for space_count in range(nSpace):
            slice_loc = '{"%s" : %d}'%(plane_header, vert_loc[space_count])
            slice_loc_dict = json.loads(slice_loc)
            arr = qoi_data.isel(Time=time_count).isel(slice_loc_dict)
            axes_name = arr.dims
            qoi_space.append(np.array(arr))
                    
        qoi_time.append(qoi_space)
       
    #print('qoi_time ', qoi_time)
    
    cmap_name = 'rainbow'
     
    for time_count in range(int(nTime)):
        time_inst_space = qoi_time[time_count]
        current_time = start_time + timedelta(seconds = time_count*dt)
        current_time_stamp = current_time.isoformat('_').split('_')[1]

        for space_count in range(nSpace):
            space = time_inst_space[space_count]
            nx = int(round(np.shape(space)[0],-2))
            ny = int(round(np.shape(space)[1],-2))
            ratio = ny/nx
            z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = vert_loc[space_count])
            z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = vert_loc[space_count] + 1)
            z = 0.5*(z1 + z2) 
            #print(ratio)
            
            plt.figure()
            cont = plt.contourf(plane_cols*DX,plane_rows*DY, space, 20, cmap=cmap_name)
            clb = plt.colorbar(cont)
            clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
            #plt.set_figheight(4.5*nTime)#figsize = (4.5*nTime)
            #plt.set_figwidth(5*nSpace*ratio)
            plt.xlabel(f"%s [m]"%axes_name[1], fontsize=14)
            plt.ylabel(f"%s [m]"%axes_name[0], fontsize=14)
            plt.tick_params(axis='x', labelsize=14)
            plt.tick_params(axis='y', labelsize=14)
            plt.title(f"%s, z = %.2f m" %(current_time_stamp, z), fontsize=14)
            plt.xlim([250*DX,350*DX])
            plt.ylim([250*DY,350*DY])
            plt.tight_layout() 
   
    
            filename = "%s_plane_%s_%s_%s_%s"%(pl, qoi, current_time_stamp, plane_header, str(vert_loc[space_count]))
            plt.savefig(os.path.join(plot_loc, filename), bbox_inches='tight')       
            
# In[]
def contourPlotTimeAvg(nc_data, qoi_data, qoi, qoi_unit, plane_header, vert_loc, plane_rows, plane_cols, pl, plot_loc, ref_time, dt):   
    DX = nc_data.DX
    DY = nc_data.DY
    
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
    start_time_stamp = start_time.isoformat('_').split('_')[1]
    
    nSpace = np.size(vert_loc)
    
    cmap_name = 'rainbow'

    time_count = 0
    for space_count in range(nSpace):
        space = qoi_data.isel(bottom_top = vert_loc[space_count]).mean(dim = 'Time')
        axes_name = space.dims
        nx = int(round(np.shape(space)[0],-2))
        ny = int(round(np.shape(space)[1],-2))
        ratio = ny/nx
        z1 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = vert_loc[space_count])
        z2 = nc_data['ZTS'].isel(south_north = 300).isel(west_east = 300).isel(Time = time_count).isel(bottom_top_stag = vert_loc[space_count] + 1)
        z = 0.5*(z1 + z2) 
        #print(ratio)
        
        plt.figure()
        cont = plt.contourf(plane_cols*DX,plane_rows*DY, space, 20, cmap=cmap_name)
        clb = plt.colorbar(cont)
        clb.ax.set_title(f"%s [%s]"%(qoi, qoi_unit), weight='bold')
        #plt.set_figheight(4.5*nTime)#figsize = (4.5*nTime)
        #plt.set_figwidth(5*nSpace*ratio)
        plt.xlabel(f"%s [m]"%axes_name[1], fontsize=14)
        plt.ylabel(f"%s [m]"%axes_name[0], fontsize=14)
        plt.tick_params(axis='x', labelsize=14)
        plt.tick_params(axis='y', labelsize=14)
        plt.title(f"time avg start = %s, z = %.2f m" %(start_time_stamp, z), fontsize=14)
        plt.xlim([250*DX,350*DX])
        plt.ylim([250*DY,350*DY])
        plt.tight_layout() 
   

        filename = "%s_plane_%s_time_avg_start_%s_%s_%s"%(pl, qoi, start_time_stamp, plane_header, str(vert_loc[space_count]))
        plt.savefig(os.path.join(plot_loc, filename), bbox_inches='tight') 
        
 
# In[]
def linePlotTimeAvg(nc_data, qoi_data, qoi, qoi_unit, we_ind, sn_ind, plot_loc, ref_time, dt):
    start_time = datetime.fromisoformat(ref_time[2]+ '_' + ref_time[3]) - timedelta(seconds = dt)
    start_time_stamp = start_time.isoformat('_').split('_')[1]
    
    time_count = 0
    z_stag = nc_data['ZTS'].isel(south_north = sn_ind).isel(west_east = we_ind).isel(Time = time_count)
    z = 0.5*(z_stag[1:] + z_stag[:-1])
    
    line_var = qoi_data.isel(south_north = sn_ind).isel(west_east = we_ind).mean(dim = 'Time')
    
    plt.figure()
    plt.plot(line_var, z)
    plt.xlabel(f"%s [%s]"%(qoi, qoi_unit), fontsize=14)
    plt.ylabel(f"Height [m]"%z, fontsize=14)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.title(f"time avg start = %s, we = %d, sn = %d m" %(start_time_stamp, we_ind, sn_ind), fontsize=14)
    #plt.xlim([250*DX,350*DX])
    plt.ylim([0, 150])
    plt.tight_layout() 


    filename = "vertline_%s_we_%s_sn_%s_time_avg_start_%s"%(qoi, str(we_ind), str(sn_ind),  start_time_stamp)
    plt.savefig(os.path.join(plot_loc, filename), bbox_inches='tight')
    
    