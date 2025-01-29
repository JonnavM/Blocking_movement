#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:47:54 2022

@author: mourik
Figures based on the output of the Celltrack Algorithm, after filtering for size and duration
"""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
from datetime import datetime
from matplotlib.colors import ListedColormap
import pandas as pd
from scipy.stats import gaussian_kde
import seaborn as sns
from matplotlib.legend_handler import HandlerTuple
import matplotlib.cm as cm

print("Modules loaded in")

dir = "Fill in directory"
fontsize=14

#Create personalized colormaps
cmap_1 = plt.cm.get_cmap('bwr', 11)
cmap_2 = plt.cm.get_cmap("YlOrRd", 10)

newcolors = cmap_2(np.linspace(0, 1, 10))
white = np.array([0, 0, 0, 0])
newcolors[:1, :] = white
cmap_2 = ListedColormap(newcolors)
print("cmaps created")

#Define function which approximates the area depending on the centre of the blocking
def gridcell(n, lat):
    R = 6371 #km
    dy = 2*np.pi*R*2.5/360
    dx = 2*np.pi*R*np.cos(np.deg2rad(lat))*2.5/360
    return n*dx*dy

def speedx(clcmassX, clcmassY):
    R = 6371 #km
    dx = clcmassX.ffill(dim="day").isel(day=-1)-clcmassX.sel(day=0)
    dx = dx.where(dx.values<180, dx.values-360)
    dx = dx.where(dx.values>-180, dx.values+360)
    dxkm = dx*2*np.pi*R*np.cos(np.deg2rad(clcmassY.sel(day=0)))*2.5/360
    v = dxkm/clcmassX.count("day")
    return v
    
print("Definitions defined")
#%% Open ERA to filter out blocks of short duration and small size. Only has to be run once. Saves files that will be loaded in later
"""
with open("/Tracking/output_data/ERA_varlat/celltrack/ERA_BI_1950_2022/tracks_all_stats.txt", "r") as f:
    ERA_area = [[num for num in line.split()] for line in f]
with open("/Tracking/output_data/ERA_varlat/celltrack/ERA_BI_1950_2022/cell_stats.txt", "r") as g:
    ERA_time = [[num for num in line.split()] for line in g]

ERA_time = np.array(ERA_time, dtype=object)[1:,:]
ERA_time[:,0]=ERA_time[:,0].astype(float)
ERA_track = xr.open_dataarray("/Blockings/per_lon/Intensity/long/block_ERA_BI_long.nc")
ERA_track_all = xr.open_dataset("/Tracking/output_data/ERA_varlat/celltrack/ERA_BI_1950_2022/tracks.nc").trackID
lon = ERA_track['lon'].data
lat = ERA_track['lat'].data

# Calculate the mean area per blocking
area = []
subarea = []
for i in range(len(ERA_area)):
    if len(ERA_area[i]) == 4: #Split when a new block is found
        block=int(ERA_area[i][1])
        area.append(subarea)
        subarea = []
    else:
        subarea.append(np.hstack([block, np.array(ERA_area[i], dtype=object)]))

long_area = []
for i in range(1, len(area)):        
    area[i] = np.vstack(area[i][1:])
    area[i] = area[i][:,:-1].astype(float)
    if len(area[i])>=4:
        long_area.append(area[i])
    else:
        continue

#Couple blocks to date of the cells
time = []
time_start = []
grid_clarea = []
clcmassX = []
clcmassY = []
BIav = []
BImax = []
a=0 #Used such that k won't start at zero everytime, but moves with the found cells
for i in range(len(long_area)):
    time_block = []
    for j in range(len(long_area[i][:,0])):
        for k in range(a, len(ERA_time[:,0])):
            if long_area[i][j,1]==ERA_time[k,0]:               
                time_block.append(datetime.strptime(ERA_time[k,15], "%Y%m%d").strftime("%Y-%m-%d"))
                a = int(long_area[i][0,1])
                break
            else:
                continue
    time.append(np.array(time_block[:], dtype=object))
    time_start.append(time_block[0])
    grid_clarea.append(long_area[i][:,4])
    clcmassX.append(long_area[i][:,7]) #was 5
    clcmassY.append(long_area[i][:,8]) #was 6
    BIav.append(long_area[i][:,14])
    BImax.append(long_area[i][:,13])
   
#Count number of blockings starting at the same day: 
count_startdate = {i:time_start.count(i) for i in time_start}    

# Set to xarray
data = np.zeros([5, len(count_startdate), len(max(long_area,key=len)), max(count_startdate.values())], dtype=float)*np.nan

c=0
for i in range(len(count_startdate)):
    for j in range(6):
        if time_start[c+j]==list(count_startdate.keys())[i]:
            for k in range(len(long_area[c+j])):
                data[0,i,k,j]=grid_clarea[c+j][k]
                data[1,i,k,j]=clcmassX[c+j][k]
                data[2,i,k,j]=clcmassY[c+j][k]
                data[3,i,k,j]=BIav[c+j][k]
                data[4,i,k,j]=BImax[c+j][k]
            if c+j==len(time_start)-1:
                break
        else:
            c+=j
            break
    
ERA_block = xr.Dataset(data_vars=dict(grid_clarea=(["time", "day", "block"], data[0,:,:,:]),
                                      clcmassX=(["time", "day", "block"], data[1,:,:,:]),
                                      clcmassY=(["time", "day", "block"], data[2,:,:,:]),
                                      BIav=(["time", "day", "block"], data[3,:,:,:]),
                                      BImax=(["time", "day", "block"], data[4,:,:,:])),
                       coords=dict(time=np.array(list(count_startdate.keys()), dtype="datetime64")),
                       attrs=dict(decription="Data on the area and location of blockings"))
ERA_block.to_netcdf("/Blockings/per_lon/Intensity/long/ERA_BI_data_weighted.nc")
#Mean area
print("Mean size of ERA blocking:", ERA_block["grid_clarea"].sel(time=slice("1951-01-01", "2014-12-31")).mean(), "+-", ERA_block["grid_clarea"].sel(time=slice("1951-01-01", "2014-12-31")).std(), "cells")
print("Mean size of ERA blocking:", gridcell(ERA_block["grid_clarea"].sel(time=slice("1951-01-01", "2014-12-31")), ERA_block["clcmassY"].sel(time=slice("1951-01-01", "2014-12-31"))).mean(), "+-", 
      gridcell(ERA_block["grid_clarea"].sel(time=slice("1951-01-01", "2014-12-31")), ERA_block["clcmassY"].sel(time=slice("1951-01-01", "2014-12-31"))).std())

#Delete ERA_time and ERA_area as they are no longer used and take up space
del(ERA_time)
del(ERA_area)

#%% Same for historical dataset
files = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]

hist_area = []
hist_time = []
hist_shape = []
dates = np.arange(19510101, 20141231, 1)

for n in range(16):
    print(n)
    with open("/Tracking/output_data/historical_varlat/celltrack/hist"+str(files[n])+"_BI_1850_2014/cell_shape.txt", "r") as h:
        d = [[num for num in line.split()] for line in h]
        hist_shape.append(np.array(d[1:]))
    with open("/Tracking/output_data/historical_varlat/celltrack/hist"+str(files[n])+"_BI_1850_2014/tracks_all_stats.txt", "r") as f:
        b = [[num for num in line.split()] for line in f]
        hist_area.append(np.array(b, dtype=object))
    with open("/Tracking/output_data/historical_varlat/celltrack/hist"+str(files[n])+"_BI_1850_2014/cell_stats.txt", "r") as g:
        c = [[num for num in line.split()] for line in g]
        append = np.array(c, dtype=object)[1:,:]
        append[:,0]=append[:,0].astype(float)
        hist_time.append(append)


long_area_hist = []
long_area_hist_16 = []
time_hist_16 = []
shape_hist_16 = []
time_start_hist_16 = []
grid_clarea_hist_16 = []
clcmassX_hist_16 = []
clcmassY_hist_16 = []
BIav_hist_16 = []
BImax_hist_16 = []
count_startdate_hist = [] 
max_count_hist = []

for n in range(16): 
    long_area_hist_prev = []
    print(n)
    # Calculate the mean area per blocking
    area2 = []
    subarea2 = []
    for i in range(len(hist_area[n])):
        if len(hist_area[n][i]) == 4:
            area2.append(subarea2)
            subarea2 = []
        else:
            subarea2.append(np.array(hist_area[n][i], dtype=object))      
    for i in range(1, len(area2)):
        area2[i] = np.vstack(area2[i][1:])
        area2[i] = area2[i][:,:-1].astype(float)
        if len(area2[i])>=4 and np.max(gridcell(area2[i][:,3], area2[i][:,5]))>=500000:
         #   long_area_hist.append(area2[i])
            long_area_hist_prev.append(area2[i])
        else:
            continue
    
    #Couple blocks to date of the cells
    time_hist = []
    shape_hist = []
    time_start_hist = []
    grid_clarea_hist = []
    clcmassX_hist = []
    clcmassY_hist = []
    BIav_hist = []
    BImax_hist = []
    
    a=0 #Used such that k won't start at zero everytime, but moves with the found cells
    for i in range(len(long_area_hist_prev)):
        time_block_hist = []
        shape_block_hist = []
        for j in range(len(long_area_hist_prev[i][:,0])):
            for k in range(a, len(hist_time[n][:,0])):
                if long_area_hist_prev[i][j,0]==hist_time[n][k,0]:
                    time_block_hist.append(datetime.strptime(hist_time[n][k,15], "%Y%m%d").strftime("%Y-%m-%d"))
                    shape_block_hist.append([hist_shape[n][k].astype(float), datetime.strptime(hist_time[n][k,15], "%Y%m%d").strftime("%Y-%m-%d")])
                    a = int(long_area_hist_prev[i][0,0])
                    break
                else:
                    continue
        time_hist.append(np.array(time_block_hist[:], dtype=object))
        shape_hist.append(shape_block_hist[:])
        time_start_hist.append(time_block_hist[0])
#        grid_clarea_hist.append(long_area_hist_prev[i][:,3])
        grid_clarea_hist.append(gridcell(long_area_hist_prev[i][:,3], long_area_hist_prev[i][:,5]))
        clcmassX_hist.append(long_area_hist_prev[i][:,6]) #was 4
        clcmassY_hist.append(long_area_hist_prev[i][:,7]) #was 5 
        BImax_hist.append(long_area_hist_prev[i][:,12])
        BIav_hist.append(long_area_hist_prev[i][:,13])
        long_area_hist.append(long_area_hist_prev[i])
 
    long_area_hist_16.append(long_area_hist_prev)
    time_hist_16.append(time_hist)
    shape_hist_16.append(shape_hist)
    time_start_hist_16.append(time_start_hist)
    grid_clarea_hist_16.append(grid_clarea_hist)
    clcmassX_hist_16.append(clcmassX_hist)
    clcmassY_hist_16.append(clcmassY_hist)
    BImax_hist_16.append(BImax_hist)
    BIav_hist_16.append(BIav_hist)
    count_startdate_hist.append({l:time_start_hist.count(l) for l in time_start_hist}) 
    max_count_hist.append(max(count_startdate_hist[n].values()))

long_area_hist_stack = np.vstack(long_area_hist)

#Delete hist_area and hist_time
del(hist_area)
del(hist_time)
del(hist_shape)

hist_clarea_list = []
hist_clcmassX_list = []
hist_clcmassY_list = []
hist_BIav_list = []
hist_BImax_list = []
hist_majlen_list = []
hist_minlen_list = []

for n in range(16):
    clarea_histx = np.zeros([len(count_startdate_hist[n]), len(max(long_area_hist, key=len)), np.max(max_count_hist)], dtype=float)*np.nan #18450,38,6
    clcmassX_histx = np.zeros([len(count_startdate_hist[n]), len(max(long_area_hist, key=len)), np.max(max_count_hist)], dtype=float)*np.nan
    clcmassY_histx = np.zeros([len(count_startdate_hist[n]), len(max(long_area_hist, key=len)), np.max(max_count_hist)], dtype=float)*np.nan
    BIav_histx = np.zeros([len(count_startdate_hist[n]), len(max(long_area_hist, key=len)), np.max(max_count_hist)], dtype=float)*np.nan
    BImax_histx = np.zeros([len(count_startdate_hist[n]), len(max(long_area_hist, key=len)), np.max(max_count_hist)], dtype=float)*np.nan
    majlen_histx = np.zeros([len(count_startdate_hist[n]), len(max(long_area_hist, key=len)), np.max(max_count_hist)], dtype=float)*np.nan
    minlen_histx = np.zeros([len(count_startdate_hist[n]), len(max(long_area_hist, key=len)), np.max(max_count_hist)], dtype=float)*np.nan

    c=0
    for i in range(len(count_startdate_hist[n])):
        for j in range(max(count_startdate_hist[n].values())+1):
            if time_start_hist_16[n][c+j]==list(count_startdate_hist[n].keys())[i]:
                for k in range(len(long_area_hist_16[n][c+j])):
                    clarea_histx[i,k,j]=grid_clarea_hist_16[n][c+j][k]
                    clcmassX_histx[i,k,j]=clcmassX_hist_16[n][c+j][k]
                    clcmassY_histx[i,k,j]=clcmassY_hist_16[n][c+j][k]
                    BIav_histx[i,k,j]=BIav_hist_16[n][c+j][k]
                    BImax_histx[i,k,j]=BImax_hist_16[n][c+j][k]
                    majlen_histx[i,k,j]=shape_hist_16[n][c+j][k][0][1]
                    minlen_histx[i,k,j]=shape_hist_16[n][c+j][k][0][2]
                if c+j==len(time_start_hist_16[n])-1:
                    break
            else:
                c+=j
                break
    hist_clarea_list.append(xr.Dataset(data_vars=dict(grid_clarea=(["time", "day", "block"], clarea_histx[:,:,:])),
                                   coords=dict(time=list(count_startdate_hist[n].keys())),
                                   attrs=dict(description="Clarea data for ensemble "+str(n)+".")))
    hist_clcmassX_list.append(xr.Dataset(data_vars=dict(clcmassX=(["time", "day", "block"], clcmassX_histx[:,:,:])),
                                   coords=dict(time=list(count_startdate_hist[n].keys())),
                                   attrs=dict(description="ClcmassX data for ensemble "+str(n)+".")))
    hist_clcmassY_list.append(xr.Dataset(data_vars=dict(clcmassY=(["time", "day", "block"], clcmassY_histx[:,:,:])),
                                   coords=dict(time=list(count_startdate_hist[n].keys())),
                                   attrs=dict(description="ClcmassY data for ensemble "+str(n)+".")))
    hist_BIav_list.append(xr.Dataset(data_vars=dict(BIav=(["time", "day", "block"], BIav_histx[:,:,:])),
                                     coords=dict(time=list(count_startdate_hist[n].keys())),
                                     attrs=dict(description="BIav data for ensemble "+str(n)+".")))
    hist_BImax_list.append(xr.Dataset(data_vars=dict(BImax=(["time", "day", "block"], BImax_histx[:,:,:])),
                                     coords=dict(time=list(count_startdate_hist[n].keys())),
                                     attrs=dict(description="BImax data for ensemble "+str(n)+".")))
    hist_majlen_list.append(xr.Dataset(data_vars=dict(majlen=(["time", "day", "block"], majlen_histx[:,:,:])),
                                       coords=dict(time=list(count_startdate_hist[n].keys())),
                                       attrs=dict(description="Major length data for ensemble "+str(n)+".")))
    hist_minlen_list.append(xr.Dataset(data_vars=dict(minlen=(["time", "day", "block"], minlen_histx[:,:,:])),
                                       coords=dict(time=list(count_startdate_hist[n].keys())),
                                       attrs=dict(description="Minor length data for ensemble "+str(n)+".")))
    
hist_clarea = xr.concat([hist_clarea_list[i] for i in range(16)], "ensemble") 
hist_clarea['time'] = pd.DatetimeIndex(hist_clarea['time'].values) 
hist_clcmassX = xr.concat([hist_clcmassX_list[i] for i in range(16)], "ensemble")
hist_clcmassX['time'] = pd.DatetimeIndex(hist_clcmassX['time'].values) 
hist_clcmassY = xr.concat([hist_clcmassY_list[i] for i in range(16)], "ensemble")    
hist_clcmassY['time'] = pd.DatetimeIndex(hist_clcmassY['time'].values)     
hist_BIav= xr.concat([hist_BIav_list[i] for i in range(16)], "ensemble")    
hist_BIav['time'] = pd.DatetimeIndex(hist_BIav['time'].values)   
hist_BImax = xr.concat([hist_BImax_list[i] for i in range(16)], "ensemble")    
hist_BImax['time'] = pd.DatetimeIndex(hist_BImax['time'].values)

hist_majlen = xr.concat([hist_majlen_list[i] for i in range(16)], "ensemble")
hist_majlen["time"] = pd.DatetimeIndex(hist_majlen["time"].values)
hist_minlen = xr.concat([hist_minlen_list[i] for i in range(16)], "ensemble")
hist_minlen["time"] = pd.DatetimeIndex(hist_minlen["time"].values)

hist_clarea.to_netcdf("/Blockings/per_lon/Intensity/long/large/hist_BI_clarea_large_weighted.nc")
hist_clcmassX.to_netcdf("/Blockings/per_lon/Intensity/long/large/hist_BI_clcmassX_large_weighted.nc")
hist_clcmassY.to_netcdf("/Blockings/per_lon/Intensity/long/large/hist_BI_clcmassY_large_weighted.nc")
hist_BIav.to_netcdf("/Blockings/per_lon/Intensity/long/large/hist_BI_BIav_large.nc")
hist_BImax.to_netcdf("/Blockings/per_lon/Intensity/long/large/hist_BI_BImax_large.nc")
hist_majlen.to_netcdf("/Blockings/per_lon/Intensity/long/large/hist_BI_majlen_large.nc")
hist_minlen.to_netcdf("/Blockings/per_lon/Intensity/long/large/hist_BI_minlen_large.nc")
"""
#%% Import data for ERA-5 and EC-Earth3bis
ERA_block = xr.open_dataset(dir+"Blockings/per_lon/Intensity/long/ERA_BI_data_weighted.nc")
ERA_speedx = speedx(ERA_block.clcmassX, ERA_block.clcmassY)
print("ERA5 data imported")

hist_clarea = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_clarea_large_weighted.nc")
hist_clcmassX = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_clcmassX_large_weighted.nc")
hist_clcmassY = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_clcmassY_large_weighted.nc")
hist_BIav = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_BIav_large.nc")
hist_BImax = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_BImax_large.nc")
hist_speedx = speedx(hist_clcmassX, hist_clcmassY)
hist_speedx = hist_speedx*1000/(24*3600)
print("ECE3p5 data imported")

#%% Blocking statistics:
##ERA5:
#Amount of blocks
print("Mean amount of blocks ERA5:", "{:.0f}".format(ERA_speedx.sel(time=slice("1951", "2014")).where(ERA_speedx.sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).mean().values), "+-", "{:.0f}".format(ERA_speedx.sel(time=slice("1951", "2014")).where(ERA_speedx.sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).std(dim=("year")).values))
print("Mean amount of blocks in winter ERA5:", "{:.0f}".format(ERA_speedx.sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951", "2014")).where(ERA_speedx.sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).mean().values), "+-", "{:.0f}".format(ERA_speedx.sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951", "2014")).where(ERA_speedx.sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).std(dim=("year")).values))
print("Mean amount of blocks in summer ERA5:", "{:.0f}".format(ERA_speedx.sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951", "2014")).where(ERA_speedx.sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).mean().values), "+-", "{:.0f}".format(ERA_speedx.sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951", "2014")).where(ERA_speedx.sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).std(dim=("year")).values))
#Size
print("Mean size ERA5:", "{:.0f}".format(gridcell(ERA_block["grid_clarea"].sel(time=slice("1951-01-01", "2014-12-31")).mean("day"), ERA_block["clcmassY"].sel(time=slice("1951-01-01", "2014-12-31")).mean("day")).mean().values), "+-", "{:.0f}".format(gridcell(ERA_block["grid_clarea"].sel(time=slice("1951-01-01", "2014-12-31")), ERA_block["clcmassY"].sel(time=slice("1951-01-01", "2014-12-31"))).std().values))
print("Mean size in winter ERA5:", "{:.0f}".format(gridcell(ERA_block["grid_clarea"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).mean("day"), ERA_block["clcmassY"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).mean("day")).mean().values), "+-", "{:.0f}".format(gridcell(ERA_block["grid_clarea"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")), ERA_block["clcmassY"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31"))).std().values))
print("Mean size in summer ERA5:", "{:.0f}".format(gridcell(ERA_block["grid_clarea"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).mean("day"), ERA_block["clcmassY"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).mean("day")).mean().values), "+-", "{:.0f}".format(gridcell(ERA_block["grid_clarea"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")), ERA_block["clcmassY"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31"))).std().values))
#Duration
print("Mean duration ERA5:", "{:.1f}".format(ERA_block["BIav"].sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(ERA_block["BIav"].sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).mean().values), "+-", "{:.1f}".format(ERA_block["BIav"].sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(ERA_block["BIav"].sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).std().values))
print("Mean duration in winter ERA5:", "{:.1f}".format(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).mean().values), "+-", "{:.1f}".format(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).std().values))
print("Mean duration in summer ERA5:", "{:.1f}".format(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).mean().values), "+-", "{:.1f}".format(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).std().values))
#BIav
print("Mean BIav ERA5:", "{:.1f}".format(ERA_block["BIav"].sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_block["BIav"].sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean BIav in winter ERA5:", "{:.1f}".format(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean BIav in summer ERA5:", "{:.1f}".format(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_block["BIav"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
#BImax
print("Mean BImax ERA5:", "{:.1f}".format(ERA_block["BImax"].sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_block["BImax"].sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean BImax in winter ERA5:", "{:.1f}".format(ERA_block["BImax"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_block["BImax"].sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean BImax in summer ERA5:", "{:.1f}".format(ERA_block["BImax"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_block["BImax"].sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
#Zonal speed
print("Mean v_x ERA5:", "{:.0f}".format(ERA_speedx.sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_speedx.sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean v_x in winter ERA5:", "{:.0f}".format(ERA_speedx.sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_speedx.sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean v_x in summer ERA5:", "{:.0f}".format(ERA_speedx.sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(ERA_speedx.sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
#Absolute velocity
print("Abs v_x ERA5:", "{:.0f}".format(abs(ERA_speedx.sel(time=slice("1951-01-01", "2014-12-31"))).mean().values), "+-", "{:.1f}".format(abs(ERA_speedx.sel(time=slice("1951-01-01", "2014-12-31"))).std().values))
print("Abs v_x in winter ERA5:", "{:.0f}".format(abs(ERA_speedx.sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31"))).mean().values), "+-", "{:.1f}".format(abs(ERA_speedx.sel(time=ERA_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31"))).std().values))
print("Abs v_x in summer ERA5:", "{:.0f}".format(abs(ERA_speedx.sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31"))).mean().values), "+-", "{:.1f}".format(abs(ERA_speedx.sel(time=ERA_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31"))).std().values))

##ECE3p5:
#Amount of blocks
print("Mean amount of blocks ECE3p5:", "{:.0f}".format(hist_speedx.sel(time=slice("1951", "2014")).where(hist_speedx.sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).mean("ensemble").mean().values), "+-", "{:.0f}".format(hist_speedx.sel(time=slice("1951", "2014")).where(hist_speedx.sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).std(dim=("ensemble", "year")).values))
print("Mean amount of blocks in winter ECE3p5:", "{:.0f}".format(hist_speedx.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951", "2014")).where(hist_speedx.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).mean("ensemble").mean().values), "+-", "{:.0f}".format(hist_speedx.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951", "2014")).where(hist_speedx.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).std(dim=("ensemble", "year")).values))
print("Mean amount of blocks in summer ECE3p5:", "{:.0f}".format(hist_speedx.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951", "2014")).where(hist_speedx.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).mean("ensemble").mean().values), "+-", "{:.0f}".format(hist_speedx.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951", "2014")).where(hist_speedx.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951", "2014")).isnull(), 1).groupby("time.year").sum(dim=("time", "block")).std(dim=("ensemble", "year")).values))
#Size
print("Mean size ECE3p5:", "{:.0f}".format(hist_clarea.sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.0f}".format(hist_clarea.sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean size in winter ECE3p5:", "{:.0f}".format(hist_clarea.sel(time=hist_clarea["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.0f}".format(hist_clarea.sel(time=hist_clarea["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean size in summer ECE3p5:", "{:.0f}".format(hist_clarea.sel(time=hist_clarea["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.0f}".format(hist_clarea.sel(time=hist_clarea["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
#Duration
print("Mean duration ECE3p5:", "{:.1f}".format(hist_BIav.sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(hist_BIav.sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).mean().values), "+-", "{:.1f}".format(hist_BIav.sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(hist_BIav.sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).std().values))
print("Mean duration in winter ECE3p5:", "{:.1f}".format(hist_BIav.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(hist_BIav.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).mean().values), "+-", "{:.1f}".format(hist_BIav.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(hist_BIav.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).std().values))
print("Mean duration in summer ECE3p5:", "{:.1f}".format(hist_BIav.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(hist_BIav.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).mean().values), "+-", "{:.1f}".format(hist_BIav.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).count("day").where(hist_BIav.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).count("day")!=0).std().values))
#BIav
print("Mean BIav ECE3p5:", "{:.1f}".format(hist_BIav.sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(hist_BIav.sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean BIav in winter ECE3p5:", "{:.1f}".format(hist_BIav.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(hist_BIav.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean BIav in summer ECE3p5:", "{:.1f}".format(hist_BIav.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(hist_BIav.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
#BImax
print("Mean BImax ECE3p5:", "{:.1f}".format(hist_BImax.sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(hist_BImax.sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean BImax in winter ECE3p5:", "{:.1f}".format(hist_BImax.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(hist_BImax.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
print("Mean BImax in summer ECE3p5:", "{:.1f}".format(hist_BImax.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).mean().values), "+-", "{:.1f}".format(hist_BImax.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31")).std().values))
#Zonal speed
print("Abs v_x ECE3p5:", "{:.0f}".format(abs(hist_speedx.sel(time=slice("1951-01-01", "2014-12-31"))).mean().values), "+-", "{:.1f}".format(abs(hist_speedx.sel(time=slice("1951-01-01", "2014-12-31"))).std().values))
print("Abs v_x in winter ECE3p5:", "{:.0f}".format(abs(hist_speedx.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31"))).mean().values), "+-", "{:.1f}".format(abs(hist_speedx.sel(time=hist_speedx["time.season"]=="DJF").sel(time=slice("1951-01-01", "2014-12-31"))).std().values))
print("Abs v_x in summer ECE3p5:", "{:.0f}".format(abs(hist_speedx.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31"))).mean().values), "+-", "{:.1f}".format(abs(hist_speedx.sel(time=hist_speedx["time.season"]=="JJA").sel(time=slice("1951-01-01", "2014-12-31"))).std().values))


#%% Spatial distribution per season
season = "DJF"
#Calculate lat and lon of blocking centers for winter
x = hist_clcmassX.sel(time=hist_clcmassX["time.season"]==season).where(hist_speedx>=hist_speedx.quantile(0.90)).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")
y = hist_clcmassY.sel(time=hist_clcmassY["time.season"]==season).where(hist_speedx>=hist_speedx.quantile(0.90)).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")

x1 = hist_clcmassX.sel(time=hist_clcmassX["time.season"]==season).where(hist_speedx<=hist_speedx.quantile(0.10)).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")
y1 = hist_clcmassY.sel(time=hist_clcmassY["time.season"]==season).where(hist_speedx<=hist_speedx.quantile(0.10)).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")

x2_cond = np.logical_and(hist_speedx<hist_speedx.quantile(0.90), hist_speedx>hist_speedx.quantile(0.10))
x2 = hist_clcmassX.sel(time=hist_clcmassX["time.season"]==season).where(x2_cond==True).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")
y2_cond = np.logical_and(hist_speedx<hist_speedx.quantile(0.90), hist_speedx>hist_speedx.quantile(0.10))
y2 = hist_clcmassY.sel(time=hist_clcmassY["time.season"]==season).where(y2_cond==True).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")

season = "JJA"
#Calculate lat and lon of blocking centers for summer
x3 = hist_clcmassX.sel(time=hist_clcmassX["time.season"]==season).where(hist_speedx>=hist_speedx.quantile(0.90)).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")
y3 = hist_clcmassY.sel(time=hist_clcmassY["time.season"]==season).where(hist_speedx>=hist_speedx.quantile(0.90)).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")

x4 = hist_clcmassX.sel(time=hist_clcmassX["time.season"]==season).where(hist_speedx<=hist_speedx.quantile(0.10)).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")
y4 = hist_clcmassY.sel(time=hist_clcmassY["time.season"]==season).where(hist_speedx<=hist_speedx.quantile(0.10)).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")

x5_cond = np.logical_and(hist_speedx<hist_speedx.quantile(0.90), hist_speedx>hist_speedx.quantile(0.10))
x5 = hist_clcmassX.sel(time=hist_clcmassX["time.season"]==season).where(x2_cond==True).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")
y5_cond = np.logical_and(hist_speedx<hist_speedx.quantile(0.90), hist_speedx>hist_speedx.quantile(0.10))
y5 = hist_clcmassY.sel(time=hist_clcmassY["time.season"]==season).where(y2_cond==True).sel(day=3).stack(z=("ensemble", "block", "time")).dropna("z")

#Plot the figure
fig = plt.figure()
#Determine place of colorbars
cbar_ax = fig.add_axes([0.12, 0.05, .24, .03])
cbar_ax2 = fig.add_axes([0.4, 0.05, 0.24, 0.03])
cbar_ax3 = fig.add_axes([0.67, 0.05, 0.24, 0.03])
#Add subplots
ax = fig.add_subplot(2,3,2, projection=ccrs.Orthographic(central_latitude=90))
fig.text(s="DJF", fontsize=fontsize, rotation="vertical", x=0.09, y=0.7)
ax.coastlines()
ax.set_title(r"$v_x \geq P90$", fontsize=fontsize)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
a = sns.kdeplot(ax=ax, x=x, y=y, cmap="Reds", shade=True, cbar=True, cbar_ax=cbar_ax2, cbar_kws={"orientation":"horizontal", "ticks":[0.6e-4, 1.2e-4, 1.8e-4, 2.4e-4], "label":"Frequency", "format":"%4.1e"}, transform=ccrs.PlateCarree(), vmax=0.00025, vmin=0.00006)#levels=10, thresh=0.00006, common_norm=True, bw_method=0.2, gridsize=100)

ax = fig.add_subplot(2,3,1, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines()
ax.set_title(r"$v_x \leq P10$", fontsize=fontsize)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
b = sns.kdeplot(ax=ax, x=x1, y=y1, cmap="Blues", shade=True, cbar=True, cbar_ax=cbar_ax, cbar_kws={"orientation":"horizontal", "ticks":[0.6e-4, 1.2e-4, 1.8e-4, 2.4e-4], "label":"Frequency", "format":"%4.1e"}, transform=ccrs.PlateCarree(), vmax=0.00025, vmin=0.00006)#levels=10, thresh=0.00006, common_norm=True, bw_method=0.2, gridsize=100)

ax = fig.add_subplot(2,3,3, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines()
ax.set_title(r"$P10<v_x<P90$", fontsize=fontsize)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
c = sns.kdeplot(ax=ax, x=x2, y=y2, cmap="Greens", shade=True, cbar=True, cbar_ax=cbar_ax3, cbar_kws={"orientation":"horizontal", "ticks":[0.6e-4, 1.2e-4, 1.8e-4, 2.4e-4], "label":"Frequency", "format":"%4.1e"}, transform=ccrs.PlateCarree(), vmax=0.00025, vmin=0.00006)#levels=10, thresh=0.00006, common_norm=True, bw_method=0.2, gridsize=100)

ax = fig.add_subplot(2,3,5, projection=ccrs.Orthographic(central_latitude=90))
fig.text(s="JJA", fontsize=fontsize, rotation="vertical", x=0.09, y=0.28)
ax.coastlines()
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
a = sns.kdeplot(ax=ax, x=x3, y=y3, cmap="Reds", shade=True, cbar=True, cbar_ax=cbar_ax2, cbar_kws={"orientation":"horizontal", "ticks":[0.6e-4, 1.2e-4, 1.8e-4, 2.4e-4], "label":"Frequency", "format":"%4.1e"}, transform=ccrs.PlateCarree(), vmax=0.00025, vmin=0.00006)#levels=10, thresh=0.00006, common_norm=True, bw_method=0.2, gridsize=100)

ax = fig.add_subplot(2,3,4, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines()
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
b = sns.kdeplot(ax=ax, x=x4, y=y4, cmap="Blues", shade=True, cbar=False, cbar_ax=cbar_ax, cbar_kws={"orientation":"horizontal", "ticks":[0.6e-4, 1.2e-4, 1.8e-4, 2.4e-4], "label":"Frequency", "format":"%4.1e"}, transform=ccrs.PlateCarree(), vmax=0.00025, vmin=0.00006)#levels=10, thresh=0.00006, common_norm=True, bw_method=0.2, gridsize=100)

ax = fig.add_subplot(2,3,6, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines()
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
c = sns.kdeplot(ax=ax, x=x5, y=y5, cmap="Greens", shade=True, cbar=True, cbar_ax=cbar_ax3, cbar_kws={"orientation":"horizontal", "ticks":[0.6e-4, 1.2e-4, 1.8e-4, 2.4e-4], "label":"Frequency", "format":"%4.1e"}, transform=ccrs.PlateCarree(), vmax=0.00025, vmin=0.00006)#levels=10, thresh=0.00006, common_norm=True, bw_method=0.2, gridsize=100)
fig.savefig(dir+"Figures_fontsize/spatialdistr_velocity_P10_P90_ng_day4_fs="+str(fontsize)+".pdf", bbox_inches='tight')
plt.show()

#%% Plot with contourf
from scipy import stats
from matplotlib.colors import ListedColormap

#Create colormaps
cmap_1 = plt.cm.get_cmap('Reds', 10)
cmap_2 = plt.cm.get_cmap("Blues", 10)
cmap_3 = plt.cm.get_cmap('Greens', 10)

newcolors1 = cmap_1(np.linspace(0, 1, 10))
newcolors2 = cmap_2(np.linspace(0, 1, 10))
newcolors3 = cmap_3(np.linspace(0, 1, 10))
white = np.array([0, 0, 0, 0])

newcolors1[:1, :] = white
newcolors2[:1, :] = white
newcolors3[:1, :] = white

cmap_1 = ListedColormap(newcolors1)
cmap_2 = ListedColormap(newcolors2)
cmap_3 = ListedColormap(newcolors3)

#Create KDE datasets
X, Y = np.mgrid[-180:180:144j, 0:90:72j]
positions = np.vstack([X.ravel(), Y.ravel()])

values = np.vstack([x, y])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)

values1 = np.vstack([x1, y1])
kernel1 = stats.gaussian_kde(values1)
Z1 = np.reshape(kernel1(positions).T, X.shape)

values2 = np.vstack([x2, y2])
kernel2 = stats.gaussian_kde(values2)
Z2 = np.reshape(kernel2(positions).T, X.shape)

values3 = np.vstack([x3, y3])
kernel3 = stats.gaussian_kde(values3)
Z3 = np.reshape(kernel3(positions).T, X.shape)

values4 = np.vstack([x4, y4])
kernel4 = stats.gaussian_kde(values4)
Z4 = np.reshape(kernel4(positions).T, X.shape)

values5 = np.vstack([x5, y5])
kernel5 = stats.gaussian_kde(values5)
Z5 = np.reshape(kernel5(positions).T, X.shape)

#Plot figure
fig = plt.figure()
cbar_ax = fig.add_axes([0.12, 0.05, .24, .03])
cbar_ax2 = fig.add_axes([0.4, 0.05, 0.24, 0.03])
cbar_ax3 = fig.add_axes([0.68, 0.05, 0.24, 0.03])

ax = fig.add_subplot(2,3,3, projection=ccrs.Orthographic(central_latitude=90))
fig.text(s="DJF", fontsize=fontsize, rotation="vertical", x=0.09, y=0.7)
ax.coastlines()
ax.set_title(r"$v_x \geq P90$", fontsize=fontsize)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
CS = ax.contourf(X, Y, Z, levels=np.linspace(0, 0.0003, 10), extend="max", cmap=cmap_1, transform=ccrs.PlateCarree())

ax = fig.add_subplot(2,3,1, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines()
ax.set_title(r"$v_x \leq P10$", fontsize=fontsize)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
CS = ax.contourf(X, Y, Z1, levels=np.linspace(0, 0.0003, 10), extend="max", cmap=cmap_2, transform=ccrs.PlateCarree())

ax = fig.add_subplot(2,3,2, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines()
ax.set_title(r"$P10<v_x<P90$", fontsize=fontsize)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
CS = ax.contourf(X, Y, Z2, levels=np.linspace(0, 0.0003, 10), extend="max", cmap=cmap_3, transform=ccrs.PlateCarree())

ax = fig.add_subplot(2,3,6, projection=ccrs.Orthographic(central_latitude=90))
fig.text(s="JJA", fontsize=fontsize, rotation="vertical", x=0.09, y=0.28)
ax.coastlines()
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
CS = ax.contourf(X, Y, Z3*10**4, levels=np.linspace(0, 0.0003e4, 10), extend="max", cmap=cmap_1, transform=ccrs.PlateCarree())
cbar = plt.colorbar(CS, ticks=[0, 1, 2, 3], orientation="horizontal", cax=cbar_ax3)
cbar.set_label(r"Frequency ($\cdot 10^{-4}$)", size=fontsize)

ax = fig.add_subplot(2,3,4, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines()
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
CS = ax.contourf(X, Y, Z4*10**4, levels=np.linspace(0, 0.0003e4, 10), extend="max", cmap=cmap_2, transform=ccrs.PlateCarree())
cbar = plt.colorbar(CS, ticks=[0, 1, 2, 3], orientation="horizontal", cax=cbar_ax)
cbar.set_label(r"Frequency ($\cdot 10^{-4}$)", size=fontsize)

ax = fig.add_subplot(2,3,5, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines()
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
CS = ax.contourf(X, Y, Z5*10**4, levels=np.linspace(0, 0.0003e4, 10), extend="max", cmap=cmap_3, transform=ccrs.PlateCarree())
cbar = plt.colorbar(CS, ticks=[0, 1, 2, 3], orientation="horizontal", cax=cbar_ax2)
cbar.set_label(r"Frequency ($\cdot 10^{-4}$)", size=fontsize)

fig.savefig(dir+"Figures_fontsize/spatialdistr_velocity_P10_P90_ng_day4_middle_fs="+str(fontsize)+".pdf", bbox_inches='tight')
plt.show()

#%% Seasonality of blocking characteristics
sizep = hist_clarea.where(hist_speedx>0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
sizen = hist_clarea.where(hist_speedx<=0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
v_xp = hist_speedx.where(hist_speedx>0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"))
v_xn = -hist_speedx.where(hist_speedx<=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"))
duration = hist_clarea.count("day")
durationp = duration.where(hist_speedx>0).where(duration!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"), skipna=True)
durationn = duration.where(duration!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"), skipna=True)
intensityp = hist_BIav.where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
intensityn = hist_BIav.where(hist_speedx<=0).where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))

index = [[sizep, sizen], [durationp, durationn], [intensityp, intensityn]]
ylabel = [r"Size (km$^2$)", "Duration (days)", "Average intensity BI (-)"]
name = ["Size", "Duration", "Intensity"]

for i in range(3):
    fig, ax = plt.subplots()
    index[i][0].plot(color="blue", alpha=0.5)
    a1, = index[i][0].rolling(dayofyear=10, center=True).mean().plot(color="blue")
    index[i][1].plot(color="dodgerblue", alpha=0.5)
    a2, = index[i][1].rolling(dayofyear=10, center=True).mean().plot(color="dodgerblue")
    ax.set_xlabel("Day of year", fontsize=fontsize)
    ax.set_ylabel(ylabel[i], color="blue")
    ax2 = ax.twinx()
    v_xp.plot(color="red", alpha=0.5)
    b1, = v_xp.rolling(dayofyear=10, center=True).mean().plot(color="red")
    v_xn.plot(color="orange", alpha=0.5)
    b2, = v_xn.rolling(dayofyear=10, center=True).mean().plot(color="orange")
    ax2.set_ylabel("Velocity (km/day)", color="red", fontsize=fontsize)
    ax2.set_ylim(150, 460)
    l = ax2.legend([(a1, b1), (a2, b2)], ["Eastward", "Westward"], handler_map={tuple: HandlerTuple(ndivide=None)}, loc="upper left", fontsize=fontsize)
    #Show seasons:
    ax.axvline(x=59, color="black", alpha=0.1)
    ax.axvline(x=151, color="black", alpha=0.1)
    ax.axvline(x=243, color="black", alpha=0.1)
    ax.axvline(x=334, color="black", alpha=0.1)
    ax.axvspan(0, 59, color="blue", alpha=0.1)
    ax.axvspan(60, 151, color="green", alpha=0.1)
    ax.axvspan(152, 243, color="red", alpha=0.1)
    ax.axvspan(244, 334, color="orange", alpha=0.1)
    ax.axvspan(334, 366, color="blue", alpha=0.1)
    ax.set_xlim(0, 366)
    ax2.text(5, 156, "DJF", fontsize=fontsize)
    ax2.text(64, 156, "MAM", fontsize=fontsize)
    ax2.text(156, 156, "JJA", fontsize=fontsize)
    ax2.text(248, 156, "SON", fontsize=fontsize)
    plt.savefig("/"+name[i]+"split_velocity_fs="+str(fontsize)+".pdf")
    plt.show()

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15,4), tight_layout=True)

for i, ax in enumerate(axes):
    index[i][0].plot(ax=ax, color="blue", alpha=0.5)
    a1, = index[i][0].rolling(dayofyear=10, center=True).mean().plot(ax=ax, color="blue")
    index[i][1].plot(ax=ax, color="dodgerblue", alpha=0.5)
    a2, = index[i][1].rolling(dayofyear=10, center=True).mean().plot(ax=ax, color="dodgerblue")
    ax.set_xlabel("Day of year", fontsize=fontsize)
    ax.set_ylabel(ylabel[i], color="blue")
    ax2 = ax.twinx()
    v_xp.plot(ax=ax2, color="red", alpha=0.5)
    b1, = v_xp.rolling(dayofyear=10, center=True).mean().plot(ax=ax2, color="red")
    v_xn.plot(ax=ax2, color="orange", alpha=0.5)
    b2, = v_xn.rolling(dayofyear=10, center=True).mean().plot(ax=ax2, color="orange")
    ax2.set_ylabel("Velocity (km/day)", color="red", fontsize=fontsize)
    ax2.set_ylim(150, 460)
    l = ax2.legend([(a1, b1), (a2, b2)], ["Eastward", "Westward"], handler_map={tuple: HandlerTuple(ndivide=None)}, loc="upper left", fontsize=fontsize)
    ax.axvline(x=59, color="black", alpha=0.1)
    ax.axvline(x=151, color="black", alpha=0.1)
    ax.axvline(x=243, color="black", alpha=0.1)
    ax.axvline(x=334, color="black", alpha=0.1)
    ax.axvspan(0, 59, color="blue", alpha=0.1)
    ax.axvspan(60, 151, color="green", alpha=0.1)
    ax.axvspan(152, 243, color="red", alpha=0.1)
    ax.axvspan(244, 334, color="orange", alpha=0.1)
    ax.axvspan(334, 366, color="blue", alpha=0.1)
    ax.set_xlim(0, 366)
    ax2.text(5, 156, "DJF", fontsize=fontsize)
    ax2.text(64, 156, "MAM", fontsize=fontsize)
    ax2.text(156, 156, "JJA", fontsize=fontsize)
    ax2.text(248, 156, "SON", fontsize=fontsize)
    ax.set_title(name[i], fontsize=fontsize)

plt.savefig(dir+"Figures_fontsize/split_velocity_fs="+str(fontsize)+".pdf")
plt.show()

#%% Splitting eastward and westward velocities

sizep = hist_clarea.where(hist_speedx>0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
sizep_ens = hist_clarea.where(hist_speedx>0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block")) #all ensembles
sizep_std = hist_clarea.where(hist_speedx>0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block")).std("ensemble") #ensemble std
#sizep_std = hist_clarea.where(hist_speedx>0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").std(dim=("time", "day", "block", "ensemble"))
sizen = hist_clarea.where(hist_speedx<=0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
sizen_ens = hist_clarea.where(hist_speedx<=0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block")) #all ensembles
sizen_std = hist_clarea.where(hist_speedx<=0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block")).std("ensemble") #ensemble std
#sizen_std = hist_clarea.where(hist_speedx<=0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").std(dim=("time", "day", "block", "ensemble"))

v_xp = hist_speedx.where(hist_speedx>0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"))
v_xp_ens = hist_speedx.where(hist_speedx>0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block")) #all ensembles
v_xp_std = hist_speedx.where(hist_speedx>0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block")).std("ensemble") #ensemble std
#v_xp_std = hist_speedx.where(hist_speedx>0, drop=True).groupby("time.dayofyear").std(dim=("time", "block", "ensemble"))
v_xn = -hist_speedx.where(hist_speedx<=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"))
v_xn_ens = -hist_speedx.where(hist_speedx<=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block")) #all ensembles
v_xn_std = -hist_speedx.where(hist_speedx<=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block")).std("ensemble") #ensemble std
#v_xn_std = -hist_speedx.where(hist_speedx<=0, drop=True).groupby("time.dayofyear").std(dim=("time", "block", "ensemble"))

duration = hist_clarea.count("day")
durationp = duration.where(hist_speedx>0).where(duration!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"), skipna=True)
durationp_ens = duration.where(hist_speedx>0).where(duration!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block"), skipna=True) #all ensembles
durationp_std = duration.where(hist_speedx>0).where(duration!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block"), skipna=True).std("ensemble") #ensemble std
#durationp_std = duration.where(hist_speedx>0).where(duration!=0, drop=True).groupby("time.dayofyear").std(dim=("time", "block", "ensemble"), skipna=True)
durationn = duration.where(duration!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"), skipna=True)
durationn_ens = duration.where(duration!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block"), skipna=True) #all ensembles
durationn_std = duration.where(duration!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "block"), skipna=True).std("ensemble") #emsemble std
#durationn_std = duration.where(duration!=0, drop=True).groupby("time.dayofyear").std(dim=("time", "block", "ensemble"), skipna=True)

intensityp = hist_BIav.where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
intensityp_ens = hist_BIav.where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block")) #all ensembles
intensityp_std = hist_BIav.where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block")).std("ensemble") #ensemble std
#intensityp_std = hist_BIav.where(hist_BIav!=0, drop=True).groupby("time.dayofyear").std(dim=("time", "day", "block", "ensemble"))
intensityn = hist_BIav.where(hist_speedx<=0).where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
intensityn_ens = hist_BIav.where(hist_speedx<=0).where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block")) #all ensembles
intensityn_std = hist_BIav.where(hist_speedx<=0).where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block")).std("ensemble")
#intensityn_std = hist_BIav.where(hist_speedx<=0).where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))

index = [[sizep, sizen], [durationp, durationn], [intensityp, intensityn], [v_xp, v_xn]]
index_std = [[sizep_std, sizen_std], [durationp_std, durationn_std], [intensityp_std, intensityn_std], [v_xp_std, v_xn_std]]
index_ens = [[sizep_ens, sizen_ens], [durationp_ens, durationn_ens], [intensityp_ens, intensityn_ens], [v_xp_ens, v_xn_ens]]
ylabel = [r"Size (km$^2$)", "Duration (days)", "Average intensity BI (-)"]
name = ["Size", "Duration", "Intensity"]

#%% Plot westward and eastward characteristics
print("Start figure split seasons")
direction = ["east", "west"]
direction_title = ["Eastward-moving blocks", "Westward-moving blocks"]
for i in range(len(direction)):
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(8,14))
    plt.subplots_adjust(hspace=0)
    plt.rcParams.update({'font.size': fontsize})
    #Hide all axes
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.spines['left'].set_visible(False)
    
    #Plot the size
    index[0][i].plot(ax=ax1, color="blue", alpha=0.5)
    a1, = index[0][i].rolling(dayofyear=15, center=True).mean().plot(ax=ax1, color="black")
    ax1.fill_between(index_std[0][i].dayofyear, index[0][i], index[0][i]+index_std[0][i], color="blue", alpha=0.2) #std
    ax1.fill_between(index_std[0][i].dayofyear, index[0][i], index[0][i]-index_std[0][i], color="blue", alpha=0.2) #std
    ax1.set_ylabel(r"Size (km$^2$)", color="blue", fontsize=fontsize)
    #ax1.set_ylim(0.78e6, 2.25e6) #original
    #ax1.set_ylim(0.7e6, 2.4e6) #v1
    ax1.set_ylim(0.7e6, 2.6e6) #v2
    ax1.set_title(direction_title[i], fontsize=fontsize)
    
    #Plot duration
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    index[1][i].plot(ax=ax2, color="orange", alpha=0.5)
    b1, = index[1][i].rolling(dayofyear=15, center=True).mean().plot(ax=ax2, color="saddlebrown")
    ax2.fill_between(index_std[1][i].dayofyear, index[1][i], index[1][i]+index_std[1][i], color="orange", alpha=0.2) #std
    ax2.fill_between(index_std[1][i].dayofyear, index[1][i], index[1][i]-index_std[1][i], color="orange", alpha=0.2) #std
    ax2.set_ylabel("Duration (days)", color="orange", fontsize=fontsize)
    #ax2.set_ylim(5.25, 6.25) #original
    ax2.set_ylim(4, 7)
    
    #Plot intensity
    index[2][i].plot(ax=ax3, color="green", alpha=0.5)
    c1, = index[2][i].rolling(dayofyear=15, center=True).mean().plot(ax=ax3, color="darkslategrey")
    ax3.fill_between(index_std[2][i].dayofyear, index[2][i], index[2][i]+index_std[2][i], color="green", alpha=0.2) #std
    ax3.fill_between(index_std[2][i].dayofyear, index[2][i], index[2][i]-index_std[2][i], color="green", alpha=0.2) #std
    ax3.set_ylabel("Average intensity BI (-)", color="green", fontsize=fontsize)
    #ax3.set_ylim(1.2, 2.3) #original
    ax3.set_ylim(1, 2.6)
    
    #Plot velocity
    (index[3][i]*1000/(24*3600)).plot(ax=ax4, color="red", alpha=0.5)
    b1, = (index[3][i]*1000/(24*3600)).rolling(dayofyear=15, center=True).mean().plot(ax=ax4, color="darkred")
    ax4.fill_between((index_std[3][i]*1000/(24*3600)).dayofyear, (index[3][i]*1000/(24*3600)), (index[3][i]*1000/(24*3600))+(index_std[3][i]*1000/(24*3600)), color="red", alpha=0.2) #std
    ax4.fill_between((index_std[3][i]*1000/(24*3600)).dayofyear, (index[3][i]*1000/(24*3600)), (index[3][i]*1000/(24*3600))-(index_std[3][i]*1000/(24*3600)), color="red", alpha=0.2) #std
    ax4.yaxis.tick_right()
    #ax4.set_ylim(200, 400) #original
    #ax4.set_ylim(180, 450) #v1
    #ax4.set_ylim(165, 470) #v2
    ax4.set_ylim(165*1000/(24*3600), 470*1000/(24*3600)) #m/s
    ax4.yaxis.set_label_position("right")
    ax4.set_ylabel("Velocity (m/s)", color="red", fontsize=fontsize)
    
    #Add seasons in the background
    axes=[ax1, ax2, ax3, ax4]
    for j in range(4):
        axes[j].axvline(x=59, color="black", alpha=0.1)
        axes[j].axvline(x=151, color="black", alpha=0.1)
        axes[j].axvline(x=243, color="black", alpha=0.1)
        axes[j].axvline(x=334, color="black", alpha=0.1)
        axes[j].axvspan(0, 59, color="blue", alpha=0.1)
        axes[j].axvspan(60, 151, color="green", alpha=0.1)
        axes[j].axvspan(152, 243, color="red", alpha=0.1)
        axes[j].axvspan(244, 334, color="orange", alpha=0.1)
        axes[j].axvspan(334, 366, color="blue", alpha=0.1)
    
    ax4.set_xlim(0, 366)
    ax4.text(5, 190*1000/(24*3600), "DJF", fontsize=fontsize) #206 for km/day
    ax4.text(64, 190*1000/(24*3600), "MAM", fontsize=fontsize)
    ax4.text(156, 190*1000/(24*3600), "JJA", fontsize=fontsize)
    ax4.text(248, 190*1000/(24*3600), "SON", fontsize=fontsize)
    ax4.set_xlabel("Day of year", fontsize=fontsize)
    
    for ax in [ax1, ax2, ax3, ax4]:
      ax.tick_params(axis='both', labelsize=fontsize)
      ax.tick_params(axis='y', which="both", labelsize=fontsize)
      ax.yaxis.offsetText.set(size=fontsize)
    plt.savefig(dir+"Figures_fontsize/"+str(direction[i])+"_velocity_fs="+str(fontsize)+"_ms.pdf")


#%% Percentiles and KDE per season for size, 
S = ["DJF", "JJA"]
size_s = []
duration_s = []
intensity_s = []
speed_s = []

region = "Pacific" #Atlantic, Pacific or NH
if region == "Atlantic":
    test = True
elif region == "Pacific":
    test = False

for i in range(len(S)):
    if region == "Atlantic" or region == "Pacific":
        cond = (hist_clcmassX.sel(time=hist_clcmassX["time.season"]==S[i], day=0)>=-60) & (hist_clcmassX.sel(time=hist_clcmassX["time.season"]==S[i], day=0)<=120) 
        size_s.append(hist_clarea.sel(time=hist_clarea["time.season"]==S[i]).where(cond==test).mean("day").stack(z=("ensemble", "block", "time")).dropna("z"))
        d = hist_BIav.sel(time=hist_BIav["time.season"]==S[i]).where(cond==test).count("day").stack(z=("ensemble", "block", "time")).dropna("z")
        duration_s.append(d[d!=0])
        intensity_s.append(hist_BIav.sel(time=hist_BIav["time.season"]==S[i]).where(cond==test).mean("day").stack(z=("ensemble", "block", "time")).dropna("z"))
        speed_s.append(hist_speedx.sel(time=hist_speedx["time.season"]==S[i]).where(cond==test).stack(z=("ensemble", "block", "time")).dropna("z"))
    else:
        size_s.append(hist_clarea.sel(time=hist_clarea["time.season"]==S[i]).mean("day").stack(z=("ensemble", "block", "time")).dropna("z"))
        d = hist_BIav.sel(time=hist_BIav["time.season"]==S[i]).count("day").stack(z=("ensemble", "block", "time")).dropna("z")
        duration_s.append(d[d!=0])
        intensity_s.append(hist_BIav.sel(time=hist_BIav["time.season"]==S[i]).mean("day").stack(z=("ensemble", "block", "time")).dropna("z"))
        speed_s.append(hist_speedx.sel(time=hist_speedx["time.season"]==S[i]).stack(z=("ensemble", "block", "time")).dropna("z"))

speed_50 = np.arange(-1960*1000/(24*3600), 2040*1000/(24*3600), 40*1000/(24*3600))
percentiles = [10, 50, 90]

size_s_matrix = np.zeros([2, 3, len(speed_50)])
duration_s_matrix = np.zeros([2, 3, len(speed_50)])
intensity_s_matrix = np.zeros([2, 3, len(speed_50)])

for s in range(len(S)):
    for p in range(len(percentiles)):
        for i in range(len(speed_50)):
            if size_s[s][np.logical_and(speed_s[s]>-2000*1000/(24*3600)+i*40*1000/(24*3600), speed_s[s]<-2000*1000/(24*3600)+(i+1)*40*1000/(24*3600))].size==0:
                size_s_matrix[s,p,i]=np.nan
            else:
                size_s_matrix[s,p,i] = np.percentile(size_s[s][np.logical_and(speed_s[s]>-2000*1000/(24*3600)+i*40*1000/(24*3600), speed_s[s]<-2000*1000/(24*3600)+(i+1)*40*1000/(24*3600))].dropna("z"), percentiles[p])
            if duration_s[s][np.logical_and(speed_s[s]>-2000*1000/(24*3600)+i*40*1000/(24*3600), speed_s[s]<-2000*1000/(24*3600)+(i+1)*40*1000/(24*3600))].size==0:
                duration_s_matrix[s,p,i]=np.nan
            else:
                duration_s_matrix[s,p,i] = np.percentile(duration_s[s][np.logical_and(speed_s[s]>-2000*1000/(24*3600)+i*40*1000/(24*3600), speed_s[s]<-2000*1000/(24*3600)+(i+1)*40*1000/(24*3600))].dropna("z"), percentiles[p])
            if intensity_s[s][np.logical_and(speed_s[s]>-2000*1000/(24*3600)+i*40*1000/(24*3600), speed_s[s]<-2000*1000/(24*3600)+(i+1)*40*1000/(24*3600))].size==0:
                intensity_s_matrix[s,p,i]=np.nan
            else:
                intensity_s_matrix[s,p,i] = np.percentile(intensity_s[s][np.logical_and(speed_s[s]>-2000*1000/(24*3600)+i*40*1000/(24*3600), speed_s[s]<-2000*1000/(24*3600)+(i+1)*40*1000/(24*3600))].dropna("z"), percentiles[p])

size_s_matrix = xr.DataArray(size_s_matrix, name="size", coords=[("season", S), ("percentiles", percentiles), ("speed", speed_50)], dims=["season", "percentiles", "speed"])        
duration_s_matrix = xr.DataArray(duration_s_matrix, name="duration", coords=[("season", S), ("percentiles", percentiles), ("speed", speed_50)], dims=["season", "percentiles", "speed"])        
intensity_s_matrix = xr.DataArray(intensity_s_matrix, name="intensity", coords=[("season", S), ("percentiles", percentiles), ("speed", speed_50)], dims=["season", "percentiles", "speed"])        

#For KDE addition
cond = (hist_clcmassX.sel(day=0)>=-60) & (hist_clcmassX.sel(day=0)<=120) 

if region == "Atlantic" or region == "Pacific":
    x = hist_speedx.sel(time=hist_speedx["time.season"]=="JJA").where(cond==test).stack(z=("ensemble", "block", "time")).dropna("z")
    y = hist_clarea.where(cond==test).sel(time=hist_clarea["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time")).dropna("z")
    y_dur = hist_BImax.where(cond==test).sel(time=hist_BImax["time.season"]=="JJA").count("day").stack(z=("ensemble", "block", "time")).dropna("z")
    y_dur = y_dur[y_dur!=0]
    y_int = hist_BIav.where(cond==test).sel(time=hist_BIav["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time")).dropna("z")
    
    xw = hist_speedx.sel(time=hist_speedx["time.season"]=="DJF").where(cond==test).stack(z=("ensemble", "block", "time")).dropna("z")
    yw = hist_clarea.where(cond==test).sel(time=hist_clarea["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time")).dropna("z")
    yw_dur = hist_BImax.where(cond==test).sel(time=hist_BImax["time.season"]=="DJF").count("day").stack(z=("ensemble", "block", "time")).dropna("z")
    yw_dur = yw_dur[yw_dur!=0]
    yw_int = hist_BIav.where(cond==test).sel(time=hist_BIav["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time")).dropna("z")
else:
    x = hist_speedx.sel(time=hist_speedx["time.season"]=="JJA").stack(z=("ensemble", "block", "time")).dropna("z")
    y = hist_clarea.sel(time=hist_clarea["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time")).dropna("z")
    y_dur = hist_BImax.sel(time=hist_BImax["time.season"]=="JJA").count("day").stack(z=("ensemble", "block", "time")).dropna("z")
    y_dur = y_dur[y_dur!=0]
    y_int = hist_BIav.sel(time=hist_BIav["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time")).dropna("z")
    
    xw = hist_speedx.sel(time=hist_speedx["time.season"]=="DJF").stack(z=("ensemble", "block", "time")).dropna("z")
    yw = hist_clarea.sel(time=hist_clarea["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time")).dropna("z")
    yw_dur = hist_BImax.sel(time=hist_BImax["time.season"]=="DJF").count("day").stack(z=("ensemble", "block", "time")).dropna("z")
    yw_dur = yw_dur[yw_dur!=0]
    yw_int = hist_BIav.sel(time=hist_BIav["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time")).dropna("z")

#%% Add all indices and DJF/JJA in one plot
from matplotlib import ticker
from matplotlib.colors import Normalize

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 


fig, axes = plt.subplots(2, 3, figsize=(15,8), sharex=True)
plt.subplots_adjust(hspace=0)
plt.rcParams.update({'font.size': fontsize})
box = axes[1,0].get_position()
box2 = axes[1,1].get_position()
box3 = axes[1,2].get_position()
cbar_ax1 = fig.add_axes([box.xmin, box.ymin-0.1, box.width, 0.03])
cbar_ax2 = fig.add_axes([box2.xmin, box2.ymin-0.1, box2.width, 0.03])
cbar_ax3 = fig.add_axes([box3.xmin, box3.ymin-0.1, box3.width, 0.03])

axes[0,0].fill_betweenx(x1=-1500*1000/(24*3600), x2=xw.quantile(0.1), y=[0,6e6], color="grey", alpha=0.1)
axes[0,0].fill_betweenx(x1=xw.quantile(0.9), x2=1500*1000/(24*3600), y=[0,6e6], color="grey", alpha=0.1)
(size_s_matrix.sel(season="DJF", percentiles=10)*1e-6).plot(ax=axes[0,0], label="10%")
(size_s_matrix.sel(season="DJF", percentiles=50)*1e-6).plot(ax=axes[0,0], label="50%")
(size_s_matrix.sel(season="DJF", percentiles=90)*1e-6).plot(ax=axes[0,0], label="90%")
a = sns.kdeplot(ax=axes[0,0], x=xw, y=yw*1e-6, cmap="Blues", fill=True, norm=Normalize(vmax=0.09))#, cbar=True)
axes[0,0].set_title("Size", fontsize=fontsize)
#axes[0,0].set_xlabel("Velocity (km/day)")
axes[0,0].set_ylabel("DJF\nSize ($\cdot 10^6$ km²)", fontsize=fontsize)
axes[0,0].legend(fontsize=fontsize)
axes[0,0].set_ylim(0.2, 6)
axes[0,0].set_xlim(-1500*1000/(24*3600), 1500*1000/(24*3600))

axes[1,0].fill_betweenx(x1=-1500*1000/(24*3600), x2=x.quantile(0.1), y=[0,6], color="grey", alpha=0.1)
axes[1,0].fill_betweenx(x1=x.quantile(0.9), x2=1500*1000/(24*3600), y=[0,6], color="grey", alpha=0.1)
axes[1,0].text(s=r"$v_x\leq P10$", x=x.quantile(0.1), y=5.5, color="grey", horizontalalignment="right")
axes[1,0].text(s=r"$v_x\geq P90$", x=x.quantile(0.9), y=5.5, color="grey")
(size_s_matrix.sel(season="JJA", percentiles=10)*1e-6).plot(ax=axes[1,0], label="10%")
(size_s_matrix.sel(season="JJA", percentiles=50)*1e-6).plot(ax=axes[1,0], label="50%")
(size_s_matrix.sel(season="JJA", percentiles=90)*1e-6).plot(ax=axes[1,0], label="90%")
a = sns.kdeplot(ax=axes[1,0], x=x, y=y*1e-6, cmap="Blues", fill=True, cbar=True, norm=Normalize(vmax=0.09), cbar_ax=cbar_ax1, cbar_kws={"orientation":"horizontal", "ticks":[0.025e-1, 0.45e-1,0.9e-1], "format":formatter, "label":"Frequency"}) # vmax=1e-9,
axes[1,0].set_title("")
axes[1,0].set_xlabel("Velocity (m/s)", fontsize=fontsize)
axes[1,0].set_ylabel("JJA\nSize ($\cdot 10^6$ km²)", fontsize=fontsize)
#axes[1,0].legend()
axes[1,0].set_ylim(0.2, 6)
axes[1,0].set_xlim(-1500*1000/(24*3600), 1500*1000/(24*3600))
    
axes[0,1].fill_betweenx(x1=-1500*1000/(24*3600), x2=xw.quantile(0.1), y=[0,15], color="grey", alpha=0.1)
axes[0,1].fill_betweenx(x1=xw.quantile(0.9), x2=1500*1000/(24*3600), y=[0,15], color="grey", alpha=0.1)
duration_s_matrix.sel(season="DJF", percentiles=10).plot(ax=axes[0,1], label="10%")
duration_s_matrix.sel(season="DJF", percentiles=50).plot(ax=axes[0,1], label="50%")
duration_s_matrix.sel(season="DJF", percentiles=90).plot(ax=axes[0,1], label="90%")
a = sns.kdeplot(ax=axes[0,1], x=xw, y=yw_dur, cmap="Blues", fill=True, norm=Normalize(vmax=0.04))
axes[0,1].set_title("Duration", fontsize=fontsize)
#axes[0,1].set_xlabel("Velocity (km/day)")
axes[0,1].set_ylabel(r"Duration (days)", fontsize=fontsize)
axes[0,1].legend(fontsize=fontsize)
axes[0,1].set_ylim(3, 15)
axes[0,1].set_xlim(-1500*1000/(24*3600), 1500*1000/(24*3600))

axes[1,1].fill_betweenx(x1=-1500*1000/(24*3600), x2=x.quantile(0.1), y=[0,15], color="grey", alpha=0.1)
axes[1,1].fill_betweenx(x1=x.quantile(0.9), x2=1500*1000/(24*3600), y=[0,15], color="grey", alpha=0.1)
duration_s_matrix.sel(season="JJA", percentiles=10).plot(ax=axes[1,1], label="10%")
duration_s_matrix.sel(season="JJA", percentiles=50).plot(ax=axes[1,1], label="50%")
duration_s_matrix.sel(season="JJA", percentiles=90).plot(ax=axes[1,1], label="90%")
a = sns.kdeplot(ax=axes[1,1], x=x, y=y_dur, cmap="Blues", fill=True, cbar=True, norm=Normalize(vmax=0.04), cbar_ax=cbar_ax2, cbar_kws={"orientation":"horizontal", "ticks":[0.1e-2, 1.1e-2, 3.8e-2], "format":formatter, "label":"Frequency"})
axes[1,1].set_title("")
axes[1,1].set_xlabel("Velocity (m/s)", fontsize=fontsize)
axes[1,1].set_ylabel(r"Duration (days)", fontsize=fontsize)
#axes[1,1].legend()
axes[1,1].set_ylim(3, 15)
axes[1,1].set_xlim(-1500*1000/(24*3600), 1500*1000/(24*3600))

axes[0,2].fill_betweenx(x1=-1500, x2=xw.quantile(0.1), y=[0,5], color="grey", alpha=0.1)
axes[0,2].fill_betweenx(x1=xw.quantile(0.9), x2=1500, y=[0,5], color="grey", alpha=0.1)
intensity_s_matrix.sel(season="DJF", percentiles=10).plot(ax=axes[0,2], label="10%")
intensity_s_matrix.sel(season="DJF", percentiles=50).plot(ax=axes[0,2], label="50%")
intensity_s_matrix.sel(season="DJF", percentiles=90).plot(ax=axes[0,2], label="90%")
a = sns.kdeplot(ax=axes[0,2], x=xw, y=yw_int, cmap="Blues", fill=True, norm=Normalize(vmax=0.07))
axes[0,2].set_title("Intensity", fontsize=fontsize)
#axes[0,2].set_xlabel("Velocity (km/day)")
axes[0,2].set_ylabel("Average BI (-)", fontsize=fontsize)
axes[0,2].legend(fontsize=fontsize)
axes[0,2].set_ylim(0, 5)
axes[0,2].set_xlim(-1500*1000/(24*3600), 1500*1000/(24*3600))

axes[1,2].fill_betweenx(x1=-1500, x2=x.quantile(0.1), y=[0,5], color="grey", alpha=0.1)
axes[1,2].fill_betweenx(x1=x.quantile(0.9), x2=1500, y=[0,5], color="grey", alpha=0.1)
intensity_s_matrix.sel(season="JJA", percentiles=10).plot(ax=axes[1,2], label="10%")
intensity_s_matrix.sel(season="JJA", percentiles=50).plot(ax=axes[1,2], label="50%")
intensity_s_matrix.sel(season="JJA", percentiles=90).plot(ax=axes[1,2], label="90%")
a = sns.kdeplot(ax=axes[1,2], x=x, y=y_int, cmap="Blues", fill=True, cbar=True, norm=Normalize(vmax=0.07), cbar_ax=cbar_ax3, cbar_kws={"orientation":"horizontal", "ticks":[0.6e-2, 3.5e-2, 6.6e-2], "format":formatter, "label":"Frequency"})
axes[1,2].set_title("")
axes[1,2].set_xlabel("Velocity (m/s)", fontsize=fontsize)
axes[1,2].set_ylabel("Average BI (-)", fontsize=fontsize)
#axes[1,2].legend()
axes[1,2].set_ylim(0, 5)
axes[1,2].set_xlim(-1500*1000/(24*3600), 1500*1000/(24*3600))
plt.savefig(dir+"Final_Figures/review/percentiles+KDE_speed_P10_P90_sharex_fs="+str(fontsize)+"_"+str(region)+"_ms.pdf", bbox_inches='tight')

#%% Make scatterplots for daily means over 16 ensembles + running mean
from sctriangulate.colors import build_custom_continuous_cmap
new_cmap = build_custom_continuous_cmap([127,127,255],[127,191,127],[255,127,127],[255,209,127],[127,127,255])

fact = 4.5 #factor increase in size

sizep = hist_clarea.where(hist_speedx>0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
sizen = hist_clarea.where(hist_speedx<=0).where(hist_clarea!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
sizep_rol = sizep.rolling(dayofyear=15, center=True).mean()
sizen_rol = sizen.rolling(dayofyear=15, center=True).mean()
v_xp = hist_speedx.where(hist_speedx>0).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"))
v_xn = -hist_speedx.where(hist_speedx<=0).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"))
v_xp_rol = v_xp.rolling(dayofyear=15, center=True).mean()
v_xn_rol = v_xn.rolling(dayofyear=15, center=True).mean()
intensityp = hist_BIav.where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
intensityn = hist_BIav.where(hist_speedx<=0).where(hist_BIav!=0, drop=True).groupby("time.dayofyear").mean(dim=("time", "day", "block", "ensemble"))
intensityp_rol = intensityp.rolling(dayofyear=15, center=True).mean()
intensityn_rol = intensityn.rolling(dayofyear=15, center=True).mean()
#%%
plt.figure(figsize=(7,4))
scatter = plt.scatter(v_xp, sizep, c=v_xp["dayofyear"], s=intensityp**fact, cmap=new_cmap, alpha=0.1)
plt.scatter(v_xp_rol, sizep_rol, c=v_xp_rol["dayofyear"], s=intensityp_rol**fact, cmap=new_cmap)
plt.scatter(v_xn, sizen, c=v_xn["dayofyear"], cmap=new_cmap, s=intensityn**fact, alpha=0.1)
plt.scatter(v_xn_rol, sizen_rol, c=v_xn_rol["dayofyear"], s=intensityn_rol**fact, cmap=new_cmap)
kw = dict(prop="sizes", num=6, color="blue", fmt="{x:.2f}",
          func=lambda s: s**(1/fact), alpha=0.3)
legend2 = plt.legend(*scatter.legend_elements(**kw), title="BI (-)", loc="upper right", bbox_to_anchor=(1.5,1.06), frameon=False, fontsize=fontsize)

plt.colorbar(label="Day of year")
plt.annotate("Westward moving", xy=(250*1000/(24*3600),1.8e6), xycoords="data", xytext=(240*1000/(24*3600), 2.3e6), arrowprops=dict(facecolor="black", width=0.5), horizontalalignment="center", verticalalignment="top")
plt.annotate("Eastward moving", xy=(330*1000/(24*3600),0.98e6), xycoords="data", xytext=(350*1000/(24*3600), 0.78e6), arrowprops=dict(facecolor="black", width=0.5), horizontalalignment="center", verticalalignment="top")
plt.xlabel("Velocity (m/s)", fontsize=fontsize)
plt.ylabel(r"Size (km$^2$)", fontsize=fontsize)
plt.ylim(0.65e6, 2.4e6)
plt.xlim(192*1000/(24*3600), 408*1000/(24*3600))
plt.savefig(dir+"Figures_fontsize/scatter_size_velocity_BI_fs="+str(fontsize)+"_ms.pdf", bbox_inches="tight")


#%% Check if timeperiod matters
# Take 10 year averages and plot velocity 
color = iter(cm.Reds(np.linspace(0.25, 1, 4)))
plt.figure()
for i in range(4):
    t1 = 1850+i*40
    t2 = 1890+i*40
    v_xp_10 = hist_speedx.where(hist_speedx>0, drop=True).sel(time=slice(str(t1)+"-01-16", str(t2)+"-01-16")).groupby("time.dayofyear").mean(dim=("time", "block", "ensemble"))
    v_xp_rol = v_xp_10.rolling(dayofyear=15, center=True).mean()
    c = next(color)
    plt.plot(v_xp_rol["dayofyear"], v_xp_rol, color=c, label=str(t1)+"-"+str(t2))
plt.xlabel("Day of year", fontsize=fontsize)
plt.ylabel("Velocity (m/s)", fontsize=fontsize)
#plt.colorbar()
plt.legend()
plt.savefig(dir+"Final_Figures/review/timeslide_velocity_fs="+str(fontsize)+"_ms.pdf")
plt.show()
