#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  12 2023

@author: mourik
Try-out for combining different temperature plots
"""
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cartopy.crs as ccrs

def speedx(clcmassX, clcmassY):
    R = 6371 #km
    dx = clcmassX.ffill(dim="day").isel(day=-1)-clcmassX.sel(day=0)
    dx = dx.where(dx.values<180, dx.values-360)
    dx = dx.where(dx.values>-180, dx.values+360)
    dxkm = dx*2*np.pi*R*np.cos(np.deg2rad(clcmassY.sel(day=0)))*2.5/360
    v = dxkm/clcmassX.count("day")
    return v

#%% Load in temperature anomalies. Aim: get a square on the north and south side of the block and save mean temperature.

tas = xr.open_dataset("/scratch/depfg/6196306/MSc_thesis/Temperature/ECEARTH3Bis/ensembles/tas/anomaly/tas_hist01_anom.nc").tas
majlen = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/per_lon/Intensity/long/large/hist_BI_majlen_large.nc")#.sel(ensemble=15)
minlen = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/per_lon/Intensity/long/large/hist_BI_minlen_large.nc")#.sel(ensemble=15)
hist_clcmassX = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/per_lon/Intensity/long/large/hist_BI_clcmassX_large_weighted.nc")#.sel(ensemble=15)
hist_clcmassY = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/per_lon/Intensity/long/large/hist_BI_clcmassY_large_weighted.nc")#.sel(ensemble=15)
hist_BIav = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/per_lon/Intensity/long/large/hist_BI_BIav_large.nc")
hist_clarea = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/per_lon/Intensity/long/large/hist_BI_clarea_large_weighted.nc")

tas["time"] = tas.indexes["time"].normalize()
hist_speedx = speedx(hist_clcmassX, hist_clcmassY)

T_block_xr = xr.open_mfdataset("/scratch/depfg/6196306/MSc_thesis/Temperature/Tblock/20x40/varlat/Tas_block_hist*_20x40_land.nc", combine="nested", concat_dim="ensemble").Temp_anom

#Try plot with density lines in the margins
season="DJF"
#perc = [0.0, 0.2, 0.4, 0.6, 0.8]
perc = np.arange(0, 1, 0.01)
print(perc)

delta_T = np.zeros([100,100])
for i in range(100):
    for j in range(100):
      v_cond = np.logical_and(hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")>=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0 + perc[i]), hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")<=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0.2 + perc[i]))
      a_cond = np.logical_and(hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))>=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(0.8 - perc[j]), hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))<=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(1.0 - perc[j]))
      delta_T[i,j]= abs(T_block_xr.sel(ensemble=slice(0,16), lat=slice(-20,19), lon=slice(-40,39)).where(np.logical_and(v_cond, a_cond))).mean(dim=("time", "block", "ensemble", "lat", "lon"))	
      print(i,j)

plt.figure()
plt.scatter(delta_T)
plt.savefig("/scratch/depfg/6196306/MSc_thesis/varlat/tryout.pdf")
