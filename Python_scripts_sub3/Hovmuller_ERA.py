#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:35:29 2024

@author: 6196306
Hövmoller diagrams for blocking velocity
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

#functions that may be necessary
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

dir = "/scratch/depfg/6196306/MSc_thesis/"

#ERA
ERA_clcmassX = xr.open_dataset(dir+"Blockings/per_lon/Intensity/long/ERA_BI_data_weighted.nc").clcmassX
ERA_clcmassY = xr.open_dataset(dir+"Blockings/per_lon/Intensity/long/ERA_BI_data_weighted.nc").clcmassY
ERA_speedx = speedx(ERA_clcmassX, ERA_clcmassY)
ERA_speedx = ERA_speedx*1000/(24*3600)



#%% Make composit for the same intervals as for the temperature-size-velocity plot
#v_cond = np.logical_and(hist_speedx_season>=hist_speedx_season.quantile(perc_vel1[i]), hist_speedx_season<=hist_speedx_season.quantile(perc_vel2[i]))

z500_ERA = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/ERA5/era5_z500.nc")
BI_ERA = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/Intensity/block_ERA_BI_1950_2022.nc")

z500_clim = z500_ERA.groupby("time.month").mean("time")
z500_anom = (z500_ERA.groupby("time.month")-z500_clim).load()

#%% Find block from figure 3
start_date = pd.Timestamp("28-05-1963")
block = 0
# Get longitude, latitude, and length for this block
lons = ERA_clcmassX.sel(block=block, time=start_date).dropna("day")
lats = ERA_clcmassY.sel(block=block, time=start_date)
length = int(lons.count().values)
end_date = start_date + pd.Timedelta(days=length)

# Extract z500 and BI data
z500_sel = z500_anom.sel(time=slice(start_date, end_date), lat=slice(lats.mean() - 5, lats.mean() + 5)).mean("lat")
BI_sel = BI_ERA.sel(time=slice(start_date, end_date), lat=slice(lats.mean() - 5, lats.mean() + 5)).mean("lat")

BI_sel_ext = BI_ERA.sel(time=slice(pd.Timestamp("27-05-1963"), pd.Timestamp("27-05-1963")+pd.Timedelta(days=length+2)), lat=slice(lats.mean()-5, lats.mean()+5)).mean("lat")
#%% Plot

fig, ax = plt.subplots(1, 2, sharey=True, figsize=(12, 4))
plt.rc('font', size=14)
z500_plot = z500_sel.plot(ax=ax[0], label="Z500 Anomaly", vmin=-240, vmax=240, cmap="coolwarm", add_colorbar=False)
cbar_z500 = fig.colorbar(z500_plot, ax=ax[0])
cbar_z500.set_label(r"$\Delta$ z500 (m)", fontsize=14)
ax[0].scatter(lons, pd.date_range(start=start_date, end=end_date, freq="D")[:-1], color="black", label="Block Path")
ax[0].contour(BI_sel_ext.lon, BI_sel_ext.time, BI_sel_ext, levels=[0], colors='black')

BI_plot = BI_sel.where(BI_sel!=0).plot(ax=ax[1], vmin=0, vmax=4, label="BI", cmap="viridis", add_colorbar=False)
cbar_BI = fig.colorbar(BI_plot, ax=ax[1])
cbar_BI.set_label("BI (-)", fontsize=14)
ax[1].scatter(lons, pd.date_range(start=start_date, end=end_date, freq="D")[:-1], color="black")
ax[1].set_ylabel(" ")
plt.tight_layout()
fig.savefig("/eejit/home/6196306/Data/MSc_thesis/Hovmuller-ERA_block0_1963-05-26_v3.pdf")
plt.show()