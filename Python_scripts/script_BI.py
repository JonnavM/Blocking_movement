#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  23 11:19:15 2024

@author: mourik

Short script to make one snapshot of the tracking mechanism including BI
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

#Download dataset:
da = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/per_lon/Intensity/block_var_ERA_BI_1950_2022.nc")
lat_min = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/LAT_min/LAT_min_lon/lat_min_lon_ERA_1950_2022.nc")
z500 = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/ERA5/era5_z500.nc")

date = "1963-05-26" 
fontsize=9
cmap = LinearSegmentedColormap.from_list('', ['white', *plt.cm.Reds(np.arange(255))])

z500_sel = z500.sel(time=date, method="nearest").sel(lon=slice(-75, 80), lat=slice(20, 90))

#Plot random day
fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
plt.rcParams.update({'font.size': fontsize})
da.sel(time=date, method="nearest").plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=6, cbar_kwargs={"label": "Blocking intensity (-)", "shrink": 0.5, "extend":"max"})
CS = ax.contour(z500_sel.lon, z500_sel.lat, z500_sel, transform=ccrs.PlateCarree(), levels=[5000, 5500, 6000], linestyles="-", colors="black", alpha=0.8, linewidths=0.8)
ax.contour(z500_sel.lon, z500_sel.lat, z500_sel, transform=ccrs.PlateCarree(), levels=[5400, 5600, 5700], linestyles="-", colors="grey", alpha=0.8, linewidths=0.8)
ax.clabel(CS, inline=True, fontsize=9)
lat_min.groupby("time.month").mean("time").sel(month=5).plot(ax=ax, transform=ccrs.PlateCarree())
ax.add_patch(mpatches.Rectangle(xy=[-18, 42], width=80, height=40,
                                facecolor='none', edgecolor='red',
                                transform=ccrs.PlateCarree(), alpha=0.5))
ax.add_patch(mpatches.Rectangle(xy=[22, 42], width=40, height=40,
                                facecolor='none', edgecolor='red',
                                transform=ccrs.PlateCarree(), linestyle="dotted", alpha=0.3))
ax.add_patch(mpatches.Rectangle(xy=[-18, 62], width=80, height=20,
                                facecolor='none', edgecolor='red',
                                transform=ccrs.PlateCarree(), linestyle="dotted", alpha=0.3))
ax.scatter(22, 62, color="black", transform=ccrs.PlateCarree())
ax.arrow(22, 62, 10, 0, head_width=2, color="black", transform=ccrs.PlateCarree())
ax.coastlines(resolution="50m", linewidth=0.5, alpha=0.8)
ax.set_extent([-75, 80, 20, 90], crs=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=True, linewidth=0.1, color='black', alpha=0.5)
#gl.xlocator = plt.FixedLocator(range(-80, 100, 40))  # Adjust xlocator
#gl.ylocator = plt.FixedLocator(range(0, 100, 20))    # Adjust ylocator
gl.top_labels = False
gl.right_labels = False
gl.xlines = False
gl.ylines = False
ax.set_title(" ")
fig.savefig("/scratch/depfg/6196306/MSc_thesis/Figures_fontsize/tracking_mechanism_v8_fs="+str(fontsize)+".pdf", bbox_inches="tight")