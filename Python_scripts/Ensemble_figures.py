#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 11:28:15 2022

@author: mourik

Time restriction of 4 days is encluded. The code for all ensemble figures in the thesis are presented here
"""
#%% Import needed modules
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
from dask.diagnostics import ProgressBar
from matplotlib.colors import ListedColormap

#%% Put files into one netcdf file and change unit of z500 of ERA5
dir = "Fill in directory"
#z500 = xr.open_mfdataset("/ERA5/*.nc").z500/9.81
#z500.to_netcdf("/ERA5/era5_z500.nc")
z500 = xr.open_dataarray(dir+"ERA5/era5_z500.nc")
print("Z500 loaded")

#%% Calculate ensemble mean for any variable you like. Last used here: blocking index historical ECE3p5

#sep_files = xr.open_mfdataset(dir+"ECEARTH3p5/historical/r*i1p5f1/ensemble/zg_day_EC-Earth3_historical_r*i1p5f1_seasmean.nc", combine="nested", concat_dim=("Ensemble"))
#ens_mean = sep_files.mean("Ensemble")
#delayed1 = ens_mean.to_netcdf(dir+"ECEARTH3p5/historical/zg_seasmean_EC-Earth3_historical_ens.nc", compute=False)
#with ProgressBar():
#    result1 = delayed1.compute()
#print("Ensemble file generated")    
    
#%% Load in the blockings for ERA and EC-Earth
block_hist_mean = xr.open_dataarray(dir+"Blockings/per_lon/long/block_hist_BI_long_seasmean_ens.nc").mean("lev") #ensemble mean
block_hist = xr.open_mfdataset(dir+"Blockings/per_lon/long/block_hist*_BI_long_seasmean.nc", combine="nested", concat_dim=("Ensemble")).trackID.mean("lev") #All loose files
block_hist_BI_mean = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/block_hist_BI_long_seasmean_ens.nc").mean("lev") #ensemble mean including intensity
block_hist_BI = xr.open_mfdataset(dir+"Blockings/per_lon/Intensity/long/block_hist*_BI_long_seasmean.nc", combine="nested", concat_dim=("Ensemble")).index.mean("lev") #All loose files including intensity

block_ERA = xr.open_dataset(dir+"Blockings/per_lon/long/block_ERA_360_long_seasmean.nc").trackID #without intensity
block_ERA_BI = xr.open_dataset(dir+"Blockings/per_lon/Intensity/long/block_ERA_BI_long_seasmean.nc").index #with intensity

z500_hist_mean = xr.open_dataarray(dir+"ECEARTH3p5/historical/zg_seasmean_EC-Earth3_historical_ens.nc")
z500_ERA = xr.open_dataset(dir+"ERA5/era5_z500_seasmean.nc").z500

print("Data loaded in")
"""
GHGN_hist_mean = xr.open_dataarray("/GHG/GHGN_hist_seasmean_ens.nc")
GHGS_hist_mean = xr.open_dataarray("/GHG/GHGS_hist_seasmean_ens.nc")
GHGN_ERA = xr.open_dataset("/GHG/GHGN_ERA_1950_2022_seasmean.nc").index
GHGS_ERA = xr.open_dataset("/GHG/GHGS_ERA_1950_2022_seasmean.nc").index
"""

# Coordinates
time_hist = block_hist_BI["time"].data
time_ERA = block_ERA_BI["time"].data
lat = block_hist_BI['lat'].data #Is the same for all of them
lon = block_hist_BI['lon'].data

#Dates in common
date = block_hist_BI.time.sel(time=slice("01-01-1951", "31-12-2014")).data

#Set fontsize
fontsize = 12

#%% Recreate figures from Davini2020: Blocking frequency
#Define new colormaps
cmap_1 = plt.cm.get_cmap('bwr', 11)
cmap_2 = plt.cm.get_cmap("YlOrRd", 10)

newcolors = cmap_2(np.linspace(0, 1, 10))
white = np.array([0, 0, 0, 0])
newcolors[:1, :] = white
cmap_2 = ListedColormap(newcolors)

#Rename datasets for shorter notation. First without blocking intensity
block_model = block_hist_mean.sel(time=slice("01-01-1951", "31-12-2014"))
block_obs = block_ERA.sel(time=slice("01-01-1951", "31-12-2014"))

# Plot yearly mean, winter and summer in one plot without blocking intensity
fig = plt.figure()
ax = fig.add_subplot(2,3,1, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines(linewidth=0.5)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
a = ax.pcolormesh(lon, lat, block_obs.groupby("time.year").mean("time").mean("year")*100, cmap=cmap_2, transform=ccrs.PlateCarree(), vmax=20)
ax.contour(lon, lat, block_obs.groupby("time.year").mean("time").mean("year")*100, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=20)
ax.set_title("Annual mean", fontsize=fontsize)
ax1 = fig.add_subplot(2,3,4, projection=ccrs.Orthographic(central_latitude=90))
ax1.coastlines(linewidth=0.5)
ax1.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
b = ax1.pcolormesh(lon, lat, block_model.groupby("time.year").mean("time").mean("year")*100-block_obs.groupby("time.year").mean("time").mean("year")*100, cmap=cmap_1, vmin=-10, vmax=10, transform=ccrs.PlateCarree())
ax1.contour(lon, lat, block_model.groupby("time.year").mean("time").mean("year")*100, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20], linewidths=0.5, transform=ccrs.PlateCarree(), colors="k")
fig.text(x=0.08, y=0.17, s="ECE3p5-ERA5", fontsize=fontsize, rotation="vertical")
fig.text(x=0.08, y=0.68, s="ERA5", fontsize=fontsize, rotation="vertical")

ax2 = fig.add_subplot(2,3,2, projection=ccrs.Orthographic(central_latitude=90))
ax2.coastlines(linewidth=0.5)
ax2.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
a = ax2.pcolormesh(lon, lat, block_obs.groupby("time.season").mean("time").sel(season="DJF")*100, cmap=cmap_2, transform=ccrs.PlateCarree(), vmax=20)
ax2.contour(lon, lat, block_obs.groupby("time.season").mean("time").sel(season="DJF")*100, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=20)
ax2.set_title("DJF", fontsize=fontsize)
ax3 = fig.add_subplot(2,3,5, projection=ccrs.Orthographic(central_latitude=90))
ax3.coastlines(linewidth=0.5)
ax3.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
b = ax3.pcolormesh(lon, lat, block_model.groupby("time.season").mean("time").sel(season="DJF")*100-block_obs.groupby("time.season").mean("time").sel(season="DJF")*100, cmap=cmap_1, vmin=-10, vmax=10, transform=ccrs.PlateCarree())
ax3.contour(lon, lat, block_model.groupby("time.season").mean("time").sel(season="DJF")*100, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20], linewidths=0.5, transform=ccrs.PlateCarree(), colors="k")
#ax3.set_title("mod-obs")

ax4 = fig.add_subplot(2,3,3, projection=ccrs.Orthographic(central_latitude=90))
ax4.coastlines(linewidth=0.5)
ax4.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
a = ax4.pcolormesh(lon, lat, block_obs.groupby("time.season").mean("time").sel(season="JJA")*100, cmap=cmap_2, transform=ccrs.PlateCarree(), vmax=20)
ax4.contour(lon, lat, block_obs.groupby("time.season").mean("time").sel(season="JJA")*100, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=20)
ax4.set_title("JJA", fontsize=fontsize)
ax5 = fig.add_subplot(2,3,6, projection=ccrs.Orthographic(central_latitude=90))
ax5.coastlines(linewidth=0.5)
ax5.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
b = ax5.pcolormesh(lon, lat, block_model.groupby("time.season").mean("time").sel(season="JJA")*100-block_obs.groupby("time.season").mean("time").sel(season="JJA")*100, cmap=cmap_1, vmin=-10, vmax=10, transform=ccrs.PlateCarree())
ax5.contour(lon, lat, block_model.groupby("time.season").mean("time").sel(season="JJA")*100, [2, 4, 6, 8, 10, 12, 14, 16, 18, 20], linewidths=0.5, transform=ccrs.PlateCarree(), colors="k")

fig.subplots_adjust(right=0.8)
cbar_ax1 = fig.add_axes([0.85, 0.14, 0.02, 0.31])
cbar_ax2 = fig.add_axes([0.85, 0.55, 0.02, 0.31])
plt.colorbar(b, extend="both", cax=cbar_ax1).set_label(label="Blocked days (%)", size=fontsize)
plt.colorbar(a, extend="max", cax=cbar_ax2).set_label(label="Blocked days (%)", size=fontsize)
plt.savefig("/Ens_diff_mod_obs_DJF_JJA_year_fs="+str(fontsize)+".pdf")
plt.show()

print("Blocking frequency plotted")

#%% Blocking intensity
block_model_BI = block_hist_BI_mean.sel(time=slice("01-01-1951", "31-12-2014"))
block_obs_BI = block_ERA_BI.sel(time=slice("01-01-1951", "31-12-2014"))

# Plot yearly mean, winter and summer in one plot with blocking intensity
fig = plt.figure()
ax = fig.add_subplot(2,3,1, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines(linewidth=0.5)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
a = ax.pcolormesh(lon, lat, block_obs_BI.groupby("time.year").mean("time").mean("year"), cmap=cmap_2, transform=ccrs.PlateCarree(), vmax=0.4)
ax.contour(lon, lat, block_obs_BI.groupby("time.year").mean("time").mean("year"), [0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=0.4)
ax.set_title("Annual mean", fontsize=fontsize)
ax1 = fig.add_subplot(2,3,4, projection=ccrs.Orthographic(central_latitude=90))
ax1.coastlines(linewidth=0.5)
ax1.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
b = ax1.pcolormesh(lon, lat, block_model_BI.groupby("time.year").mean("time").mean("year")-block_obs_BI.groupby("time.year").mean("time").mean("year"), cmap=cmap_1, vmin=-0.2, vmax=0.2, transform=ccrs.PlateCarree())
ax1.contour(lon, lat, block_model_BI.groupby("time.year").mean("time").mean("year"), [0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4], linewidths=0.5, transform=ccrs.PlateCarree(), colors="k")
fig.text(x=0.08, y=0.17, s="ECE3p5-ERA5", fontsize=fontsize, rotation="vertical")
fig.text(x=0.08, y=0.68, s="ERA5", fontsize=fontsize, rotation="vertical")

ax2 = fig.add_subplot(2,3,2, projection=ccrs.Orthographic(central_latitude=90))
ax2.coastlines(linewidth=0.5)
ax2.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
a = ax2.pcolormesh(lon, lat, block_obs_BI.groupby("time.season").mean("time").sel(season="DJF"), cmap=cmap_2, transform=ccrs.PlateCarree(), vmax=0.4)
ax2.contour(lon, lat, block_obs_BI.groupby("time.season").mean("time").sel(season="DJF"), [0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=0.4)
ax2.set_title("DJF", fontsize = fontsize)
ax3 = fig.add_subplot(2,3,5, projection=ccrs.Orthographic(central_latitude=90))
ax3.coastlines(linewidth=0.5)
ax3.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
b = ax3.pcolormesh(lon, lat, block_model_BI.groupby("time.season").mean("time").sel(season="DJF")-block_obs_BI.groupby("time.season").mean("time").sel(season="DJF"), cmap=cmap_1, vmin=-0.2, vmax=0.2, transform=ccrs.PlateCarree())
ax3.contour(lon, lat, block_model_BI.groupby("time.season").mean("time").sel(season="DJF"), [0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4], linewidths=0.5, transform=ccrs.PlateCarree(), colors="k")
#ax3.set_title("mod-obs")

ax4 = fig.add_subplot(2,3,3, projection=ccrs.Orthographic(central_latitude=90))
ax4.coastlines(linewidth=0.5)
ax4.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
a = ax4.pcolormesh(lon, lat, block_obs_BI.groupby("time.season").mean("time").sel(season="JJA"), cmap=cmap_2, transform=ccrs.PlateCarree(), vmax=0.4)
ax4.contour(lon, lat, block_obs_BI.groupby("time.season").mean("time").sel(season="JJA"), [0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=0.4)
ax4.set_title("JJA", fontsize=fontsize)
ax5 = fig.add_subplot(2,3,6, projection=ccrs.Orthographic(central_latitude=90))
ax5.coastlines(linewidth=0.5)
ax5.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
b = ax5.pcolormesh(lon, lat, block_model_BI.groupby("time.season").mean("time").sel(season="JJA")-block_obs_BI.groupby("time.season").mean("time").sel(season="JJA"), cmap=cmap_1, vmin=-0.2, vmax=0.2, transform=ccrs.PlateCarree())
ax5.contour(lon, lat, block_model_BI.groupby("time.season").mean("time").sel(season="JJA"), [0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4], linewidths=0.5, transform=ccrs.PlateCarree(), colors="k")

fig.subplots_adjust(right=0.8)
cbar_ax1 = fig.add_axes([0.85, 0.14, 0.02, 0.31])
cbar_ax2 = fig.add_axes([0.85, 0.55, 0.02, 0.31])
plt.colorbar(b, extend="both", cax=cbar_ax1).set_label( label="Blocked days (%)", size=fontsize)
plt.colorbar(a, extend="max", cax=cbar_ax2).set_label(label="Blocked days (%)", size=fontsize)
plt.savefig("Ens_BI_diff_mod_obs_DJF_JJA_year_fs="+str(fontsize)+".pdf")
plt.show()

print("Blocking intensity plotted")

#%% Intramodel variation
cmap_3 = plt.cm.get_cmap("Blues", 10)
newcolors = cmap_3(np.linspace(0, 1, 10))
white = np.array([0, 0, 0, 0])
newcolors[:1, :] = white
cmap_3 = ListedColormap(newcolors)

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(2,2,2, projection=ccrs.Orthographic(central_latitude=90))
ax.coastlines(linewidth=0.5)
a = ax.pcolormesh(lon, lat, block_hist_BI.groupby("time.season").mean("time").sel(season="DJF").std("Ensemble"), cmap=cmap_3, transform=ccrs.PlateCarree(), vmax=0.03)
ax.contour(lon, lat, block_hist_BI.groupby("time.season").mean("time").sel(season="DJF").std("Ensemble"), [0.003, 0.006, 0.009, 0.012, 0.015, 0.018, 0.021, 0.024, 0.027, 0.03], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=0.03)
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax.set_title("Intensity", fontsize=fontsize)

ax1 = fig.add_subplot(2,2,4, projection=ccrs.Orthographic(central_latitude=90))
ax1.coastlines(linewidth=0.5)
b = ax1.pcolormesh(lon, lat, block_hist_BI.groupby("time.season").mean("time").sel(season="JJA").std("Ensemble"), cmap=cmap_3, transform=ccrs.PlateCarree(), vmax=0.03)
ax1.contour(lon, lat, block_hist_BI.groupby("time.season").mean("time").sel(season="JJA").std("Ensemble"), [0.003, 0.006, 0.009, 0.012, 0.015, 0.018, 0.021, 0.024, 0.027, 0.03], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=0.03)
ax1.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())

ax2 = fig.add_subplot(2,2,1, projection=ccrs.Orthographic(central_latitude=90))
ax2.coastlines(linewidth=0.5)
a = ax2.pcolormesh(lon, lat, block_hist.groupby("time.season").mean("time").sel(season="DJF").std("Ensemble")*100, cmap=cmap_3, transform=ccrs.PlateCarree(), vmax=1.5)
ax2.contour(lon, lat, block_hist.groupby("time.season").mean("time").sel(season="DJF").std("Ensemble")*100, [0.15, 0.3, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.5], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=1.5)
ax2.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax2.set_title("Frequency", fontsize=fontsize)
ax2.text(-0.07, 0.55, "DJF", va="bottom", ha="center", rotation="vertical", rotation_mode="anchor", transform=ax2.transAxes, fontsize=fontsize)

ax3 = fig.add_subplot(2,2,3, projection=ccrs.Orthographic(central_latitude=90))
ax3.coastlines(linewidth=0.5)
c = ax3.pcolormesh(lon, lat, block_hist.groupby("time.season").mean("time").sel(season="JJA").std("Ensemble")*100, cmap=cmap_3, transform=ccrs.PlateCarree(), vmax=1.5)
ax3.contour(lon, lat, block_hist.groupby("time.season").mean("time").sel(season="JJA").std("Ensemble")*100, [0.15, 0.3, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.5], linewidths=0.5, transform=ccrs.PlateCarree(), colors= "k", vmax=1.5)
ax3.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax3.text(-0.07, 0.55, "JJA", va="bottom", ha="center", rotation="vertical", rotation_mode="anchor", transform=ax3.transAxes, fontsize=fontsize)

box1 = ax1.get_position()
box2 = ax3.get_position()

cbar_ax1 = fig.add_axes([box1.xmin, box1.ymin-0.1, box1.width, 0.03])
cbar_ax2 = fig.add_axes([box2.xmin, box2.ymin-0.1, box2.width, 0.03])
plt.colorbar(b, cax = cbar_ax1, orientation="horizontal", extend="max").set_label(label=r"Blocking intensity $\sigma$ (-)", size=fontsize)
plt.colorbar(c, cax = cbar_ax2, orientation="horizontal", extend="max").set_label(label=r"Blocking frequency $\sigma$ (%)", size=fontsize)
plt.savefig("Ens_4_std_mod_DJF_JJA_fs="+str(fontsize)+".pdf", bbox_inches="tight")
plt.show()

print("Intramodel variation plotted")

#%% Plot z500 field: In first figure absolute values of Z500 and in the second figure their gradients
fig = plt.figure()
#plt.suptitle("z500_hist-z500_ERA, 1951-2014")
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Orthographic(central_latitude=90))
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax.coastlines(linewidth=0.5)
a = ax.pcolormesh(lon, lat[14:31], z500_hist_mean.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).mean("time") - z500_ERA.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).mean("time"), transform=ccrs.PlateCarree(), cmap=cmap_1, vmin=-50, vmax=50)
plt.colorbar(a, location="bottom", extend="both").set_label(label=r"$\Delta$ zg [m]", size=fontsize)
ax.set_title("Year", fontsize=fontsize)
ax1 = fig.add_subplot(1, 3, 2, projection=ccrs.Orthographic(central_latitude=90))
ax1.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax1.coastlines(linewidth=0.5)
a = ax1.pcolormesh(lon, lat[14:31], z500_hist_mean.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).groupby("time.season").mean("time").sel(season="DJF") - z500_ERA.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).groupby("time.season").mean("time").sel(season="DJF"), transform=ccrs.PlateCarree(), cmap=cmap_1, vmin=-50, vmax=50)
plt.colorbar(a, location="bottom", extend="both").set_label(label=r"$\Delta$ zg [m]", size=fontsize)
ax1.set_title("DJF", fontsize=fontsize)
ax2 = fig.add_subplot(1, 3, 3, projection=ccrs.Orthographic(central_latitude=90))
ax2.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax2.coastlines(linewidth=0.5)
a = ax2.pcolormesh(lon, lat[14:31], z500_hist_mean.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).groupby("time.season").mean("time").sel(season="JJA") - z500_ERA.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).groupby("time.season").mean("time").sel(season="JJA"), transform=ccrs.PlateCarree(), cmap=cmap_1, vmin=-50, vmax=50)
plt.colorbar(a, location="bottom", extend="both").set_label(label=r"$\Delta$ zg [m]", size=fontsize)
ax2.set_title("JJA", fontsize=fontsize)
plt.tight_layout()
plt.savefig("Diff_z500hist_era_fs="+str(fontsize)+".pdf", bbox_inches="tight")
plt.show()

#Show gradient
grad_z500_year = np.gradient(z500_hist_mean.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).mean("time") - z500_ERA.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).mean("time"))
grad_z500_DJF = np.gradient(z500_hist_mean.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).groupby("time.season").mean("time").sel(season="DJF") - z500_ERA.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).groupby("time.season").mean("time").sel(season="DJF"))
grad_z500_JJA = np.gradient(z500_hist_mean.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).groupby("time.season").mean("time").sel(season="JJA") - z500_ERA.sel(time=slice("01-01-1951", "31-12-2014"), lat=slice(35,75)).groupby("time.season").mean("time").sel(season="JJA"))

fig = plt.figure()
#plt.suptitle("gradient z500_hist-z500_ERA, 1951-2014")
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Orthographic(central_latitude=90))
ax.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax.coastlines(linewidth=0.5)
a = ax.pcolormesh(lon, lat[14:31], grad_z500_year[0], transform=ccrs.PlateCarree(), cmap=cmap_1, vmin=-20, vmax=20)
plt.colorbar(a, location="bottom", extend="both").set_label(label=r"$\Delta$ zg [m/deg]", size=fontsize)
ax.set_title("Year", fontsize=fontsize)
ax1 = fig.add_subplot(1, 3, 2, projection=ccrs.Orthographic(central_latitude=90))
ax1.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax1.coastlines(linewidth=0.5)
a = ax1.pcolormesh(lon, lat[14:31], grad_z500_DJF[0], transform=ccrs.PlateCarree(), cmap=cmap_1, vmin=-20, vmax=20)
plt.colorbar(a, location="bottom", extend="both").set_label(label=r"$\Delta$ zg [m/deg]", size=fontsize)
ax1.set_title("DJF", fontsize=fontsize)
ax2 = fig.add_subplot(1, 3, 3, projection=ccrs.Orthographic(central_latitude=90))
ax2.set_extent([180, -180, 35, 75], crs=ccrs.PlateCarree())
ax2.coastlines(linewidth=0.5)
a = ax2.pcolormesh(lon, lat[14:31], grad_z500_JJA[0], transform=ccrs.PlateCarree(), cmap=cmap_1, vmin=-20, vmax=20)
plt.colorbar(a, location="bottom", extend="both").set_label(label=r"$\Delta$ zg [m/deg]", size=fontsize)
ax2.set_title("JJA", fontsize=fontsize)
plt.tight_layout()
plt.savefig("Grad_z500_hist_era_fs="+str(fontsize)+".pdf", bbox_inches="tight")
plt.show()

print("Z500 fields plotted")
