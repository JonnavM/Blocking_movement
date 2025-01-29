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

dir = "Fill in directory"

#ECE3p5
hist_clarea = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_clarea_large_weighted.nc")
hist_clcmassX = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_clcmassX_large_weighted.nc")
hist_clcmassY = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_clcmassY_large_weighted.nc")
hist_BIav = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_BIav_large.nc")
hist_BImax = xr.open_dataarray(dir+"Blockings/per_lon/Intensity/long/large/hist_BI_BImax_large.nc")
hist_speedx = speedx(hist_clcmassX, hist_clcmassY)
hist_speedx = hist_speedx*1000/(24*3600)

#%% Make composit for the same intervals as for the temperature-size-velocity plot
#v_cond = np.logical_and(hist_speedx_season>=hist_speedx_season.quantile(perc_vel1[i]), hist_speedx_season<=hist_speedx_season.quantile(perc_vel2[i]))


z500_ECE = xr.open_mfdataset("/ECEARTH3p5/historical/r*i1p5f1/ensemble/zg_day_EC-Earth3_historical_r*i1p5f1.nc", combine="nested", concat_dim="ensemble").zg.persist()
BI_ECE = xr.open_mfdataset("/Blockings/Intensity/block_hist*_BI_1850_2014.nc", combine="nested", concat_dim="ensemble").zg.load()

z500_clim = z500_ECE.groupby("time.month").mean("time")
z500_anom = (z500_ECE.groupby("time.month")-z500_clim).load() #Dit duurt te lang. Via cdo uitrekenen en inladen

#%% 
season="DJF"

perc_vel1 = [0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95]
perc_vel2 = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.00]

v_cond = np.logical_and(hist_speedx.sel(time=hist_clcmassX["time.season"]==season)>=hist_speedx.sel(time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(perc_vel1[4]), 
                        hist_speedx.sel(time=hist_clcmassX["time.season"]==season)<=hist_speedx.sel(time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(perc_vel2[4]))

# Iterate over ensembles and blocks. This leads to separate figures for each block
for ensemble in range(v_cond.sizes["ensemble"]):
    for block in range(v_cond.sizes["block"]):
        # Get all time points where the condition is True
        valid_times = v_cond.sel(ensemble=ensemble, block=block).where(v_cond.sel(ensemble=ensemble, block=block), drop=True)

        for start_date in valid_times["time"].values:
            start_date = pd.Timestamp(start_date)
            
            # Get longitude, latitude, and length for this block
            lons = hist_clcmassX.sel(ensemble=ensemble, block=block, time=start_date).dropna("day")
            lats = hist_clcmassY.sel(ensemble=ensemble, block=block, time=start_date)
            length = int(lons.count().values)
            end_date = start_date + pd.Timedelta(days=length)

            # Extract z500 and BI data
            z500_sel = z500_anom.sel(ensemble=ensemble, time=slice(start_date, end_date), lat=slice(lats.mean() - 5, lats.mean() + 5)).mean("lat")
            BI_sel = BI_ECE.sel(ensemble=ensemble, time=slice(start_date, end_date), lat=slice(lats.mean() - 5, lats.mean() + 5)).mean("lat")

            # Plot
            
            fig, ax = plt.subplots(1, 2, sharey=True, figsize=(12, 6))
            z500_plot = z500_sel.plot(ax=ax[0], label="Z500 Anomaly", vmin=-240, vmax=240, cmap="coolwarm", levels=20, add_colorbar=False)
            cbar_z500 = fig.colorbar(z500_plot, ax=ax[0])
            cbar_z500.set_label(r"$\Delta$ z500 (m)")
            ax[0].scatter(lons, pd.date_range(start=start_date, end=end_date, freq="D")[:-1], color="black", label="Block Path")
            ax[0].legend()
            ax[0].set_title(f"Z500 Anomaly (Ensemble {ensemble})")

            BI_plot = BI_sel.plot(ax=ax[1], vmin=0, vmax=3, levels=20, label="BI", cmap="viridis", add_colorbar=False)
            cbar_BI = fig.colorbar(BI_plot, ax=ax[1])
            cbar_BI.set_label("BI (-)")
            ax[1].scatter(lons, pd.date_range(start=start_date, end=end_date, freq="D")[:-1], color="black")
            ax[1].set_title(f"Block Intensity (Ensemble {ensemble})")
            ax[1].set_ylabel(" ")
            plt.tight_layout()
            #fig.savefig("Hovmuller-P0-005_ens0_block0_1986-12-19.pdf")
            plt.show()

#%% Now as a composite
season="DJF"

#perc_vel1 = [0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95]
#perc_vel2 = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.00]
perc_vel1 = [0.0, 0.1, 0.9]
perc_vel2 = [0.1, 0.9, 1.00]
for i in range(len(perc_vel1)):
    v_cond = np.logical_and(hist_speedx.sel(time=hist_clcmassX["time.season"]==season)>=hist_speedx.sel(time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(perc_vel1[i]), 
                            hist_speedx.sel(time=hist_clcmassX["time.season"]==season)<=hist_speedx.sel(time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(perc_vel2[i]))
    
    # Initialize lists to store composite data
    composite_z500 = []
    composite_BI = []
    day_counts = []
    lons_blocks = []
    
    original_lon_grid = np.arange(-180, 180, 2.5)
    
    def map_to_grid(lon_values, grid):
        """Map lon_values to the nearest values in the grid."""
        return grid[np.abs(grid[:, None] - lon_values).argmin(axis=0)]
    
    
    # Iterate over ensembles and blocks
    for ensemble in range(v_cond.sizes["ensemble"]):
        for block in range(v_cond.sizes["block"]):
            # Get all time points where the condition is True
            valid_times = v_cond.sel(ensemble=ensemble, block=block).where(v_cond.sel(ensemble=ensemble, block=block), drop=True)
    
            for start_date in valid_times["time"].values:
                start_date = pd.Timestamp(start_date)
    
                # Get longitude, latitude, and length for this block
                lons = hist_clcmassX.sel(ensemble=ensemble, block=block, time=start_date).dropna("day")
                lats = hist_clcmassY.sel(ensemble=ensemble, block=block, time=start_date)
                length = int(lons.count().values)
                end_date = start_date + pd.Timedelta(days=length)
    
                # Extract z500 and BI data
                z500_sel = z500_anom.sel(ensemble=ensemble, time=slice(start_date, end_date), lat=slice(lats.mean() - 5, lats.mean() + 5)).mean("lat")
                BI_sel = BI_ECE.sel(ensemble=ensemble, time=slice(start_date, end_date), lat=slice(lats.mean() - 5, lats.mean() + 5)).mean("lat")
    
                # Center longitudes
                block_center_lon = lons.values[0]
                lon_anomaly = (lons - block_center_lon + 360) % 360
                lon_anomaly = lon_anomaly.where(lon_anomaly <= 180, lon_anomaly - 360)
                lon_anomaly_da = xr.DataArray(lon_anomaly.values, coords=[('time', lon_anomaly["day"].values)], name=lon_anomaly.name)
                lons_blocks.append(lon_anomaly_da)
                
                # Shift longitudes in z500_sel to center at block_center_lon
                shifted_lon = (z500_sel.lon - block_center_lon + 360) % 360  # Ensure wrapping to [0, 360]
                shifted_lon = shifted_lon.where(shifted_lon <= 180, shifted_lon - 360)  # Shift back to [-180, 180]
                adjusted_lon = map_to_grid(shifted_lon.values, original_lon_grid)
    
                # Resample to daily indices
                days = np.arange(len(z500_sel["time"]))
                #composite_z500.append(z500_sel.assign_coords(lon=adjusted_lon, time=days))
                #composite_BI.append(BI_sel.assign_coords(lon=adjusted_lon, time=days))
                # Group and average data by the new longitudes to align to the grid
                if (adjusted_lon == original_lon_grid).all():
                    z500_block = z500_sel.assign_coords(lon=adjusted_lon, time=days)
                    BI_block = BI_sel.assign_coords(lon=adjusted_lon, time=days)
                else:
                    z500_block = z500_sel.assign_coords(lon=adjusted_lon, time=days).groupby("lon").mean()
                    BI_block = BI_sel.assign_coords(lon=adjusted_lon, time=days).groupby("lon").mean()
    
                composite_z500.append(z500_block)
                composite_BI.append(BI_block)
    
                day_counts.append(days)
    
    # Combine the aligned blocks into final composites
    composite_z500_all = xr.concat(composite_z500, dim="block").mean(dim=("block"))
    composite_BI_all = xr.concat(composite_BI, dim="block").mean("block")
    lons_all = xr.concat(lons_blocks, dim="block").mean("block")
    
    # Plot the composite figure
    #titles = ["P0-0.05 of velocity in JJA", "P0.05-0.1 of velocity in JJA", "P0.1-0.2 of velocity in JJA", "P0.2-0.4 of velocity in JJA", 
    #          "P0.4-0.6 of velocity in JJA", "P0.6-0.8 of velocity in JJA", "P0.8-0.9 of velocity in JJA", "P0.9-0.95 of velocity in JJA", "P0.95-1 of velocity in JJA"]
    #percentiles = ["0-005_JJA", "005-01_JJA", "01-02_JJA", "02-04_JJA", "04-06_JJA", "06-08_JJA", "08-09_JJA", "09-095_JJA", "095-1_JJA"]
    titles = ["P0-0.1 of velocity in DJF", "P0.1-0.9 of velocity in DJF", "P0.9-1 of velocity in DJF"]
    percentiles = ["0-01_DJF", "01-09_DJF", "09-1_DJF"]
    fig, ax = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    fig.suptitle(titles[i])
    
    # Plot Z500 anomaly composite
    z500_plot = composite_z500_all.plot.contourf(
        ax=ax[0], x="lon", y="time", levels=20, cmap="coolwarm", extend="both", vmin=-240, vmax=240, add_colorbar=False
    )
    cbar_z500 = fig.colorbar(z500_plot, ax=ax[0])
    cbar_z500.set_label(r"$\Delta$ z500 (m)")
    
    ax[0].scatter(lons_all.values[:-1], lons_all["time"][:-1], marker="o", color="black")
    ax[0].set_xlabel(r"$\Delta$ Longitude (°)")
    ax[0].set_ylabel("Day within Block")
    
    # Plot Block Intensity (BI) composite
    BI_plot = composite_BI_all.plot.contourf(
        ax=ax[1], x="lon", y="time", levels=20, cmap="viridis", extend="both", vmin=0, vmax=3.5, add_colorbar=False
    )
    cbar_BI = fig.colorbar(BI_plot, ax=ax[1])
    cbar_BI.set_label("BI (-)")
    
    ax[1].scatter(lons_all.values[:-1], lons_all["time"][:-1], marker="o", color="black")
    ax[1].set_xlabel(r"$\Delta$ Longitude (°)")
    ax[1].set_ylabel(" ")
    
    plt.tight_layout()
    plt.show()
    fig.savefig("/eejit/home/6196306/Data/MSc_thesis/Hovmuller-"+str(percentiles[i])+".pdf")
    
    #With some adjustments
    fig, ax = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    fig.suptitle(titles[i])
    
    # Plot Z500 anomaly composite
    z500_plot = composite_z500_all.plot.contourf(
        ax=ax[0], x="lon", y="time", levels=20, cmap="coolwarm", extend="both", vmin=-240, vmax=240, add_colorbar=False
    )
    cbar_z500 = fig.colorbar(z500_plot, ax=ax[0])
    cbar_z500.set_label(r"$\Delta$ z500 (m)")
    
    ax[0].scatter(lons_all.values[:-1], lons_all["time"][:-1], marker="o", color="black")
    ax[0].set_xlabel(r"$\Delta$ Longitude (°)")
    ax[0].set_ylabel("Day within Block")
    ax[0].set_ylim(0,20)
    # Plot Block Intensity (BI) composite
    BI_plot = composite_BI_all.plot.contourf(
        ax=ax[1], x="lon", y="time", levels=20, cmap="viridis", extend="both", vmin=0, vmax=3, add_colorbar=False
    )
    cbar_BI = fig.colorbar(BI_plot, ax=ax[1])
    cbar_BI.set_label("BI (-)")
    
    ax[1].scatter(lons_all.values[:-1], lons_all["time"][:-1], marker="o", color="black")
    ax[1].set_xlabel(r"$\Delta$ Longitude (°)")
    ax[1].set_ylabel(" ")
    ax[1].set_ylim(0,20)
    
    plt.tight_layout()
    plt.show()
    fig.savefig("Hovmuller-"+str(percentiles[i])+"_v2.pdf")
    
 
