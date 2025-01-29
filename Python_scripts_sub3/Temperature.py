#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 14:49:36 2022

@author: mourik
Template for all figures surrounding temperatures
"""
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cartopy.crs as ccrs
from matplotlib.colors import Normalize


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
hist_clarea_dmean = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Blockings/per_lon/Intensity/long/large/hist_BI_clarea_large_weighted_daymean.nc")

tas["time"] = tas.indexes["time"].normalize()
hist_speedx = speedx(hist_clcmassX, hist_clcmassY)
hist_speedx = hist_speedx*1000/(24*3600)

#Import temp field around blocks, either with or without oceans. This can only be done if you already produced these files in "save temperature field for different blocks"
T_block_xr = xr.open_mfdataset("/scratch/depfg/6196306/MSc_thesis/Temperature/Tblock/20x40/varlat/Tas_block_hist*_20x40_land.nc", combine="nested", concat_dim="ensemble").Temp_anom
z500_block_xr = xr.open_mfdataset("/scratch/depfg/6196306/MSc_thesis/ECEARTH3p5/historical/z500block/z500_block_hist*_20x40.nc", combine="nested", concat_dim="ensemble").z500
#Import different quadrants, either with or without oceans. This can only be done if you already produced these files in "Create quadrants"
UL = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Temperature/UL_varlat_land.nc")
UR = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Temperature/UR_varlat_land.nc")
LL = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Temperature/LL_varlat_land.nc")
LR = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Temperature/LR_varlat_land.nc")


#%% Create quadrants
#UL = np.zeros((len(hist_BIav["ensemble"]), len(hist_BIav["time"]), len(hist_BIav["day"]), len(hist_BIav["block"])))
#UR = UL
#LL = UL
#LR = UL

#UL = np.load("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/UL_varlat_land.npy")
#UR = np.load("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/UR_varlat_land.npy")
#LL = np.load("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/LL_varlat_land.npy")
#LR = np.load("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/LR_varlat_land.npy")


#i=15
#tas["time"] = tas.indexes["time"].normalize()
#for t in range(len(majlen["time"])):
#    for b in range(len(majlen["block"])):
#        for d in range(len(majlen["day"])):
#            if majlen.sel(time=majlen["time"][t], day=majlen["day"][d], block=majlen["block"][b]).notnull():
#                Time = pd.to_datetime(str(majlen.indexes["time"][t])).strftime("%Y-%m-%d")
#                day = majlen["day"][d]
#                block = majlen["block"][b]
#                X = hist_clcmassX.sel(time=Time, day=day, block=block)
#                Y = hist_clcmassY.sel(time=Time, day=day, block=block)
#                del_x = 40 #majlen.sel(time=Time, day=day, block=block)
#                del_y = 20 #minlen.sel(time=Time, day=day, block=block)
#                UL[i, t, d, b] = tas.sel(time=Time, lat=slice(Y, Y+del_y), lon=slice(X-del_x, X)).mean().values
#                UR[i, t, d, b] = tas.sel(time=Time, lat=slice(Y, Y+del_y), lon=slice(X, X+del_x)).mean().values
#                LL[i, t, d, b] = tas.sel(time=Time, lat=slice(Y-del_y, Y), lon=slice(X-del_x, X)).mean().values
#                LR[i, t, d, b] = tas.sel(time=Time, lat=slice(Y-del_y, Y), lon=slice(X, X+del_x)).mean().values
#            else:
#                break

#Save outputs to work on later
#np.save("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/UL_varlat_land.npy", UL)
#np.save("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/UR_varlat_land.npy", UR)
#np.save("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/LL_varlat_land.npy", LL)
#np.save("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/LR_varlat_land.npy", LR)

# Change outputs to xarray datasets
#UL = xr.Dataset(data_vars=dict(upperleft=(["ensemble", "time", "day", "block"], UL[:,:,:,:])), coords=dict(time=hist_clcmassX["time"])).upperleft
#UR = xr.Dataset(data_vars=dict(upperright=(["ensemble", "time", "day", "block"], UR[:,:,:,:])), coords=dict(time=hist_clcmassX["time"])).upperright
#LL = xr.Dataset(data_vars=dict(lowerleft=(["ensemble", "time", "day", "block"], LL[:,:,:,:])), coords=dict(time=hist_clcmassX["time"])).lowerleft
#LR = xr.Dataset(data_vars=dict(lowerright=(["ensemble", "time", "day", "block"], LR[:,:,:,:])), coords=dict(time=hist_clcmassX["time"])).lowerright

#UL = UL.where(UL!=0, np.nan)
#UR = UR.where(UR!=0, np.nan)
#LL = LL.where(LL!=0, np.nan)
#LR = LR.where(LR!=0, np.nan)

#Save as xarray files
#UL.to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/UL_varlat_land.nc")
#UR.to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/UR_varlat_land.nc")
#LL.to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/LL_varlat_land.nc")
#LR.to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/LR_varlat_land.nc")


#%% Save temperature field for different blocks 20x40
ens_str = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
ens = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
for i in range(len(ens)):
    #tas = xr.open_dataset("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/ECEARTH3Bis/ensembles/tas/anomaly/land/tas_hist"+str(ens_str[i])+"_anom_land.nc").tas
    z500 = xr.open_dataset("/scratch/depfg/6196306/MSc_thesis/ECEARTH3p5/historical/r"+str(ens[i]+1)+"i1p5f1/ensemble/zg_day_EC-Earth3_historical_r"+str(ens[i]+1)+"i1p5f1.nc").zg
    #tas["time"] = tas.indexes["time"].normalize()
    z500["time"] = z500.indexes["time"].normalize()
    T_block_shape = np.zeros([59653, 5, 17, 33])*np.nan
    for t in range(len(majlen["time"])):
        for b in range(len(majlen["block"])):
            if majlen.sel(ensemble=ens[i], time=majlen["time"][t], day=majlen["day"][0], block=majlen["block"][b]).notnull():
                Time = pd.to_datetime(str(majlen.indexes["time"][t])).strftime("%Y-%m-%d")
                block = majlen["block"][b]
                T_1day = []
                for d in range(len(majlen["day"])):
                    X = hist_clcmassX.sel(ensemble=ens[i], time=Time, block=block, day=d)
                    Y = hist_clcmassY.sel(ensemble=ens[i], time=Time, block=block, day=d)
                    del_x = 40 
                    del_y = 20 
                    X_min = X-del_x
                    X_plus = X+del_x
                    lon_min = X_min if X_min>=-180 else X_min+360
                    lon_max = X_plus if X_plus<=180 else X_plus-360
                    if lon_max>lon_min:
                        #T = tas.sel(time=Time, lat=slice(Y-del_y, Y+del_y), lon=slice(lon_min, lon_max))
                        T = z500.sel(time=Time, lat=slice(Y-del_y, Y+del_y), lon=slice(lon_min, lon_max))
                    elif lon_max<lon_min:
                        #condition1 = tas.lon.where(tas.lon<lon_max, drop=True)
                        condition1 = z500.lon.where(z500.lon<lon_max, drop=True)
                        #condition2 = tas.lon.where(tas.lon>lon_min, drop=True)
                        condition2 = z500.lon.where(z500.lon>lon_min, drop=True)
                        condition = np.concatenate([condition2.values, condition1.values])
                        #T = tas.sel(time=Time, lat=slice(Y-del_y, Y+del_y), lon=condition)    
                        T = z500.sel(time=Time, lat=slice(Y-del_y, Y+del_y), lon=condition)
                    if len(T["lat"])==16:
                        T["lat"] = np.arange(-20, 20, 2.5)
                    elif len(T["lat"])==17:
                        T["lat"] = np.arange(-20, 22.5, 2.5)
                    elif len(T["lat"])==15:
                        T["lat"] = np.arange(-17.5, 20, 2.5)
                    if len(T["lon"])==32:
                        T["lon"] = np.arange(-40, 40, 2.5)
                    elif len(T["lon"])==33:
                        T["lon"] = np.arange(-40, 42.5, 2.5)
                    elif len(T["lon"])==31:
                        T["lon"] = np.arange(-37.5, 40, 2.5)
                    T_1day.append(T)
                    if X.isnull()==True:
                        break
                    else:
                        continue
                T_daymean = xr.concat(T_1day, "day").mean("day", skipna=True)
                for lat in range(len(T_daymean["lat"])):
                    for lon in range(len(T_daymean["lon"])):
                        T_block_shape[t, b, lat, lon] = T_daymean.sel(lat=T_daymean["lat"][lat], lon=T_daymean["lon"][lon])     
            else:
                break
    #T_block_xr = xr.Dataset(data_vars=dict(Temp_anom=(["time", "block", "lat", "lon"], T_block_shape)), coords=dict(time=hist_clcmassX.sel(ensemble=0)["time"], lon=np.arange(-40, 42.5, 2.5), lat=np.arange(-20, 22.5, 2.5))).Temp_anom
    #T_block_xr.to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/Tblock/20x40/varlat/Tas_block_hist"+str(ens_str[i])+"_20x40_land.nc")
    z500_block_xr = xr.Dataset(data_vars=dict(z500=(["time", "block", "lat", "lon"], T_block_shape)), coords=dict(time=hist_clcmassX.sel(ensemble=0)["time"], lon=np.arange(-40, 42.5, 2.5), lat=np.arange(-20, 22.5, 2.5))).z500
    #z500_block_xr.to_netcdf("/scratch/depfg/6196306/MSc_thesis/ECEARTH3p5/historical/z500block/z500_block_hist"+str(ens_str[i])+"_20x40.nc")

#%% Percentile and KDE plot for temperature against velocity per season for the quadrants
title = ["Upperleft", "Upperright", "Lowerleft", "Lowerright"]
comp = [UL, UR, LL, LR]

#Winter values
temp_UL_sw = UL.sel(ensemble=slice(0,15), time=UL["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")
temp_UR_sw = UR.sel(ensemble=slice(0,15), time=UR["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")
temp_LL_sw = LL.sel(ensemble=slice(0,15), time=LL["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")
temp_LR_sw = LR.sel(ensemble=slice(0,15), time=LR["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")

speed_w = hist_speedx.sel(ensemble=slice(0,15), time=hist_speedx["time.season"]=="DJF").stack(z=("ensemble", "block", "time"))#.dropna("z")

temp_UL_matrix_sw = np.zeros([3, 100]) #Save 10, 50 and 90 percentiles
temp_UR_matrix_sw = np.zeros([3, 100]) #Save 10, 50 and 90 percentiles
temp_LL_matrix_sw = np.zeros([3, 100]) #Save 10, 50 and 90 percentiles
temp_LR_matrix_sw = np.zeros([3, 100]) #Save 10, 50 and 90 percentiles

#summer values
temp_UL_ss = UL.sel(ensemble=slice(0,15), time=UL["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")
temp_UR_ss = UR.sel(ensemble=slice(0,15), time=UR["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")
temp_LL_ss = LL.sel(ensemble=slice(0,15), time=LL["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")
temp_LR_ss = LR.sel(ensemble=slice(0,15), time=LR["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")

speed_s = hist_speedx.sel(ensemble=slice(0,15), time=hist_speedx["time.season"]=="JJA").stack(z=("ensemble", "block", "time"))#.dropna("z")

#Create a matrix to save the percentile values in
temp_UL_matrix_ss = np.zeros([3, 100]) #Save 10, 50 and 90 percentiles
temp_UR_matrix_ss = np.zeros([3, 100]) #Save 10, 50 and 90 percentiles
temp_LL_matrix_ss = np.zeros([3, 100]) #Save 10, 50 and 90 percentiles
temp_LR_matrix_ss = np.zeros([3, 100]) #Save 10, 50 and 90 percentiles

#Define the steps in velocity
step = 40*1000/(24*3600)
speed_40 = np.arange(-1960*1000/(24*3600), 2040*1000/(24*3600), step)
percentiles = [10, 50, 90]
temp_mat_w = [temp_UL_matrix_sw, temp_UR_matrix_sw, temp_LL_matrix_sw, temp_LR_matrix_sw]
temp_mat_s = [temp_UL_matrix_ss, temp_UR_matrix_ss, temp_LL_matrix_ss, temp_LR_matrix_ss]

temp_sw = [temp_UL_sw, temp_UR_sw, temp_LL_sw, temp_LR_sw]
temp_ss = [temp_UL_ss, temp_UR_ss, temp_LL_ss, temp_LR_ss]

#Loop through the percentiles and velocity bins
#For winter
for p in range(len(percentiles)):
    for q in range(4):
        for i in range(100):
            if temp_sw[q][np.logical_and(speed_w>-2000*1000/(24*3600)+i*step, speed_w<-2000*1000/(24*3600)+(i+1)*step)].size==0:
                temp_mat_w[q][p,i]=np.nan
            else:
                temp_mat_w[q][p, i] = np.percentile(temp_sw[q][np.logical_and(speed_w>-2000*1000/(24*3600)+i*step, speed_w<-2000*1000/(24*3600)+(i+1)*step)].dropna("z"), percentiles[p])
#And summer
for p in range(len(percentiles)):
    for q in range(4):
        for i in range(100):
            if temp_ss[q][np.logical_and(speed_s>-2000*1000/(24*3600)+i*step, speed_s<-2000*1000/(24*3600)+(i+1)*step)].size==0:
                temp_mat_s[q][p,i]=np.nan
            else:
                temp_mat_s[q][p, i] = np.percentile(temp_ss[q][np.logical_and(speed_s>-2000*1000/(24*3600)+i*step, speed_s<-2000*1000/(24*3600)+(i+1)*step)].dropna("z"), percentiles[p])

#Change the matrices to dataarrays                
temp_UL_matrix_sw = xr.DataArray(temp_UL_matrix_sw, name="temperature", coords=[("percentiles", percentiles), ("speed", speed_40)], dims=["percentiles", "speed"])        
temp_UR_matrix_sw = xr.DataArray(temp_UR_matrix_sw, name="temperature", coords=[("percentiles", percentiles), ("speed", speed_40)], dims=["percentiles", "speed"])        
temp_LL_matrix_sw = xr.DataArray(temp_LL_matrix_sw, name="temperature", coords=[("percentiles", percentiles), ("speed", speed_40)], dims=["percentiles", "speed"])        
temp_LR_matrix_sw = xr.DataArray(temp_LR_matrix_sw, name="temperature", coords=[("percentiles", percentiles), ("speed", speed_40)], dims=["percentiles", "speed"])        

temp_UL_matrix_ss = xr.DataArray(temp_UL_matrix_ss, name="temperature", coords=[("percentiles", percentiles), ("speed", speed_40)], dims=["percentiles", "speed"])        
temp_UR_matrix_ss = xr.DataArray(temp_UR_matrix_ss, name="temperature", coords=[("percentiles", percentiles), ("speed", speed_40)], dims=["percentiles", "speed"])        
temp_LL_matrix_ss = xr.DataArray(temp_LL_matrix_ss, name="temperature", coords=[("percentiles", percentiles), ("speed", speed_40)], dims=["percentiles", "speed"])        
temp_LR_matrix_ss = xr.DataArray(temp_LR_matrix_ss, name="temperature", coords=[("percentiles", percentiles), ("speed", speed_40)], dims=["percentiles", "speed"])

temp_matrix_sw = [temp_UL_matrix_sw, temp_UR_matrix_sw, temp_LL_matrix_sw, temp_LR_matrix_sw]
temp_matrix_ss = [temp_UL_matrix_ss, temp_UR_matrix_ss, temp_LL_matrix_ss, temp_LR_matrix_ss]

xw = hist_speedx.sel(ensemble=slice(0,15), time=hist_speedx["time.season"]=="DJF").stack(z=("ensemble", "block", "time"))#.dropna("z")
#yw = LR.sel(ensemble=slice(0,15), time=LR["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")

xs = hist_speedx.sel(ensemble=slice(0,15), time=hist_speedx["time.season"]=="JJA").stack(z=("ensemble", "block", "time"))#.dropna("z")
#ys = LR.sel(ensemble=slice(0,15), time=LR["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")

from matplotlib import ticker
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 

#Plot percentiles for size
for i in range(4):
    yw = comp[i].sel(ensemble=slice(0,15), time=comp[i]["time.season"]=="DJF").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")
    ys = comp[i].sel(ensemble=slice(0,15), time=comp[i]["time.season"]=="JJA").mean("day").stack(z=("ensemble", "block", "time"))#.dropna("z")

    fig, ax = plt.subplots(1,2, figsize=(10,5), sharey=True)
    plt.subplots_adjust(wspace=0.12)
    cbar_ax = fig.add_axes([0.92, 0.124, .03, .76])
    ax[0].set_title("DJF")
    ax[0].set_xlabel("Velocity (m/s)")
    ax[0].set_ylabel(r"$\Delta T$ (K)")
    ax[0].set_ylim(-10, 10)
    ax[0].set_xlim(-1500*1000/(24*3600), 1500*1000/(24*3600))
    ax[0].fill_betweenx(x1=-1500*1000/(24*3600), x2=xw.quantile(0.1), y=[-10,10], color="grey", alpha=0.1)
    ax[0].fill_betweenx(x1=xw.quantile(0.9), x2=1500*1000/(24*3600), y=[-10,10], color="grey", alpha=0.1)
    ax[0].plot(speed_40, temp_matrix_sw[i].sel(percentiles=10), label="10%")
    ax[0].plot(speed_40, temp_matrix_sw[i].sel(percentiles=50), label="50%")
    ax[0].plot(speed_40, temp_matrix_sw[i].sel(percentiles=90), label="90%")
    a = sns.kdeplot(ax=ax[0], x=xw, y=yw, cmap="Blues", fill=True, cbar=False, levels=10, norm=Normalize(vmax=0.03, vmin=0.00126))
    ax[0].hlines(y=0, xmin=-2000*1000/(24*3600), xmax=2000*1000/(24*3600), color="black", alpha=0.5, linestyle="--")
    
    ax[1].set_title("JJA")
    ax[1].set_xlabel("Velocity (m/s)")
    ax[1].set_ylim(-10, 10)
    ax[1].set_xlim(-1500*1000/(24*3600), 1500*1000/(24*3600))
    ax[1].fill_betweenx(x1=-1500*1000/(24*3600), x2=xs.quantile(0.1), y=[-10,10], color="grey", alpha=0.1)
    ax[1].fill_betweenx(x1=xs.quantile(0.9), x2=1500*1000/(24*3600), y=[-10,10], color="grey", alpha=0.1)
    ax[1].text(s=r"$v_x\leq P10$", x=xs.quantile(0.1), y=-9, color="grey", horizontalalignment="right")
    ax[1].text(s=r"$v_x\geq P90$", x=xs.quantile(0.9), y=-9, color="grey")
    ax[1].plot(speed_40, temp_matrix_ss[i].sel(percentiles=10), label="10%")
    ax[1].plot(speed_40, temp_matrix_ss[i].sel(percentiles=50), label="50%")
    ax[1].plot(speed_40, temp_matrix_ss[i].sel(percentiles=90), label="90%")
    a = sns.kdeplot(ax=ax[1], x=xs, y=ys, cmap="Blues", fill=True, cbar=True, cbar_ax=cbar_ax, cbar_kws={"label":"Frequency", "format":formatter}, levels=10, norm=Normalize(vmax=0.03, vmin=0.00126))
    ax[1].legend()
    ax[1].hlines(y=0, xmin=-2000*1000/(24*3600), xmax=2000*1000/(24*3600), color="black", alpha=0.5, linestyle="--")
    fig.savefig("/scratch/depfg/6196306/MSc_thesis/Final_Figures/review/percentiles_temp_speed_P10_P90_"+str(title[i])+"_ms.pdf", bbox_inches='tight')
    plt.show()

#%% Single plot for composite: mean values
season="JJA"
#Both middle for velocity as size
v_cond = np.logical_and(hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")>=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0.4), hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")<=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0.6))
a_cond = np.logical_and(hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))>=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(0.4), hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))<=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(0.6))
plt.figure()
T_block_xr.sel(ensemble=slice(0,16), lat=slice(-20,19), lon=slice(-40,39)).where(np.logical_and(v_cond==True, a_cond==True)).mean(dim=("time", "block", "ensemble")).plot(cmap="RdBu_r", vmin=-2.5, vmax=2.5, cbar_kwargs={"label":r"$\Delta T$ (K)"})
plt.title(season + r", $0.4 \leq $P$ \leq 0.6$")
plt.xlabel(r"Relative longitude ($^\circ$)")
plt.ylabel("Relative latitude ($^\circ$)")
plt.savefig("/usr/people/mourik/Documents/Python/Figures/varlat/composit_vel_size_"+str(season)+"_P40-60_20x40_land.pdf")
plt.show()

#%% Composite figure of the temperature for different velocities and sizes

season="JJA"
perc = [0.0, 0.2, 0.4, 0.6, 0.8]
fig = plt.figure()
gs = fig.add_gridspec(5,5, hspace=0, wspace=0)
ax = gs.subplots(sharex=True, sharey=True)
fig.suptitle(season)

for i in range(5):
    for j in range(5):
        v_cond = np.logical_and(hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")>=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0 + perc[i]), hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")<=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0.2 + perc[i]))
        a_cond = np.logical_and(hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))>=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(0.8 - perc[j]), hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))<=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(1.0 - perc[j]))
        im = ax[j,i].pcolormesh(T_block_xr.sel(ensemble=slice(0,16), lat=slice(-20,19), lon=slice(-40,39)).where(np.logical_and(v_cond, a_cond)).mean(dim=("time", "block", "ensemble")), vmin=-2.5, vmax=2.5, cmap="RdBu_r")
        print(i,j)
for ax in ax.flatten():
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.76])
fig.colorbar(im, cax=cbar_ax, label=r"$\Delta T$ (K)")

fig.supylabel(r"Percentile size")
fig.supxlabel("Percentile velocity")
fig.text(0.12 , 0.08, "0%             20%           40%             60%             80%           100%")
fig.text(0.10 , 0.1, "0%         20%          40%          60%          80%          100%", rotation="vertical")
#fig.savefig("/scratch/depfg/6196306/MSc_thesis/varlat/composit_size_vel_"+str(season)+"_20x40_land_v2.pdf")
fig.savefig("/scratch/depfg/6196306/MSc_thesis/Final_Figures/review/composit_size_vel_"+str(season)+"_20x40_land_v3.pdf")

#%% Plot relative to mean state
fig = plt.figure()
gs = fig.add_gridspec(5,5, hspace=0, wspace=0)
ax = gs.subplots(sharex=True, sharey=True)
fig.suptitle(season)

v_mean = np.logical_and(hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")>=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0 + perc[2]), hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")<=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0.2 + perc[2]))
a_mean = np.logical_and(hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))>=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(0 + perc[2]), hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))<=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(0.2 + perc[2]))
block_mean = T_block_xr.sel(ensemble=slice(0,16), lat=slice(-20,19), lon=slice(-40,39)).where(np.logical_and(v_mean, a_mean)).mean(dim=("time", "block", "ensemble"))

for i in range(5):
    for j in range(5):
        if i==2 and j==2:
            im = ax[i,j].pcolormesh(T_block_xr.sel(ensemble=slice(0,16), lat=slice(-20,19), lon=slice(-40,39)).where(np.logical_and(v_cond, a_cond)).mean(dim=("time", "block", "ensemble")), vmin=-2.5, vmax=2.5, cmap="RdBu_r")
        else:
            v_cond = np.logical_and(hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")>=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0 + perc[i]), hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")<=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(0.2 + perc[i]))
            a_cond = np.logical_and(hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))>=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(0.8 - perc[j]), hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))<=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(1.0 - perc[j]))
            im = ax[j,i].pcolormesh(T_block_xr.sel(ensemble=slice(0,16), lat=slice(-20,19), lon=slice(-40,39)).where(np.logical_and(v_cond, a_cond)).mean(dim=("time", "block", "ensemble"))-block_mean, vmin=-2.5, vmax=2.5, cmap="RdBu_r")
        print(i,j)
for ax in ax.flatten():
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.76])
fig.colorbar(im, cax=cbar_ax, label=r"$\Delta T$ (K)")

fig.supylabel(r"Percentile size")
fig.supxlabel("Percentile velocity")
fig.text(0.12 , 0.08, "0%             20%           40%             60%             80%           100%")
fig.text(0.10 , 0.1, "0%         20%          40%          60%          80%          100%", rotation="vertical")
#fig.savefig("/scratch/depfg/6196306/MSc_thesis/varlat/composit_size_vel_"+str(season)+"_20x40_relative_land_v2.pdf")
fig.savefig("/scratch/depfg/6196306/MSc_thesis/Final_Figures/review/composit_size_vel_"+str(season)+"_20x40_relative_land_v3.pdf")

#%% Same but with different percentiles
season="JJA"
perc_vel1 = [0.0, 0.05, 0.1, 0.9, 0.95]
perc_vel2 = [0.05, 0.1, 0.9, 0.95, 1.00]
perc_size = [0.0, 0.2, 0.4, 0.6, 0.8]
fig = plt.figure()
gs = fig.add_gridspec(5,5, hspace=0, wspace=0)
ax = gs.subplots(sharex=True, sharey=True)
fig.suptitle(season)

for i in range(5):
    for j in range(5):
        v_cond = np.logical_and(hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")>=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(perc_vel1[i]), hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble")<=hist_speedx.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean("ensemble").quantile(perc_vel2[i]))
        a_cond = np.logical_and(hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))>=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(0.8 - perc_size[j]), hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble"))<=hist_clarea.sel(ensemble=slice(0,16), time=hist_clcmassX["time.season"]==season).mean(dim=("day", "ensemble")).quantile(1.0 - perc_size[j]))
        im = ax[j,i].pcolormesh(T_block_xr.sel(ensemble=slice(0,16), lat=slice(-20,19), lon=slice(-40,39)).where(np.logical_and(v_cond, a_cond)).mean(dim=("time", "block", "ensemble")), vmin=-2.5, vmax=2.5, cmap="RdBu_r")
        print(i,j)
for ax in ax.flatten():
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.76])
fig.colorbar(im, cax=cbar_ax, label=r"$\Delta T$ (K)")

fig.supylabel(r"Percentile size")
fig.supxlabel("Percentile velocity")
fig.text(0.12 , 0.08, "0%             5%           10%             90%             95%           100%")
fig.text(0.10 , 0.1, "0%         20%         40%         60%         80%         100%", rotation="vertical")
fig.savefig("/scratch/depfg/6196306/MSc_thesis/Final_Figures/review/composit_size_vel_"+str(season)+"_20x40_land_v3.pdf")

#%% Also try with z500 background
# 52:27-1:04:10 = ~12 minutes
# 1:08:00-
season="JJA"
perc_vel1 = [0.0, 0.05, 0.1, 0.9, 0.95]
perc_vel2 = [0.05, 0.1, 0.9, 0.95, 1.00]
perc_size = [0.0, 0.2, 0.4, 0.6, 0.8]

hist_speedx_season = hist_speedx.sel(time=hist_clcmassX["time.season"]==season)#.mean("ensemble")
hist_clarea_season = hist_clarea_dmean.sel(time=hist_clcmassX["time.season"]==season)#.mean("ensemble")
T_block_sel = T_block_xr.sel(lat=slice(-20,19), lon=slice(-40,39), time=hist_clcmassX["time.season"]==season)
z500_block_sel = z500_block_xr.sel(lat=slice(-20,19), lon=slice(-40,39), time=hist_clcmassX["time.season"]==season)

for i in range(5):
    for j in range(5):
        v_cond = np.logical_and(hist_speedx_season>=hist_speedx_season.quantile(perc_vel1[i]), hist_speedx_season<=hist_speedx_season.quantile(perc_vel2[i]))
        print("v")
        a_cond = np.logical_and(hist_clarea_season>=hist_clarea_season.quantile(0.8 - perc_size[j]), hist_clarea_season<=hist_clarea_season.quantile(1.0 - perc_size[j]))
        print("a")
        T_block_save = T_block_sel.where(np.logical_and(v_cond, a_cond)).mean(dim=("time", "block", "ensemble"))
        T_block_save.to_netcdf("/scratch/depfg/6196306/MSc_thesis/Temperature/Tblock/20x40/varlat/T_block_5-10-90-95_mean_"+str(season)+"_"+str(j)+"_"+str(i)+".nc")
        print("T")
        z500_block_save = z500_block_sel.where(np.logical_and(v_cond, a_cond)).mean(dim=("time", "block", "ensemble"))
        z500_block_save.to_netcdf("/scratch/depfg/6196306/MSc_thesis/ECEARTH3p5/historical/z500block/z500_block_5-10-90-95_mean_"+str(season)+"_"+str(j)+"_"+str(i)+".nc")
        print(i,j)

#%% Plot figure
fig = plt.figure()
gs = fig.add_gridspec(5,5, hspace=0, wspace=0)
ax = gs.subplots(sharex=True, sharey=True)
fig.suptitle(season)

for i in range(5):
    for j in range(5):
        T_block_ij = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/Temperature/Tblock/20x40/varlat/T_block_5-10-90-95_mean_"+str(season)+"_"+str(j)+"_"+str(i)+".nc")
        im = ax[j,i].pcolormesh(T_block_ij, vmin=-2.5, vmax=2.5, cmap="RdBu_r")
        z500_block_ij = xr.open_dataarray("/scratch/depfg/6196306/MSc_thesis/ECEARTH3p5/historical/z500block/z500_block_5-10-90-95_mean_"+str(season)+"_"+str(j)+"_"+str(i)+".nc")
        im2 = ax[j,i].contour(z500_block_ij, cmap="Greys", levels=4, linewidths=1)
    
for ax in ax.flatten():
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.76])
fig.colorbar(im, cax=cbar_ax, label=r"$\Delta T$ (K)")
fig.supylabel(r"Percentile size")
fig.supxlabel("Percentile velocity")
fig.text(0.12 , 0.08, "0%             5%            10%            90%            95%          100%")
fig.text(0.10 , 0.1, "0%         20%       40%       60%       80%      100%", rotation="vertical")

fig.savefig("/scratch/depfg/6196306/MSc_thesis/Final_Figures/review/composit_size_vel_"+str(season)+"_20x40_land_z500_v3.pdf")


#%% Compare ERA5 to ECE3Bis for temperature
tas_ERA = xr.open_dataset("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/ERA5/era5_t2m_zonalmean.nc", decode_times=False).t2m
time = pd.date_range("1950-01", periods=879, freq="M")
tas_ERA["time"]=time
tas_ECE = xr.open_mfdataset("/net/pc190625/nobackup_1/users/mourik/DATA/Temperature/ECEARTH3Bis/ensembles/tas/tas_hist*.nc", combine="nested", concat_dim="ensemble").tas

tas_ERA_lat = tas_ERA.sel(time=slice("1951", "2014")).mean(dim=("time", "lon")).sel(lat=slice(90, 0))
tas_ECE_lat = tas_ECE.sel(time=slice("1951", "2014")).mean(dim=("ensemble", "time", "lon"))
tas_ECE_lat_std = tas_ECE.sel(time=slice("1951", "2014")).mean(dim=("time", "lon")).std("ensemble").load()
tas_ECE_lat.load()

tas_ECE_year = tas_ECE.sel(time=slice("1950", "2014")).groupby("time.year").mean(dim=("ensemble", "lon", "lat", "time")).load()
tas_ECE_std = tas_ECE.sel(time=slice("1950", "2014")).mean(dim=("lon", "lat")).groupby("time.year").mean("time").mean("year").std("ensemble")


diff_lat = tas_ECE_lat - tas_ERA_lat
lat = diff_lat["lat"]

plt.figure()
plt.plot(diff_lat, lat, c="blue")
plt.fill_betweenx(lat, x1 = diff_lat-tas_ECE_lat_std, x2 = diff_lat+tas_ECE_lat_std, color="blue", alpha=0.3)
plt.xlabel(r"$\Delta T$ (K)")
plt.ylabel(r"Latitude ($^\circ$)")
plt.vlines(x=0, ymin=0, ymax=90, colors="black")
plt.xlim(-3, 3)
plt.savefig("/usr/people/mourik/Documents/Python/Figures/varlat/NH_lat_temp.pdf")
plt.show()

plt.figure()
plt.plot(np.arange(1950, 2015, 1), tas_ERA.sel(lat=slice(90,0)).mean(dim=("lat", "lon")).groupby("time.year").mean("time").sel(year=slice("1950", "2014")), label="ERA5", c="black")
plt.plot(np.arange(1950, 2015, 1), tas_ECE_year, label="ECE3p5", c="blue")
plt.fill_between(np.arange(1950, 2015, 1), y1=tas_ECE_year-tas_ECE_std, y2=tas_ECE_year+tas_ECE_std, color="blue", alpha=0.3)
plt.xlabel("Year")
plt.ylabel("T (K)")
plt.legend()
plt.savefig("/usr/people/mourik/Documents/Python/Figures/varlat/NH_annual_mean_temp.pdf")
plt.show()

coev_ECE = tas_ECE_year.sel(year=slice("1980", "2010")).polyfit(dim="year", deg=1).polyfit_coefficients
coev_ERA = tas_ERA.sel(lat=slice(90,0)).mean(dim=("lat", "lon")).groupby("time.year").mean("time").sel(year=slice("1980", "2010")).polyfit(dim="year", deg=1).polyfit_coefficients

#%% Compare ERA5 to ECE3p5 for Z500
z500_ECE = xr.open_mfdataset("/net/pc190625/nobackup/users/mourik/DATA/ECEARTH3Bis/remapbil2_5x2_5/northernhemisphere/historical/r*i1p5f1/ensemble/zg_day_EC-Earth3_historical_r*i1p5f1_seasmean.nc", combine="nested", concat_dim="ensemble").zg
z500_ERA = xr.open_dataset("/net/pc190625/nobackup/users/mourik/DATA/ERA5/era5_z500_seasmean.nc").z500

z500_ECE_year = z500_ECE.sel(time=slice("1950", "2014")).groupby("time.year").mean(dim=("time", "ensemble", "lat", "lon"))
z500_ERA_year = z500_ERA.sel(time=slice("1950", "2014")).groupby("time.year").mean(dim=("time", "lat", "lon"))
z500_ECE_std = z500_ECE.sel(time=slice("1950", "2014")).mean(dim=("lon", "lat")).groupby("time.year").mean("time").mean("year").std("ensemble")

plt.figure()
plt.plot(np.arange(1950, 2015, 1), z500_ERA_year, label="ERA5", c="black")
plt.plot(np.arange(1950, 2015, 1), z500_ECE_year, label="ECE3p5", c="blue")
plt.fill_between(np.arange(1950, 2015, 1), y1=z500_ECE_year-z500_ECE_std, y2=z500_ECE_year+z500_ECE_std, color="blue", alpha=0.3)
plt.xlabel("Year")
plt.ylabel("Z500 (m)")
plt.legend()
plt.savefig("/usr/people/mourik/Documents/Python/Figures/varlat/NH_annual_mean_z500.pdf")
plt.show()

coev_ECE = z500_ECE_year.sel(year=slice("1980", "2010")).polyfit(dim="year", deg=1).polyfit_coefficients
coev_ERA = z500_ERA.sel(lat=slice(90,0)).mean(dim=("lat", "lon")).groupby("time.year").mean("time").sel(year=slice("1980", "2010")).polyfit(dim="year", deg=1).polyfit_coefficients


#%% histogram with all velocities
hist_speedx_ms = hist_speedx*1000/(3600*24)
fig, ax = plt.subplots(1,2, sharey=True, sharex=True, figsize=(12,6)) 
hist_speedx_ms.sel(time=hist_speedx_ms["time.season"]=="DJF").plot.hist(density=True, bins=30, ax=ax[0])
ax[0].axvline(hist_speedx_ms.sel(time=hist_speedx["time.season"]=="DJF").quantile(0.9), color="black")
ax[0].axvline(hist_speedx_ms.sel(time=hist_speedx["time.season"]=="DJF").quantile(0.1), color="black")
ax[0].axvline(hist_speedx_ms.sel(time=hist_speedx["time.season"]=="DJF").quantile(0.95), color="grey")
ax[0].axvline(hist_speedx_ms.sel(time=hist_speedx["time.season"]=="DJF").quantile(0.05), color="grey")
ax[0].axvline(10, color="red")
ax[0].set_title("DJF")
ax[0].set_ylabel("Fraction")
ax[0].set_xlim(-35, 35)

hist_speedx_ms.sel(time=hist_speedx_ms["time.season"]=="JJA").plot.hist(density=True, bins=30, ax=ax[1])
ax[1].axvline(hist_speedx_ms.sel(time=hist_speedx["time.season"]=="JJA").quantile(0.9), color="black")
ax[1].axvline(hist_speedx_ms.sel(time=hist_speedx["time.season"]=="JJA").quantile(0.1), color="black")
ax[1].axvline(hist_speedx_ms.sel(time=hist_speedx["time.season"]=="JJA").quantile(0.95), color="grey")
ax[1].axvline(hist_speedx_ms.sel(time=hist_speedx["time.season"]=="JJA").quantile(0.05), color="grey")
ax[1].axvline(10, color="red")
ax[1].set_title("JJA")
fig.supxlabel("Velocity (m/s)", fontsize=14)