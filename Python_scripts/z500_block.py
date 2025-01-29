import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def speedx(clcmassX, clcmassY):
    R = 6371 #km
    dx = clcmassX.ffill(dim="day").isel(day=-1)-clcmassX.sel(day=0)
    dx = dx.where(dx.values<180, dx.values-360)
    dx = dx.where(dx.values>-180, dx.values+360)
    dxkm = dx*2*np.pi*R*np.cos(np.deg2rad(clcmassY.sel(day=0)))*2.5/360
    v = dxkm/clcmassX.count("day")
    return v

majlen = xr.open_dataarray("/Blockings/per_lon/Intensity/long/large/hist_BI_majlen_large.nc")#.sel(ensemble=15)
minlen = xr.open_dataarray("/Blockings/per_lon/Intensity/long/large/hist_BI_minlen_large.nc")#.sel(ensemble=15)
hist_clcmassX = xr.open_dataarray("/Blockings/per_lon/Intensity/long/large/hist_BI_clcmassX_large_weighted.nc")#.sel(ensemble=15)
hist_clcmassY = xr.open_dataarray("/Blockings/per_lon/Intensity/long/large/hist_BI_clcmassY_large_weighted.nc")#.sel(ensemble=15)
hist_BIav = xr.open_dataarray("/Blockings/per_lon/Intensity/long/large/hist_BI_BIav_large.nc")
hist_clarea = xr.open_dataarray("/Blockings/per_lon/Intensity/long/large/hist_BI_clarea_large_weighted.nc")

hist_speedx = speedx(hist_clcmassX, hist_clcmassY)

ens_str = ["02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
ens = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
for i in range(len(ens)):
    print(ens[i])
    #tas = xr.open_dataset("/Temperature/ECEARTH3Bis/ensembles/tas/anomaly/land/tas_hist"+str(ens_str[i])+"_anom_land.nc").tas
    z500 = xr.open_dataset("/ECEARTH3p5/historical/r"+str(ens[i]+1)+"i1p5f1/ensemble/zg_day_EC-Earth3_historical_r"+str(ens[i]+1)+"i1p5f1.nc").zg
    #tas["time"] = tas.indexes["time"].normalize()
    z500["time"] = z500.indexes["time"].normalize()
    T_block_shape = np.zeros([59653, 5, 17, 33])*np.nan
    for t in range(len(majlen["time"])):
        print(t)
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
    #T_block_xr.to_netcdf("/Temperature/Tblock/20x40/varlat/Tas_block_hist"+str(ens_str[i])+"_20x40_land.nc")
    z500_block_xr = xr.Dataset(data_vars=dict(z500=(["time", "block", "lat", "lon"], T_block_shape)), coords=dict(time=hist_clcmassX.sel(ensemble=0)["time"], lon=np.arange(-40, 42.5, 2.5), lat=np.arange(-20, 22.5, 2.5))).z500
    z500_block_xr.to_netcdf("/ECEARTH3p5/historical/z500block/z500_block_hist"+str(ens_str[i])+"_20x40.nc")

