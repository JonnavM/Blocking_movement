#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 06:35:38 2022

@author: mourik
This is the basecode for finding atmospheric blockings in datasets. 
It results in a file containing LAT_min, a file containing the blockings and a file 
containing GHGS and GHGN.
"""
#%% Import needed modules
import numpy as np
import xarray as xr

#%% Import data and calculate [Z500]
" The first part of the code is used to import the data and get an hemisphere and time mean geopotential height, [Z500]."
" This [Z500] will later be used as a minimum latitude from which to calculate the blockings."

#Choose which dataset to use: "ERA5" or "ECE3p5" historical
var = "ERA5"

if var == "ECE3p5":
    members = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
elif var == "ERA5":
    members = [1]

data = []
z500 = []
avg_z500 = []

for i in members:
    if var == "ECE3p5":
        data.append(xr.open_mfdataset("/net/pc190625/nobackup/users/mourik/DATA/ECEARTH3Bis/remapbil2_5x2_5/northernhemisphere/historical/r"+str(i)+"i1p5f1/day/zg/*.nc"))
        z500.append(data[i-members[0]].zg.mean("plev")) 
    elif var == "ERA5":
        data.append(xr.open_dataset("/net/pc190625/nobackup/users/mourik/DATA/ERA5/regrid/era5_z500_daily_-180-180E_0-90N.nc"))
        z500.append(data[i-members[0]].z500/9.81)
        
time_var = z500[0]['time'].data
lat = z500[0]['lat'].data #Is the same for all of them
lon = z500[0]['lon'].data

# Calculate time mean over previous 15 days, excluding the day itself, and area mean over NH, to get the hemisphere-wide mean
# Average over the whole northern hemisphere. --> Weighted average
gw = np.cos(np.deg2rad(z500[0].lat))
rol = 15 #number of days over which to average

for i in range(len(members)):
    avg_z500.append(z500[i].weighted(gw).mean(dim=("lat", "lon")).rolling(time=rol).mean())

time_var = time_var[rol:]

#%% Calculate lat_min for determining subtropical belt
" In this part of the code, the minimum latitude is extracted by comparing [Z500] to Z500."
res = lon[2]-lon[1] #Resolution of the dataset. Assumming 

def lat_min(z500_var, avg_z500_var, time_var):
    "Function which calculates the minimum latitude as one value per day. z500_var = z500[i], avg_z500_var=avg_z500[i], time_var=time_var"
    avg_z500_var = avg_z500_var.sel(time=time_var)
    z500_var = z500_var.sel(time=time_var)
    condition = (z500_var<=avg_z500_var)
    box = condition.cumsum(dim="lat")
    LAT_min = box.where(box!=0, np.nan).idxmin(dim="lat").mean("lon")-res
    return LAT_min

#Make new definition for saving a different LAT_min per longitude
def lat_min_lon(z500_var, avg_z500_var, time_var):
    "Function which calculates the minimum latitude per longitude. z500_var = z500[i], avg_z500_var=avg_z500[i], time_var=time_var"
    avg_z500_var = avg_z500_var.sel(time=time_var)
    z500_var = z500_var.sel(time=time_var)
    condition = (z500_var<=avg_z500_var)
    box = condition.cumsum(dim="lat")
    LAT_min = box.where(box!=0, np.nan).idxmin(dim="lat")-res
    return LAT_min

#The minimum latitudes are saved in these lists:
m_LAT_min = []
m_LAT_min_lon = []

#Loop through the ensemble members to fill the lists
for i in range(len(members)):
    print("member", members[i])
    #m_LAT_min.append(lat_min(z500[i], avg_z500[i], time_var)) #Use this one if you want one value of LATmin per day
    m_LAT_min_lon.append(lat_min_lon(z500[i], avg_z500[i], time_var)) #Use this one if you want one value per longitude per day
    print("Saving")
    if var == "ECE3p5":
        m_LAT_min_lon[i].to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/LAT_min/LAT_min_lon/lat_min_lon_hist"+str(members[i])+"_1850_2014.nc")
    elif var == "ERA5":
        #m_LAT_min_lon[i].to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/LAT_min/LAT_min_lon/lat_min_lon_ERA_1950_2022.nc")
        m_LAT_min_lon[i].to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/LAT_min/LAT_min_lon/lat_min_lon_ERA_2023_mei.nc")
#%% Make a matrix filled with the values of the minimum latitude:
LAT_min = []
for i in members:
    if var =="ECE3p5":
        LAT_min.append(xr.open_dataarray("/net/pc190625/nobackup_1/users/mourik/DATA/LAT_min/LAT_min_lon/lat_min_lon_hist"+str(i)+"_1850_2014.nc"))
    elif var == "ERA5":
        LAT_min.append(xr.open_dataarray("/net/pc190625/nobackup_1/users/mourik/DATA/LAT_min/LAT_min_lon/lat_min_lon_ERA_1950_2022.nc"))

def matrix_latmin(LAT_min, time):
    matrix_2d = np.zeros([len(time), len(lon)])
    matrix_2d[:,:]=LAT_min
    matrix_3d = np.repeat(matrix_2d[:, np.newaxis, :], repeats=len(lat), axis=1)
    matrix_latmin = xr.DataArray(matrix_3d, coords=[("time", time), ("lat", lat), ("lon", lon)], dims=["time", "lat", "lon"])
    return matrix_latmin

#%% Finds blockings
" In this part of the code, the blockings are filtered out of the data. This is done by defining the geopotential height gradients"
" and requiring GHGS>0 and GHGN<0. GHGS2 is also defined here, but will not be used unless we want to distinguish between "
" different kinds of blockings."

#Make an array with shifted latitude, both for plus and for minus dphi degrees
dphi = 15 #Delta phi, typical latitudinal extension of blockings according to Sousa et al., 2021.

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

phi_plus = np.concatenate([lat[find_nearest(lat, dphi):], lat[:find_nearest(lat, dphi)]])
phi_min = np.concatenate([lat[find_nearest(lat, 90-dphi):], lat[:find_nearest(lat, 90-dphi)]])
phi_2min = np.concatenate([lat[find_nearest(lat, 90-2*dphi):], lat[:find_nearest(lat, 90-2*dphi)]])

# Define geopotential height gradients 
def Z500(var, labda, phi, time):
    "Puts Z500 in the same form as the geopotetial gradients. Makes it easier to work with. var=z500[i]"
    return var.sel(time=time, lat=phi, lon=labda)
 
def GHGS(var, lat_var, labda, phi, time):
    "Occurence of equatorward easterlies"
    a = (Z500(var, labda, phi, time).values - Z500(var, labda, phi_min, time).values)/dphi
    b = xr.DataArray(a, coords=[("time", time), ("lat", phi), ("lon", labda)], dims=["time", "lat", "lon"])
    result = b.where(b.lat>=matrix_latmin(lat_var, time))
    return result
        
def GHGN(var, lat_var, labda, phi, time):
    "Occurence of poleward westerlies"
    a = (Z500(var, labda, phi_plus, time).values - Z500(var, labda, phi, time).values)/dphi
    b = xr.DataArray(a, coords=[("time", time), ("lat", phi), ("lon", labda)], dims=["time", "lat", "lon"])
    result = b.where(np.logical_and(b.lat>=matrix_latmin(lat_var, time), b.lat<=75))
    return result

def GHGS2(var, lat_var, labda, phi, time):
    "For distinction between Rex and Omega"
    a = (Z500(var, labda, phi_min, time).values - Z500(var, labda, phi_2min, time).values)/dphi
    b = xr.DataArray(a, coords=[("time", time), ("lat", phi), ("lon", labda)], dims=["time", "lat", "lon"])
    result = b.where(b.lat>=matrix_latmin(lat_var, time))
    return result


#%% Save GHGN and GHGS to files. This step is only necessary if you want to study this data separately. 
GHGN_hist = []
GHGS_hist = []
for i in range(len(members)):
    GHGN_hist.append(GHGN(z500[i], LAT_min[i], lon, lat, time_var))
    GHGS_hist.append(GHGS(z500[i], LAT_min[i], lon, lat, time_var))
    GHGN_hist[i].to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/GHG/GHGN_585_"+str(members[i])+"_2015_2166.nc")
    GHGS_hist[i].to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/GHG/GHGS_585_"+str(members[i])+"_2015_2166.nc")

#%% Local spatial filtering
"Keep only gridpoints embedded in the high pressure systems. GHG>20m per degree is defined as strong wind flow and is therefore neglected in the analysis."
dphi_z = lon[2]-lon[1]
dphi_m = lat[2]-lat[1]
phi_z_plus = np.concatenate([lat[1:], lat[:1]])
phi_z_min = np.concatenate([lat[-1:], lat[:-1]])
labda_plus = np.concatenate([lon[1:], lon[:1]])
labda_min = np.concatenate([lon[-1:], lon[:-1]])
                            
def GHG_zonal(var, labda, phi, time):
    "Calculate zonal gradient in geopotential height"
    a = (Z500(var, labda, phi_z_plus, time).values - Z500(var, labda, phi_z_min, time).values)/(2*dphi_z)
    b = xr.DataArray(a, coords=[("time", time), ("lat", phi), ("lon", labda)], dims=["time", "lat", "lon"])
    return b

def GHG_meridional(var, labda, phi, time):
    "Calculate meridional gradient in geopotential height"
    a = (Z500(var, labda_plus, phi, time).values - Z500(var, labda_min, phi, time).values)/(2*dphi_m)
    b = xr.DataArray(a, coords=[("time", time), ("lat", phi), ("lon", labda)], dims=["time", "lat", "lon"])
    return b

def GHGu20(var, lat_var, labda, phi, time):
    "Calculate geopotential height gradient and check where the geostrophic wind magnitude is smaller than 20 m per degree"
    a = np.sqrt(GHG_zonal(var, labda, phi, time).values**2 + GHG_meridional(var, labda, phi, time).values**2)
    b = xr.DataArray(a, coords=[("time", time), ("lat", phi), ("lon", labda)], dims=["time", "lat", "lon"])
    result = b.where(np.logical_and(b.lat>=matrix_latmin(lat_var, time), b.lat<=75))
    return result

#Different ways of saving blockings:
    #1) Only with 0's and 1's
    #2) With sqrt(GHGS**2+GHGN**2)
    #3) With Wiedenmann's method BI=100[(MZ/RC)-1]. This is the one that we will use in this study.

#1)
#def block(var, lat_var, labda, phi, time):
#    return np.logical_and(np.logical_and(GHGS(var, lat_var, labda, phi, time)>0, GHGN(var, lat_var, labda, phi, time)<0), GHGu20(var, lat_var, labda, phi, time)<20)*1

#2)
#def block(var, lat_var, labda, phi, time):
#    GHG_sqrt = np.sqrt(GHGS(var, lat_var, labda, phi, time)**2+GHGN(var, lat_var, labda, phi, time)**2)
#    condition = np.logical_and(np.logical_and(GHGS(var, lat_var, labda, phi, time)>0, GHGN(var, lat_var, labda, phi, time)<0), GHGu20(var, lat_var, labda, phi, time)<20)
#    return GHG_sqrt.where(condition==True, 0)

#3) 
lon_step = int(60/res) #Find datapoint 60 degrees upstrean
end_lon = int(len(lon)+lon_step)
def block(var, lat_var, labda, phi, time):
    "This function returns the locations and intensities of blocks on a labdaxphi grid for the given time"
    MZ = Z500(var, labda, phi, time) #Initial z500 datapoint that will be compared to its surroundings
    zu = MZ.pad(lon=[lon_step,lon_step], mode="wrap").rolling(lon=lon_step).min() #Upstream z500 value
    zd = zu.shift(lon=-lon_step) #Downstream z500 value
    zu = zu.isel(lon=slice(lon_step, end_lon))
    zd = zd.isel(lon=slice(lon_step, end_lon))
    RC = 0.5*((zu+MZ)/2+(zd+MZ)/2)
    BI = 100*((MZ/RC)-1)
    condition = np.logical_and(np.logical_and(GHGS(var, lat_var, labda, phi, time)>0, GHGN(var, lat_var, labda, phi, time)<0), GHGu20(var, lat_var, labda, phi, time)<20)
    return BI.where(condition==True, 0)

#Save blocking index
block_hist = []
for i in range(len(members)):
    print(i)
    block_hist.append(block(z500[i], LAT_min[i], lon, lat, time_var))   #Choose which type of block to use
    if var == "ECE3p5":
        block_hist[i].to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/Blockings/per_lon/Intensity/block_var_hist"+str(members[i])+"_BI_1850_2014.nc")
    elif var == "ERA5":
        #block_hist[i].to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/Blockings/per_lon/Intensity/block_var_ERA_BI_1950_2022.nc")
        block_hist[i].to_netcdf("/net/pc190625/nobackup_1/users/mourik/DATA/Blockings/per_lon/Intensity/block_var_ERA_BI_2023_mei.nc")
