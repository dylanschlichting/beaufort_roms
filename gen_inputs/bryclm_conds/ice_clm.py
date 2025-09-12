###################### Climatology Conditions for Ice ###############################
# The purpose of this script is to make climatology file for ice 
# for the Beaufort shelf for 2020 so that we can run with the sea ice model.
#
# Notes:
# - This has everything set to zero for now but will eventually be replaced with
#   real data 
#################################################################################


# Load in packages
#%matplotlib widget # widget not currently working but instead prevents plots from showing
import matplotlib.pyplot as plt
#import ipywidgets as widgets
import numpy as np
import xarray as xr
#import xesmf as xe
import pandas as pd
#import ESMF
import math
from netCDF4 import Dataset
from datetime import datetime, timedelta
from cftime import num2date, date2num


# Not sure how much of this we need....
# Load in the ROMS grid vertical coordinates
#grid_vertical = xr.open_dataset('/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Final_bryclm_conds/ROMS_grid_depth_hpluszeta_2020_003.nc', drop_variables='z_w')  #UPDATE PATH 
#salt_clm = xr.open_dataset('/pl/active/moriarty_lab/BriannaU/Paper1/2020_version/Inputs/Boundary_clm_conds/Attempt001/salt_clm_001.nc')
#salt_clm = xr.open_dataset('/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Final_bryclm_conds/Attempt001/salt_clm_001.nc')

# Load in the model grid
grid = xr.open_dataset('/global/homes/b/bundzis/Projects/Beaufort_ROMS_2020_test_nosed/Include/KakAKgrd_shelf_big010_smooth006_thin_sponge.nc')

# pull out the angle to rotate the currents to match the grid's u,v
#phi = grid_vertical.angle[0,0].values # radians 

# Read in the dimensions
#time_len = len(grid_vertical.time)
eta_rho_len = len(grid.eta_rho) # 206
xi_rho_len = len(grid.xi_rho) # 608
eta_u_len = len(grid.eta_u) # 206
xi_u_len = len(grid.xi_u) # 607
eta_v_len = len(grid.eta_v) # 205
xi_v_len = len(grid.xi_v) # 608

# Define other dimension lengths
# eta rho
Mp = len(grid.eta_rho)

# xi rho
Lp = len(grid.xi_rho)

# for u/v points 
Lm, Mm = (Lp-2), (Mp-2) #number/dimension of cells
L,  M  = Lm+1, Mm+1 #number/dimension of psi points

# srho
#N = len(salt_clm.s_rho)

# latitude
lat_len = len(grid.lat_rho)

# longitude
lon_len = len(grid.lon_rho)


# Read in the sea ice data
hycom_ice_data = xr.open_dataset('/pscratch/sd/b/bundzis/External_data/HYCOM_data/hycom_daily_ncss/hycom_ice_2019_2024.nc')
hycom_ice_data

# time
# Make a datetime array to use
# OLD - open water season 
# time_akdt = np.arange(datetime(2020,7,1,hour=1,minute=0,second=0), datetime(2020,11,2,hour=4,minute=0, second=0),timedelta(hours=3))
# time_akdt_dt = pd.to_datetime(time_akdt)
# time_tmp_len = len(time_akdt_dt)
# datetime1 = time_akdt_dt
# NEW 
# Full time - likely needs to be changed since ice data is probably not daily...but I guess for our purposes we 
# only want daily for nudging sooooo fill gaps and make daily, then use that here
time_akdt = np.arange(datetime(2019,1,1,hour=1,minute=0,second=0), datetime(2024,9,6,hour=1,minute=0, second=0),timedelta(days=1))
time_akdt_dt = pd.to_datetime(time_akdt)
time_tmp_len = len(time_akdt_dt)
datetime1 = time_akdt_dt

# Convert all time to seconds since 2000-01-01 (really 1999-12-31 so it starts at beginning
# of year=hour 0)
time_tmp = ((datetime1[:] - datetime(1999,12,31)).total_seconds() - 86400)

# copy the actual time values - seconds since 2000-01-01
#time_tmp = grid_vertical.time.values


# name the variables to fill the netcdf 
# latitude
lat_tmp = grid.lat_rho.values

# longitude
lon_tmp = grid.lon_rho.values

# Make a z_rho to be used in the netcdfs (0 is deepest, values are positive) # this took ~7 minutes
#z_rho_pos = salt_clm.z_rho


# Make zeros that are the correct size to file the netcdf 
# (planning on using data first but this will be used for dye_03)
ice_zeros_clm = np.zeros((time_tmp_len, eta_rho_len, xi_rho_len))
ice_zeros_clm_u = np.zeros((time_tmp_len, eta_u_len, xi_u_len))
ice_zeros_clm_v = np.zeros((time_tmp_len, eta_v_len, xi_v_len))


# Set up the netcdf for the climatology 
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# If doing real values for the first two
#vert_dye_clm = '/scratch/alpine/brun1463/ROMS_scratch/Beaufort_ROMS_Perlmutter_Stuff_scratch/Final_bryclm_conds/Attempt001/passive_dye_clm_zeros_002.nc' 
# If doing all zeros 
ice_clm = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_test_sed_scratch/Model_Inputs/Bry_Clm_Conds/Attempt001/ice_clm_zeros_001.nc'   #UPDATE PATH

#create file to write to
nc = Dataset(ice_clm, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006.nc',
                      type = 'ROMS ice climatolgoy zeros',
                      history = 'Created by Brianna Undzis',
                      Conventions = 'CF',
                      Institution = 'University of Colorado Boulder',
                      date = str(datetime.today()))

#create dictionary for model
d = {}
d = global_defaults

for att, value in d.items():
    setattr(nc, att, value)

# Create dimensions
nc.createDimension('xi_rho',  Lp)   # RHO
nc.createDimension('eta_rho', Mp)
nc.createDimension('xi_u', L) # U
nc.createDimension('xi_v', Lp)
nc.createDimension('eta_u', Mp) # V
nc.createDimension('eta_v', M)
#nc.createDimension('s_rho', N)
#nc.createDimension('s_w',   (N+1))
#nc1.createDimension('salt_time', None)
nc.createDimension('ocean_time', None)
#nc1.createDimension('time_HYCOM', None)
nc.createDimension('one',     1)

# Create variables
# --------------------
# Coordinate Variables
# --------------------
# xi rho
xi_rho = nc.createVariable('xi_rho', 'd', ('xi_rho',), zlib=True)
xi_rho.long_name = 'xi coordinate of RHO-points'
xi_rho.standard_name = 'projection_xi_coordinate'
xi_rho.units = 'meter'
xi_rho_tmp = np.arange(0, Lp)
xi_rho[:] = xi_rho_tmp[:]

# eta rho
eta_rho = nc.createVariable('eta_rho', 'd', ('eta_rho',), zlib=True)
eta_rho.long_name = 'eta coordinate of RHO-points'
eta_rho.standard_name = 'projection_eta_coordinate'
eta_rho.units = 'meter'
eta_rho_tmp = np.arange(0, Mp)
eta_rho[:] = eta_rho_tmp[:]

# # salt_time (in seconds)
# salt_time_g = nc2.createVariable('salt_time', None, ('salt_time'), zlib=True)
# salt_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
# salt_time_g.units = 'second'
# salt_time_g.field = 'time, scalar, series'
# salt_time_g[:] = time_tmp[:]

# ocean_time (in seconds)
ocean_time_g = nc.createVariable('ocean_time', None, ('ocean_time'), zlib=True)
ocean_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
ocean_time_g.units = 'second'
ocean_time_g.field = 'time, scalar, series'
ocean_time_g[:] = time_tmp[:]
    
# # time (HYCOM version)
# time_g2 = nc2.createVariable('time_HYCOM', None, ('time_HYCOM'), zlib=True)
# time_g2.long_name = '3-hour time steps' 
# time_g2.units = 'datetime'
# time_g2.field = 'time_HYCOM, scalar, series'
# time_g2[:] = time_tmp2[:]


# --------------------
# Climatology Ice
# --------------------

# Sea Ice
# Ice area fraction (Aice) 
Aice_clm_g = nc.createVariable('Aice', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
Aice_clm_g.standard_name = 'sea_ice_area_fraction'
Aice_clm_g.long_name = 'fraction of cell covered by ice'
Aice_clm_g.units = 'nondimensional'
#Aice_clm_g[:,:,:] = ice_zeros_clm

# Ice thickness (ice_thickness) 
ice_thickness_clm_g = nc.createVariable('ice_thickness', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
ice_thickness_clm_g.standard_name = 'sea_ice_thickness'
ice_thickness_clm_g.long_name = 'average ice thickness in cell'
ice_thickness_clm_g.units = 'meter'
#ice_thickness_clm_g[:,:,:] = ice_zeros_clm

# Meltpond thickness (meltpond_thickness)
# meltpond_thickness_g = nc.createVariable('meltpond_thickness', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
# meltpond_thickness_g.standard_name = 'melt_pond_water_thickness_on_sea_ice'
# meltpond_thickness_g.long_name = 'surface melt water thickness on ice'
# meltpond_thickness_g.units = 'meter'
# meltpond_thickness_g[:,:,:] = meltpond_thickness_tmp 

# Ice age (ice_age) 
ice_age_clm_g = nc.createVariable('ice_age', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
ice_age_clm_g.standard_name = 'age_of_sea_ice'
ice_age_clm_g.long_name = 'age of sea ice'
ice_age_clm_g.units = 'second'
ice_age_clm_g[:,:,:] = ice_zeros_clm

# Snow thickness (snow_thickness) 
snow_thickness_clm_g = nc.createVariable('snow_thickness', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
snow_thickness_clm_g.standard_name = 'snowfall_thickness_above_sea_ice'
snow_thickness_clm_g.long_name = 'thickness of snow cover'
snow_thickness_clm_g.units = 'meter'
snow_thickness_clm_g[:,:,:] = ice_zeros_clm  

# # Sea ice temperature (sea_ice_temperature)
# tice_g = nc.createVariable('Tice', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
# tice_g.standard_name = 'sea_ice_temperature'
# tice_g.long_name = 'interior ice temperature'
# tice_g.units = 'Celcius'
# tice_g[:,:] = ice_temp_tmp

# # Under ice temperature (under_ice_temp)
# under_ice_temp_g = nc.createVariable('under_ice_temp', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
# under_ice_temp_g.standard_name = 'temperature_of_molecular_sub_layer_under_sea_ice'
# under_ice_temp_g.long_name = 'temperature of molecular sub-layer under ice'
# under_ice_temp_g.units = 'Celcius'
# under_ice_temp_g[:,:] = under_ice_temp_tmp

# # Under ice salinity (under_ice_salt)
# under_ice_salt_g = nc.createVariable('under_ice_salt', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
# under_ice_salt_g.standard_name = 'salinity_of_molecular_sub_layer_under_sea_ice'
# under_ice_salt_g.long_name = 'salinity of molecular sub-layer under ice'
# under_ice_salt_g.units = 'Celcius'
# under_ice_salt_g[:,:] = under_ice_salt_tmp

# Sea ice surface temperature  (ice_sst)
sea_ice_surface_temperature_clm_g = nc.createVariable('ice_sst', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
sea_ice_surface_temperature_clm_g.standard_name = 'sea_ice_surface_temperature'
sea_ice_surface_temperature_clm_g.long_name = 'temperature of ice/snow surface'
sea_ice_surface_temperature_clm_g.units = 'Celcius'
sea_ice_surface_temperature_clm_g[:,:,:] = ice_zeros_clm

# Ice u-velocity (Uice) 
Uice_clm_g = nc.createVariable('Uice', 'f8', ('ocean_time', 'eta_u', 'xi_u'), zlib=True)
Uice_clm_g.standard_name = 'sea_ice_x_velocity'
Uice_clm_g.long_name = 'u-component of ice velocity'
Uice_clm_g.units = 'meter second-1'
Uice_clm_g[:,:,:] = ice_zeros_clm_u

# Ice v-velocity (Vice)
Vice_clm_g = nc.createVariable('Vice', 'f8', ('ocean_time', 'eta_v', 'xi_v'), zlib=True)
Vice_clm_g.standard_name = 'sea_ice_y_velocity'
Vice_clm_g.long_name = 'v-component of ice velocity'
Vice_clm_g.units = 'meter second-1'
Vice_clm_g[:,:,:] = ice_zeros_clm_v


# ---------------------------------------- End netCDF set up ----------------------------------------

# Close the netcdf
nc.close()
