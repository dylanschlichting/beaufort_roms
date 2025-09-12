###################### Climatology Conditions for Passive Tracers #############
# The purpose of this script is to make climatology file for passive 
# tracers for the Beaufort shelf for 2020 so that we can use DVD to look at 
# numerical mixing of salinity (and eventually temperature, too).
#
# Notes:
# - dye_01 = salinity, dye_02 = salinity^2, dye_03 = numerical mixing of salinity 
# - This will first be done where dye_01 is set equal to salinity from the real
#   salt_clm file, dye_02 will be set equal to salinity^2 from the real saly_clm file, 
#   and dye_03 will be set to 0
#   - This instead pushes everything to zero
# - This has been updated to add in dyes for temperature numerical mixing, too
#################################################################################


# Load in packages
#%matplotlib widget # widget not currently working but instead prevents plots from showing
import matplotlib.pyplot as plt
#import ipywidgets as widgets
import numpy as np
import xarray as xr
import xesmf as xe
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
#salt_clm = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/salt_clm_2019_2024_20vert_001.nc')
salt_clm = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/salt_clm_2019_2024_30vert_001.nc')

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
N = len(salt_clm.s_rho)

# latitude
lat_len = len(grid.lat_rho)

# longitude
lon_len = len(grid.lon_rho)

# time
# OLD
# # Make a datetime array to use
# time_akdt = np.arange(datetime(2020,7,1,hour=1,minute=0,second=0), datetime(2020,11,2,hour=4,minute=0, second=0),timedelta(hours=3))
# time_akdt_dt = pd.to_datetime(time_akdt)
# time_tmp_len = len(time_akdt_dt)
# datetime1 = time_akdt_dt
# NEW
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
z_rho_pos = salt_clm.z_rho


# Make zeros that are the correct size to file the netcdf 
# (planning on using data first but this will be used for dye_03)
dye_zeros_clm = np.zeros((time_tmp_len, N, eta_rho_len, xi_rho_len))



# Set up the netcdf for the climatology 
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# If doing real values for the first two
#vert_dye_clm = '/scratch/alpine/brun1463/ROMS_scratch/Beaufort_ROMS_Perlmutter_Stuff_scratch/Final_bryclm_conds/Attempt001/passive_dye_clm_zeros_002.nc' 
# If doing all zeros 
#vert_dye_clm = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/passive_dye_clm_zeros_20vert_001.nc'   #UPDATE PATH
vert_dye_clm = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/passive_dye_clm_zeros_30vert_001.nc'   #UPDATE PATH

#create file to write to
nc1 = Dataset(vert_dye_clm, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006_thin_sponge.nc',
                      type = 'ROMS grid passive tracer climatolgoy',
                      history = 'Created by Brianna Undzis',
                      Conventions = 'CF',
                      Institution = 'University of Colorado Boulder',
                      date = str(datetime.today()))

#create dictionary for model
d = {}
d = global_defaults

for att, value in d.items():
    setattr(nc1, att, value)

# Create dimensions
nc1.createDimension('xi_rho',  Lp)   # RHO
nc1.createDimension('eta_rho', Mp)
nc1.createDimension('s_rho', N)
nc1.createDimension('s_w',   (N+1))
#nc1.createDimension('salt_time', None)
nc1.createDimension('ocean_time', None)
nc1.createDimension('dye_01_time', None)
nc1.createDimension('dye_02_time', None)
nc1.createDimension('dye_03_time', None)
nc1.createDimension('dye_04_time', None)
nc1.createDimension('dye_05_time', None)
nc1.createDimension('dye_06_time', None)
#nc1.createDimension('time_HYCOM', None)
nc1.createDimension('one',     1)

# Create variables
# --------------------
# Coordinate Variables
# --------------------
# xi rho
xi_rho = nc1.createVariable('xi_rho', 'd', ('xi_rho',), zlib=True)
xi_rho.long_name = 'xi coordinate of RHO-points'
xi_rho.standard_name = 'projection_xi_coordinate'
xi_rho.units = 'meter'
xi_rho_tmp = np.arange(0, Lp)
xi_rho[:] = xi_rho_tmp[:]

# eta rho
eta_rho = nc1.createVariable('eta_rho', 'd', ('eta_rho',), zlib=True)
eta_rho.long_name = 'eta coordinate of RHO-points'
eta_rho.standard_name = 'projection_eta_coordinate'
eta_rho.units = 'meter'
eta_rho_tmp = np.arange(0, Mp)
eta_rho[:] = eta_rho_tmp[:]

# s rho
s_rho = nc1.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
s_rho.long_name = 's coordinate of RHO-points'
s_rho.standard_name = 'projection_s_coordinate'
s_rho.units = 'meter'
s_rho_tmp = np.arange(0, N)
s_rho[:] = s_rho_tmp[:]

# # salt_time (in seconds)
# salt_time_g = nc2.createVariable('salt_time', None, ('salt_time'), zlib=True)
# salt_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
# salt_time_g.units = 'second'
# salt_time_g.field = 'time, scalar, series'
# salt_time_g[:] = time_tmp[:]

# ocean_time (in seconds)
ocean_time_g = nc1.createVariable('ocean_time', None, ('ocean_time'), zlib=True)
ocean_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
ocean_time_g.units = 'second'
ocean_time_g.field = 'time, scalar, series'
ocean_time_g[:] = time_tmp[:]

# dye_01_time (in seconds)
dye_01_time_g = nc1.createVariable('dye_01_time', None, ('dye_01_time'), zlib=True)
dye_01_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
dye_01_time_g.units = 'second'
dye_01_time_g.field = 'time, scalar, series'
dye_01_time_g[:] = time_tmp[:]

# dye_02_time (in seconds)
dye_02_time_g = nc1.createVariable('dye_02_time', None, ('dye_02_time'), zlib=True)
dye_02_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
dye_02_time_g.units = 'second'
dye_02_time_g.field = 'time, scalar, series'
dye_02_time_g[:] = time_tmp[:]

# dye_03_time (in seconds)
dye_03_time_g = nc1.createVariable('dye_03_time', None, ('dye_03_time'), zlib=True)
dye_03_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
dye_03_time_g.units = 'second'
dye_03_time_g.field = 'time, scalar, series'
dye_03_time_g[:] = time_tmp[:]

# dye_04_time (in seconds)
dye_04_time_g = nc1.createVariable('dye_04_time', None, ('dye_04_time'), zlib=True)
dye_04_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
dye_04_time_g.units = 'second'
dye_04_time_g.field = 'time, scalar, series'
dye_04_time_g[:] = time_tmp[:]

# dye_05_time (in seconds)
dye_05_time_g = nc1.createVariable('dye_05_time', None, ('dye_05_time'), zlib=True)
dye_05_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
dye_05_time_g.units = 'second'
dye_05_time_g.field = 'time, scalar, series'
dye_05_time_g[:] = time_tmp[:]

# dye_06_time (in seconds)
dye_06_time_g = nc1.createVariable('dye_06_time', None, ('dye_06_time'), zlib=True)
dye_06_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
dye_06_time_g.units = 'second'
dye_06_time_g.field = 'time, scalar, series'
dye_06_time_g[:] = time_tmp[:]
    
# # time (HYCOM version)
# time_g2 = nc2.createVariable('time_HYCOM', None, ('time_HYCOM'), zlib=True)
# time_g2.long_name = '3-hour time steps' 
# time_g2.units = 'datetime'
# time_g2.field = 'time_HYCOM, scalar, series'
# time_g2[:] = time_tmp2[:]


# --------------------
# Boundary Passive Tracers
# --------------------

# dye_01 (salinity)
dye_01_g = nc1.createVariable('dye_01', 'f8', ('dye_01_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True)
dye_01_g.long_name = 'dye_01 concentration'
dye_01_g.units = 'kilogram meter-3'
#dye_01_g[:,:,:] = salt_clm.salt[:,:,:,:] # real
dye_01_g[:,:,:] = dye_zeros_clm[:,:,:] # zeros

# dye_02 (salinity^2)
# dye_02
dye_02_g = nc1.createVariable('dye_02', 'f8', ('dye_02_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True)
dye_02_g.long_name = 'dye_02 concentration'
dye_02_g.units = 'kilogram meter-3'
#dye_02_g[:,:,:] = salt_clm.salt[:,:,:,:]**2 # real
dye_02_g[:,:,:] = dye_zeros_clm[:,:,:] # zeros

# dye_03 (numerical mixing of salinity)
dye_03_g = nc1.createVariable('dye_03', 'f8', ('dye_03_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True)
dye_03_g.long_name = 'dye_03 concentration'
dye_03_g.units = 'kilogram meter-3'
dye_03_g[:,:,:] = dye_zeros_clm[:,:,:] # zeros

# dye_04 (temperature)
dye_04_g = nc1.createVariable('dye_04', 'f8', ('dye_04_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True)
dye_04_g.long_name = 'dye_04 concentration'
dye_04_g.units = 'kilogram meter-3'
#dye_04_g[:,:,:] = salt_clm.salt[:,:,:,:] # real
dye_04_g[:,:,:] = dye_zeros_clm[:,:,:] # zeros

# dye_05 (temperature^2)
dye_05_g = nc1.createVariable('dye_05', 'f8', ('dye_05_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True)
dye_05_g.long_name = 'dye_05 concentration'
dye_05_g.units = 'kilogram meter-3'
#dye_05_g[:,:,:] = salt_clm.salt[:,:,:,:]**2 # real
dye_05_g[:,:,:] = dye_zeros_clm[:,:,:] # zeros

# dye_06 (numerical mixing of temperature)
dye_06_g = nc1.createVariable('dye_06', 'f8', ('dye_06_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True)
dye_06_g.long_name = 'dye_06 concentration'
dye_06_g.units = 'kilogram meter-3'
dye_06_g[:,:,:] = dye_zeros_clm[:,:,:] # zeros


# ---------------------------------------- End netCDF set up ----------------------------------------

# Close the netcdf
nc1.close()
