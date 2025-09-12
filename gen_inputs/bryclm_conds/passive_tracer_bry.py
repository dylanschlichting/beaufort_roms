###################### Boundary Conditions for Passive Tracers #############
# The purpose of this script is to make boundary conditions for passive 
# tracers for the Beaufort shelf for 2020 so that we can use DVD to look at 
# numerical mixing of salinity (and eventually temperature, too).
#
# Notes:
# - dye_01 = salinity, dye_02 = salinity^2, dye_03 = numerical mixing of salinity 
# - This will first be done where dye_01 is set equal to salinity from the real
#   salt_bry file, dye_02 will be set equal to salinity^2 from the real saly_bry file, 
#   and dye_03 will be set to 0
# - This will also be set up so that they are all set to 0 and make another 
#   forcing file using that
#   - We are using the version that sets everything to zero
# - This has been updated to add values for temperature numerical mixing 
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
#salt_bry = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/salt_bry_2019_2024_20vert_001.nc')
salt_bry = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/salt_bry_2019_2024_30vert_001.nc')

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
N = len(salt_bry.s_rho)

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
z_rho_pos = salt_bry.z_rho


# Make zeros that are the correct size to fille the netcdf 
dye_zeros_west = np.zeros((time_tmp_len, N, eta_rho_len))

dye_zeros_north = np.zeros((time_tmp_len, N, xi_rho_len)) 

dye_zeros_east = np.zeros((time_tmp_len, N, eta_rho_len)) 


# Set up the netcdf for the climatology 
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# If doing real values for the first two
#vert_dye_bry = '/scratch/alpine/brun1463/ROMS_scratch/Beaufort_ROMS_Perlmutter_Stuff_scratch/Final_bryclm_conds/Attempt001/passive_dye_bry_salt_salt2_zerofor03_001.nc' 
# If doing all zeros 
#vert_dye_bry = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/passive_dye_bry_zeros_20vert_001.nc'   #UPDATE PATH
vert_dye_bry = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/passive_dye_bry_zeros_30vert_001.nc'   #UPDATE PATH

#create file to write to
nc1 = Dataset(vert_dye_bry, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006_thin_sponge.nc',
                      type = 'ROMS grid passive tracer boundary',
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
nc1.createDimension('dye_time', None)
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

# dye_time (in seconds)
dye_time_g = nc1.createVariable('dye_time', None, ('dye_time'), zlib=True)
dye_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
dye_time_g.units = 'second'
dye_time_g.field = 'time, scalar, series'
dye_time_g[:] = time_tmp[:]
    
# # time (HYCOM version)
# time_g2 = nc2.createVariable('time_HYCOM', None, ('time_HYCOM'), zlib=True)
# time_g2.long_name = '3-hour time steps' 
# time_g2.units = 'datetime'
# time_g2.field = 'time_HYCOM, scalar, series'
# time_g2[:] = time_tmp2[:]


# --------------------
# Vertical variables
# --------------------

# z_r (depths of rho points)
z_r_g = nc1.createVariable('z_rho', 'd', ('dye_time', 's_rho', 'eta_rho', 'xi_rho'), zlib = True)
z_r_g.long_name = 'Z coordinate of rho-points'
z_r_g.units = 'meter'
z_r_g[:,:,:,:] = z_rho_pos[:,:,:,:]


# --------------------
# Boundary Passive Tracers
# --------------------

# dye_01
# dye_west_01
dye_west_01_g = nc1.createVariable('dye_west_01', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_west_01_g.long_name = 'dye_01 at western boundary'
dye_west_01_g.units = 'kilogram meter-3'
#dye_west_01_g[:,:,:] = salt_bry.salt_west[:,:,:] # real
dye_west_01_g[:,:,:] = dye_zeros_west[:,:,:] # zeros


# dye_north_01
dye_north_01_g = nc1.createVariable('dye_north_01', 'f8', ('dye_time', 's_rho', 'xi_rho'), zlib=True)
dye_north_01_g.long_name = 'dye_01 at northern boundary'
dye_north_01_g.units = 'kilogram meter-3'
#dye_north_01_g[:,:,:] = salt_bry.salt_north[:,:,:] # real
dye_north_01_g[:,:,:] = dye_zeros_north[:,:,:] # zeros

# dye_east_01
dye_east_01_g = nc1.createVariable('dye_east_01', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_east_01_g.long_name = 'dye_01 at eastern boundary'
dye_east_01_g.units = 'kilogram meter-3'
#dye_east_01_g[:,:,:] = salt_bry.salt_east[:,:,:] # real
dye_east_01_g[:,:,:] = dye_zeros_east[:,:,:] # zeros


# dye_02
# dye_west_02
dye_west_02_g = nc1.createVariable('dye_west_02', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_west_02_g.long_name = 'dye_02 at western boundary'
dye_west_02_g.units = 'kilogram meter-3'
#dye_west_02_g[:,:,:] = salt_bry.salt_west[:,:,:]**2 # real
dye_west_02_g[:,:,:] = dye_zeros_west[:,:,:] # zeros

# dye_north_02
dye_north_02_g = nc1.createVariable('dye_north_02', 'f8', ('dye_time', 's_rho', 'xi_rho'), zlib=True)
dye_north_02_g.long_name = 'dye_02 at northern boundary'
dye_north_02_g.units = 'kilogram meter-3'
#dye_north_02_g[:,:,:] = salt_bry.salt_north[:,:,:]**2 # real
dye_north_02_g[:,:,:] = dye_zeros_north[:,:,:] # zeros

# dye_east_02
dye_east_02_g = nc1.createVariable('dye_east_02', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_east_02_g.long_name = 'dye_02 at eastern boundary'
dye_east_02_g.units = 'kilogram meter-3'
#dye_east_02_g[:,:,:] = salt_bry.salt_east[:,:,:]**2 # real
dye_east_02_g[:,:,:] = dye_zeros_east[:,:,:] # zeros


# dye_03
# dye_west_03
dye_west_03_g = nc1.createVariable('dye_west_03', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_west_03_g.long_name = 'dye_03 at western boundary'
dye_west_03_g.units = 'kilogram meter-3'
dye_west_03_g[:,:,:] = dye_zeros_west[:,:,:] # zeros

# dye_north_03
dye_north_03_g = nc1.createVariable('dye_north_03', 'f8', ('dye_time', 's_rho', 'xi_rho'), zlib=True)
dye_north_03_g.long_name = 'dye_03 at northern boundary'
dye_north_03_g.units = 'kilogram meter-3'
dye_north_03_g[:,:,:] = dye_zeros_north[:,:,:] # zeros

# dye_east_03
dye_east_03_g = nc1.createVariable('dye_east_03', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_east_03_g.long_name = 'dye_03 at eastern boundary'
dye_east_03_g.units = 'kilogram meter-3'
dye_east_03_g[:,:,:] = dye_zeros_east[:,:,:] # zeros

# dye_04
# dye_west_04
dye_west_04_g = nc1.createVariable('dye_west_04', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_west_04_g.long_name = 'dye_04 at western boundary'
dye_west_04_g.units = 'kilogram meter-3'
dye_west_04_g[:,:,:] = dye_zeros_west[:,:,:] # zeros

# dye_north_04
dye_north_04_g = nc1.createVariable('dye_north_04', 'f8', ('dye_time', 's_rho', 'xi_rho'), zlib=True)
dye_north_04_g.long_name = 'dye_04 at northern boundary'
dye_north_04_g.units = 'kilogram meter-3'
dye_north_04_g[:,:,:] = dye_zeros_north[:,:,:] # zeros

# dye_east_04
dye_east_04_g = nc1.createVariable('dye_east_04', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_east_04_g.long_name = 'dye_04 at eastern boundary'
dye_east_04_g.units = 'kilogram meter-3'
dye_east_04_g[:,:,:] = dye_zeros_east[:,:,:] # zeros

# dye_05
# dye_west_05
dye_west_05_g = nc1.createVariable('dye_west_05', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_west_05_g.long_name = 'dye_05 at western boundary'
dye_west_05_g.units = 'kilogram meter-3'
dye_west_05_g[:,:,:] = dye_zeros_west[:,:,:] # zeros

# dye_north_05
dye_north_05_g = nc1.createVariable('dye_north_05', 'f8', ('dye_time', 's_rho', 'xi_rho'), zlib=True)
dye_north_05_g.long_name = 'dye_05 at northern boundary'
dye_north_05_g.units = 'kilogram meter-3'
dye_north_05_g[:,:,:] = dye_zeros_north[:,:,:] # zeros

# dye_east_05
dye_east_05_g = nc1.createVariable('dye_east_05', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_east_05_g.long_name = 'dye_05 at eastern boundary'
dye_east_05_g.units = 'kilogram meter-3'
dye_east_05_g[:,:,:] = dye_zeros_east[:,:,:] # zeros

# dye_06
# dye_west_06
dye_west_06_g = nc1.createVariable('dye_west_06', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_west_06_g.long_name = 'dye_06 at western boundary'
dye_west_06_g.units = 'kilogram meter-3'
dye_west_06_g[:,:,:] = dye_zeros_west[:,:,:] # zeros

# dye_north_06
dye_north_06_g = nc1.createVariable('dye_north_06', 'f8', ('dye_time', 's_rho', 'xi_rho'), zlib=True)
dye_north_06_g.long_name = 'dye_06 at northern boundary'
dye_north_06_g.units = 'kilogram meter-3'
dye_north_06_g[:,:,:] = dye_zeros_north[:,:,:] # zeros

# dye_east_06
dye_east_06_g = nc1.createVariable('dye_east_06', 'f8', ('dye_time', 's_rho', 'eta_rho'), zlib=True)
dye_east_06_g.long_name = 'dye_06 at eastern boundary'
dye_east_06_g.units = 'kilogram meter-3'
dye_east_06_g[:,:,:] = dye_zeros_east[:,:,:] # zeros

# ---------------------------------------- End netCDF set up ----------------------------------------

# Close the netcdf
nc1.close()


