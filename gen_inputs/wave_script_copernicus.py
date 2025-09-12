################ Wave Forcing - Copernicus ############################
# The purpose of this script is to look 
# at the Copernicus wave data for 2019 - 2024
# to see if it would be good to use for the wave 
# forcing. If it is, this script will then make 
# the wave forcing file for ROMS.
#
# Notes:
# - This has been adapted from wave_script_copernicus.py on 
#   Alpine in Beaufort_Shelf_Rivers_Alpine_002/Scripts/Waves
#   but ignores the bottom wave parameters since we are not doing
#   sediment here
#######################################################################


# Load in the packages 
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd
#import ESMF
import math
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from cftime import num2date, date2num


# Load in the ROMS grid 
# List variables to drop 
#drop_vars = ['z_w', 'z_rho']
grid = xr.open_dataset('/global/homes/b/bundzis/Projects/Beaufort_ROMS_2020_test_nosed/Include/KakAKgrd_shelf_big010_smooth006_thin_sponge.nc')  # UPDATE PATH
#'/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/ROMS_grid_depth_hpluszeta_011.nc', drop_variables=drop_vars)


# Set a bunch of dimensions 

# Pull out the angle to rotate the currents to match the grid's u,v
phi = grid.angle[0,0].values # radians 

# Read in the dimensions
# rho
eta_rho_len = len(grid.eta_rho) # 206
xi_rho_len = len(grid.xi_rho) # 608
# u
eta_u_len = len(grid.eta_u) # 206
xi_u_len = len(grid.xi_u) # 607
# v
eta_v_len = len(grid.eta_v) #
xi_v_len = len(grid.xi_v) # 

# Define other dimension lengths
# eta rho
Mp = len(grid.eta_rho)
# xi rho
Lp = len(grid.xi_rho)

# latitude
lat_u_len = len(grid.lat_u)
lat_v_len = len(grid.lat_v)

# longitude
lon_u_len = len(grid.lon_u)
lon_v_len = len(grid.lon_v)


# Load in the Copernicus wave data
#wave_data = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Waves/ww3_global_WaveWatch_III_Global_Wave_Model_best.nc')  # UPDATE PATH
wave_data = xr.open_dataset('/pscratch/sd/b/bundzis/External_data/Copernicus/Waves/ARCTIC_MULTIYEAR_WAV_002_013/cmems_mod_arc_wav_my_3km_PT1H-i_202012/copernicus_waves_2019_2024_beaufort.nc')  # UPDATE PATH

print(wave_data.time[0:5].values, flush=True)
print(wave_data.time[-5:-1].values, flush=True)


# Make a time array of the times we want
# ERA5 data spans 2017-01-01 hour 1 to 2024-12-31 hour 23 (hourly)
# Make a datetime array of the corresponding datetime values in AKDT to use for the netcdf
#time_utc = np.arange(datetime(2017,7,1,hour=1,minute=0,second=0), datetime(2024,12,31,hour=23,minute=0, second=0),timedelta(hours=1)) # in UTC, ERA5 time
time_akdt = np.arange(datetime(2019,1,1,hour=0,minute=0,second=0), datetime(2024,6,30,hour=23,minute=0, second=0),timedelta(hours=1)) # in AKT (daylight savings...March - November...)

time_akdt_dt = pd.to_datetime(time_akdt)
print(time_akdt_dt[0:7], flush=True)
print(time_akdt_dt[-9:-1], flush=True)
print(time_akdt_dt[-1], flush=True)
print(len(time_akdt_dt), flush=True)

time_len = len(time_akdt_dt)
datetime1 = time_akdt_dt


# Convert all time to seconds since 2000-01-01 (really 1999-12-31 so it starts at beginning
# of year=hour 0)
time_tmp = ((datetime1[:] - datetime(1999,12,31)).total_seconds() - 86400)

# name the variables to fill the netcdf 
# latitude
lat_tmp = grid.lat_rho.values

# longitude
lon_tmp = grid.lon_rho.values

# ww3 number of lats
cop_lat_len = len(wave_data.lat.values)
#print('era5_lat_len', era5_lat_len)

# era5 number of lons
cop_lon_len = len(wave_data.lon.values)
#print('era5_lon_len', era5_lon_len)


# Set up the netcdf for the forcing
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# OG
wave_frc = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/wave_forcing_file_kaktovik_shelf_cop_2019_2024_001.nc'  # UPDATE PATH
# TEMP
#wave_frc = '/projects/brun1463/ROMS/Kakak3_Alpine/Include/wave_forcing_file_kaktovik_shelf_era5_2020_data001.nc'  # UPDATE PATH

#create file to write to
nc1 = Dataset(wave_frc, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006_thin_sponge.nc',
                      type = 'ROMS grid wave forcing file',
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
nc1.createDimension('wave_time', None)
nc1.createDimension('one',     1)


# Create variables # this took several minutes (~9)
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

# wave_time (in seconds)
wave_time_g = nc1.createVariable('wave_time', None, ('wave_time'), zlib=True)
wave_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
wave_time_g.units = 'second'
wave_time_g.field = 'time, scalar, series'
wave_time_g[:] = time_tmp[:]
    
# =============================================================================
# # time (HYCOM version)
# time_g2 = nc1.createVariable('time_HYCOM', None, ('time_HYCOM'), zlib=True)
# time_g2.long_name = '3-hour time steps' 
# time_g2.units = 'datetime'
# time_g2.field = 'time_HYCOM, scalar, series'
# time_g2[:] = time_tmp2[:]
# =============================================================================

# --------------------
# Waves
# --------------------
# ************* copying and editing HYCOM2ROMS_salt_bryclm.py *****************8
# Significant wave height (Hwave)
swh_interp_g = nc1.createVariable('Hwave', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
swh_interp_g.long_name = 'wind-induced significant wave height'
swh_interp_g.units = 'meter' 

# REAL NAME OF WAVE DATA FROM ERA5
# Mean wind-induced wave direction (Dwave)
dir_interp_g = nc1.createVariable('Dwave', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
dir_interp_g.long_name = 'wind-induced wave direction - mean'
dir_interp_g.units = 'degrees' # **CHECK THE CONVENTION

# *** FAKE NAME OF WAVE DATA FROM ERA5 ***
# Mean wind-induced wave direction (Dwave)
# Renaming as peak to test something in ROMS
pdir_interp_g = nc1.createVariable('Dwavep', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
pdir_interp_g.long_name = 'wind-induced wave direction - peak'
pdir_interp_g.units = 'degrees' # **CHECK THE CONVENTION - I think this is right
# ****************************************

# Peak wind-induced surface wave period (Pwave_top)
pwavet_interp_g = nc1.createVariable('Pwave_top', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
pwavet_interp_g.long_name = 'wind-induced peak surface wave period'
pwavet_interp_g.units = 'second' 

# # Bottom orbital velocity 
# ubr_g = nc1.createVariable('Uwave_rms', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
# ubr_g.long_name = 'wind-induced bottom orbital velocity'
# ubr_g.units = 'meter second-1'

# # Bottom wave period 
# pwaveb_g = nc1.createVariable('Pwave_bot', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
# pwaveb_g.long_name = 'wind-induced bottom wave period'
# pwaveb_g.units = 'second'

# ------------------------------- End netCDF file setup ---------------------------


# Set the input and output grids, and sepcify the lat/lon
# Since we are looking at waves, we will use lon_rho and lat_rho as the primary lat/lon for the grid 
# Input grid (era5)
ds_in_cop = wave_data.copy() 
#ds_in_era5['lon_360'] = ds_in_hycom.lon.values
ds_in_cop['lon'] = (('rlat', 'rlon'), ds_in_cop.lon.values)
ds_in_cop['lat'] = (('rlat', 'rlon'), ds_in_cop.lat.values)

# Output grid (ROMS rho grid)
#ds_out_rho = grid_vertical
ds_out_rho = grid.copy()
ds_out_rho['lat'] = (('eta_rho', 'xi_rho'), ds_out_rho.lat_rho.values)
ds_out_rho['lon'] = (('eta_rho', 'xi_rho'), ds_out_rho.lon_rho.values)

# Add masks 
# ex: ds["mask"] = xr.where(~np.isnan(ds["zeta"].isel(ocean_time=0)), 1, 0)
# Input grid (copernicus)
# this is only a surface mask - which is what we want 
ds_in_cop_mask = xr.where(~np.isnan(ds_in_cop['VHM0'][0,:,:].values), 1, 0) 
print('got through nan mask', flush=True)
ds_in_cop['mask'] = (('lat', 'lon'), ds_in_cop_mask)

# Output grid (ROMS rho grid)
ds_out_rho['mask'] = (('eta_rho', 'xi_rho'), ds_out_rho.mask_rho.values)

# Regrid from era5 grid to rho grid with the masks included and extrapolation used 
regridder_cop2rho = xe.Regridder(ds_in_cop, ds_out_rho, method="bilinear", extrap_method='inverse_dist') #extrap_method="nearest_s2d"
regridder_cop2rho

# Save the weights - only need to do this once
fn_cop2rho = regridder_cop2rho.to_netcdf('regrid_cop2rho_weights.nc')
#print(fn_era52rho)

# Now use the regridder/weights to regrid the significant wave height
dr_cop2rho_Thgt = wave_data['VHM0'][:,:,:].copy() # total significant wave height 
dr_out_cop2rho_Thgt = regridder_cop2rho(dr_cop2rho_Thgt) 
dr_out_cop2rho_Thgt

# Now use the regridder/weights to regrid the mean wave direction
dr_cop2rho_Tdir = wave_data['VPED'][:,:,:].copy() # peak wave direction 
dr_out_cop2rho_Tdir = regridder_cop2rho(dr_cop2rho_Tdir) 
dr_out_cop2rho_Tdir

# Now use the regridder/weights to regrid the peak wave period
dr_cop2rho_Tper = wave_data['VTPK'][:,:,:].copy() # total peak wave period 
dr_out_cop2rho_Tper = regridder_cop2rho(dr_cop2rho_Tper) 
dr_out_cop2rho_Tper 

# Use fill.f90 to fill the nans in the array
import fill

# Define a function to call to do the filling, taken from model2roms
def laplacefilter(field, threshold, toxi, toeta):
    undef = 2.0e+35 
    tx = 0.9 * undef
    critx = 0.01
    cor = 1.6
    mxs = 10

    field = np.where(abs(field) > threshold, undef, field)

    field = fill.extrapolate.fill(int(1), int(toxi),
                                int(1), int(toeta),
                                float(tx), float(critx), float(cor), float(mxs),
                                np.asarray(field, order='F'),
                                int(toxi),
                                int(toeta))
    return field

# Loop through time to fill in the nans
# Make some variables first
toxi = xi_rho_len
toeta = eta_rho_len

# Make an array to hold the new data without nans
#print('got here 1')
#dr_out_era52rho_zeta_nonan = np.empty((time_len, eta_rho_len, xi_rho_len))
# Significant wave height
dr_out_cop2rho_Thgt_nonan = np.empty((time_len, eta_rho_len, xi_rho_len))

# Mean wave direction
dr_out_cop2rho_Tdir_nonan = np.empty((time_len, eta_rho_len, xi_rho_len)) 

# Peak wave period 
dr_out_cop2rho_Tper_nonan = np.empty((time_len, eta_rho_len, xi_rho_len))
#print('got here 2')

# Make a copy of the OG array to work with
#print('got here 3')
# Significant wave height
dr_out_cop2rho_Thgt_cp1 = dr_out_cop2rho_Thgt.copy() 

# Mean wave direction
dr_out_cop2rho_Tdir_cp1 = dr_out_cop2rho_Tdir.copy() 

# Peak wave period 
dr_out_cop2rho_Tper_cp1 = dr_out_cop2rho_Tper.copy() 
#print('got here 4')


# Loop through depth to replace all the nans with real values 
# Loop through time
for t in range(time_len):
    # Print the time 
    print('t: ', t, flush=True)

    # Pull out the horizontal 'field' for that time
    # Significant wave height
    #field = dr_out_hycom2rho_zeta_cp1[t,:,:]
    field1 = dr_out_cop2rho_Thgt_cp1[t,:,:]
    
    # Mean wave direction 
    field2 = dr_out_cop2rho_Tdir_cp1[t,:,:]
    
    # Peak wave period
    field3 = dr_out_cop2rho_Tper_cp1[t,:,:]


    # Use the Laplace Filter to get rid of nans
    #field = laplacefilter(field, 1000, toxi, toeta)
    field1 = laplacefilter(field1, 1000, toxi, toeta) # Thgt
    field2 = laplacefilter(field2, 1000, toxi, toeta) # Tdir
    field3 = laplacefilter(field3, 1000, toxi, toeta) # Tper
    

    # Multiply by the rho mask 
    #field = field * ds_out_rho.mask.values
    field1 = field1 * ds_out_rho.mask.values # Thgt
    field2 = field2 * ds_out_rho.mask.values # Tdir
    field3 = field3 * ds_out_rho.mask.values # Tper

    # Check to see if there are any nans
    #print('nans: ', np.where(np.isnan(field)))
    #print('nanmin: ', np.nanmin(field))
    #print('nanmax: ', np.nanmax(field))
    #input('press enter to continue...')

    # Save this field to a new array
    #dr_out_hycom2rho_zeta_nonan[t,:,:] = field
    dr_out_cop2rho_Thgt_nonan[t,:,:] = field1 # Thgt
    dr_out_cop2rho_Tdir_nonan[t,:,:] = field2 # Tdir
    dr_out_cop2rho_Tper_nonan[t,:,:] = field3 # Tper

    # Save this to the output array and to the netcdf
    #zeta_interp_g[:,:,:] = dr_out_hycom2rho_zeta_nonan[:,:,:]
    swh_interp_g[t,:,:] = dr_out_cop2rho_Thgt_nonan[t,:,:] # Thgt
    dir_interp_g[t,:,:] = dr_out_cop2rho_Tdir_nonan[t,:,:] # Tdir
    pdir_interp_g[t,:,:] = dr_out_cop2rho_Tdir_nonan[t,:,:] # Tdir
    pwavet_interp_g[t,:,:] = dr_out_cop2rho_Tper_nonan[t,:,:] # Tper

    # Force save to the netcdfs
    nc1.sync()


# Close the netcdfs
nc1.close()


