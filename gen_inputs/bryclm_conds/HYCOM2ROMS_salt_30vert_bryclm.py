########################### xESMF Interpolate Salinity from HYCOM to ROMS ###########################
# The purpose of this script is to use the wonderful xesmf package to interpolate salinity
# from HYCOM data to the ROMS grid. This script then saves these values to the bry and clm files.

###### IMPORTANT NOTES ######
# - HYCOM vertical grid has 40 vertical layers, with depth[0] = 0 (water surface)
#   and depth[-1] = 5000.0 (bottom/beyond bottom/seabed/deepest depth)
# - ROMS vertical grid has 20 vertical layers, with depth[0] = seabed/deepest depth
#   and depth[-1] = shallow/water surface
# - Remember that HYCOM lat/lon indices are not the same as ROMS grid lat/lon indices
#   - Important to keep in mind when comparing results
# - The model2roms workflow is to do the xesmf regridding then run the intepolation.f90
#   after to make up for the different vertical resolutions 
# - To run on summit, use a highmem node with 200GB RAM for a run time of 2 hours
#   - It should only take ~1 hour but it's safest to give it more in case

##############################

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


# Load in the HYCOM salinity data
hycom_salt_all = xr.open_dataset('/pscratch/sd/b/bundzis/External_data/HYCOM_data/bryclm_data/Salt/hycom_salt_2019_2024_daily_nogap.nc')
# only 2019 - 2020 for now
#hycom_salt_all = xr.open_dataset('/pscratch/sd/b/bundzis/External_data/HYCOM_data/bryclm_data/Salt/hycom_salt_2019_2020_daily_nogap.nc')

# Cut out just the time period we care about for the model - the open water season
# JK we want all time so just use it all and comment this part out
# salt
#hycom_salt_all = hycom_salt_all.sel(time=slice('2019-07-01','2019-11-01 21:00:00')) # UTC
#hycom_salt_all = hycom_salt_all.sel(time=slice('2020-07-01 09:00:00', '2020-11-02 09:00:00')) # this is 2019-07-01 01:00:00 - 2019-11-02 01:00:00 in AKDT
#print('hycom time length: ', len(hycom_salt_all.time.values))
#print(hycom_salt_all.time[-6:-1].values)
#input('press enter to continue...')

# Grid's lat/lon is in different convention than HYCOM lat/lon 
# need to make HYCOM match grid's lat lon convention 
hycom_salt_all['lon_180'] = -(360 - hycom_salt_all.lon.values)

# Load in the ROMS grid vertical coordinates
# OG grid with 20 vertical layers
#grid_vertical = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/ROMS_grid_depth_hpluszeta_2019_2024_20vert_001.nc', drop_variables='z_w')  #UPDATE PATH 
# New grid with 30 vertical layers
grid_vertical = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/ROMS_grid_depth_hpluszeta_2019_2024_30vert_001.nc', drop_variables='z_w')  #UPDATE PATH 


# # TEMPORARY
# # Subset grid_vertical to just the times that match the salt data since we 
# # don't have all of the salt data yet
# # Make time a datetime 
# reference_time = pd.Timestamp('2000-01-01')
# grid_vertical['time'] = reference_time + pd.to_timedelta(grid_vertical.time.values, unit='s')
# print('grid_vertical.time: ', grid_vertical.time, flush=True)
# grid_vertical = grid_vertical.sel(time=slice('2019-01-01', '2020-12-31'))
# print('new grid_vertical time: ', grid_vertical.time, flush=True)


# pull out the angle to rotate the currents to match the grid's u,v
phi = grid_vertical.angle[0,0].values # radians 

# Read in the dimensions
time_len = len(grid_vertical.time)
eta_rho_len = len(grid_vertical.eta_rho) # 206
xi_rho_len = len(grid_vertical.xi_rho) # 608
eta_u_len = len(grid_vertical.eta_u) # 206
xi_u_len = len(grid_vertical.xi_u) # 607
eta_v_len = len(grid_vertical.eta_v) # 205
xi_v_len = len(grid_vertical.xi_v) # 608

# Define other dimension lengths
# eta rho
Mp = len(grid_vertical.eta_rho)

# xi rho
Lp = len(grid_vertical.xi_rho)

# for u/v points 
Lm, Mm = (Lp-2), (Mp-2) #number/dimension of cells
L,  M  = Lm+1, Mm+1 #number/dimension of psi points

# srho
N = len(grid_vertical.s_rho)

# latitude
lat_len = len(grid_vertical.lat_rho)

# longitude
lon_len = len(grid_vertical.lon_rho)

# time
# Make a datetime array to use
# OLD
# time_akdt = np.arange(datetime(2020,7,1,hour=1,minute=0,second=0), datetime(2020,11,2,hour=4,minute=0, second=0),timedelta(hours=3))
# time_akdt_dt = pd.to_datetime(time_akdt)
# time_tmp_len = len(time_akdt_dt)
# datetime1 = time_akdt_dt
# NEW - might need to shorten this to match however long the salt data are
# Full time
time_akdt = np.arange(datetime(2019,1,1,hour=1,minute=0,second=0), datetime(2024,9,6,hour=1,minute=0, second=0),timedelta(days=1))
time_akdt_dt = pd.to_datetime(time_akdt)
time_tmp_len = len(time_akdt_dt)
datetime1 = time_akdt_dt
# Shortened to match 2019 - 2020 data
# time_akdt = np.arange(datetime(2019,1,1,hour=1,minute=0,second=0), datetime(2021,1,1,hour=1,minute=0, second=0),timedelta(days=1))
# time_akdt_dt = pd.to_datetime(time_akdt)
# time_tmp_len = len(time_akdt_dt)
# datetime1 = time_akdt_dt

# Print time lengths to compare 
print('grid_vertical time len: ', time_len, flush=True)
print('grid_vertical times[0:5]: ', grid_vertical.time[0:5].values, flush=True)
print('grid_vertical time[-1]: ', grid_vertical.time[-1].values, flush=True)
print('akdt_time len: ', time_tmp_len, flush=True)
print('time_akdt time[0:5]: ' , time_akdt_dt[0:5], flush=True)
print('time_akdt time[-1]: ' , time_akdt_dt[-1], flush=True)

# Convert all time to seconds since 2000-01-01 (really 1999-12-31 so it starts at beginning
# of year=hour 0)
time_tmp = ((datetime1[:] - datetime(1999,12,31)).total_seconds() - 86400)

# copy this from grd_vertical since we already did that there
#time2_tmp_len = len(grid_vertical.time_HYCOM)

# copy the actual time values - seconds since 2000-01-01
#time_tmp = grid_vertical.time.values

# HYCOM time 
#time_tmp2 = grid_vertical.time_HYCOM.values

# name the variables to fill the netcdf 
# latitude
lat_tmp = grid_vertical.lat_rho.values

# longitude
lon_tmp = grid_vertical.lon_rho.values

# HYCOM number of vertical layers
N_hycom = len(hycom_salt_all.depth)
#print(N_hycom)

# HYCOM number of lats
hycom_lat_len = len(hycom_salt_all.lat.values)
#print('hycom_lat_len', hycom_lat_len)

# HYCOM number of lons
hycom_lon_len = len(hycom_salt_all.lon.values)
#print('hycom_lon_len', hycom_lon_len)

# Make a z_rho to be used in the netcdfs (0 is deepest, values are positive) # this took ~7 minutes
z_rho_pos = grid_vertical.z_rho.values * (-1)

# Set up the netcdf for the climatology 
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# OG grid with 20 vertical layers
#vert_salt_clm = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/salt_clm_2019_2024_20vert_001.nc'   #UPDATE PATH
# New grid with 30 vertical layers
vert_salt_clm = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/salt_clm_2019_2024_30vert_001.nc'   #UPDATE PATH

#create file to write to
nc1 = Dataset(vert_salt_clm, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006_thin_sponge.nc',
                      type = 'ROMS grid vertically interpolated HYCOM salinity climatology',
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
nc1.createDimension('salt_time', None)
nc1.createDimension('time_HYCOM', None)
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

# s rho
s_rho = nc1.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
s_rho.long_name = 's coordinate of RHO-points'
s_rho.standard_name = 'projection_s_coordinate'
s_rho.units = 'meter'
s_rho_tmp = np.arange(0, N)
s_rho[:] = s_rho_tmp[:]

# salt_time (in seconds)
salt_time_g = nc1.createVariable('salt_time', None, ('salt_time'), zlib=True)
salt_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
salt_time_g.units = 'second'
salt_time_g.field = 'time, scalar, series'
salt_time_g[:] = time_tmp[:]
    
# # time (HYCOM version)
# time_g2 = nc1.createVariable('time_HYCOM', None, ('time_HYCOM'), zlib=True)
# time_g2.long_name = '3-hour time steps' 
# time_g2.units = 'datetime'
# time_g2.field = 'time_HYCOM, scalar, series'
# time_g2[:] = time_tmp2[:]


# --------------------
# Vertical variables
# --------------------

# z_r (depths of rho points)
z_r_g = nc1.createVariable('z_rho', 'd', ('salt_time', 's_rho', 'eta_rho', 'xi_rho'), zlib = True)
z_r_g.long_name = 'Z coordinate of rho-points'
z_r_g.units = 'meter'
z_r_g[:,:,:,:] = z_rho_pos[:,:,:,:]

# --------------------
# Salinity
# --------------------

# salinity
salt_interp_g = nc1.createVariable('salt', 'f8', ('salt_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True)
salt_interp_g.long_name = 'salinity'
salt_interp_g.units = 'PSS' #** Is this the same as PSU?


# Set up the netcdf for the boundary 
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# OG grid with 20 vertical layers
#vert_salt_bry = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/salt_bry_2019_2024_20vert_001.nc'  #UPDATE PATH
# New grid with 30 vertical layers
vert_salt_bry = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/salt_bry_2019_2024_30vert_001.nc'  #UPDATE PATH

#create file to write to
nc2 = Dataset(vert_salt_bry, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006_thin_sponge.nc',
                      type = 'ROMS grid vertically interpolated HYCOM salinity boundary',
                      history = 'Created by Brianna Undzis',
                      Conventions = 'CF',
                      Institution = 'University of Colorado Boulder',
                      date = str(datetime.today()))
    
#create dictionary for model
d = {}
d = global_defaults

for att, value in d.items():
    setattr(nc2, att, value)

# Create dimensions
nc2.createDimension('xi_rho',  Lp)   # RHO
nc2.createDimension('eta_rho', Mp)
nc2.createDimension('s_rho', N)
nc2.createDimension('s_w',   (N+1))
nc2.createDimension('salt_time', None)
nc2.createDimension('time_HYCOM', None)
nc2.createDimension('one',     1)

# Create variables
# --------------------
# Coordinate Variables
# --------------------
# xi rho
xi_rho = nc2.createVariable('xi_rho', 'd', ('xi_rho',), zlib=True)
xi_rho.long_name = 'xi coordinate of RHO-points'
xi_rho.standard_name = 'projection_xi_coordinate'
xi_rho.units = 'meter'
xi_rho_tmp = np.arange(0, Lp)
xi_rho[:] = xi_rho_tmp[:]

# eta rho
eta_rho = nc2.createVariable('eta_rho', 'd', ('eta_rho',), zlib=True)
eta_rho.long_name = 'eta coordinate of RHO-points'
eta_rho.standard_name = 'projection_eta_coordinate'
eta_rho.units = 'meter'
eta_rho_tmp = np.arange(0, Mp)
eta_rho[:] = eta_rho_tmp[:]

# s rho
s_rho = nc2.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
s_rho.long_name = 's coordinate of RHO-points'
s_rho.standard_name = 'projection_s_coordinate'
s_rho.units = 'meter'
s_rho_tmp = np.arange(0, N)
s_rho[:] = s_rho_tmp[:]

# salt_time (in seconds)
salt_time_g = nc2.createVariable('salt_time', None, ('salt_time'), zlib=True)
salt_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
salt_time_g.units = 'second'
salt_time_g.field = 'time, scalar, series'
salt_time_g[:] = time_tmp[:]
    
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
z_r_g = nc2.createVariable('z_rho', 'd', ('salt_time', 's_rho', 'eta_rho', 'xi_rho'), zlib = True)
z_r_g.long_name = 'Z coordinate of rho-points'
z_r_g.units = 'meter'
z_r_g[:,:,:,:] = z_rho_pos[:,:,:,:]

# --------------------
# Boundary Salinity
# --------------------

# salinity
# salt_west
salt_west_g = nc2.createVariable('salt_west', 'f8', ('salt_time', 's_rho', 'eta_rho'), zlib=True)
salt_west_g.long_name = 'salinity at western boundary'
salt_west_g.units = 'PSS'

# salt_north
salt_north_g = nc2.createVariable('salt_north', 'f8', ('salt_time', 's_rho', 'xi_rho'), zlib=True)
salt_north_g.long_name = 'salinity at northern boundary'
salt_north_g.units = 'PSS'

# salt_east
salt_east_g = nc2.createVariable('salt_east', 'f8', ('salt_time', 's_rho', 'eta_rho'), zlib=True)
salt_east_g.long_name = 'salinity at eastern boundary'
salt_east_g.units = 'PSS'

# ---------------------------------------- End netCDF set up ----------------------------------------

# Delete the z_rho used in the netcdf to save memory
del(z_rho_pos)


# Set the input and output grids, and sepcify the lat/lon
# Since we are looking at salt for now, we will use lon_rho and lat_rho as the primary lat/lon for the grid 
# Input grid (HYCOM)
ds_in_hycom = hycom_salt_all.copy() # need to use lon_180 for this grid 
ds_in_hycom['lon_360'] = ds_in_hycom.lon.values
ds_in_hycom['lon'] = ds_in_hycom.lon_180.values

# Output grid (ROMS rho grid, but keeps HYCOM vertical levels)
#ds_out_rho = grid_vertical
ds_out_rho = grid_vertical.copy()
ds_out_rho['lat'] = (('eta_rho', 'xi_rho'), ds_out_rho.lat_rho.values)
ds_out_rho['lon'] = (('eta_rho', 'xi_rho'), ds_out_rho.lon_rho.values)

# Add masks 
# ex: ds["mask"] = xr.where(~np.isnan(ds["zeta"].isel(ocean_time=0)), 1, 0)
# Input grid (HYCOM)
# this is only a surface mask - which is what we want 
ds_in_hycom_mask = xr.where(~np.isnan(ds_in_hycom['salinity'][0,0,:,:].values), 1, 0) 
ds_in_hycom['mask'] = (('lat', 'lon'), ds_in_hycom_mask)

# Output grid (ROMS rho grid)
ds_out_rho['mask'] = (('eta_rho', 'xi_rho'), ds_out_rho.mask_rho.values)

# Regrid from HYCOM grid to rho grid with the masks included and extrapolation used 
regridder_hycom2rho = xe.Regridder(ds_in_hycom, ds_out_rho, method="bilinear", extrap_method='nearest_s2d') #extrap_method="nearest_s2d"
regridder_hycom2rho

# Save the weights - only need to do this once
fn_hycom2rho = regridder_hycom2rho.to_netcdf('regrid_hycom2rho_weights.nc')
#print(fn_hycom2u)

# Now use the regridder/weights to regrid the salinity  
dr_hycom2rho_salt = hycom_salt_all['salinity'].copy()
dr_out_hycom2rho_salt = regridder_hycom2rho(dr_hycom2rho_salt) 
dr_out_hycom2rho_salt

# Prep the data to be input for interpolation.f90
# Depths must be negative so multiply them by -1
# Save HYCOM depths as arrays so we can clear the HYCOM data
# HYCOM
# Memory-conscious new way
hycom_depth = np.asarray(hycom_salt_all.depth.values)
#print(hycom_depth)
hycom_depth_interp = hycom_depth * (-1) 
#print(hycom_depth_interp) 

# ROMS grid
# grid_vertical.z_rho shoould already be negative
#print(grid_vertical.z_rho[0,:,200,200].values)
#print(grid_vertical.z_rho[100,:,160,450].values)
# grid_vertical['z_rho'][:,:,:,:] = grid_vertical.z_rho[:,:,:,:].values * (-1) # this takes ~2 minutes
# print(grid_vertical.z_rho[0,:,200,200].values)
# print(grid_vertical.z_rho[100,:,160,450].values)

# Now that we are done with the HYCOM data, let's remove it from memory!
# delete salt files
del(hycom_salt_all)

# Input data must be arranged with deepest value at highest index
# HYCOM already does this
z_rho_flip = np.flip(grid_vertical.z_rho[:,:,:,:].values, axis=1) # better for memory?

# z_rho and bathymetry are already on rho points so there is
# no need to regrid them 

# This works to import interpolation.f90 as a python package!
# and this works even though it is in a different location 
# because we specify the pathname below
# from numpy import f2py
# with open('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Bryclm_conds/interpolation.f90') as sourcefile:
#     sourcecode = sourcefile.read()
# f2py.compile(sourcecode, modulename='vertInterp', extension='.f90')
import vertInterp

# Use fill.f90 to fill the nans in the array
# Import fill.f90 from model2roms to see how to 
# use this/if it can be used to get rid of nans 
# from numpy import f2py
# with open('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Bryclm_conds/fill.f90') as sourcefile2:
#     sourcecode2 = sourcefile2.read()
# f2py.compile(sourcecode2, modulename='fill', extension='.f90')
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

# Loop through depth levels to fill in the nans
# Make some variables first
toxi = xi_rho_len
toeta = eta_rho_len

# Make an array to hold the new data without nans
#print('got here 1')
dr_out_hycom2rho_salt_nonan = np.empty((time_len, N_hycom, eta_rho_len, xi_rho_len))
#print('got here 2')

# Make a copy of the OG array to work with
#print('got here 3')
dr_out_hycom2rho_salt_cp1 = dr_out_hycom2rho_salt.copy() # this takes ~2 minutes
#print('got here 4')

# OLD SERIAL
# Loop through depth to replace all the nans with real values 
# Loop through time
for t in range(time_len):
    # Print the time 
    print('t: ', t)
    
    # Loop through depth
    for k in range(N_hycom):
        # Print the level we are on
        #print(k)

        # Pull out the horizontal 'field' for that level
        field = dr_out_hycom2rho_salt_cp1[t,k,:,:]

        # Use the Laplace Filter to get rid of nans
        field = laplacefilter(field, 1000, toxi, toeta)

        # Multiply by the rho mask 
        field = field * ds_out_rho.mask.values

        # Check to see if there are any nans
        #print('nans: ', np.where(np.isnan(field)))
        #print('nanmin: ', np.nanmin(field))
        #print('nanmax: ', np.nanmax(field))
        #input('press enter to continue...')

        # Save this field to a new array
        dr_out_hycom2rho_salt_nonan[t,k,:,:] = field


# NEW PARALLEL
# from joblib import Parallel, delayed

# # Make the above loop a function so that it can be parallelized
# # Define the function for one (t, k) pair
# def replace_nan_1time(t, k):

#     # Pull out the horizontal 'field' for that time and level
#     field = dr_out_hycom2rho_salt_cp1[t, k, :, :]

#     # Use the Laplace Filter to get rid of nans
#     field = laplacefilter(field, 1000, toxi, toeta)

#     # Multiply by the u mask 
#     field = field * ds_out_rho.mask.values

#     # Return these values
#     return (t, k, field)
    
# # Call the function in parallel
# results = Parallel(n_jobs=120)(
#     delayed(replace_nan_1time)(t, k) 
#     for t in range(time_len) 
#     for k in range(N_hycom)
# )

# # Save this data to the array
# for t, k, field in results:
#     dr_out_hycom2rho_salt_nonan[t, k, :, :] = field



# Now use the interpolation.f90 functions to vertically interpolate
# from HYCOM vertical grid to ROMS vertical grid

# Make an array to hold the output (3D, only give it one time for now)
salt_interp_vert_nonan = np.empty((N, eta_rho_len, xi_rho_len))
# For time loop, will need to save to 4d array so set that up here
salt_roms = np.empty((time_len, N, eta_rho_len, xi_rho_len))

# OLD SERIAL
# Now call the function!
# Loop through time
for tt in range(time_len):
    # Print the time we are on
    print('tt: ', tt)
    
    # Interpolate this time step
    salt_interp_vert_nonan = vertInterp.interpolation.dovertinter(salt_interp_vert_nonan, dr_out_hycom2rho_salt_nonan[tt,:,:,:], 
                                                                   grid_vertical.h[:,:].values, z_rho_flip[tt,:,:,:], 
                                                                   hycom_depth_interp) 
    
    # Now flip the axis 
    salt_interp_vert_nonan_flip = np.flip(salt_interp_vert_nonan, axis=0)
    
    # Save this to the output array and to the netcdf
    # to array
    salt_roms[tt,:,:,:] = salt_interp_vert_nonan_flip
    
    # to climatology
    salt_interp_g[tt,:,:,:] = salt_interp_vert_nonan_flip
    
    # to boundary
    salt_west_g[tt,:,:] = salt_interp_vert_nonan_flip[:,:,0]
    salt_north_g[tt,:,:] = salt_interp_vert_nonan_flip[:,-1,:]
    salt_east_g[tt,:,:] = salt_interp_vert_nonan_flip[:,:,-1]
    
    # Force save to the netcdfs
    nc1.sync()
    nc2.sync()


# # NEW PARALLEL
# from tqdm import tqdm

# print('started parallal vertical interpolation', flush=True)

# def vert_interp_1time(tt):
#     print(f"Processing time step: {tt}")

#     # Get the data for this time step to be given to the function
#     input_field = dr_out_hycom2rho_salt_nonan[tt, :, :, :]
#     grid_vert_h = grid_vertical.h[:, :]
#     zrho_flip = z_rho_flip[tt, :, :, :]
    
#     # Interpolate this time step
#     interp = vertInterp.interpolation.dovertinter(
#         np.empty((N, eta_v_len, xi_v_len)), 
#         input_field, 
#         grid_vert_h, 
#         zrho_flip, 
#         hycom_depth_interp
#     )
    
#     # Now flip the axis to be ROMS convention
#     interp_flip = np.flip(interp, axis=0)

#     # Return all variabels and slices for writing to netcdf
#     return tt, interp_flip, interp_flip[:, :, 0], interp_flip[:, -1, :], interp_flip[:, :, -1]

# # Call this function in parallel
# results2 = Parallel(n_jobs=120)(
#     delayed(vert_interp_1time)(tt) 
#     for tt in tqdm(range(time_len), desc="Interpolating time steps")
# )

# # Now save the results to the netcdfs
# for tt, full_field, west_slice, north_slice, east_slice in results2:
#     # In-memory array
#     salt_roms[tt, :, :, :] = full_field

#     # Climatology NetCDF
#     salt_interp_g[tt, :, :, :] = full_field

#     # Boundary NetCDFs
#     salt_west_g[tt, :, :] = west_slice
#     salt_north_g[tt, :, :] = north_slice
#     salt_east_g[tt, :, :] = east_slice

# # Sync once at the end
# nc1.sync()
# nc2.sync()




# Close the netcdfs
nc1.close()
nc2.close()







