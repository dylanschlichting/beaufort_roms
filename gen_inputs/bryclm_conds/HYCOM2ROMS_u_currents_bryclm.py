################################ xESMF Interpolate u Currents from HYCOM to ROMS ###################################
# The purpose of this script is to vertically interpolate the u currents from
# HYCOM onto the ROMS grid to make the boundary and climatology forcing files.

###### IMPORTANT NOTES ######
# - HYCOM vertical grid has 40 vertical layers, with depth[0] = 0 (water surface)
#   and depth[-1] = 5000.0 (bottom/beyond bottom/seabed/deepest depth)
# - ROMS vertical grid has 20 vertical layers, with depth[0] = seabed/deepest depth
#   and depth[-1] = shallow/water surface
# - Remember that HYCOM lat/lon indices are not the same as ROMS grid lat/lon indices
#   - Important to keep in mind when comparing results
# - The model2roms workflow is to do the xesmf regridding, then run fill.f90 to get rid of nans,
#   then run intepolation.f90 after to make up for the different vertical resolutions 
# - To run on summit, use a highmem node with 400GB RAM for a run time of 3 hours
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


# NEW
# Load in the HYCOM currents data
# u currents 
hycom_u_all = xr.open_dataset('/pscratch/sd/b/bundzis/External_data/HYCOM_data/bryclm_data/U/hycom_u_2019_2024_daily_nogap.nc')  
# v currents 
hycom_v_all = xr.open_dataset('/pscratch/sd/b/bundzis/External_data/HYCOM_data/bryclm_data/V/hycom_v_2019_2024_daily_nogap.nc')  

print('loaded in the HYCOM u and v data', flush=True)

# Cut out just the time period we care about for the model - the open water season
# JK we want all time so just use it all and comment this part out
# u currents 
# hycom_u_all = hycom_u_all.sel(time=slice('2019-07-01', '2019-11-01 21:00:00')) # UTC
#hycom_u_all = hycom_u_all.sel(time=slice('2020-07-01 09:00:00', '2020-11-02 09:00:00')) # this is 2019-07-01 01:00:00 - 2019-11-02 01:00:00 in AKDT

# v currents 
#hycom_v_all = hycom_v_all.sel(time=slice('2019-07-01', '2019-11-01 21:00:00')) # UTC
#hycom_v_all = hycom_v_all.sel(time=slice('2020-07-01 09:00:00', '2020-11-02 09:00:00')) # this is 2019-07-01 01:00:00 - 2019-11-02 01:00:00 in AKDT

#print('hycom time length: ', len(hycom_u_all.time.values))
#print(hycom_u_all.time[-6:-1].values)

# Grid's lat/lon is in different convention than HYCOM lat/lon 
# need to make HYCOM match grid's lat lon convention 
hycom_u_all['lon_180'] = -(360 - hycom_u_all.lon.values)
hycom_v_all['lon_180'] = -(360 - hycom_v_all.lon.values)


# Load in the ROMS grid vertical coordinates
# OG grid with 20 vertical layers
grid_vertical = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/ROMS_grid_depth_hpluszeta_2019_2024_20vert_001.nc', drop_variables='z_w')  #UPDATE PATH 
# New grid with 30 vertical layers
#grid_vertical = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/ROMS_grid_depth_hpluszeta_2019_2024_30vert_001.nc', drop_variables='z_w')  #UPDATE PATH 

print('loaded in the vertical grid', flush=True)

# Pull out the angle to rotate the currents to match the grid's u,v
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
# --- AKDT version --- 
# Make a datetime array to use
# OLD
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

# copy time from grid_vertical since we did this there
time2_tmp_len = len(grid_vertical.time_HYCOM)


# name the variables to fill the netcdf 
# latitude
lat_tmp = grid_vertical.lat_rho.values

# longitude
lon_tmp = grid_vertical.lon_rho.values

# HYCOM number of vertical layers
N_hycom = len(hycom_u_all.depth)
#print(N_hycom)

# HYCOM number of lats
hycom_lat_len = len(hycom_u_all.lat.values)
#print('hycom_lat_len', hycom_lat_len)

# HYCOM number of lons
hycom_lon_len = len(hycom_u_all.lon.values)
#print('hycom_lon_len', hycom_lon_len)


# Make a z_rho to be used in the netcdfs (0 is deepest, values are positive) # this took ~7 minutes
z_rho_pos = grid_vertical.z_rho.values * (-1)

print('Setting up netcdfs', flush=True)

# Set up the netcdf for the climatology 
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# OG grid with 20 vertical layers
vert_u_currents_clm = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/u_currents_clm_2019_2024_20vert_002.nc'  #UPDATE PATH
# New grid with 30 vertical layers
#vert_u_currents_clm = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/u_currents_clm_2019_2024_30vert_001.nc'  #UPDATE PATH

#create file to write to
nc1 = Dataset(vert_u_currents_clm, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006_thin_sponge.nc',
                      type = 'ROMS grid vertically interpolated HYCOM u currents climatology',
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
nc1.createDimension('xi_u',    L)    # u
nc1.createDimension('eta_u',   Mp)
nc1.createDimension('v3d_time', None)
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

# xi u
xi_u = nc1.createVariable('xi_u', 'd', ('xi_u',), zlib = True)
xi_u.long_name = 'xi coordinate of U-points'
xi_u.standard_name = 'projection_xi_coordinate'
xi_u.units = 'meter'
xi_u[:]=np.arange(0,L)

# eta u
eta_u = nc1.createVariable('eta_u', 'd', ('eta_u',), zlib = True)
eta_u.long_name = 'eta coordinate of U-points'
eta_u.standard_name = 'projection_eta_coordinate'
eta_u.units = 'meter'
eta_u[:]=np.arange(0,Mp)

# s rho
s_rho = nc1.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
s_rho.long_name = 's coordinate of RHO-points'
s_rho.standard_name = 'projection_s_coordinate'
s_rho.units = 'meter'
s_rho_tmp = np.arange(0, N)
s_rho[:] = s_rho_tmp[:]

# v3d_time (in seconds)
v3d_time_g = nc1.createVariable('v3d_time', None, ('v3d_time'), zlib=True)
v3d_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
v3d_time_g.units = 'second'
v3d_time_g.field = 'time, scalar, series'
v3d_time_g[:] = time_tmp[:]
    
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
z_r_g = nc1.createVariable('z_rho', 'd', ('v3d_time', 's_rho', 'eta_rho', 'xi_rho'), zlib = True)
z_r_g.long_name = 'Z coordinate of rho-points'
z_r_g.units = 'meter'
z_r_g[:,:,:,:] = z_rho_pos[:,:,:,:]

# --------------------
# Boundary Currents
# --------------------

# u current
u_interp_g = nc1.createVariable('u', 'f8', ('v3d_time', 's_rho', 'eta_u', 'xi_u'), zlib=True)
u_interp_g.long_name = 'water u momentum'
u_interp_g.units = 'meter per second'


# Set up the netcdf for the boundary 
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# OG grid with 20 vertical layers
vert_u_currents_bry = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/u_currents_bry_2019_2024_20vert_002.nc'  #UPDATE PATH
# New grid with 30 vertical layers
#vert_u_currents_bry = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/Attempt001/u_currents_bry_2019_2024_30vert_001.nc'  #UPDATE PATH

#create file to write to
nc2 = Dataset(vert_u_currents_bry, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006_thin_sponge.nc',
                      type = 'ROMS grid vertically interpolated HYCOM u currents boundary',
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
nc2.createDimension('xi_u',    L)    # u
nc2.createDimension('eta_u',   Mp)
nc2.createDimension('v3d_time', None)
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

# xi u
xi_u = nc2.createVariable('xi_u', 'd', ('xi_u',), zlib = True)
xi_u.long_name = 'xi coordinate of U-points'
xi_u.standard_name = 'projection_xi_coordinate'
xi_u.units = 'meter'
xi_u[:]=np.arange(0,L)

# eta u
eta_u = nc2.createVariable('eta_u', 'd', ('eta_u',), zlib = True)
eta_u.long_name = 'eta coordinate of U-points'
eta_u.standard_name = 'projection_eta_coordinate'
eta_u.units = 'meter'
eta_u[:]=np.arange(0,Mp)

# s rho
s_rho = nc2.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
s_rho.long_name = 's coordinate of RHO-points'
s_rho.standard_name = 'projection_s_coordinate'
s_rho.units = 'meter'
s_rho_tmp = np.arange(0, N)
s_rho[:] = s_rho_tmp[:]

# v3d_time (in seconds)
v3d_time_g = nc2.createVariable('v3d_time', None, ('v3d_time'), zlib=True)
v3d_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
v3d_time_g.units = 'second'
v3d_time_g.field = 'time, scalar, series'
v3d_time_g[:] = time_tmp[:]
    
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
z_r_g = nc2.createVariable('z_rho', 'd', ('v3d_time', 's_rho', 'eta_rho', 'xi_rho'), zlib = True)
z_r_g.long_name = 'Z coordinate of rho-points'
z_r_g.units = 'meter'
z_r_g[:,:,:,:] = z_rho_pos[:,:,:,:]

# --------------------
# Boundary Currents
# --------------------

# u current
# u_west
u_west_g = nc2.createVariable('u_west', 'f8', ('v3d_time', 's_rho', 'eta_u'), zlib=True)
u_west_g.long_name = 'water u momentum western boundary'
u_west_g.units = 'meter per second'

# u_north
u_north_g = nc2.createVariable('u_north', 'f8', ('v3d_time', 's_rho', 'xi_u'), zlib=True)
u_north_g.long_name = 'water u momentum northern boundary'
u_north_g.units = 'meter per second'

# u_east
u_east_g = nc2.createVariable('u_east', 'f8', ('v3d_time', 's_rho', 'eta_u'), zlib=True)
u_east_g.long_name = 'water u momentum eastern boundary'
u_east_g.units = 'meter per second'

# ------------------------------------------- End netcdf set up -------------------------------------------

print('done with netcdf setup', flush=True)

# Delete the z_rho used in the netcdf to save memory
del(z_rho_pos)

# Now before we can regird and such, we need to rotate the HYCOM data
# Define a function to do this rotation 
def rotate2(point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    
    This equation can be verified here: https://www.myroms.org/forum/viewtopic.php?f=3&t=295 
    where we are rotating from lon,lat to the same lon,lat (rotating HYCOM data about an angle)
    """
    px, py = point

    qx = math.cos(angle) * (px) - math.sin(angle) * (py)
    qy = math.sin(angle) * (px) + math.cos(angle) * (py)
    return qx, qy

# Apply this rotation to all u and v in HYCOM
# Declare empty arrays to fill with the rotated values 
water_u_rot = np.empty((time_len,N_hycom,hycom_lat_len,hycom_lon_len))
water_v_rot = np.empty((time_len,N_hycom,hycom_lat_len,hycom_lon_len))

print('rotating u and v', flush=True)

# Rotate the currents and save the values to the arrays; this takes ~1 minute
water_u_rot[:,:,:,:], water_v_rot[:,:,:,:] = rotate2((hycom_u_all.water_u[:,:,:,:].values, hycom_v_all.water_v[:,:,:,:].values), phi)

# Check to see if this rotation worked 
# check u
#print('Unrotated u: ', hycom_u_all.water_u[0,10,80,200].values)
#print('Expected rotated u: ', math.cos(phi) * (hycom_u_all.water_u[0,10,80,200].values) - math.sin(phi) * (hycom_v_all.water_v[0,10,80,200].values))
#print('Actual Rotated u: ', water_u_rot[0,10,80,200])

# check v
#print('Unrotated v: ', hycom_v_all.water_v[0,10,80,200].values)
#print('Expected rotated v: ', math.sin(phi) * (hycom_u_all.water_u[0,10,80,200].values) + math.cos(phi) * (hycom_v_all.water_v[0,10,80,200].values))
#print('Actual Rotated v: ', water_v_rot[0,10,80,200])

# This rotation seemed to work correctly and was super fast!! So yay! 
# Now make these rotated currents part of the HYCOM data set
# ex: ds['nmap'] = (('y', 'x'), nmap)
hycom_u_all['water_u_rot'] = (('time', 'depth', 'lat', 'lon'), water_u_rot)
hycom_v_all['water_v_rot'] = (('time', 'depth', 'lat', 'lon'), water_v_rot)

# Set the input and output grids, and sepcify the lat/lon
# Since we are looking at u for now, we will use lon_u and lat_u as the primary lat/lon for the grid 
# Input grid (HYCOM)
ds_in_hycom = hycom_u_all.copy() # need to use lon_180 for this grid 
ds_in_hycom['lon_360'] = ds_in_hycom.lon.values
ds_in_hycom['lon'] = ds_in_hycom.lon_180.values

# Output grid (ROMS u grid, but keeps HYCOM vertical levels)
ds_out_u = grid_vertical.copy()
ds_out_u['lat'] = (('eta_u', 'xi_u'), ds_out_u.lat_u.values)
ds_out_u['lon'] = (('eta_u', 'xi_u'), ds_out_u.lon_u.values)

# Add masks 
# ex: ds["mask"] = xr.where(~np.isnan(ds["zeta"].isel(ocean_time=0)), 1, 0)
# Input grid (HYCOM)
# this is only a surface mask - which is what we want
ds_in_hycom_mask = xr.where(~np.isnan(ds_in_hycom['water_u_rot'][0,0,:,:].values), 1, 0) 
ds_in_hycom['mask'] = (('lat', 'lon'), ds_in_hycom_mask)

# Output grid (ROMS u)
ds_out_u['mask'] = (('eta_u', 'xi_u'), ds_out_u.mask_u.values)

# Regrid from HYCOM grid to u grid with the masks included and extrapolation used 
regridder_hycom2u = xe.Regridder(ds_in_hycom, ds_out_u, method="bilinear", extrap_method='nearest_s2d') #extrap_method="nearest_s2d"
regridder_hycom2u

# Save the weights - only need to do this once
#fn_hycom2u = regridder_hycom2u.to_netcdf('regrid_hycom2u_weights.nc')
#print(fn_hycom2u)

print('regridding u', flush=True)

# Now use the regridder/weights to regrid the pre-rotated water u 
dr_hycom2u_u = hycom_u_all['water_u_rot'].copy()
dr_out_hycom2u_u = regridder_hycom2u(dr_hycom2u_u) # this takes ~10 minutes
dr_out_hycom2u_u

# Prep the data to be input for interpolation.f90
# Depths must be negative so multiply them by -1
# Save HYCOM depths as arrays so we can clear the HYCOM data
# HYCOM
# Memory-conscious new way
hycom_depth = np.asarray(hycom_u_all.depth.values)
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
# delete u files
del(hycom_u_all)

# delete v files
del(hycom_v_all)

# Input data must be arranged with deepest value at highest index
# HYCOM already does this
# Do this for ROMS grid (z_rho)
z_rho_flip = np.flip(grid_vertical.z_rho[:,:,:,:].values, axis=1) # better for memory?

# z_rho has rho dimensions instead of u dimensions...
# So we need to make z_rho the correct shape 
# same with bathymetry 
# Let's do this using the xesmf regridding 
# ds_out_u is good to keep the same - this is the u grid
# Need our input grid - rho grid
ds_in_rho = grid_vertical.copy()
ds_in_rho['lat'] = (('eta_rho', 'xi_rho'), ds_in_rho.lat_rho.values)
ds_in_rho['lon'] = (('eta_rho', 'xi_rho'), ds_in_rho.lon_rho.values)

# Add masks 
# ex: ds["mask"] = xr.where(~np.isnan(ds["zeta"].isel(ocean_time=0)), 1, 0)
# Input grid (ROMS rho)
ds_in_rho['mask'] = (('eta_rho', 'xi_rho'), ds_in_rho.mask_rho.values)

# Regrid from ROMS rho grid to u grid with the masks included and extrapolation used 
regridder_rho2u = xe.Regridder(ds_in_rho, ds_out_u, "bilinear", extrap_method='nearest_s2d') #extrap_method="nearest_s2d"
regridder_rho2u
#print(regridder_rho2u.shape)

# Save the weights - only need to do this once
#fn_rho2u = regridder_rho2u.to_netcdf('regrid_rho2u_weights.nc')

# Now use the regridder/weights to regrid the bathymetry from rho to u points 
dr_rho2u_h = grid_vertical['h'].copy()
dr_out_rho2u_h = regridder_rho2u(dr_rho2u_h) # this took a few seconds 
dr_out_rho2u_h

# Now use the regridder/weights to regrid z_rho from rho to u points 
#dr_rho2u_zrho_flip = np.copy(z_rho_flip)
dr_rho2u_zrho_flip = z_rho_flip
dr_out_rho2u_zrho_flip = regridder_rho2u(dr_rho2u_zrho_flip) # this took ~ 4 minutes ish 
dr_out_rho2u_zrho_flip

# Now use the regridder/weights to regrid the hpluszeta bathymetry from rho to u points 
dr_rho2u_hpluszeta = grid_vertical['total_depth'].copy()
dr_out_rho2u_hpluszeta = regridder_rho2u(dr_rho2u_hpluszeta) # this took a few seconds 
dr_out_rho2u_hpluszeta

# Now I think we can delete grid_vertical to save memory space
del(grid_vertical)

print('importing fortran functions', flush=True)

# This works to import interpolation.f90 as a python package!
# and this works even though it is in a different location 
# because we specify the path below
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
toxi = xi_u_len
toeta = eta_u_len

# Make an array to hold the new data without nans
#print('got here 1')
dr_out_hycom2u_u_nonan = np.empty((time_len, N_hycom, eta_u_len, xi_u_len))
#print('got here 2')

# Make a copy of the OG array to work with
#print('got here 3')
dr_out_hycom2u_u_cp1 = dr_out_hycom2u_u.copy() # this takes ~2 minutes
#print('got here 4')

# OLD SERIAL
# Loop through time & depth to replace all the nans with real values 
# Loop through time
for t in range(time_len):
    # Print the time 
    print('t: ', t)
    
    # Loop through depth
    for k in range(N_hycom):
        # Print the level we are on
        #print(k)

        # Pull out the horizontal 'field' for that level
        field = dr_out_hycom2u_u_cp1[t,k,:,:]

        # Use the Laplace Filter to get rid of nans
        field = laplacefilter(field, 1000, toxi, toeta)

        # Multiply by the u mask 
        field = field * ds_out_u.mask.values

        # Check to see if there are any nans
        #print('nans: ', np.where(np.isnan(field)))
        #print('nanmin: ', np.nanmin(field))
        #print('nanmax: ', np.nanmax(field))
        #input('press enter to continue...')

        # Save this field to a new array
        dr_out_hycom2u_u_nonan[t,k,:,:] = field


# NEW PARALLEL
# from joblib import Parallel, delayed

# print('started parallal Laplace filter', flush=True)

# # Make the above loop a function so that it can be parallelized
# # Define the function for one (t, k) pair
# def replace_nan_1time(t, k):

#     # Pull out the horizontal 'field' for that time and level
#     field = dr_out_hycom2u_u_cp1[t, k, :, :]

#     # Use the Laplace Filter to get rid of nans
#     field = laplacefilter(field, 1000, toxi, toeta)

#     # Multiply by the u mask 
#     field = field * ds_out_u.mask.values

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
#     dr_out_hycom2u_u_nonan[t, k, :, :] = field



# Now use the interpolation.f90 functions to vertically interpolate
# from HYCOM vertical grid to ROMS vertical grid

# Make an array to hold the output 
u_rot_interp_vert_nonan = np.empty((N, eta_u_len, xi_u_len))
# For time loop, will need to save to 4d array so set that up here
u_roms = np.empty((time_len, N, eta_u_len, xi_u_len))

# OLD SERIAL
# Now call the function!
# Loop through time
for tt in range(time_len):
    # Print the time we are on
    print('tt: ', tt)
    
    # Interpolate this time step
    u_rot_interp_vert_nonan = vertInterp.interpolation.dovertinter(u_rot_interp_vert_nonan, dr_out_hycom2u_u_nonan[tt,:,:,:], 
                                                                   dr_out_rho2u_h[:,:], dr_out_rho2u_zrho_flip[tt,:,:,:], 
                                                                   hycom_depth_interp) 
    
    # Now flip the axis to be ROMS convention
    u_rot_interp_vert_nonan_flip = np.flip(u_rot_interp_vert_nonan, axis=0)
    
    ############### TEMPORARY TO TRY TO GET RID OF MASK PROBLEM ################### 
    # This works so we are leaving it in
    # Multiply by the mask one more time before saving 
    # Create an array to hold the ouput
    u_rot_interp_vert_nonan_flip_m = np.empty((N,eta_u_len,xi_u_len))
    
    # Loop through vertical levels
    for kk in range(N):
        # Multiply this layer by the mask 
        u_rot_interp_vert_nonan_flip_m[kk,:,:] = u_rot_interp_vert_nonan_flip[kk,:,:] * ds_out_u.mask.values
        
    ########### END TEMPORARY ##################
    
    # Save this to the output array and to the netcdf
    # to array
    u_roms[tt,:,:,:] = u_rot_interp_vert_nonan_flip_m
    
    # to climatology netcdf
    u_interp_g[tt,:,:,:] = u_rot_interp_vert_nonan_flip_m
    
    # to boundary netcdf
    u_west_g[tt,:,:] = u_rot_interp_vert_nonan_flip_m[:,:,0]
    u_north_g[tt,:,:] = u_rot_interp_vert_nonan_flip_m[:,-1,:]
    u_east_g[tt,:,:] = u_rot_interp_vert_nonan_flip_m[:,:,-1]
    
    # Force save to the netcdfs
    nc1.sync()
    nc2.sync()


# # # NEW PARALLEL
# from tqdm import tqdm

# print('started parallal vertical interpolation', flush=True)

# def vert_interp_1time(tt):
#     print(f"Processing time step: {tt}")

#     # Get the data for this time step to be given to the function
#     input_field = dr_out_hycom2u_u_nonan[tt, :, :, :]
#     rho2u_h = dr_out_rho2u_h[:, :]
#     rho2u_zrho = dr_out_rho2u_zrho_flip[tt, :, :, :]
    
#     # Interpolate this time step
#     interp = vertInterp.interpolation.dovertinter(
#         np.empty((N, eta_u_len, xi_u_len)), 
#         input_field, 
#         rho2u_h, 
#         rho2u_zrho, 
#         hycom_depth_interp
#     )
    
#     # Now flip the axis to be ROMS convention
#     interp_flip = np.flip(interp, axis=0)
    
#     # Multiply by the mask one more time before saving
#     masked = np.empty((N, eta_u_len, xi_u_len))
#     mask = ds_out_u.mask.values
#     for kk in range(N):
#         masked[kk, :, :] = interp_flip[kk, :, :] * mask

#     # Return all variabels and slices for writing to netcdf
#     return tt, masked, masked[:, :, 0], masked[:, -1, :], masked[:, :, -1]

# # Call this function in parallel
# results2 = Parallel(n_jobs=120)(
#     delayed(vert_interp_1time)(tt) 
#     for tt in tqdm(range(time_len), desc="Interpolating time steps")
# )

# # Now save the results to the netcdfs
# for tt, full_field, west_slice, north_slice, east_slice in results2:
#     # In-memory array
#     u_roms[tt, :, :, :] = full_field

#     # Climatology NetCDF
#     u_interp_g[tt, :, :, :] = full_field

#     # Boundary NetCDFs
#     u_west_g[tt, :, :] = west_slice
#     u_north_g[tt, :, :] = north_slice
#     u_east_g[tt, :, :] = east_slice

# # Sync once at the end
# nc1.sync()
# nc2.sync()


# Close the netcdfs
nc1.close()
nc2.close()





