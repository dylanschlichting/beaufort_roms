#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 12:14:14 2022

@author: brun1463
"""

############# Find ROMS Vertical Deths ###########
# Find the depths of the ROMS model grid
# and play with different settings to get different vertical profiles
# 
#  
# This code has two sections commented out - the 
# first one is used to make the netcdf file so 
# if you are ready to make the netcdf file, 
# uncomment this section; the second is used to 
# print plots of the different vertical settings 
# so you can see how they look before
# committing to them.
##################################################


# Load in packages
import numpy as np
import xarray as xr
from depths import *
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import xarray as xr
#import gridmap.gridmap.projection
from projection import PolarStereographic
from datetime import datetime, timedelta
from cftime import num2date, date2num

# Read in the model data 
model_grid = xr.open_dataset('/global/homes/b/bundzis/Projects/Beaufort_ROMS_2020_test_nosed/Include/KakAKgrd_shelf_big010_smooth006_thin_sponge.nc') # UPDATE PATH
#print(model_grid.head())

# Read in the interpolated zeta
zeta_interp = xr.open_dataset('/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/zeta_clm_2019_2024_001.nc') #UPDATE PATH
#print(zeta_interp.head())

# Add the model bathymetry to the interpolated zeta over all time to get 
# the true depths over time 
# make an empty array to fill with the new values
h_plus_zeta = np.empty((len(zeta_interp.zeta_time), len(zeta_interp.eta_rho), len(zeta_interp.xi_rho))) #(time, eta, xi)
for t in range(len(zeta_interp.zeta_time)):
    h_plus_zeta[t,:,:] = zeta_interp.zeta[t,:,:].values + model_grid.h[:,:].values

# check to see if it worked 
#print('model depth: ', model_grid.h[92,100].values)
#print('interpolated zeta: ', zeta_interp.zeta[103,92,100].values)
#print('h plus zeta: ', h_plus_zeta[103,92,100])


# Assign values to the vertical parameteres
Vtransform = 2 #    OG: 2
Vstretching = 4 #   OG: 3 
theta_s = 5 # 1, 0.65, OG: 1
theta_b = 3 # 1, 3, OG: 3
hc = 5
#hc = 10
N = 30   # number of vertical layers, OG: 20
#igrid = 5

# read in the dimensions
time_len = len(zeta_interp.zeta_time)
eta_len = len(model_grid.eta_rho) # 206?
xi_len = len(model_grid.xi_rho) # 608?

# check dimension shapes
print('time: ', time_len)
print('eta_len: ', eta_len)
print('xi_len: ', xi_len)

# ----- Set up the netCDF file -----
# As the interpolation loops through time, we want to save the output to  
# the netcdf after each time step so that it is not lost if the script 
# crashes/times out. Thus, we will set up the netcdf file before all 
# the interpolation stuff.

#define the dimension lengths
# eta rho
Mp = len(model_grid.eta_rho)

# xi rho
Lp = len(model_grid.xi_rho)

# other dimensions
Lm, Mm = (Lp-2), (Mp-2) #number/dimension of cells
L,  M  = Lm+1, Mm+1 #number/dimension of psi points

# latitude
lat_len = len(model_grid.lat_rho)

# longitude
lon_len = len(model_grid.lon_rho)

# time
# for this, I think we need to make the actualy time steps in terms of the reference time 
#make time, stepping in 3 hours in terms of seconds
time_tmp_len = np.copy(time_len)
time2_tmp_len = len(zeta_interp.time_HYCOM)

#seconds since 2000-01-01 00:00:00 (hour:minute:second)
#seconds_since_init = 615254400.0

# copy time in seconds for from file 
time_tmp = zeta_interp.zeta_time

#or in hours since 2000-01-01 00:00:00
#hours_since_init = 170904
#time_tmp = np.linspace(hours_since_init, hours_since_init+int(time_tmp_len), 
                        #int(time_tmp_len)), endpoint=True)
# HYCOM time 
time_tmp2 = zeta_interp.time_HYCOM

#OR
#time_tmp = UNLIMITED

# name the variables to fill the netcdf
# latitude
lat_tmp = model_grid.lat_rho.values

# longitude
lon_tmp = model_grid.lon_rho.values



### ------ Section 1: Uncomment this section if you want to make the netcdf file ----- ###
#-----Create the netCDF file-----
#name of file I am writing to
grid_depths = '/pscratch/sd/b/bundzis/Beaufort_ROMS_2020_dvd_myroms_ice_scratch/Forcing_files/Bryclm/ROMS_grid_depth_hpluszeta_2019_2024_30vert_001.nc' # UPDATE PATH

# Create a variable for grid mapping
gridmap_varname = 'grid_mapping'

#create file to write to
nc = Dataset(grid_depths, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006_thin_sponge.nc',
                      type = 'ROMS grid depths using HYCOM zeta and IBCAO bathymetry',
                      history = 'Created by Brianna Undzis',
                      Conventions = 'CF',
                      Institution = 'University of Colorado Boulder',
                      date = str(datetime.today()))
    
#create dictionary for model
d = {}
d = global_defaults

for att, value in d.items():
    setattr(nc, att, value)

#Create dimensions
nc.createDimension('xi_rho',  Lp)   # RHO
nc.createDimension('eta_rho', Mp)
nc.createDimension('s_rho', N)
nc.createDimension('s_w', N+1)
nc.createDimension('lat_rho', lat_len)
nc.createDimension('lon_rho', lon_len)

nc.createDimension('xi_u',    L)    # u
nc.createDimension('eta_u',   Mp)

nc.createDimension('xi_v',    Lp)   # v
nc.createDimension('eta_v',   M)

nc.createDimension('xi_psi',  L)    # PSI
nc.createDimension('eta_psi', M)

nc.createDimension('time', None)
nc.createDimension('time_HYCOM', None)
nc.createDimension('one',     1)
    
# Create variables
# --------------------
# Coordinate Variables
# --------------------
# xi rho
xi_rho = nc.createVariable('xi_rho', 'd', ('xi_rho',), zlib=True)
xi_rho.long_name = 'X coordinate of RHO-points'
xi_rho.standard_name = 'projection_x_coordinate'
xi_rho.units = 'meter'
xi_rho_tmp = np.arange(0, Lp)
xi_rho[:] = xi_rho_tmp[:]

# eta rho
eta_rho = nc.createVariable('eta_rho', 'd', ('eta_rho',), zlib=True)
eta_rho.long_name = 'Y coordinate of RHO-points'
eta_rho.standard_name = 'projection_y_coordinate'
eta_rho.units = 'meter'
eta_rho_tmp = np.arange(0, Mp)
eta_rho[:] = eta_rho_tmp[:]

# xi u
xi_u = nc.createVariable('xi_u', 'd', ('xi_u',), zlib = True)
xi_u.long_name = 'X coordinate of U-points'
xi_u.standard_name = 'projection_x_coordinate'
xi_u.units = 'meter'
xi_u[:]=np.arange(0,L)

# eta u
eta_u = nc.createVariable('eta_u', 'd', ('eta_u',), zlib = True)
eta_u.long_name = 'Y coordinate of U-points'
eta_u.standard_name = 'projection_y_coordinate'
eta_u.units = 'meter'
eta_u[:]=np.arange(0,Mp)

# xi v
xi_v = nc.createVariable('xi_v', 'd', ('xi_v',), zlib = True)
xi_v.long_name = 'X coordinate of V-points'
xi_v.standard_name = 'projection_x_coordinate'
xi_v.units = 'meter'
xi_v[:] = np.arange(0,Lp)

# eta v
eta_v = nc.createVariable('eta_v', 'd', ('eta_v',), zlib = True)
eta_v.long_name = 'Y coordinate of V-points'
eta_v.standard_name = 'projection_y_coordinate'
eta_v.units = 'meter'
eta_v[:] = np.arange(0,M)

# xi psi
xi_psi = nc.createVariable('xi_psi', 'd', ('xi_psi',), zlib = True)
xi_psi.long_name = 'X coordinate of PSI-points'
xi_psi.standard_name = 'projection_x_coordinate'
xi_psi[:] = np.arange(0,L)

# eta psi
eta_psi = nc.createVariable('eta_psi', 'd', ('eta_psi',), zlib = True)
eta_psi.long_name = 'Y coordinate of PSI-points'
eta_psi.standard_name = 'projection_y_coordinate'
eta_psi[:] = np.arange(0,M)

# s rho
s_rho = nc.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
s_rho.long_name = 'Y coordinate of RHO-points'
s_rho.standard_name = 'projection_y_coordinate'
s_rho.units = 'meter'
s_rho_tmp = np.arange(0, N)
s_rho[:] = s_rho_tmp[:]

# time (in seconds)
time_g = nc.createVariable('time', None, ('time'), zlib=True)
time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
time_g.units = 'second'
time_g.field = 'time, scalar, series'
time_g[:] = time_tmp[:]
    
# time (HYCOM version)
time_g2 = nc.createVariable('time_HYCOM', None, ('time'), zlib=True)
time_g2.long_name = '3-hour time steps' 
time_g2.units = 'datetime'
time_g2.field = 'time_HYCOM, scalar, series'
time_g2[:] = time_tmp2[:]

# --------------------
# Geographic variables
# --------------------
# -- rho points --
# x_rho
x_rho_g = nc.createVariable('x_rho', 'd', ('eta_rho', 'xi_rho'), zlib = True)
x_rho_g.long_name = 'X coordinate of RHO-points'
x_rho_g.standard_name = 'X'
x_rho_g.units = 'meter'
x_rho_g[:,:] = model_grid.x_rho.values

# y_rho
y_rho_g = nc.createVariable('y_rho', 'd', ('eta_rho', 'xi_rho'), zlib = True)
y_rho_g.long_name = 'Y coordinate of RHO-points'
y_rho_g.standard_name = 'Y'
y_rho_g.units = 'meter'
y_rho_g[:,:]= model_grid.y_rho.values

# lon rho
lon_rho_g = nc.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'), zlib=True)
lon_rho_g.long_name = 'longitude of RHO-points'
lon_rho_g.standard_name = 'longitude'
lon_rho_g.units = 'degrees_east'
lon_rho_g[:,:] = lon_tmp[:,:]

# lat rho
lat_rho_g = nc.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'), zlib=True)
lat_rho_g.long_name = 'latitude of RHO-points'
lat_rho_g.standard_name = 'latitude'
lat_rho_g.units = 'degrees_north'
lat_rho_g[:,:] = lat_tmp[:,:]

# -- u points --
# x_u
x_u_g = nc.createVariable('x_u', 'd', ('eta_u', 'xi_u'), zlib = True)
x_u_g.long_name = 'X coordinate of U-points'
x_u_g.standard_name = 'X'
x_u_g.units = 'meter'
x_u_g[:,:]= model_grid.x_u.values
    
# y_u
y_u_g = nc.createVariable('y_u', 'd', ('eta_u', 'xi_u'), zlib = True) 
y_u_g.long_name = 'Y coordinate of U-points'
y_u_g.standard_name = 'Y'
y_u_g.units = 'meter'
y_u_g[:,:]= model_grid.y_u.values

# lon_u
lon_u_g = nc.createVariable('lon_u', 'd', ('eta_u', 'xi_u'), zlib=True)
lon_u_g.long_name = 'longitude of U-points'
lon_u_g.standard_name = 'longitude'
lon_u_g.units = 'degrees_east'
lon_u_g[:,:] = model_grid.lon_u.values

# lat_u
lat_u_g = nc.createVariable('lat_u', 'd', ('eta_u', 'xi_u'), zlib=True)
lat_u_g.long_name = 'latitude of U-points'
lat_u_g.standard_name = 'latitude'
lat_u_g.units = 'degrees_north'
lat_u_g[:,:] = model_grid.lat_u.values

# -- v points --
# x_v
x_v_g = nc.createVariable('x_v', 'd', ('eta_v', 'xi_v'), zlib = True)
x_v_g.long_name = 'X coordinate of V-points'
x_v_g.standard_name = 'X'
x_v_g.units = 'meter'
x_v_g[:,:]= model_grid.x_v.values
    
# y_v
y_v_g = nc.createVariable('y_v', 'd', ('eta_v', 'xi_v'), zlib = True)
y_v_g.long_name = 'Y coordinate of V-points'
y_v_g.standard_name = 'Y'
y_v_g.units = 'meter'
y_v_g[:,:]= model_grid.y_v.values

# lon_v
lon_v_g = nc.createVariable('lon_v', 'd', ('eta_v', 'xi_v'), zlib=True)
lon_v_g.long_name = 'longitude of V-points'
lon_v_g.standard_name = 'longitude'
lon_v_g.units = 'degrees_east'
lon_v_g[:,:] = model_grid.lon_v.values

# lat_v
lat_v_g = nc.createVariable('lat_v', 'd', ('eta_v', 'xi_v'), zlib=True)
lat_v_g.long_name = 'latitude of V-points'
lat_v_g.standard_name = 'latitude'
lat_v_g.units = 'degrees_north'
lat_v_g[:,:] = model_grid.lat_v.values

# -- psi points --
# x_psi
x_psi_g = nc.createVariable('x_psi', 'd', ('eta_psi', 'xi_psi'), zlib = True)
x_psi_g.long_name = 'X coordinate of PSI-points'
x_psi_g.standard_name = 'X'
x_psi_g.units = 'meter'
x_psi_g[:,:]= model_grid.x_psi.values
    
# y_psi
y_psi_g = nc.createVariable('y_psi', 'd', ('eta_psi', 'xi_psi'), zlib = True)
y_psi_g.long_name = 'Y coordinate of PSI-points'
y_psi_g.standard_name = 'Y'
y_psi_g.units = 'meter'
y_psi_g[:,:]= model_grid.y_psi.values

# lon_psi
lon_psi_g = nc.createVariable('lon_psi', 'd', ('eta_psi', 'xi_psi'), zlib=True)
lon_psi_g.long_name = 'longitude of PSI-points'
lon_psi_g.standard_name = 'longitude'
lon_psi_g.units = 'degrees_east'
lon_psi_g[:,:] = model_grid.lon_psi.values

# lat_psi
lat_psi_g = nc.createVariable('lat_psi', 'd', ('eta_psi', 'xi_psi'), zlib=True)
lat_psi_g.long_name = 'latitude of PSI-points'
lat_psi_g.standard_name = 'latitude'
lat_psi_g.units = 'degrees_north'
lat_psi_g[:,:] = model_grid.lat_psi.values

# ----------
# Bathymetry
# ----------
# bathymetry
h = nc.createVariable('h', 'd', ('eta_rho', 'xi_rho'), zlib=True)
h.long_name = "Final bathymetry at RHO-points"
h.standard_name = "sea_floor_depth"
h.units = "meter"
h.field = "bath, scalar"
h.coordinates = "lon_rho lat_rho"
h.grid_mapping = gridmap_varname 
h[:,:]= model_grid.h.values

# ------------
# Metric terms
# ------------
pm_g = nc.createVariable('pm', 'd', ('eta_rho', 'xi_rho'), zlib = True)
pm_g.long_name = 'curvilinear coordinate metric in XI'
pm_g.units = 'meter-1'
pm_g.field = 'pm, scalar'
pm_g.coordinates = 'lon_rho lat_rho'
pm_g.grid_mapping = gridmap_varname
pm_g[:,:]= model_grid.pm.values
    
pn_g = nc.createVariable('pn', 'd', ('eta_rho', 'xi_rho'), zlib = True)
pn_g.long_name = 'curvilinear coordinate metric in ETA'
pn_g.units = 'meter-1'
pn_g.field = 'pn, scalar'
pn_g.coordinates = 'lon_rho lat_rho'
pn_g.grid_mapping = gridmap_varname
pn_g[:,:]= model_grid.pn.values

dndx_g = nc.createVariable('dndx', 'd', ('eta_rho', 'xi_rho'), zlib = True)
dndx_g.long_name = 'xi derivative of inverse metric factor pn'
dndx_g.units = 'meter'
dndx_g.field = 'dndx, scalar'
dndx_g.coordinates = 'lon_rho lat_rho'
dndx_g.grid_mapping = gridmap_varname
dndx_g[:,:]= model_grid.dndx.values
    
dmde_g = nc.createVariable('dmde', 'd', ('eta_rho', 'xi_rho'), zlib = True)
dmde_g.long_name = 'eta derivative of inverse metric factor pm'
dmde_g.units = 'meter'
dmde_g.field = 'dmde, scalar'
dmde_g.coordinates = 'lon_rho lat_rho'
dmde_g.grid_mapping = gridmap_varname
dmde_g[:,:]= model_grid.dmde.values

# ---------
# Spherical
# ---------
spherical = nc.createVariable('spherical', 'c', ())
spherical.long_name = 'grid type logical switch'
spherical.flag_values = 'T, F'
spherical.flag_meanings = 'spherical Cartesian'
#spherical[:] = 0 # OG
spherical[:] = 'T' # NEW

# --------
# Coriolis
# --------
f_g = nc.createVariable('f', 'd', ('eta_rho', 'xi_rho'), zlib = True)
f_g.long_name = 'Coriolis parameter at RHO-points'
f_g.standard_name = 'coriolis_parameter'
f_g.units = 'second-1'
f_g.field = 'Coriolos, scalar'
f_g.coordinates = 'lon_rho lat_rho'
f_g.grid_mapping = gridmap_varname
f_g[:,:] = model_grid.f.values

# -----
# Angle
# -----
angle_g = nc.createVariable('angle', 'd', ('eta_rho', 'xi_rho'), zlib = True)
angle_g.long_name = 'angle between xi axis and east'
angle_g.standard_name = 'angle_of_rotation_from_east_to_x'
angle_g.units = 'radian'
angle_g.coordinates = 'lon_rho lat_rho'
angle_g.grid_mapping = gridmap_varname
angle_g[:,:] = model_grid.angle.values

# -----
# Masks
# -----
mask_rho_g = nc.createVariable('mask_rho', 'd', ('eta_rho', 'xi_rho'), zlib = True)
mask_rho_g.long_name = 'mask on RHO-points'
mask_rho_g.option_0 = 'land'
mask_rho_g.option_1 = 'water'
mask_rho_g.coordinates = 'lon_rho lat_rho'
mask_rho_g.grid_mapping = gridmap_varname
mask_rho_g[:,:]= model_grid.mask_rho.values

mask_psi_g = nc.createVariable('mask_psi', 'd', ('eta_psi', 'xi_psi'), zlib = True)
mask_psi_g.long_name = 'mask on PSI-points'
mask_psi_g.option_0 = 'land'
mask_psi_g.option_1 = 'water'
mask_psi_g.coordinates = 'lon_psi lat_psi'
mask_psi_g.grid_mapping = gridmap_varname
mask_psi_g[:,:]= model_grid.mask_psi.values 

mask_u_g = nc.createVariable('mask_u', 'd', ('eta_u', 'xi_u'), zlib = True)
mask_u_g.long_name = 'mask on U-points'
mask_u_g.option_0 = 'land'
mask_u_g.option_1 = 'water'
mask_u_g.coordinates = 'lon_u lat_u'
mask_u_g.grid_mapping = gridmap_varname
mask_u_g[:,:]= model_grid.mask_u.values

mask_v_g = nc.createVariable('mask_v', 'd', ('eta_v', 'xi_v'), zlib = True)
mask_v_g.long_name = 'mask on V-points'
mask_v_g.option_0 = 'land'
mask_v_g.option_1 = 'water'
mask_v_g.coordinates = 'lon_v lat_v'
mask_v_g.grid_mapping = gridmap_varname
mask_v_g[:,:]= model_grid.mask_v.values

# --------
# Grid map
# --------
gridmap_varname_g = nc.createVariable(gridmap_varname, 'i', ())
gridmap_varname_g.long_name = 'grid mapping'

#nc2 = nc.variables[gridmap_varname]
#gmap = gridmap.gridmap.projection.PolarStereographic(0, 0, 500, (90-(360-((model_grid.angle[0,0].values/(2*np.pi))*360.0)))) # (x grid coordinate of north pole, y grid coordinate of north pole, grid resolution (m), angle of y-axis (deg))
gmap = PolarStereographic(0, 0, 500, (90-(360-((model_grid.angle[0,0].values/(2*np.pi))*360.0)))) # (x grid coordinate of north pole, y grid coordinate of north pole, grid resolution (m), angle of y-axis (deg))


d2 = gmap.CFmapping_dict()
for att in d2:
    #setattr(gridmap_varname_g, att, d2[att]) # makes lat lon not work in Panoply | don't use
    setattr(nc, att, d2[att]) # works in Panoply but puts all this info in global attributes | use
gridmap_varname_g.proj4string = gmap.proj4string


# --------
# Grid depth
# --------
# grid depth of rho points, z_rho
z_rho_g = nc.createVariable('z_rho', 'f8', ('time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True)
z_rho_g.long_name = 'ROMS dpeth of rho points'
z_rho_g.units = 'meter'
z_rho_g.coordinates = 'lon lat'
z_rho_g.time = 'time'

# total grid depth 
total_depth_g = nc.createVariable('total_depth', 'f8', ('time', 'eta_rho', 'xi_rho'), zlib=True)
total_depth_g.long_name = 'bathymetry plus surface elevation'
total_depth_g.units = 'm'
total_depth_g.coordinates = 'lon lat'
total_depth_g.time = 'time'
total_depth_g[:,:,:] = h_plus_zeta[:,:,:]

# grid depth of w points, z_w
z_w_g = nc.createVariable('z_w', 'f8', ('time', 's_w', 'eta_rho', 'xi_rho'), zlib=True)
z_w_g.long_name = 'ROMS dpeth of w points'
z_w_g.units = 'meter'
z_w_g.coordinates = 'lon lat'
z_w_g.time = 'time'



# ----- call the function -----
# make an array to save the values to
# This will be useful when we have a time varying component 
# whole grid
#zrho_all = np.empty((time_len, eta_len, xi_len, N))
zrho_all = np.empty((time_len, N, eta_len, xi_len))
zw_all = np.empty((time_len, N+1, eta_len, xi_len))

# pull out the grid's bathymetry
h_tmp = model_grid.h[:,:].values

# Do a fancy ChatGPT version that parallelizes this 
# NEW
# ----- Find depths over time -----
from joblib import Parallel, delayed

def compute_depth_for_1timestep(tt):
    print(f'Processing time step {tt}')
    
    # read in zeta for this time
    zeta_tmp = zeta_interp.zeta[tt,:,:].values

    # flatten the bathymetry and zeta arrays
    h_tmp_flat = h_tmp.flatten()
    zeta_tmp_flat = zeta_tmp.flatten()

    # find C = something with Vstretching
    # If OG 20 vertical layer grid
    #C = get_Vstretching_3(theta_s=theta_s, theta_b=theta_b)
    # If new 30 vertical layer grid 
    C = get_Vstretching_4(theta_s=theta_s, theta_b=theta_b)

    # Calculate depths
    depths = get_depths(Vtransform, C, h_tmp_flat, hc)

    # Get s_rho
    srho = get_srho(N)

    # Find depths of rho points
    zr = depths(srho).reshape((N, eta_len, xi_len))

    # Find zrho
    h_tmp2 = np.copy(h_tmp_flat)
    zrho2 = get_zrho(Vtransform, Vstretching, N, theta_s, theta_b, h_tmp2, hc, zeta_tmp_flat)
    zrho2 = zrho2.reshape((N, eta_len, xi_len))
    zrho3 = np.copy(zrho2)

    # Find z_w
    zw1 = get_zw(Vtransform, Vstretching, N+1, theta_s, theta_b, h_tmp2, hc, zeta_tmp_flat)
    zw1 = zw1.reshape((N+1, eta_len, xi_len))

    return tt, zrho3, zw1

results = Parallel(n_jobs=120)(delayed(compute_depth_for_1timestep)(tt) for tt in range(time_len))

for tt, zrho3, zw1 in results:
    zrho_all[tt,:,:,:] = zrho3
    z_rho_g[tt,:,:,:] = zrho3
    zw_all[tt,:,:,:] = zw1
    z_w_g[tt,:,:,:] = zw1

# Now sync netcdf once
nc.sync()


# OLD
# # ----- Find depths over time -----
# # start looping through time
# for tt in range(time_len):
# #for tt in range(3):
#     print(tt)
#     # read in the bathymetry at this time (h_plus_zeta)
#     #h_tmp = h_plus_zeta[tt,:,:] #OG when we used hpluszeta as h
    
#     # read in zeta for this time
#     zeta_tmp = zeta_interp.zeta[tt,:,:].values
    
#     # flatten the bathymetry and zeta arrays
#     h_tmp = h_tmp.flatten()
#     zeta_tmp = zeta_tmp.flatten()
    
#     # checkpoint
#     #print('1')
    
#     # find C = something with Vstretching
#     C = get_Vstretching_3(theta_s=theta_s, theta_b=theta_b)
    
#     # checkpoint
#     #print('2')
    
#     # calculate the depths
#     depths = get_depths(Vtransform, C, h_tmp, hc)
    
#     # checkpoint
#     #print('3')
    
#     # find s_rho (evenly spaced)
#     srho = get_srho(N)
    
#     # checkpoint 
#     #print('4')
    
#     # find the depths of rho points
#     zr = depths(srho)
    
#     # checkpoint
#     #print('5')
    
#     # unflatten zeta
#     zr = zr.reshape((N, eta_len, xi_len))
#     #print(zr.shape)
    
#     # find zrho
#     h_tmp2 = np.copy(h_tmp)
#     zrho2 = get_zrho(Vtransform, Vstretching, N, theta_s, theta_b, h_tmp2, hc, zeta_tmp) 
#     #print('zrho2: ', zrho2.shape)
#     zrho2 = zrho2.reshape((N, eta_len, xi_len))
#     #print(zrho2.shape)
#     # move some axes 
#     zrho3 = np.copy(zrho2)
#     #zrho3 = np.moveaxis(zrho3, 0, -1)
#     #print(zrho3.shape)
    
#     # find z_w 
#     zw1 = get_zw(Vtransform, Vstretching, N+1, theta_s, theta_b, h_tmp2, hc, zeta_tmp)
#     #print('zw1: ', zw1.shape)
#     zw1 = zw1.reshape((N+1, len(model_grid.eta_rho), len(model_grid.xi_rho)))
#     #print(zw1.shape)
    
#     # Now save this zrho into the array or just to the netcdf
#     # z_rho
#     zrho_all[tt,:,:,:] = zrho3
#     z_rho_g[tt,:,:,:] = zrho3
#     # z_w
#     zw_all[tt,:,:,:] = zw1
#     z_w_g[tt,:,:,:] = zw1
    
#     # force sync to netcdf file in case the script crashes 
#     nc.sync()
    
#     # slight pause to see if this worked
#     #input('Press enter to continue...')

    
    
    
# close the netcdf
nc.close()
print('closed netcdf')

### ----- End Section 1 ----- ###




### ----- Section 2: Uncomment this section if you want to make plots of the settings ----- ###
# =============================================================================
# # ----------- OG we know the code below works  so let's try making the time-looped version above ------------
# # read in the bathymetry 
# #h = np.copy(h_plus_zeta)
# h = np.asarray(model_grid.h)
# print(h.shape)
# # flatten the bathymetry
# h = h.flatten()
# print(h.shape)
# 
# # read in zeta
# #zeta = zeta_interp.surf_el
# 
# # flatten zeta
# 
# 
# 
# # checkpoint
# print('1')
# 
# # find C = something with Vstretching
# C = get_Vstretching_3(theta_s=theta_s, theta_b=theta_b)
# 
# # checkpoint
# print('2')
# 
# # calculate the depths 
# depths = get_depths(Vtransform, C, h, hc)
# 
# # checkpoint
# print('3')
# 
# # find s_rho (evenly spaced)
# srho = get_srho(N)
# 
# # checkpoint
# print('4')
# 
# # find the depths of rho points 
# zr = depths(srho)
# 
# # checkpoint
# print('5')
# 
# # unflatten zeta
# zr = zr.reshape((N, len(model_grid.eta_rho), len(model_grid.xi_rho)))
# print(zr.shape)
# 
# # play with get_zrho
# # I think this is the first time we need zeta
# h12 = np.copy(h)
# #zeta12 = np.asarray(model_grid.zeta[0,:,:])
# #zeta12 = zeta12.flatten()
# 
# # checkpoint
# print('6')
# 
# # find z_rho
# zrho12 = get_zrho(Vtransform, Vstretching, N, theta_s, theta_b, h12, hc) # add zeta here when we have it
# print('zrh2o: ', zrho12.shape)
# zrho12 = zrho12.reshape((N, len(model_grid.eta_rho), len(model_grid.xi_rho)))
# print(zrho12.shape)
# 
# # find z_w
# zw1 = get_zw(Vtransform, Vstretching, N+1, theta_s, theta_b, h12, hc)
# print('zw1: ', zw1.shape)
# zw1 = zw1.reshape((N+1, len(model_grid.eta_rho), len(model_grid.xi_rho)))
# print(zw1.shape)
# 
# # plot the results for different depths/sections of the grid 
# fig1, ax1 = plt.subplots()
# ax1.plot(zrho12[:,2,300], color='r', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Max Bathymetry: ~ 5 meters')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# fig2, ax2 = plt.subplots()
# ax2.plot(zrho12[:,64,300], color='orange', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Max Bathymetry: ~ 15 meters')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# fig3, ax3 = plt.subplots()
# ax3.plot(zrho12[:,75,300], color='gold', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Max Bathymetry: ~ 25 meters')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# fig4, ax4 = plt.subplots()
# ax4.plot(zrho12[:,101,300], color='green', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Max Bathymetry: ~ 40 meters')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# fig4, ax4 = plt.subplots()
# ax4.plot(zrho12[:,160,300], color='blue', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Max Bathymetry: ~ 100 meters')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# fig5, ax5 = plt.subplots()
# ax5.plot(zrho12[:,177,300], color='purple', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Max Bathymetry: ~ 1000 meters')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# fig6, ax6 = plt.subplots()
# ax6.plot(zrho12[:,205,300], color='magenta', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Max Bathymetry: ~ 2000 meters')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# # I feel like these are weird dimensions, is there a way to flip all of this at the end?
# # want (time, eta, xi, N) or (eta, xi, N) = (206, 608, 20)
# zrho12_switch = np.copy(zrho12)
# zrho12_switch = np.moveaxis(zrho12_switch, 0, -1)
# print(zrho12_switch.shape)
# 
# # Plot to see if  it still worked the same
# fig7, ax7 = plt.subplots()
# ax7.plot(zrho12_switch[2,300,:], color='r', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Max Bathymetry: ~ 5 meters, swapped')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# # From the plot, we can see that this works!
# 
# # Plot z_w and z_rho together to make sure it looks right
# fig8, ax8 = plt.subplots()
# ax8.plot(zrho12[:,205,300], color='darkmagenta', marker='x', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# ax8.plot(zw1[:,205,300], color='limegreen', marker='+', label=('Vtransform=' + str(Vtransform) + ', Vstretching=' + str(Vstretching) +', theta_s=' + str(theta_s) +', theta_b=' + str(theta_b)))
# plt.title('Comparison of z_w and z_rho')
# plt.xlabel('Vertical Layer')
# plt.ylabel('Depth (m)')
# plt.legend(bbox_to_anchor=(1.04,1), loc='upper left')
# 
# # -------------------------------------------------------------------------------------------------------------------
# =============================================================================
### -----  End Section 2 ----- ###

