#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals

"""
sum plots
Savannah Byron
This code helps me select a variable and plot the average value in time
"""

############################## model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
datafile_ice_fmt = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.00{year:02d}-{month:02d}-01.nc'
datafile_ocean_fmt = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpaso.hist.am.timeSeriesStatsMonthly.00{year:02d}-{month:02d}-01.nc'

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os

#define variables
deg2rad = 3.14159 / 180.0
nCells = 465044
months = range(1, 13)  # Adjusted to include all 12 months
years = range(20, 22)  # Adjusted to include all

# Choose a specific time point (e.g., year 20 and month 1)
specific_year = 20
specific_month = 1

# Initialize arrays for different factors
BC1_profile_NH = []
BC1_profile_SH = []
BC2_profile_NH = []
BC2_profile_SH = []
Dust_profile_NH = []
Dust_profile_SH = []

# Format the filename for the specific month
datafile_ice = datafile_ice_fmt.format(year=specific_year, month=specific_month)
filepath_ice = os.path.join(runDir, datafile_ice)
datafile_ocean = datafile_ocean_fmt.format(year=specific_year, month=specific_month)
filepath_ocean = os.path.join(runDir, datafile_ocean)

# Check if the file exists
if os.path.exists(filepath_ice) & os.path.exists(filepath_ocean):
    print(f'Processing file: {filepath_ice}')
    print(f'Processing file: {filepath_ocean}')
    
    # Read the corresponding output file
    output_filepath_ice = os.path.join(runDir, datafile_ice)
    output_filepath_ocean = os.path.join(runDir, datafile_ocean)
    
    if os.path.exists(output_filepath_ice) & os.path.exists(output_filepath_ocean):
        # Open the dataset
        output_ice = xr.open_dataset(output_filepath_ice)
        output_ocean = xr.open_dataset(output_filepath_ocean)
        
        # Get variables from the current output dataset
        varNZ = output_ice.variables['timeMonthly_avg_verticalAerosolsIceCell']
        varIV = output_ice.variables['timeMonthly_avg_iceVolumeCell']
        varAC = output_ice.variables['timeMonthly_avg_iceAreaCell']
        varFW = output_ice.variables['timeMonthly_avg_oceanFreshWaterFlux']
        
        # Get mesh data
        mesh = xr.open_dataset(meshFileName)
        latCell = mesh.variables['latCell']
        nCells = mesh.sizes['nCells']
        latCell = mesh.variables['latCell']  # in radians
        lonCell = mesh.variables['lonCell']
        xCell = mesh.variables['xCell']  # in meters
        yCell = mesh.variables['yCell']
        areaCell = mesh.variables['areaCell']
        
        # Get NH indices
        NHindices = np.where(latCell > 50 * deg2rad)[0]
        SHindices = np.where(latCell < -50 * deg2rad)[0]
        
        # Extract the vertical profiles for different factors
        BC1_profile_NH = np.mean(varNZ.isel(nCells=NHindices).isel(nzAerosolsIceLayers=slice(0, 7)), axis=(0, 1))
        BC1_profile_SH = np.mean(varNZ.isel(nCells=SHindices).isel(nzAerosolsIceLayers=slice(0, 7)), axis=(0, 1))
        
        BC2_profile_NH = np.mean(varNZ.isel(nCells=NHindices).isel(nzAerosolsIceLayers=slice(8, 16)), axis=(0, 1))
        BC2_profile_SH = np.mean(varNZ.isel(nCells=SHindices).isel(nzAerosolsIceLayers=slice(8, 16)), axis=(0, 1))
        
        Dust_profile_NH = np.mean(varNZ.isel(nCells=NHindices).isel(nzAerosolsIceLayers=slice(17, 24)), axis=(0, 1))
        Dust_profile_SH = np.mean(varNZ.isel(nCells=SHindices).isel(nzAerosolsIceLayers=slice(17, 24)), axis=(0, 1))
        
        # Close datasets to free up resources
        output_ice.close()
        mesh.close()
        
    else:
        print(f'Output file not found: {output_filepath_ice}')
else:
    print(f'File not found: {filepath_ice}')

# Plot the vertical profiles
fig, axs = plt.subplots(1, 2, figsize=(15, 8), sharey=True)

# Plot NH data
axs[0].plot(BC1_profile_NH, range(0, 7), label='Black Carbon 1 NH', color='blue')
axs[0].plot(BC2_profile_NH, range(8, 16), label='Black Carbon 2 NH', color='green')
axs[0].plot(Dust_profile_NH, range(17, 24), label='Dust NH', color='red')
axs[0].set_title('Northern Hemisphere Vertical Profile')
axs[0].set_xlabel('Values')
axs[0].set_ylabel('Vertical Levels')
axs[0].invert_yaxis()
axs[0].legend()

# Plot SH data
axs[1].plot(BC1_profile_SH, range(0, 7), label='Black Carbon 1 SH', color='blue')
axs[1].plot(BC2_profile_SH, range(8, 16), label='Black Carbon 2 SH', color='green')
axs[1].plot(Dust_profile_SH, range(17, 24), label='Dust SH', color='red')
axs[1].set_title('Southern Hemisphere Vertical Profile')
axs[1].set_xlabel('Values')
axs[1].invert_yaxis()
axs[1].legend()

plt.tight_layout()
plt.savefig('NH_SH_Products_VerticalProfile_SpecificTime.png')
plt.show()

