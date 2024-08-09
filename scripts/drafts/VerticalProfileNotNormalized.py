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

# Initialize arrays dynamically based on variable choice
timeNZ_product_NH = []
timeNZ_product_SH = []
timeFW_product_NH = []
timeFW_product_SH = []
timeIV_product_NH = []
timeIV_product_SH = []
timeNZS_product_NH = []
timeNZS_product_SH = []

# Initialize arrays for different factors
timeBC1_product_NH = []
timeBC1_product_SH = []
timeBC2_product_NH = []
timeBC2_product_SH = []
timeDust_product_NH = []
timeDust_product_SH = []

for year in years:
    for month in months:
        # Format the filename for the current month
        datafile_ice = datafile_ice_fmt.format(year=year, month=month)
        filepath_ice = os.path.join(runDir, datafile_ice)
        datafile_ocean = datafile_ocean_fmt.format(year=year, month=month)
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
                
                # Sum for different factors
                BC1_product_NH = np.sum(varNZ[:, :, 0:7].isel(nCells=NHindices), axis=(1, 2))
                timeBC1_product_NH.append([year + month / 12, np.sum(BC1_product_NH)])
                
                BC1_product_SH = np.sum(varNZ[:, :, 0:7].isel(nCells=SHindices), axis=(1, 2))
                timeBC1_product_SH.append([year + month / 12, np.sum(BC1_product_SH)])
                
                BC2_product_NH = np.sum(varNZ[:, :, 8:16].isel(nCells=NHindices), axis=(1, 2))
                timeBC2_product_NH.append([year + month / 12, np.sum(BC2_product_NH)])
                
                BC2_product_SH = np.sum(varNZ[:, :, 8:16].isel(nCells=SHindices), axis=(1, 2))
                timeBC2_product_SH.append([year + month / 12, np.sum(BC2_product_SH)])
                
                Dust_product_NH = np.sum(varNZ[:, :, 17:24].isel(nCells=NHindices), axis=(1, 2))
                timeDust_product_NH.append([year + month / 12, np.sum(Dust_product_NH)])
                
                Dust_product_SH = np.sum(varNZ[:, :, 17:24].isel(nCells=SHindices), axis=(1, 2))
                timeDust_product_SH.append([year + month / 12, np.sum(Dust_product_SH)])
                
                # Close datasets to free up resources
                output_ice.close()
                mesh.close()
                
            else:
                print(f'Output file not found: {output_filepath_ice}')
        else:
            print(f'File not found: {filepath_ice}')

# Convert lists to numpy arrays
timeBC1_product_NH = np.array(timeBC1_product_NH)
timeBC1_product_SH = np.array(timeBC1_product_SH)
timeBC2_product_NH = np.array(timeBC2_product_NH)
timeBC2_product_SH = np.array(timeBC2_product_SH)
timeDust_product_NH = np.array(timeDust_product_NH)
timeDust_product_SH = np.array(timeDust_product_SH)

# Define years and months based on your data
start_year = min(years)
end_year = max(years)

fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True, sharey=True)

# Create ticks for x-axis
start_month = 1
end_month = 12
ticks = [(y + m / 12) for y in range(start_year, end_year + 1) for m in range(start_month, end_month + 1)]
tick_labels = [f'{year}-{month:02d}' for year in range(start_year, end_year + 1) for month in range(start_month, end_month + 1)]

# Plot NH data
axs[0].plot(timeBC1_product_NH[:, 0], timeBC1_product_NH[:, 1], label='Black Carbon 1 NH')
axs[0].plot(timeBC2_product_NH[:, 0], timeBC2_product_NH[:, 1], label='Black Carbon 2 NH')
axs[0].plot(timeDust_product_NH[:, 0], timeDust_product_NH[:, 1], label='Dust NH')
axs[0].set_title('Northern Hemisphere Vertical Profile')
axs[0].legend()

# Plot SH data
axs[1].plot(timeBC1_product_SH[:, 0], timeBC1_product_SH[:, 1], label='Black Carbon 1 SH')
axs[1].plot(timeBC2_product_SH[:, 0], timeBC2_product_SH[:, 1], label='Black Carbon 2 SH')
axs[1].plot(timeDust_product_SH[:, 0], timeDust_product_SH[:, 1], label='Dust SH')
axs[1].set_title('Southern Hemisphere Vertical Profile')
axs[1].legend()

# Add common labels
plt.xlabel('Time')
fig.text(0.04, 0.5, 'Values', va='center', rotation='vertical')

plt.tight_layout()
plt.savefig('VertNotNormalNH_SH_Products_VerticalProfile.png')
plt.show()

