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

# Define variable choices
var1_choice = 'NZ'  # Change to select the first variable (options: 'NZ', 'FW', 'VC', 'AC')
var2_choice = 'FW'  # Change to select the second variable

# Initialize arrays based on variable choices
time_var1_product_NH = []
time_var1_product_SH = []
time_var2_product_NH = []
time_var2_product_SH = []

deg2rad = 3.14159 / 180.0
months = range(1, 13)  # Include all 12 months
years = range(20, 22)  # Adjusted to include all

for year in years:
    for month in months:
        # Format the filename for the current month
        datafile_ice = datafile_ice_fmt.format(year=year, month=month)
        filepath_ice = os.path.join(runDir, datafile_ice)
        
        datafile_ocean = datafile_ocean_fmt.format(year=year, month=month)
        filepath_ocean = os.path.join(runDir, datafile_ocean)

        # Check if the file exists 
        if os.path.exists(filepath_ice) and os.path.exists(filepath_ocean):
            print(f'Processing file: {filepath_ice}')
            print(f'Processing file: {filepath_ocean}')

            # Open the dataset
            output_ice = xr.open_dataset(filepath_ice)
            output_ocean = xr.open_dataset(filepath_ocean)

            # Choose variables based on the selected choices
            if var1_choice == 'NZ':
                var1 = output_ice.variables['timeMonthly_avg_verticalAerosolsIceCell']
            elif var1_choice == 'FW':
                var1 = output_ice.variables['timeMonthly_avg_oceanFreshWaterFlux']
            elif var1_choice == 'VC':
                var1 = output_ice.variables['timeMonthly_avg_iceVolumeCell']
            elif var1_choice == 'AC':
                var1 = output_ice.variables['timeMonthly_avg_iceAreaCell']

            if var2_choice == 'NZ':
                var2 = output_ice.variables['timeMonthly_avg_verticalAerosolsIceCell']
            elif var2_choice == 'FW':
                var2 = output_ice.variables['timeMonthly_avg_oceanFreshWaterFlux']
            elif var2_choice == 'VC':
                var2 = output_ice.variables['timeMonthly_avg_iceVolumeCell']
            elif var2_choice == 'AC':
                var2 = output_ice.variables['timeMonthly_avg_iceAreaCell']

            # Get mesh data 
            mesh = xr.open_dataset(meshFileName)
            latCell = mesh.variables['latCell']
            nCells = mesh.sizes['nCells']
            areaCell = mesh.variables['areaCell']

            # Get NH indices
            NHindices = np.where(latCell > 60 * deg2rad)[0]
            SHindices = np.where(latCell < -60 * deg2rad)[0]

            # Sum for var1 (NZAerosols) and var2 (Ocean Freshwater Flux)
            var1_product_NH = np.sum(var1.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(0, 1))
            time_var1_product_NH.append([year + month / 12, np.sum(var1_product_NH)])

            var1_product_SH = np.sum(var1.isel(nCells=SHindices) * areaCell.isel(nCells=SHindices), axis=(0, 1))
            time_var1_product_SH.append([year + month / 12, np.sum(var1_product_SH)])

            var2_product_NH = np.sum(var2.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(0, 1))
            time_var2_product_NH.append([year + month / 12, np.sum(var2_product_NH)])

            var2_product_SH = np.sum(var2.isel(nCells=SHindices) * areaCell.isel(nCells=SHindices), axis=(0, 1))
            time_var2_product_SH.append([year + month / 12, np.sum(var2_product_SH)])

            # Close datasets to free up resources
            output_ice.close()
            mesh.close()
        else:
            print(f'File not found: {filepath_ice}')

# Convert lists to numpy arrays
time_var1_product_NH = np.array(time_var1_product_NH)
time_var1_product_SH = np.array(time_var1_product_SH)
time_var2_product_NH = np.array(time_var2_product_NH)
time_var2_product_SH = np.array(time_var2_product_SH)

# Define years and months based on your data
start_year = min(years)
end_year = max(years)

# Create ticks for x-axis
ticks = [(y + m / 12) for y in range(start_year, end_year + 1) for m in range(1, 13)]
tick_labels = [f'{year}-{month:02d}' for year in range(start_year, end_year + 1) for month in range(1, 13)]

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 18))  # Two rows, one column

# Plot for var1 (NZAerosols) on the primary y-axis
line1 = ax1.plot(time_var1_product_NH[:, 0], time_var1_product_NH[:, 1], linestyle='-', color='b', label=f'NH - {var1_choice}')
ax1.set_ylabel(f'Sum of {var1_choice}', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create a twin axis for var2 (Freshwater Flux)
ax1_twin = ax1.twinx()
line2 = ax1_twin.plot(time_var2_product_NH[:, 0], time_var2_product_NH[:, 1], linestyle='-', color='g', label=f'NH - {var2_choice}')
ax1_twin.set_ylabel(f'Sum of {var2_choice}', color='g')
ax1_twin.tick_params(axis='y', labelcolor='g')

# Set x-axis limits and ticks
ax1.set_xticks(ticks)
ax1.set_xticklabels(tick_labels, rotation=45)
ax1.set_xlim(start_year + 1 / 24, end_year + 12 / 24)
ax1.set_xlabel('Time (Year-Month)')

# Title and grid
ax1.set_title(f'Sum of {var1_choice} and {var2_choice} Over Time (NH)')
ax1.grid(True)

# Show legend
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper right')

# Plot for var1 (FWAerosols) on the primary y-axis
line3 = ax2.plot(time_var1_product_SH[:, 0], time_var1_product_SH[:, 1], linestyle='-', color='b', label=f'SH - {var1_choice}')
ax2.set_ylabel(f'Sum of {var1_choice}', color='b')
ax2.tick_params(axis='y', labelcolor='b')

# Create a twin axis for var2 (Freshwater Flux)
ax2_twin = ax2.twinx()
line4 = ax2_twin.plot(time_var2_product_SH[:, 0], time_var2_product_SH[:, 1], linestyle='-', color='g', label=f'SH - {var2_choice}')
ax2_twin.set_ylabel(f'Sum of {var2_choice}', color='g')
ax2_twin.tick_params(axis='y', labelcolor='g')

# Set x-axis limits and ticks
ax2.set_xticks(ticks)
ax2.set_xticklabels(tick_labels, rotation=45)
ax2.set_xlim(start_year + 1 / 24, end_year + 12 / 24)
ax2.set_xlabel('Time (Year-Month)')

# Title and grid
ax2.set_title(f'Sum of {var1_choice} and {var2_choice} Over Time (SH)')
ax2.grid(True)

# Show legend
lines = line3 + line4
labels = [l.get_label() for l in lines]
ax2.legend(lines, labels, loc='upper right')

# Save the figure
plt.savefig(f'SeparateSeasonalCycle{start_year}-{end_year}sum{var1_choice}and{var2_choice}_time_series_SH.png')
plt.show()
plt.close(fig)

