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

# Define constants and variables
deg2rad = 3.14159 / 180.0
months = range(1, 13)  # Adjusted to include all 12 months
years = range(20, 22)  # Adjusted to include all

# Initialize arrays to store summed values over time
sum_timefirstvar_product_NH = []
sum_timefirstvar_product_SH = []
sum_timesecondvar_product_NH = []
sum_timesecondvar_product_SH = []

for year in years:
    for month in months:
        # Format the filename for the current month
        datafile_ice = datafile_ice_fmt.format(year=year, month=month)
        filepath_ice = os.path.join(runDir, datafile_ice)
        datafile_ocean = datafile_ocean_fmt.format(year=year, month=month)
        filepath_ocean = os.path.join(runDir, datafile_ocean)

        # Check if the file exists
        if os.path.exists(filepath_ice) and os.path.exists(filepath_ocean):
            print(f'Processing files: {filepath_ice} and {filepath_ocean}')

            # Read the datasets
            output_ice = xr.open_dataset(filepath_ice)
            output_ocean = xr.open_dataset(filepath_ocean)

            # Get ice volume and area data
            var2 = output_ice['timeMonthly_avg_iceVolumeCell']
            var3 = output_ice['timeMonthly_avg_iceAreaCell']

            # Get mesh data
            mesh = xr.open_dataset(meshFileName)
            latCell = mesh['latCell']
            areaCell = mesh['areaCell']

            # Get NH and SH indices
            NHindices = np.where(latCell > 50 * deg2rad)[0]
            SHindices = np.where(latCell < -50 * deg2rad)[0]

            # Calculate sums for NH and SH regions for iceVolumeCell
            sum_NH_iceVolume = np.sum(var2.isel(nCells=NHindices), axis=1)
            sum_SH_iceVolume = np.sum(var2.isel(nCells=SHindices), axis=1)
            
            # Calculate sums for NH and SH regions for iceAreaCell
            sum_NH_iceArea = np.sum(var3.isel(nCells=NHindices), axis=1)
            sum_SH_iceArea = np.sum(var3.isel(nCells=SHindices), axis=1)

            sum_timefirstvar_product_NH.append(sum_NH_iceVolume)
            sum_timefirstvar_product_SH.append(sum_SH_iceVolume)
            sum_timesecondvar_product_NH.append(sum_NH_iceArea)
            sum_timesecondvar_product_SH.append(sum_SH_iceArea)

            # Close datasets
            output_ice.close()
            output_ocean.close()
            mesh.close()

        else:
            print(f'Files not found: {filepath_ice} or {filepath_ocean}')

# Debug print statements to check the structure
print(f'sum_timefirstvar_product_NH: {sum_timefirstvar_product_NH}')
print(f'sum_timefirstvar_product_SH: {sum_timefirstvar_product_SH}')
print(f'sum_timesecondvar_product_NH: {sum_timesecondvar_product_NH}')
print(f'sum_timesecondvar_product_SH: {sum_timesecondvar_product_SH}')

# Convert lists to numpy arrays
sum_timefirstvar_product_NH = np.array(sum_timefirstvar_product_NH)
sum_timefirstvar_product_SH = np.array(sum_timefirstvar_product_SH)
sum_timesecondvar_product_NH = np.array(sum_timesecondvar_product_NH)
sum_timesecondvar_product_SH = np.array(sum_timesecondvar_product_SH)


# Print shapes for further debugging
print(f'sum_timefirstvar_product_NH shape: {sum_timefirstvar_product_NH.shape}')
print(f'sum_timefirstvar_product_SH shape: {sum_timefirstvar_product_SH.shape}')
print(f'sum_timesecondvar_product_NH shape: {sum_timesecondvar_product_NH.shape}')
print(f'sum_timesecondvar_product_SH shape: {sum_timesecondvar_product_SH.shape}')

# Define variables names
varNam2 = 'iceVolumeCell in kg'
varNam3 = 'iceAreaCell'

# Define years and months based on your data
start_year = min(years)
end_year = max(years)

# Create ticks for x-axis
start_month = 1
end_month = 12
ticks = [(y + m / 12) for y in range(start_year, end_year + 1) for m in range(start_month, end_month + 1)]
tick_labels = [f'{year}-{month:02d}' for year in range(start_year, end_year + 1) for month in range(start_month, end_month + 1)]

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 18))  # Two rows, one column

# Plot iceVolumeCell on the primary y-axis for NH
line1 = ax1.plot(sum_timefirstvar_product_NH[:, 0], linestyle='-', color='b', label='NH - iceVolumeCell')
ax1.set_ylabel(f'Sum of {varNam2}', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create a twin axis for iceAreaCell
ax1_twin = ax1.twinx()
line2 = ax1_twin.plot(sum_timesecondvar_product_NH[:, 0], linestyle='-', color='g', label='NH - iceAreaCell')
ax1_twin.set_ylabel(f'Sum of {varNam3}', color='g')
ax1_twin.tick_params(axis='y', labelcolor='g')

# Plot iceVolumeCell on the primary y-axis for SH
line3 = ax2.plot(sum_timefirstvar_product_SH[:, 0], linestyle='-', color='b', label='SH - iceVolumeCell')
ax2.set_ylabel(f'Sum of {varNam2}', color='b')
ax2.tick_params(axis='y', labelcolor='b')

# Create a twin axis for iceAreaCell
ax2_twin = ax2.twinx()
line4 = ax2_twin.plot(sum_timesecondvar_product_SH[:, 0], linestyle='-', color='g', label='SH - iceAreaCell')
ax2_twin.set_ylabel(f'Sum of {varNam3}', color='g')
ax2_twin.tick_params(axis='y', labelcolor='g')

# Set x-axis limits and ticks (assuming you have defined 'ticks' and 'tick_labels')
ax2.set_xticks(ticks)
ax2.set_xticklabels(tick_labels, rotation=45)
ax2.set_xlim(start_year + start_month / 12, end_year + end_month / 12)
ax2.set_xlabel('Time (Year-Month)')

# Title and grid
ax1.set_title(f'Sum of {varNam2} and {varNam3} Over Time (NH)')
ax1.grid(True)
ax2.set_title(f'Sum of {varNam2} and {varNam3} Over Time (SH)')
ax2.grid(True)

# Show legend
lines = line1 + line2 + line3 + line4
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper right')

# Adjust layout
plt.tight_layout()



# Set x-axis limits and ticks
ax2.set_xticks(ticks)
ax2.set_xticklabels(tick_labels, rotation=45)
ax2.set_xlim(start_year + start_month / 12, end_year + end_month / 12)
ax2.set_xlabel('Time (Year-Month)')

# Title and grid
ax2.set_title(f'Sum of {varNam2} and {varNam3} Over Time (SH)')
ax2.grid(True)

# Show legend
lines = line3 + line4
labels = [l.get_label() for l in lines]
ax2.legend(lines, labels, loc='upper right')

# Save the figure
plt.savefig(f'SeparateSeasonalCycle{start_year}-{end_year}sum{varNam2}and{varNam3}_time_series_SH.png')
plt.show()
plt.close(fig)

