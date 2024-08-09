#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals

"""
Sum plots
Savannah Byron
This code helps me select variables to plot and visualize average values over time.
"""

# Model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
datafile_ice_fmt = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.00{year:02d}-{month:02d}-01.nc'
datafile_ocean_fmt = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpaso.hist.am.timeSeriesStatsMonthly.00{year:02d}-{month:02d}-01.nc'

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os

# Define constants
deg2rad = 3.14159 / 180.0
months = range(1, 13)
years = range(20, 22)

# Define available variables
available_variables = {
    'NZ': 'timeMonthly_avg_verticalAerosolsIceCell',
    'FW': 'timeMonthly_avg_oceanFreshWaterFlux',
    'VC': 'timeMonthly_avg_iceVolumeCategory',
    'AC': 'timeMonthly_avg_iceAreaCell'
}

# User input for variable selection
print("Available variables to plot:")
for key in available_variables:
    print(f"- {key}: {available_variables[key]}")
selected_keys = input("Enter the keys of the variables you want to plot (comma-separated, e.g., NZ,FW): ").split(",")

# Strip whitespace from keys
selected_keys = [key.strip() for key in selected_keys]

# Initialize arrays for selected variables
time_data = {key: [] for key in selected_keys if key in available_variables}

# Check if all selected variables are available
for key in selected_keys:
    if key not in available_variables:
        print(f"Warning: Variable '{key}' is not available. Skipping.")

for year in years:
    for month in months:
        datafile_ice = datafile_ice_fmt.format(year=year, month=month)
        filepath_ice = os.path.join(runDir, datafile_ice)

        if os.path.exists(filepath_ice) and os.path.exists(filepath_ocean):
            print(f'Processing file: {filepath_ice}')
            print(f'Processing file: {filepath_ocean}')


            # Get mesh data
            mesh = xr.open_dataset(meshFileName)
            latCell = mesh.variables['latCell']
            NHindices = np.where(latCell > 60 * deg2rad)[0]
            SHindices = np.where(latCell < -60 * deg2rad)[0]

            # Process selected variables
            for key in selected_keys:
                if key == 'NZ':
                    varNZ = output_ice[available_variables['NZ']]
                    varAC = output_ice[available_variables['AC']]
                    
                    NZAER_product_NH = np.sum(varNZ.isel(nCells=NHindices) * varAC.isel(nCells=NHindices) * mesh.variables['areaCell'].isel(nCells=NHindices), axis=(1, 2))
                    time_data['NZ'].append([year + month / 12, np.sum(NZAER_product_NH)])

                    NZAER_product_SH = np.sum(varNZ.isel(nCells=SHindices) * varAC.isel(nCells=SHindices) * mesh.variables['areaCell'].isel(nCells=SHindices), axis=(1, 2))
                    time_data['NZ'].append([year + month / 12, np.sum(NZAER_product_SH)])

                elif key == 'FW':
                    varFW = output_ocean[available_variables['FW']]
                    varAC = output_ice[available_variables['AC']]
                    
                    oceanFreshWaterFlux_product_NH = np.sum(varFW.isel(nCells=NHindices) * varAC.isel(nCells=NHindices), axis=(0, 1))
                    time_data['FW'].append([year + month / 12, np.sum(oceanFreshWaterFlux_product_NH)])

                    oceanFreshWaterFlux_product_SH = np.sum(varFW.isel(nCells=SHindices) * varAC.isel(nCells=SHindices), axis=(0, 1))
                    time_data['FW'].append([year + month / 12, np.sum(oceanFreshWaterFlux_product_SH)])

                else:
                    # For other variables like VC and AC, just take the mean
                    var = output_ice[available_variables[key]]
                    mean_value_NH = np.mean(var.isel(nCells=NHindices))
                    mean_value_SH = np.mean(var.isel(nCells=SHindices))
                    time_data[key].append([year + month / 12, mean_value_NH])
                    time_data[key].append([year + month / 12, mean_value_SH])

            # Close datasets
            output_ice.close()
            output_ocean.close()
            mesh.close()

# Convert lists to numpy arrays
for key in time_data.keys():
    time_data[key] = np.array(time_data[key])

# Define years and months for ticks
start_year = min(years)
end_year = max(years)
ticks = [(y + m / 12) for y in range(start_year, end_year + 1) for m in range(1, 13)]
tick_labels = [f'{year}-{month:02d}' for year in range(start_year, end_year + 1) for month in range(1, 13)]

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 18))  # Two rows, one column

# Plot for Northern Hemisphere
for key in selected_keys:
    ax1.plot(time_data[key][:len(months), 0], time_data[key][:len(months), 1], label=f'{key} (NH)')

ax1.set_ylabel('Values (NH)')
ax1.tick_params(axis='y')
ax1.set_xticks(ticks)
ax1.set_xticklabels(tick_labels, rotation=45)
ax1.set_xlim(start_year + 1/24, end_year + 12/24)
ax1.set_title('Selected Variables Over Time (NH)')
ax1.grid(True)
ax1.legend()

# Plot for Southern Hemisphere
for key in selected_keys:
    ax2.plot(time_data[key][len(months):, 0], time_data[key][len(months):, 1], label=f'{key} (SH)')

ax2.set_ylabel('Values (SH)')
ax2.tick_params(axis='y')
ax2.set_xticks(ticks)
ax2.set_xticklabels(tick_labels, rotation=45)
ax2.set_xlim(start_year + 1/24, end_year + 12/24)
ax2.set_title('Selected Variables Over Time (SH)')
ax2.grid(True)
ax2.legend()

# Save the figure
plt.savefig(f'SeparateSeasonalCycle{start_year}-{end_year}_selected_variables.png')
plt.show()
plt.close(fig)

