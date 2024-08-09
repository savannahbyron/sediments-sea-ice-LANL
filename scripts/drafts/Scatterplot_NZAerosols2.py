#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import xarray
import matplotlib.pyplot as plt
import os
import calendar

# Model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
datafile = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.0030'

# Target months
target_months = [1, 4, 7, 10]

# Convert degrees to radians
deg2rad = np.pi / 180.0

regions = ['global','NSA','RS','SVA']
# you can put this in a loop later:
for region in regions:

    #Define Lat long bounds in rads
    #North Slope of Alaska (NSA)
    if region=='NSA':
        lat_min_rad = 68 * deg2rad
        lat_max_rad = 71 * deg2rad
        lon_min_rad = (-165+360) * deg2rad
        lon_max_rad = (-145+360) * deg2rad
    #RussianShelf (RS)
    elif region=='RS':
        lat_min_rad = 68 * deg2rad
        lat_max_rad = 81 * deg2rad
        lon_min_rad = 20 * deg2rad
        lon_max_rad = 180 * deg2rad
    #Near Svalbard (SVA)
    elif region=='SVA':
        lat_min_rad = 74 * deg2rad
        lat_max_rad = 80 * deg2rad
        lon_min_rad = 10 * deg2rad
        lon_max_rad = 30 * deg2rad

# Ask user for confirmation to run all regions or just global
run_all_regions = input("Do you want to run all regions (NSA, Russian Shelf, Near Svalbard) [y/n]? ")

# Check user input
if run_all_regions.lower() == 'y':
    regions_to_plot = regions  # Plot all regions
elif run_all_regions.lower() == 'n':
    regions_to_plot = ['global']  # Plot only global view
else:
    print("Invalid input. Running global view by default.")
    regions_to_plot = ['global']  # Default to global view if input is invalid


# Create a folder to save plots
output_folder = 'seaice_plots'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Iterate over regions
for region_key in regions_to_plot:
    # Create a figure for each region
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))

    # Read mesh file
    mesh = xarray.open_dataset(meshFileName)
    latCell = mesh.variables['latCell'][:] 
    lonCell = mesh.variables['lonCell'][:]
    xCell = mesh.variables['xCell'][:]
    yCell = mesh.variables['yCell'][:]

    # Iterate over target months
    for idx, month in enumerate(target_months):
        # Generate the output file name
        outputFileName = f'{datafile}-{month:02d}-01.nc'

        # Read the corresponding output file
        output = xarray.open_dataset(os.path.join(runDir, outputFileName))

        # Extract variable for plotting
        varName = 'timeMonthly_avg_verticalAerosolsIceCell'
        var1 = output.variables[varName][0, :, 0]

        if region_key == 'global':
          valid_indices = np.arange(len(latCell))  # Use all indices
        elif region_key == 'NSA':
          valid_indices = np.where((latCell >= 68 * deg2rad) & (latCell <= 71 * deg2rad) & 
          (lonCell >= (-165 + 360) * deg2rad) & (lonCell <= (-145 + 360) * deg2rad))[0]
        elif region_key == 'RS':
          valid_indices = np.where((latCell >= 68 * deg2rad) & (latCell <= 81 * deg2rad) & 
          (lonCell >= 20 * deg2rad) & (lonCell <= 180 * deg2rad))[0]
        elif region_key == 'SVA':
          valid_indices = np.where((latCell >= 74 * deg2rad) & (latCell <= 80 * deg2rad) &
          (lonCell >= 10 * deg2rad) & (lonCell <= 30 * deg2rad))[0]

# Debug prints
        print(f'len(valid_indices): {len(valid_indices)}')

        # Plot on the corresponding subplot
        ax = axs[idx // 2, idx % 2]
        # Debug prints
       # print(f'len(lon_indices): {len(lon_indices)}')
       # print(f'len(lat_indices): {len(lat_indices)}')

        scatter = ax.scatter(lonCell[valid_indices], latCell[valid_indices], c=var1[valid_indices], cmap='bwr', s=5.4)

        # Set title and labels
        month_name = calendar.month_name[month]
        ax.set_title(f'{month_name} - {region_key}')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

        # Add color bar to the last subplot
        if idx == len(target_months) - 1:
            cbar = fig.colorbar(scatter, ax=ax, orientation='vertical')
            cbar.set_label('Ice Area')

    # Adjust layout and save the figure for the current region
    fig.tight_layout()
    plot_filename = f'Lat_Long_seaice_Output_Year30_{region_key.replace(" ", "_")}.png'
    plot_path = os.path.join(output_folder, plot_filename)
    plt.savefig(plot_path)

    # Close the figure to free up memory
    plt.close(fig)

print(f'Plots saved in {output_folder}')

