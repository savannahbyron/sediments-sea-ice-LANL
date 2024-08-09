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

# Define Lat long bounds in rads
# North Slope of Alaska (NSA)
lat_min_rad_NSA = 68 * deg2rad
lat_max_rad_NSA = 71 * deg2rad
lon_min_rad_NSA = (-165 + 360) * deg2rad
lon_max_rad_NSA = (-145 + 360) * deg2rad

# RussianShelf (RS)
lat_min_rad_RS = 68 * deg2rad
lat_max_rad_RS = 81 * deg2rad
lon_min_rad_RS = 20 * deg2rad
lon_max_rad_RS = 180 * deg2rad

# Near Svalbard (SVA)
lat_min_rad_SVA = 74 * deg2rad
lat_max_rad_SVA = 80 * deg2rad
lon_min_rad_SVA = 10 * deg2rad
lon_max_rad_SVA = 30 * deg2rad

# Ask user for confirmation to run all regions or just global
run_all_regions = input("Do you want to run all regions (NSA, Russian Shelf, Near Svalbard) [y/n]? ")

# Check user input
if run_all_regions.lower() == 'y':
    regions_to_plot = ['NSA', 'RS', 'SVA']  # Plot all regions
elif run_all_regions.lower() == 'n':
    regions_to_plot = ['Everywhere']  # Plot only global view
else:
    print("Invalid input. Running global view by default.")
    regions_to_plot = ['Everywhere']  # Default to global view if input is invalid

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

        if region_key == 'Everywhere':
            lat_indices = slice(None)  # Use all latitude indices
            lon_indices = slice(None)  # Use all longitude indices
        else:
            if region_key == 'NSA':
                lat_min = lat_min_rad_NSA
                lat_max = lat_max_rad_NSA
                lon_min = lon_min_rad_NSA
                lon_max = lon_max_rad_NSA
            elif region_key == 'RS':
                lat_min = lat_min_rad_RS
                lat_max = lat_max_rad_RS
                lon_min = lon_min_rad_RS
                lon_max = lon_max_rad_RS
            elif region_key == 'SVA':
                lat_min = lat_min_rad_SVA
                lat_max = lat_max_rad_SVA
                lon_min = lon_min_rad_SVA
                lon_max = lon_max_rad_SVA

            lat_indices = np.where((latCell >= lat_min) & (latCell <= lat_max))[0]
            lon_indices = np.where((lonCell >= lon_min) & (lonCell <= lon_max))[0]

        # Plot on the corresponding subplot
        ax = axs[idx // 2, idx % 2]
        scatter = ax.scatter(lonCell[lon_indices], latCell[lat_indices], c=var1[lat_indices], cmap='bwr', s=0.4)

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

