#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import xarray
import matplotlib.pyplot as plt
import os
import calendar
from mpl_toolkits.basemap import Basemap

# Model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
datafile = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.0030'

# Target months
target_months = [1, 4, 7, 10]

# Convert degrees to radians
deg2rad = np.pi / 180.0

regions = ['global', 'NSA', 'RS', 'SVA', 'NH', 'SH']

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

# Prompt the user for their choice
plot_choice = input("Do you want to plot by month or by year? (enter 'month' or 'year'): ").strip().lower()

# Create a folder to save plots
output_folder = 'seaice_plots_July9th'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

if plot_choice == 'month':
    # Prompt the user for the year
    year = input("Enter the year you want to plot (e.g., 0030): ")

    # Adjust the datafile path based on the specified year
    datafile = f'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.{year}'

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
            elif region_key == "NH":
                valid_indices = np.where(latCell > 60 * deg2rad)[0]
            elif region_key == "SH":
                valid_indices = np.where(latCell < -60 * deg2rad)[0]


            # Plot on the corresponding subplot
            ax = axs[idx // 2, idx % 2]
             # Draw map boundaries
            m.drawcoastlines()
            m.drawcountries()
            scatter = ax.scatter(lonCell[valid_indices], latCell[valid_indices], c=var1[valid_indices], cmap='bwr', s=10)

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
        plot_filename = f'Lat_Long_seaice_Output_{year}_{region_key.replace(" ", "_")}.png'
        plot_path = os.path.join(output_folder, plot_filename)
        plt.savefig(plot_path)

        # Close the figure to free up memory
        plt.close(fig)

elif plot_choice == 'year':
    # Prompt the user for the years
    years = [input(f"Enter year {i + 1} (e.g., 0030): ") for i in range(4)]

    # Adjust the datafile path based on the specified years
    datafiles = [f'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.{year}' for year in years]

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

        # Iterate over target years
        for idx, datafile in enumerate(datafiles):
            # Generate the output file name
            outputFileName = f'{datafile}-01-01.nc'

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
            elif region_key == "NH":
                valid_indices = np.where(latCell > 60 * deg2rad)[0]
            elif region_key == "SH":
                valid_indices = np.where(latCell < -60 * deg2rad)[0]


            # Plot on the corresponding subplot
            ax = axs[idx // 2, idx % 2]
            scatter = ax.scatter(lonCell[valid_indices], latCell[valid_indices], c=var1[valid_indices], cmap='bwr', s=10)

            # Set title and labels
            year_name = years[idx]
            ax.set_title(f'{year_name} - {region_key}')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

            # Add color bar to the last subplot
            if idx == len(years) - 1:
                cbar = fig.colorbar(scatter, ax=ax, orientation='vertical')
                cbar.set_label('Ice Area')

        # Adjust layout and save the figure for the current region
        fig.tight_layout()
        plot_filename = f'Lat_Long_seaice_Output_{region_key.replace(" ", "_")}_years.png'
        plot_path = os.path.join(output_folder, plot_filename)
        plt.savefig(plot_path)

        # Close the figure to free up memory
        plt.close(fig)

else:
    print("Invalid choice. Please run the script again and choose either 'month' or 'year'.")

print(f'Plots saved in {output_folder}')

