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
timeRFW_product_NH = []
timeRFW_product_SH = []
timeNZS_product_NH = []
timeNZS_product_SH = []

for year in years:
    for month in months:
        # Format the filename for the current month
        datafile_ice = datafile_ice_fmt.format(year=year, month=month)
        filepath_ice = os.path.join(runDir, datafile_ice)
        # same thing for ocean
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
                varRFW = output_ocean.variables['timeMonthly_avg_riverRunoffFlux']
                
                # Get mesh data
                mesh = xr.open_dataset(meshFileName)
                latCell = mesh.variables['latCell']
                lonCell = mesh.variables['lonCell']
                xCell = mesh.variables['xCell']  # in meters
                yCell = mesh.variables['yCell']
                areaCell = mesh.variables['areaCell']

                # Get NH and SH indices
                NHindices = np.where(latCell > 50 * deg2rad)[0]
                SHindices = np.where(latCell < -50 * deg2rad)[0]

                # Process NH data
                NZ_layers_NH = varNZ.isel(nCells=NHindices)
                RFW_product_NH = varRFW.isel(nCells=NHindices)
                IV_product_NH = varIV.isel(nCells=NHindices)
                AC_product_NH = varAC.isel(nCells=NHindices)
                areaCell_NH = areaCell.isel(nCells=NHindices)

                for layer in range(varNZ.shape[2]):
                    NZAER_product_NH = np.sum(NZ_layers_NH[:, :, layer] * IV_product_NH * AC_product_NH * areaCell_NH, axis=(1))
                    timeNZ_product_NH.append([year + month / 12, layer, np.sum(NZAER_product_NH)])

                RiverrunoffFreshWaterFlux_product_NH = np.sum(RFW_product_NH * AC_product_NH * areaCell_NH, axis=(0, 1))
                timeRFW_product_NH.append([year + month / 12, np.sum(RiverrunoffFreshWaterFlux_product_NH)])

                NZS_product_NH = NZ_layers_NH.sum(dim='nCells')
                total_volume_NH = IV_product_NH.where(IV_product_NH != 0).sum(dim='nCells')
                NZAER_ratio_NH = (NZS_product_NH / total_volume_NH).squeeze()

                # Debugging: print shapes and values
                #print(f'Year: {year}, Month: {month}, NZAER_ratio_NH shape: {NZAER_ratio_NH.shape}')
                #print(NZAER_ratio_NH)

                for layer in range(NZAER_ratio_NH.shape[0]):
                    timeNZS_product_NH.append([year + month / 12, layer, NZAER_ratio_NH[layer].values])

                # Process SH data
                NZ_layers_SH = varNZ.isel(nCells=SHindices)
                RFW_product_SH = varRFW.isel(nCells=SHindices)
                IV_product_SH = varIV.isel(nCells=SHindices)
                AC_product_SH = varAC.isel(nCells=SHindices)
                areaCell_SH = areaCell.isel(nCells=SHindices)

                for layer in range(varNZ.shape[2]):
                    NZAER_product_SH = np.sum(NZ_layers_SH[:, :, layer] * IV_product_SH * AC_product_SH * areaCell_SH, axis=(1))
                    timeNZ_product_SH.append([year + month / 12, layer, np.sum(NZAER_product_SH)])

                RiverrunoffFreshWaterFlux_product_SH = np.sum(RFW_product_SH * AC_product_SH * areaCell_SH, axis=(0, 1))
                timeRFW_product_SH.append([year + month / 12, np.sum(RiverrunoffFreshWaterFlux_product_SH)])

                NZS_product_SH = NZ_layers_SH.sum(dim='nCells')
                total_volume_SH = IV_product_SH.where(IV_product_SH != 0).sum(dim='nCells')
                NZAER_ratio_SH = (NZS_product_SH / total_volume_SH).squeeze()

                # Debugging: print shapes and values
                #print(f'Year: {year}, Month: {month}, NZAER_ratio_SH shape: {NZAER_ratio_SH.shape}')
                #print(NZAER_ratio_SH)

                for layer in range(NZAER_ratio_SH.shape[0]):
                    timeNZS_product_SH.append([year + month / 12, layer, NZAER_ratio_SH[layer].values])

                # Close datasets to free up resources
                output_ice.close()
                mesh.close()

            else:
                print(f'Output file not found: {output_filepath_ice}')
        else:
            print(f'File not found: {filepath_ice}')

# Convert lists to numpy arrays
timeNZ_product_NH = np.array(timeNZ_product_NH)
timeNZ_product_SH = np.array(timeNZ_product_SH)
timeRFW_product_NH = np.array(timeRFW_product_NH)
timeRFW_product_SH = np.array(timeRFW_product_SH)
timeNZS_product_NH = np.array(timeNZS_product_NH)
timeNZS_product_SH = np.array(timeNZS_product_SH)

# Define years and months based on your data
start_year = min(years)
end_year = max(years)
ticks = [(y + m / 12) for y in range(start_year, end_year + 1) for m in range(1, 13)]
tick_labels = [f'{year}-{month:02d}' for year in range(start_year, end_year + 1) for month in range(1, 13)]

# Plotting function
def plot_data(time_data, rfw_data, title, filename, varNZ_shape):
    fig, axs = plt.subplots(3, 1, figsize=(15, 18))
    
    # NZ Aerosols plot
    for layer in range(varNZ_shape[2]):
        layer_data = time_data[time_data[:, 1] == layer]
        axs[0].plot(layer_data[:, 0], layer_data[:, 2], label=f'Layer {layer}')
    axs[0].set_title(f'{title} - NZ Aerosols', fontsize=16)
    axs[0].legend(loc='upper right', fontsize=12)
    axs[0].set_xticks(ticks)
    axs[0].set_xticklabels(tick_labels, rotation=45, fontsize=12)
    axs[0].set_xlim(start_year + 1 / 12, end_year + 12 / 12)
    axs[0].set_xlabel('Time (Year-Month)', fontsize=14)

    # River Runoff Freshwater plot
    axs[1].plot(rfw_data[:, 0], rfw_data[:, 1], label='River Runoff Freshwater Flux')
    axs[1].set_title(f'{title} - River Runoff Freshwater Flux', fontsize=16)
    axs[1].legend(loc='upper right', fontsize=12)
    axs[1].set_xticks(ticks)
    axs[1].set_xticklabels(tick_labels, rotation=45, fontsize=12)
    axs[1].set_xlim(start_year + 1 / 12, end_year + 12 / 12)
    axs[1].set_xlabel('Time (Year-Month)', fontsize=14)

    # NZ Aerosol Concentration plot
    for layer in range(varNZ_shape[2]):
        layer_data = time_data[time_data[:, 1] == layer]
        axs[2].plot(layer_data[:, 0], layer_data[:, 2], label=f'Layer {layer}')
    axs[2].set_title(f'{title} - NZ Aerosol Concentration', fontsize=16)
    axs[2].legend(loc='upper right', fontsize=12)
    axs[2].set_xticks(ticks)
    axs[2].set_xticklabels(tick_labels, rotation=45, fontsize=12)
    axs[2].set_xlim(start_year + 1 / 12, end_year + 12 / 12)
    axs[2].set_xlabel('Time (Year-Month)', fontsize=14)

    plt.tight_layout()
    plt.savefig(filename)
    plt.show()

# Plot data for NH
plot_data(timeNZ_product_NH, timeRFW_product_NH, 'Northern Hemisphere', 'NH_Products.png', varNZ.shape)

# Plot data for SH
plot_data(timeNZ_product_SH, timeRFW_product_SH, 'Southern Hemisphere', 'SH_Products.png', varNZ.shape)

