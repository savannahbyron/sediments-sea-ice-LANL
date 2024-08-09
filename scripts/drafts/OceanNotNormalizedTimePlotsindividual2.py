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
timeIV_product_NH = []
timeIV_product_SH = []
timeNZS_product_NH = []
timeNZS_product_SH = []

for year in years:
    for month in months:
        # Format the filename for the current month
        datafile_ice = datafile_ice_fmt.format(year=year, month=month)
        filepath_ice = os.path.join(runDir, datafile_ice)
        # Same thing for ocean
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
                #print(varNZ.dims)
                #print(varIV.dims)
                varNameNZ = 'verticalAerosolsIce in kg'
                varNameIV = 'iceVolumeCell in kg'
                varNameAC = 'iceAreaCell'
                varNameRFW = 'RiverrunoffFreshwaterFlux in kg s-1'

                # Get mesh data
                mesh = xr.open_dataset(meshFileName)
                latCell = mesh.variables['latCell']
                #print('read: '+meshFileName)
                nCells = mesh.sizes['nCells']
                latCell = mesh.variables['latCell'] # in radians
                lonCell = mesh.variables['lonCell']
                xCell = mesh.variables['xCell'] # in meters
                yCell = mesh.variables['yCell']
                areaCell = mesh.variables['areaCell']

                # Get NH indices
                NHindices = np.where(latCell > 60 * deg2rad)[0]
                SHindices = np.where(latCell < -60 * deg2rad)[0]

                # Sum for NZAerosols
                NZAER_product_NH = np.sum(varNZ.isel(nCells=NHindices) * varIV.isel(nCells=NHindices) * varAC.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(1, 2))
                timeNZ_product_NH.append([year + month / 12, np.sum(NZAER_product_NH)])

                NZAER_product_SH = np.sum(varNZ[:, SHindices, :] * varIV[:, SHindices] * varAC.isel(nCells=SHindices) * areaCell.isel(nCells=SHindices), axis=(1, 2))
                timeNZ_product_SH.append([year + month / 12, np.sum(NZAER_product_SH)])

                #RiverRunoff Freshwater
                RiverrunoffFreshWaterFlux_product_NH = np.sum(varRFW.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(0, 1))
                timeRFW_product_NH.append([year + month / 12, np.sum(RiverrunoffFreshWaterFlux_product_NH)])

                RiverrunoffFreshWaterFlux_product_SH = np.sum(varRFW.isel(nCells=SHindices) * areaCell.isel(nCells=SHindices), axis=(0, 1))
                timeRFW_product_SH.append([year + month / 12, np.sum(RiverrunoffFreshWaterFlux_product_SH)])

                #ICEVOLUME
                IV_product_NH = np.sum(varIV.isel(nCells=NHindices), axis=(0,1))
                timeIV_product_NH.append([year + month / 12, np.sum(IV_product_NH)])
                IV_product_SH = np.sum(varIV.isel(nCells=SHindices), axis=(0,1))
                timeIV_product_SH.append([year + month / 12, np.sum(IV_product_SH)])

        #NZAEROSOL CONCENTRATION

                # ICEVOLUME
                IV_product_NH = varIV.isel(nCells=NHindices)
                IV_product_SH = varIV.isel(nCells=SHindices)

                # Total volume of each cell, ignoring cells where volume is zero
                total_volume_NH = IV_product_NH.where(IV_product_NH != 0)
                total_volume_SH = IV_product_SH.where(IV_product_SH != 0)

                # NZAEROSOL CONCENTRATION
                NZS_product_NH = varNZ.isel(nCells=NHindices).sum(dim= 'nCells')
                NZS_product_SH = varNZ.isel(nCells=SHindices).sum(dim= 'nCells')

                # Calculate the ratio, ignoring cells where volume is zero
                NZAER_ratio_NH = NZS_product_NH / total_volume_NH.sum(dim='nCells')
                NZAER_ratio_SH = NZS_product_SH / total_volume_SH.sum(dim='nCells')

                # Append the results
                timeNZS_product_NH.append([year + month / 12, np.sum(NZAER_ratio_NH)])
                timeNZS_product_SH.append([year + month / 12, np.sum(NZAER_ratio_SH)])

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
timeIV_product_NH = np.array(timeIV_product_NH)
timeIV_product_SH = np.array(timeIV_product_SH)
timeNZS_product_NH = np.array(timeNZS_product_NH)
timeNZS_product_SH = np.array(timeNZS_product_SH)

# Define years and months based on your data
start_year = min(years)
end_year = max(years)
start_month = 1
end_month = 12

# Create ticks for x-axis with month names
month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
ticks = [(y + m / 12) for y in range(start_year, end_year + 1) for m in range(start_month, end_month + 1)]
tick_labels = [month_names[(month-1)%12] for year in range(start_year, end_year + 1) for month in range(start_month, end_month + 1)]

variables = {
    'Total Sediment': (timeNZ_product_NH, timeNZ_product_SH, 'blue', 'Vertical Aerosols Ice (NZ) Product'),
    'River Freshwater Flux': (timeRFW_product_NH, timeRFW_product_SH, 'green', 'RiverRunoff Freshwater Flux (RFW) Product'),
    'Ice Volume': (timeIV_product_NH, timeIV_product_SH, 'red', 'Ice Volume (IV) Product'),
    'Sediment Concentration': (timeNZS_product_NH, timeNZS_product_SH, 'purple', 'NZ Aerosol Concentration (NZS) Product')
}

for var_name, (time_data_NH, time_data_SH, color, title) in variables.items():
    plt.figure(figsize=(18, 10))
    plt.plot(time_data_NH[:, 0], time_data_NH[:, 1], label=f'{var_name} NH', linewidth=2, color=color)
    plt.plot(time_data_SH[:, 0], time_data_SH[:, 1], label=f'{var_name} SH', linewidth=2, linestyle='--', color=color)
    plt.title(title, fontsize=24)
    plt.xticks(ticks, tick_labels, rotation=45, fontsize=16)
    plt.xlim(start_year + start_month / 12, end_year + end_month / 12)
    plt.xlabel('Time (Month)', fontsize=20)
    plt.ylabel('Values', fontsize=20)
#    plt.legend(fontsize=16)
#plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{var_name}_NH_SH.png')
    plt.close()

