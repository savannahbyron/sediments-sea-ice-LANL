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
 

for year in years:
    for month in months:
        # Format the filename for the current month
        datafile_ice = datafile_ice_fmt.format(year=year, month=month)
        filepath_ice = os.path.join(runDir, datafile_ice)
#same thing for ocean
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
                #print(varNZ.dims)
                #print(varIV.dims)
                varNameNZ = 'verticalAerosolsIce in kg'
                varNameIV = 'iceVolumeCell in kg'
                varNameAC = 'iceAreaCell'
                varNameFW = 'oceanFreshwaterFlux in kg s-1'

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
     
                # Overlay the plot with timeMonthly_avg_oceanFreshWaterFlux * timeMonthly_avg_iceAreaCell
                oceanFreshWaterFlux_product_NH = np.sum(varFW.isel(nCells=NHindices) * varAC.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(0, 1)) 
                timeFW_product_NH.append([year + month / 12, np.sum(oceanFreshWaterFlux_product_NH)])

                oceanFreshWaterFlux_product_SH = np.sum(varFW.isel(nCells=SHindices) * varAC.isel(nCells=SHindices) * areaCell.isel(nCells=SHindices), axis=(0, 1)) 
                timeFW_product_SH.append([year + month / 12, np.sum(oceanFreshWaterFlux_product_SH)])

                #ICEVOLUME
                IV_product_NH = np.sum(varIV.isel(nCells=NHindices), axis=(0,1))
                timeIV_product_NH.append([year + month / 12, np.sum(IV_product_NH)])
                IV_product_SH = np.sum(varIV.isel(nCells=SHindices), axis=(0,1))
                timeIV_product_SH.append([year + month / 12, np.sum(IV_product_SH)])
                
                #NZAEROSOL CONCENTRATION
                NZS_product_NH = np.sum(varNZ.isel(nCells=NHindices), axis=(0,1))
                timeNZS_product_NH.append([year + month / 12, np.sum(NZS_product_NH)])
                NZS_product_SH = np.sum(varNZ.isel(nCells=SHindices), axis=(0,1))
                timeNZS_product_SH.append([year + month / 12, np.sum(NZS_product_SH)])


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
timeFW_product_NH = np.array(timeFW_product_NH)
timeFW_product_SH = np.array(timeFW_product_SH)
timeIV_product_NH = np.array(timeIV_product_NH)
timeIV_product_SH = np.array(timeIV_product_SH)
timeNZS_product_NH = np.array(timeNZS_product_NH)
timeNZS_product_SH = np.array(timeNZS_product_SH)

# Define years and months based on your data
start_year = min(years)
end_year = max(years)

fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True, sharey=True)

# Create ticks for x-axis
start_month = 1
end_month = 12
ticks = [(y + m / 12) for y in range(start_year, end_year + 1) for m in range(start_month, end_month + 1)]
tick_labels = [f'{year}-{month:02d}' for year in range(start_year, end_year + 1) for month in range(start_month, end_month + 1)]
# Set x-axis limits and ticks
#ax1.set_xticks(ticks)
#ax1.set_xticklabels(tick_labels, rotation=45)
#ax1.set_xlim(start_year + start_month / 12, end_year + end_month / 12)
#ax1.set_xlabel('Time (Year-Month)')

# Combine all NH and SH data to determine the y-axis limits
NH_data = np.concatenate([timeNZ_product_NH, timeFW_product_NH, timeIV_product_NH, timeNZS_product_NH])
SH_data = np.concatenate([timeNZ_product_SH, timeFW_product_SH, timeIV_product_SH, timeNZS_product_SH])


min_value = min(np.min(NH_data), np.min(SH_data))
max_value = max(np.max(NH_data), np.max(SH_data))

# Plot NH data
axs[0].plot(timeNZ_product_NH, label='NZ Product NH')
axs[0].plot(timeFW_product_NH, label='FW Product NH')
axs[0].plot(timeIV_product_NH, label='IV Product NH')
axs[0].plot(timeNZS_product_NH, label='NZS Product NH')
axs[0].set_ylim([min_value, max_value])
axs[0].set_title('Northern Hemisphere Products')
axs[0].legend()

# Plot SH data
axs[1].plot(timeNZ_product_SH, label='NZ Product SH')
axs[1].plot(timeFW_product_SH, label='FW Product SH')
axs[1].plot(timeIV_product_SH, label='IV Product SH')
axs[1].plot(timeNZS_product_SH, label='NZS Product SH')
axs[1].set_ylim([min_value, max_value])
axs[1].set_title('Southern Hemisphere Products')
axs[1].legend()

# Add common labels
plt.xlabel('Time')
fig.text(0.04, 0.5, 'Normalized Values', va='center', rotation='vertical')

plt.tight_layout()
plt.savefig('NH_SH_Products2.png')
plt.show()
