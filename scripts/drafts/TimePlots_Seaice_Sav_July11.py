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
                varName = 'timeMonthly_avg_verticalAerosolsIceCell'
                var1 = output_ice.variables['timeMonthly_avg_verticalAerosolsIceCell']
                var2 = output_ice.variables['timeMonthly_avg_iceVolumeCell']
                var3 = output_ice.variables['timeMonthly_avg_iceAreaCell']
                var4 = output_ice.variables['timeMonthly_avg_oceanFreshWaterFlux']
                #print(var1.dims)
                #print(var2.dims)
                varNam1 = 'verticalAerosolsIce in kg'
                varNam2 = 'iceVolumeCell in kg'
                varNam3 = 'iceAreaCell'
                varNam4 = 'oceanFreshwaterFlux in kg s-1'

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
                NZAER_product_NH = np.sum(var1.isel(nCells=NHindices) * var2.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(1, 2))
                #NZAER_product_NH = np.sum(var1.isel(nCells=NHindices) * var2.isel(nCells=NHindices) * var3.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(1, 2))
                timeNZ_product_NH.append([year + month / 12, np.sum(NZAER_product_NH)])

                NZAER_product_SH = np.sum(var1[:, SHindices, :] * var2[:, SHindices] * areaCell.isel(nCells=SHindices), axis=(1, 2))
                #NZAER_product_SH = np.sum(var1[:, SHindices, :] * var2[:, SHindices] * var3.isel(nCells=SHindices) * areaCell.isel(nCells=SHindices), axis=(1, 2))
                timeNZ_product_SH.append([year + month / 12, np.sum(NZAER_product_SH)])
     
                # Overlay the plot with timeMonthly_avg_oceanFreshWaterFlux * timeMonthly_avg_iceAreaCell
                oceanFreshWaterFlux_product_NH = np.sum(var4.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(0, 1)) 
                #oceanFreshWaterFlux_product_NH = np.sum(var4.isel(nCells=NHindices) * var3.isel(nCells=NHindices) * areaCell.isel(nCells=NHindices), axis=(0, 1)) 
                timeFW_product_NH.append([year + month / 12, np.sum(oceanFreshWaterFlux_product_NH)])

                oceanFreshWaterFlux_product_SH = np.sum(var4.isel(nCells=SHindices) * areaCell.isel(nCells=SHindices), axis=(0, 1)) 
                #oceanFreshWaterFlux_product_SH = np.sum(var4.isel(nCells=SHindices) * var3.isel(nCells=SHindices) * areaCell.isel(nCells=SHindices), axis=(0, 1)) 
                timeFW_product_SH.append([year + month / 12, np.sum(oceanFreshWaterFlux_product_SH)])

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

# Plotting
fig, ax = plt.subplots(figsize=(15, 9))
ax.plot(timeNZ_product_NH[:, 0], timeNZ_product_NH[:, 1],linestyle='-', color='b', label='NH -NZAerosols')
ax.plot(timeNZ_product_SH[:, 0], timeNZ_product_SH[:, 1],linestyle='--', color='b', label='SH - NZAerosols')

ax2 = ax.twinx()

ax2.plot(timeFW_product_NH[:,0], timeFW_product_NH[: ,1], linestyle='-', color='g', label='NH - oceanFreshWaterFlux ')
ax2.plot(timeFW_product_SH[:,0], timeFW_product_SH[: ,1], linestyle='--', color='g', label='SH - oceanFreshWaterFlux ', linewidth =2)

# Determine x-axis limits and ticks 
start_year = min(years)
start_month = 1
end_year = max(years)
end_month = 12
ticks = [(y + m / 12) for y in range(start_year, end_year + 1) for m in range(start_month, 13)]
tick_labels = [f'{year}-{month:02d}' for year in range(start_year, end_year + 1) for month in range(start_month, 13)]
ax.set_xticks(ticks)
ax.set_xticklabels(tick_labels, rotation=45)
ax.set_xlim(start_year + start_month / 12, end_year + end_month / 12)

ax.set_xlabel('Time (Year-Month)')  # Changed
ax.set_ylabel('Sum of NZAerosol')
ax2.set_ylabel('Sum of oceanFreshWaterFlux')
ax.set_title(f'Sum of over Time - with Volume')
ax.grid(True)
ax.legend(loc='upper left')
ax2.legend(loc=0)
ax.tick_params(axis='y', labelcolor='blue')
ax2.tick_params(axis='y', labelcolor='green')

# Save the figure
plt.savefig(f'Overlaidyears{start_year}-{end_year}sum_NH_SH_time_series.png') 
plt.show()
plt.close(fig)

