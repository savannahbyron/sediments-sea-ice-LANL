#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals

"""
sum plot
Savannah Byron
This code helps me select a variable and plot the average value in time, i am currently wokring on adding different location places, clearer xaxis labels
"""
  
############################## model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
datafile_pattern = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.00{year:02d}-{month:02d}-01.nc'

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os

#define variables
deg2rad = 3.14159 / 180.0
nCells = 465044
months = range(1, 13)  # Adjusted to include all 12 months
years = range(20, 22)  # Adjusted to include all 

# Choose variable to process
var_choice = input("Choose variable to process (1 for iceAreaCell, 2 for verticalAerosolsIceCell, 3 for iceVolumeCell): ")
if var_choice == '1':
    varName = 'timeMonthly_avg_iceAreaCell'#Fraction of grid cell covered in sea ice
    varNam2 = 'iceAreaCell in m'
elif var_choice == '2':
    varName = 'timeMonthly_avg_verticalAerosolsIceCell'  #units are kg m-3
    varNam2 = 'verticalAerosolsIceCell in kg m-2'
elif var_choice == '3':
    varName = 'timeMonthly_avg_iceVolumeCell' #units are m
    varNam2 = 'iceVolumeCell'
else:
    print("Invalid choice. Exiting.")
    exit()

# Initialize arrays dynamically based on variable choice
time_NZaerosols_NH = []#this i want to change to be a function of the above, so the name changes but we shall see
time_NZaerosols_SH = []
time_product_NH = []
time_product_SH = []

 
for year in years:
    for month in months:
        # Format the filename for the current month
        datafile = datafile_pattern.format(year=year, month=month)
        filepath = os.path.join(runDir, datafile)
    
        # Check if the file exists 
        if os.path.exists(filepath):
            print(f'Processing file: {filepath}')
    
            # Read the corresponding output file
            output_filepath = os.path.join(runDir, datafile)
    
            if os.path.exists(output_filepath):
                # Open the dataset
                output = xr.open_dataset(output_filepath)
    
                # Get variables from the current output dataset
                varName = 'timeMonthly_avg_verticalAerosolsIceCell'
                var1 = output.variables[varName]
                var2 = output.variables['timeMonthly_avg_iceVolumeCell']
                #print(var1.dims)
                #print(var2.dims)

                # Get mesh data 
                mesh = xr.open_dataset(meshFileName)
                latCell = mesh.variables['latCell']
                #print('read: '+meshFileName)
                nCells = mesh.sizes['nCells']
                latCell = mesh.variables['latCell'] # in radians
                lonCell = mesh.variables['lonCell']
                xCell = mesh.variables['xCell'] # in meters
                yCell = mesh.variables['yCell']

                # Get NH indices
                NHindices = np.where(latCell > 60 * deg2rad)[0]
                SHindices = np.where(latCell < -60 * deg2rad)[0]

                # Sum over NH indices
                var_product_NH = np.sum(var1[:, NHindices, :] * var2[:, NHindices], axis=(1, 2))
                time_product_NH.append([year + month / 12, np.sum(var_product_NH)])

                # Sum over SH indices
                var_product_SH = np.sum(var1[:, SHindices, :] * var2[:, SHindices], axis=(1, 2))
                time_product_SH.append([year + month / 12, np.sum(var_product_SH)])
    
              # Close datasets to free up resources
                output.close()
                mesh.close()
    
            else:
                print(f'Output file not found: {output_filepath}')
    
        else:
            print(f'File not found: {filepath}')

# Convert lists to numpy arrays
time_product_NH = np.array(time_product_NH)
time_product_SH = np.array(time_product_SH)
# Plotting
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(time_product_NH[:, 0], time_product_NH[:, 1], marker='o', linestyle='-', color='b', label='NH')
ax.plot(time_product_SH[:, 0], time_product_SH[:, 1], marker='o', linestyle='-', color='g', label='SH')

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
ax.set_ylabel(f'Sum of {varNam2}')
ax.set_title(f'Sum of {varNam2} over Time - with Volume')
ax.grid(True)
ax.legend()

# Save the figure
plt.savefig(f'years{start_year}-{end_year}sum_of{varName}_NH_SH_label_time_series.png')
plt.show()
plt.close(fig)

