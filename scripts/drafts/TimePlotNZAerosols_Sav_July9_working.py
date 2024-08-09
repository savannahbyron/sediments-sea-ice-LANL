#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals

"""
sum plot
Savannah Byron
July 2024 I was working with this with mark on July 2, 2024, this works by printing the sum of NZaerosols in loop forms 
"""
  
############################## model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
datafile_pattern = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.00{year:02d}-{month:02d}-01.nc'

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import pickle

#define variables
deg2rad = 3.14159 / 180.0
nCells = 465044
months = range(1, 13)  # Adjusted to include all 12 months
years = range(10, 41)  # Adjusted to include all 

# Choose variable to process
var_choice = input("Choose variable to process (1 for iceAreaCell, 2 for verticalAerosolsIceCell, 3 for iceVolumeCell): ")
if var_choice == '1':
    varName = 'timeMonthly_avg_iceAreaCell'
elif var_choice == '2':
    varName = 'timeMonthly_avg_verticalAerosolsIceCell'
elif var_choice == '3':
    varName = 'timeMonthly_avg_iceVolumeCell'
else:
    print("Invalid choice. Exiting.")
    exit()

# Initialize arrays dynamically based on variable choice
time_NZaerosols_NH = []
time_NZaerosols_SH = []

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
    
                # Get mesh data 
                mesh = xr.open_dataset(meshFileName)
                latCell = mesh.variables['latCell']
    
                # Get NH indices
                NHindices = np.where(latCell > 60 * deg2rad)[0]
                SHindices = np.where(latCell < -60 * deg2rad)[0]

                # Sum over NH indices
                var1D_sum_NH = np.sum(var1[:, NHindices, :], axis=(1, 2))
                time_NZaerosols_NH.append([year + month / 12, np.sum(var1D_sum_NH)])

                # Sum over SH indices
                var1D_sum_SH = np.sum(var1[:, SHindices, :], axis=(1, 2))
                time_NZaerosols_SH.append([year + month / 12, np.sum(var1D_sum_SH)])
    
              # Close datasets to free up resources
                output.close()
                mesh.close()
    
            else:
                print(f'Output file not found: {output_filepath}')
    
        else:
            print(f'File not found: {filepath}')

# Convert lists to numpy arrays
time_NZaerosols_NH = np.array(time_NZaerosols_NH)
time_NZaerosols_SH = np.array(time_NZaerosols_SH)
# Plotting
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(time_NZaerosols_NH[:, 0], time_NZaerosols_NH[:, 1], marker='o', linestyle='-', color='b', label=f'{hemisphere_label} NH')
ax.plot(time_NZaerosols_SH[:, 0], time_NZaerosols_SH[:, 1], marker='o', linestyle='-', color='g', label=f'{hemisphere_label} SH')
ax.set_xlabel('Time-Years')
ax.set_ylabel(f'Sum of {varName}')
ax.set_title(f'Sum of {varName} over Time - {hemisphere_label}')
ax.grid(True)
ax.legend()

# Save the figure
plt.savefig(f'sum_{hemisphere_label}_time_series.png')
plt.show()
plt.close(fig)

