#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
sum plot
Savannah Byron
July 2024
"""

############################## model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
datafile_pattern = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.0030-{month:02d}-01.nc'

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import calendar

# Convert degrees to radians
deg2rad = 3.14159 / 180.0


# Define the size of nCells based on your data
nCells = 465044

# Initialize an array to accumulate the sum
sum_over_months = np.zeros((12, nCells))  # 

# Define the months you want to iterate through (e.g., 1 to 12 for all months in a year)
months = range(1, 12)

# Create a figure for the cumulative sum plot
fig, ax = plt.subplots(figsize=(10, 6))

# Iterate over the target months
for month in months:
    # Format the filename for the current month
    datafile = datafile_pattern.format(month=month)
    filepath = os.path.join(runDir, datafile)

    # Check if the file exists before trying to open it
    if os.path.exists(filepath):
        print(f'Processing file: {filepath}')

        # Read the corresponding output file
        outputFileName = f'{datafile}'
        output_filepath = os.path.join(runDir, outputFileName)

        if os.path.exists(output_filepath):
            # Open the dataset
            output = xr.open_dataset(output_filepath)

            # Get variables from the current output dataset
            varName = 'timeMonthly_avg_verticalAerosolsIceCell'
            var1 = output.variables[varName]
            print(f'var1 shape: {var1.shape}')  # Add this line to check shape DEBUG

            # Get mesh data (assuming this is constant across months)
            mesh = xr.open_dataset(meshFileName)
            latCell = mesh.variables['latCell']
            lonCell = mesh.variables['lonCell']
            xCell = mesh.variables['xCell']
            yCell = mesh.variables['yCell']

            # Get NH indices
            NHindices = np.where(latCell > 60 * deg2rad)[0]

            # Sum over specific range
            print(f'var1 shape: {var1[:, NHindices, 0:7].shape}')  # Check the shape of sliced var
            var1D_sum = np.sum(var1[:, NHindices, 0:7], axis = (0,1,2))
            print(f'var1D_sum shape: {var1D_sum.shape}')  # Add this line to check shape DEBUG
            print(f'var1D_sum value: {var1D_sum}')  # Check the value of var1D_sum

            # Accumulate the sum over months
            sum_over_months += var1D_sum

            # Close datasets to free up resources
            output.close()
            mesh.close()

        else:
            print(f'Output file not found: {output_filepath}')

    else:
        print(f'File not found: {filepath}')

# Plot the cumulative sum over months
ax.plot(sum_over_months, marker='o', linestyle='-', color='b')

# Set labels and title
ax.set_xlabel('Cell Index')
ax.set_ylabel('Sum of NZaerosol')
ax.set_title('Cumulative Sum of NZaerosol over Months')

# Save the figure
plt.savefig('cumulative_sum_NZaerosol.png')
plt.close(fig)  # Close the figure to free up memory

