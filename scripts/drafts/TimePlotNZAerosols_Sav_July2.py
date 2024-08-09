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

# Convert degrees to radians
deg2rad = 3.14159 / 180.0

# Define the size of nCells 
nCells = 465044

# Make  array for sum over months
pickle_filename = 'sum_over_months.pkl'
if os.path.exists(pickle_filename):
    with open(pickle_filename, 'rb') as f:
        sum_over_months = pickle.load(f)
    print('Loaded sum_over_months from pickle file.')
else:
# Initialize an array to store the sum of NZaerosol over months and cells
    sum_over_months = np.zeros((len(time_array), nCells))

time_array = []
for year in range(10, 41):  # Adjusted to 10 to 40 (30 years)
    for month in range(1, 13):
        time_array.append(f'{year}-{month:02d}')

# Define the months to iterate through (e.g., 1 to 12 for all months in a year)
months = range(1, 13)  # Adjusted to include all 12 months

# Define the months to iterate through (e.g., 1 to 12 for all months in a year)
years = range(10, 41)  # Adjusted to include all 


# Create a figure for the cumulative sum plot
fig, ax = plt.subplots(figsize=(10, 6))

for year in years:

# Iterate over the target months
    for month in months:
        # Format the filename for the current month
        datafile = datafile_pattern.format(year=year, month=month)
        filepath = os.path.join(runDir, datafile)
    
        # Check if the file exists 
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
    
                # Get mesh data 
                mesh = xr.open_dataset(meshFileName)
                latCell = mesh.variables['latCell']
    
                # Get NH indices
                NHindices = np.where(latCell > 60 * deg2rad)[0]
    
                # Sum over specific range (time steps and NH indices)
                var1D_sum = np.sum(var1[:, NHindices,:], axis=(1, 2))  # Sum over NHindices and can do first 7 time steps by doing NHindices, 0:7
    
                # Accumulate the sum over months
                sum_over_months[month - 1] = var1D_sum  # Store in the correct month index
                with open(pickle_filename, 'wb') as f:
                    pickle.dump(sum_over_months, f)
                print('Saved sum_over_months to pickle file')


                # Close datasets to free up resources
                output.close()
                mesh.close()
    
            else:
                print(f'Output file not found: {output_filepath}')
    
        else:
            print(f'File not found: {filepath}')

# Plot the cumulative sum over months as a time series
time_months = np.arange(1, 12)  # Months for x-axis
time_years = np.arange(10, 40)  # years for x-axis

total_months = (40 - 10 + 1) * 12
xaxis = np.arange(0, total_months) / 12
#xaxis = np.arange(0, np.size(sum_over_months[:, NHindices],0))/12 
for NHindices in range(nCells):
    ax.plot(xaxis, sum_over_months[:, NHindices], marker='o', linestyle='-', label=f'Cell Index {NHindices}')

# Set labels and title
ax.set_xlabel('Month and Year')
ax.set_ylabel('Sum of NZaerosol')
ax.set_title('Sum of NZaerosol over Months and Years')

# Set the x-axis ticks and labels
ax.set_xticks(np.arange(0, len(years) * 12, 12))
ax.set_xticklabels([f'Month {month} Year {year}' for year in years for month in range(1, 13)], rotation=90)

# Save the figure
plt.savefig('sum_NZaerosol.png')
plt.close(fig)  # Close the figure to free up memory
