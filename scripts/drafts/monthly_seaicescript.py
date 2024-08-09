#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
simple scatter plot
Savannah Byron
May 2024
"""

############################## model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
outputfilename = 'dustFromOceanOn.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.0012'
#import needed directories/programs
import numpy as np
import xarray
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

target_months = [1, 4, 7, 10]

# Iterate over target months
for month in target_months:
    # Generate the output file name
    outputFileName = f'{outputfilename}-{month:02d}-01.nc'
    # Read the corresponding output file
    output = xarray.open_dataset(os.path.join(runDir, outputFileName))

    if month == 1:
        subplot_index = (0, 0)
    elif month == 4:
        subplot_index = (0, 1)
    elif month == 7:
        subplot_index = (1, 0)
    elif month == 10:
        subplot_index = (1, 1)

#convert degrees to radians 
    deg2rad = 3.14159/180.0
    rad2deg = 180.0/3.14159
# Create a figure with two subplots
    fig, axs = plt.subplots(2, 2, figsize=(14, 9))
#define variables
    print('read: '+meshFileName)
    mesh = xarray.open_dataset(meshFileName)
    nCells = mesh.sizes['nCells']
    latCell = mesh.variables['latCell'] # in radians
    lonCell = mesh.variables['lonCell']
    xCell = mesh.variables['xCell'] # in meters
    yCell = mesh.variables['yCell']

    print('read: '+runDir+outputFileName)
    output = xarray.open_dataset(runDir+outputFileName)
    varName = 'timeMonthly_avg_iceAreaCell'
    var1 = output.variables[varName]

# Add var name to the center of the plot
    fig.text(0.5, 0.5, varName, ha='center', va='top', fontsize=16)

    NHindices = np.where(latCell>60*deg2rad)
# reduce the variable to 1D so we can use these NHindices
    var1D = var1[0,:]

# Plot Arctic data in the first subplott
    scatter1=axs[0,0].scatter(xCell[NHindices], yCell[NHindices], c=var1D[NHindices], cmap='bwr', s=0.4)
    axs[0,0].axis('off')

# Plot Antarctic data in the second subplot
    scatter2=axs[0,1].scatter(yCell[NHindices], xCell[NHindices], c=var1D[NHindices], cmap='bwr', s=0.4, zorder=0)
    axs[0,1].axis('off')

# Plot Arctic data in the first subplott
    scatter3=axs[1,0].scatter(xCell[NHindices], yCell[NHindices], c=var1D[NHindices], cmap='bwr', s=0.4)
    axs[1,0].axis('off')

# Plot Antarctic data in the second subplot
    scatter4=axs[1,1].scatter(yCell[NHindices], xCell[NHindices], c=var1D[NHindices], cmap='bwr', s=0.4)
    axs[1,1].axis('off')

#Set Lat Long Bounds
#NSA
    axs[0,0].set_xlim([68,71])
    axs[0,0].set_ylim([-165, -145])
#Russian Shelf
    axs[0,1].set_xlim([68,81])
    axs[0,1].set_ylim([20,180])
#Near Svalbard
    axs[1,0].set_xlim([74,80])
    axs[1,0].set_ylim([10,30])
#Total
   # axs[1,1].set_xlim([desired_min_lat, desired_max_lat])
    #axs[1,1].set_ylim([desired_min_lon, desired_max_lon])

# name graph
    month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    month_name = month_names[month - 1]
    axs[0,0].set_title(f'{month_name} Arctic Sea Ice')
    axs[0,1].set_title(f'{month_name}')
    axs[1,0].set_title(f'{month_name}')
    axs[1,1].set_title(f'{month_name}')

# Add color bars
    cbar1 = fig.colorbar(scatter1, ax=axs[0, 0], orientation='vertical')
    cbar1.set_label('Ice Area')
    cbar2 = fig.colorbar(scatter2, ax=axs[0, 1], orientation='vertical')
    cbar2.set_label('Ice Area')
    cbar3 = fig.colorbar(scatter3, ax=axs[1, 0], orientation='vertical')
    cbar3.set_label('Ice Area')
    cbar4 = fig.colorbar(scatter4, ax=axs[1, 1], orientation='vertical')
    cbar4.set_label('Ice Area')

plt.savefig('Lat_Long_seaice_Output.png')
