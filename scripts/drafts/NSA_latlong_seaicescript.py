#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
This is a code that Nicole requested to find the points for riverrunofflux for the russian shelf, NSA, and svalbard
Savannah Byron
June 2024
"""

############################## model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
outputFileName = 'dustFromOceanOn.icepack.GMPAS.IcoswISC30E3r5.pm.mpaso.hist.am.timeSeriesStatsMonthly.0012-01-01.nc'
#import needed directories/programs
import numpy as np
import xarray
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

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

#Define Lat long bounds in rads
#North Slope of Alaska (NSA)
lat_min_rad_NSA = 68 * deg2rad
lat_max_rad_NSA = 71 * deg2rad
lon_min_rad_NSA = (-165+360) * deg2rad
lon_max_rad_NSA = (-145+360) * deg2rad
#RussianShelf (RS)
lat_min_rad_RS = 68 * deg2rad
lat_max_rad_RS = 81 * deg2rad
lon_min_rad_RS = 20 * deg2rad
lon_max_rad_RS = 180 * deg2rad
#Near Svalbard (SVA)
lat_min_rad_SVA = 74 * deg2rad
lat_max_rad_SVA = 80 * deg2rad
lon_min_rad_SVA = 10 * deg2rad
lon_max_rad_SVA = 30 * deg2rad

print('read: '+runDir+outputFileName)
output = xarray.open_dataset(runDir+outputFileName)
varName = 'timeMonthly_avg_riverRunoffFlux'
var1 = output.variables[varName]

# Add var name to the center of the plot
fig.text(0.5, 0.5, varName, ha='center', va='top', fontsize=16)

#define indices
NHindices = np.where(latCell>60*deg2rad)
NHindicesNSA = np.where((latCell>lat_min_rad_NSA) & (latCell<lat_max_rad_NSA) & (lonCell > lon_min_rad_NSA) & (lonCell < lon_max_rad_NSA))
NHindicesRS = np.where((latCell>lat_min_rad_RS) & (latCell<lat_max_rad_RS) & (lonCell > lon_min_rad_RS) & (lonCell < lon_max_rad_RS))
NHindicesSVA = np.where((latCell>lat_min_rad_SVA) & (latCell<lat_max_rad_SVA) & (lonCell > lon_min_rad_SVA) & (lonCell < lon_max_rad_SVA))




# reduce the variable to 1D so we can use these NHindices
var1D = var1[0,:]

# Plot Arctic data in the first subplott
scatter1=axs[0,0].scatter(xCell[NHindicesNSA], yCell[NHindicesNSA], c=var1D[NHindicesNSA], cmap='bwr', s=25)#increase s size to see easier
axs[0,0].axis('off')

# Plot Antarctic data in the second subplot
scatter2=axs[0,1].scatter(xCell[NHindicesRS], yCell[NHindicesRS], c=var1D[NHindicesRS], cmap='bwr', s=5, zorder=0)
axs[0,1].axis('off')

# Plot Arctic data in the first subplott
scatter3=axs[1,0].scatter(xCell[NHindicesSVA], yCell[NHindicesSVA], c=var1D[NHindicesSVA], cmap='bwr', s=20)
axs[1,0].axis('off')

# Plot Antarctic data in the second subplot
scatter4=axs[1,1].scatter(xCell[NHindices], yCell[NHindices], c=var1D[NHindices], cmap='bwr', s=0.3)
axs[1,1].axis('off')
# name graph
#month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
#month_name = month_names[month - 1]
axs[0,0].set_title('Arctic Sea Ice')
# axs[0,1].set_title(f'{month_name}')
# axs[1,0].set_title(f'{month_name}')
# axs[1,1].set_title(f'{month_name}')

# Add color bars
#cbar1 = fig.colorbar(scatter1, ax=axs[0, 0], orientation='vertical')
#cbar1.set_label('Ice Area')
cbar2 = fig.colorbar(scatter2, ax=axs[0, 1], orientation='vertical')
cbar2.set_label('Ice Area')
cbar3 = fig.colorbar(scatter3, ax=axs[1, 0], orientation='vertical')
cbar3.set_label('Ice Area')
cbar4 = fig.colorbar(scatter4, ax=axs[1, 1], orientation='vertical')
cbar4.set_label('Ice Area')


# Create a folder to save plots
plot_dir = '/global/u1/s/sbyron/python/plots'
os.makedirs(plot_dir, exist_ok=True)  # Create the folder if it doesn't exist

# Load output data
print('read: ' + os.path.join(runDir, outputFileName))
output = xarray.open_dataset(os.path.join(runDir, outputFileName))

# Save the figure
plt.savefig(os.path.join(plot_dir, 'Lat_Long_seaice_Output2.png'))
