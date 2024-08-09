#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
scatter plot looking at variables, in this case river runoff, some changes that i noticed were how to set the plot title in the center top, messing with color bars, getting the right mesh
Savannah Byron
June 2024, 
"""

############################## model files, run dirs
# running in /lustre/scratch5/mpeterse/runs/240126_wind_stress_test/ocean/global_ocean/QU240/WOA23/performance_test/forward
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
#outputFileName = 'dustFromOceanOn.icepack.GMPAS.IcoswISC30E3r5.pm.mpaso.hist.am.timeSeriesStatsMonthly.0012-01-01.nc'

import numpy as np
import xarray
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import subprocess

# Define the output directory for the frames
frames_dir = 'monthly_frames'
os.makedirs(frames_dir, exist_ok=True)

#Iterate through months
for month in range(1,12):
#update outputfile name with current month
  outputFileName = f'dustFromOceanOn.icepack.GMPAS.IcoswISC30E3r5.pm.mpaso.hist.am.timeSeriesStatsMonthly.0012-{"{:02d}".format(month)}-01.nc'

  deg2rad = 3.14159/180.0
  rad2deg = 180.0/3.14159
# Create a figure with two subplots
  fig, axs = plt.subplots(1, 2, figsize=(14, 6))
# fig = plt.figure(figsize=(20,12))

  print('read: '+meshFileName)
  mesh = xarray.open_dataset(meshFileName)
  nCells = mesh.sizes['nCells']
  latCell = mesh.variables['latCell'] # in radians
  lonCell = mesh.variables['lonCell']
  xCell = mesh.variables['xCell'] # in meters
  yCell = mesh.variables['yCell']

  print('read: '+runDir+outputFileName)
  output = xarray.open_dataset(runDir+outputFileName)
  varName = 'timeMonthly_avg_riverRunoffFlux'
  var = output.variables[varName]

# Split the variable name by underscores and take the last part
  var_title = varName.split('_')[-1]


  NHindices = np.where(latCell>60*deg2rad)
# reduce the variable to 1D so we can use these NHindices
  var1D = var[0,:]

#Add Title to center of plot
  fig.text(0.5, 0.85, var_title, ha='center', va='top', fontsize=16)


# Plot Arctic data in the first subplott
  scatter1=axs[0].scatter(xCell[NHindices], yCell[NHindices], c=var1D[NHindices], cmap='bwr', s=0.4)
  axs[0].set_title('Arctic')
  axs[0].axis('off')

  SHindices = np.where(latCell<-60*deg2rad)
# reduce the variable to 1D so we can use these NHindices
  var1DSH = var[0,:][SHindices]


# Plot Antarctic data in the second subplot
  scatter2=axs[1].scatter(yCell[SHindices], xCell[SHindices], c=var1D[SHindices], cmap='bwr', s=0.4)
  axs[1].set_title('Antarctic')

  plt.axis('off')

#set colorbars
  cbar1 = fig.colorbar(scatter1, ax=axs[0], orientation='vertical')
#cbar1.set_label(f'{var_title}')

  cbar2 = fig.colorbar(scatter2, ax=axs[1], orientation='vertical')
#cbar1.set_label(f'{var_title}')

 # Save the plot with the current month in the filename
  plt.savefig(os.path.join(frames_dir, f'{var_title}_12year_{month:02d}.png'))
  plt.close()  # Close the figure to release memory

 # plt.savefig(f'{var_title}_12year_Output.png')

# Use ffmpeg to create a movie from the frames
output_movie = 'monthly_movie.mp4'
subprocess.run(['ffmpeg', '-framerate', '1', '-i', f'{frames_dir}/{var_title}_12year_%02d.png', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', output_movie])

