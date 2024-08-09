#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals

"""
sum plots
Savannah Byron
plotting using cartopy
"""

############################## model files, run dirs
runDir = '/global/cfs/cdirs/m3958/njeffery/E3SMv3/dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm/run/'
meshFileName = '/global/cfs/projectdirs/e3sm/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc'
datafile = 'dustFromOceanOn2.icepack.GMPAS.IcoswISC30E3r5.pm.mpassi.hist.am.timeSeriesStatsMonthly.0020-01-01.nc'

import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import numpy as np

from cartopy import config
import cartopy.crs as ccrs

#define variables
deg2rad = 3.14159 / 180.0

fname_ice = os.path.join(runDir, datafile)
fname_mesh = os.path.join(meshFileName)

dataset_ice = netcdf_dataset(fname_ice)
dataset_mesh = netcdf_dataset(fname_mesh)

areaCell = dataset_ice.variables['timeMonthly_avg_iceAreaCell'][0, :]
latCell = dataset_mesh.variables['latCell'][:]
lonCell = dataset_mesh.variables['lonCell'][:]
latCelldeg = np.degrees(latCell)
lonCelldeg = np.degrees(lonCell)


ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(0, 90))
dot_size=10

sc = ax.scatter(lonCelldeg, latCelldeg,
                               c=areaCell, cmap='bwr',
                               s=dot_size, transform=ccrs.PlateCarree(),
                               )
ax.coastlines()

plt.savefig('Cartopy_Sample')
