#####################################################################################################################
#####################################################################################################################
# Get familiar with TROPOMI NO2 L2 offline products using a sample file

# Data manual (2019) https://sentinel.esa.int/documents/247904/2474726/Sentinel-5P-Level-2-Product-User-Manual-Nitrogen-Dioxide
# NASA webinar (2019) https://www.youtube.com/watch?v=-yOInEUJTYM
# Example codes by NASA: https://appliedsciences.nasa.gov/join-mission/training/english/high-resolution-no2-monitoring-space-tropomi

import os
import glob
import numpy  as np
import pandas as pd
from netCDF4 import Dataset
#####################################################################################################################
# set working directory and import a sample TROPOMI NO2 file
os.chdir("/rds/projects/s/shiz-shi-aphh/TROPOMI_NO2_0_RAW")
test_TROPOMI_file = "S5P_OFFL_L2__NO2____20190801T021405_20190801T035535_09317_01_010302_20190807T033743.nc"

# open the TROPOMI file
test = Dataset(test_TROPOMI_file, "r", format="NETCDF4")

# "root group" only describes the data
print(test)

# variables are stored within "PRODUCT" and its subgroups. "METADATA" groups store some retrieval details
# subgroups are always listed at the bottom unless there is no more
# can also use data user guide or NASA Panoply to find the paths to data fields
print(test.groups['PRODUCT'])                                   
print(test.groups['PRODUCT']['SUPPORT_DATA'])                     
print(test.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']) 
print(test.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']) 
print(test.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'])
#####################################################################################################################
# read some relevant variables (reading both values and attributes)
lat = test.groups['PRODUCT'].variables['latitude']
lat_bounds = test.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds']
lon = test.groups['PRODUCT'].variables['longitude']
lon_bounds = test.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds']
Flag = test.groups['PRODUCT'].variables['qa_value']
Pre = test.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision']
wind_east = test.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind']
wind_north = test.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind']
NO2 = test.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column']

# fill value and unit conversion factor for NO2 VCD (and for "precision")
NO2_fill_value = NO2._FillValue
NO2_unit_convert = NO2.multiplication_factor_to_convert_to_molecules_percm2

# dimensions for variables (e.g. longitude_bounds(time, scanline, ground_pixel, corner))
Y_axis = test.groups['PRODUCT'].variables['scanline'] # along-track dimension index
X_axis = test.groups['PRODUCT'].variables['ground_pixel'] # across-track dimension index

# use examples to understand lon_bounds and lat_bounds
print(lon.shape)
print(lon_bounds.shape)
print(lon[0][:][1][0:2])
print(lon_bounds[0][:][1][0:2])
#####################################################################################################################
# extract numeric values only from variables
# fill values for NO2 VCD and precesion are masked by default, use "np.array()" will remove the mask
lat_value  = test.groups['PRODUCT'].variables['latitude'][0][:][:]
lon_value  = test.groups['PRODUCT'].variables['longitude'][0][:][:]
Flag_value = test.groups['PRODUCT'].variables['qa_value'][0][:][:]
Pre_value  = np.array(test.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision'][0][:][:])
NO2_value  = np.array(test.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:]) 

# replace fill values with "NA"
NO2_value[NO2_value==NO2_fill_value]=np.nan

# the NO2 value and precision value should always have the same unit
NO2_value = NO2_value*NO2_unit_convert
NO2_value = NO2_value/1e15
Pre_value = Pre_value*NO2_unit_convert
Pre_value = Pre_value/1e15

# check NO2 VCD data range
print("min NO2 VCD:",np.nanmin(NO2_value))
print("max NO VCD:",np.nanmax(NO2_value))
print("mean NO VCD:",np.nanmean(NO2_value))
print("median NO VCD:",np.nanmedian(NO2_value))
#####################################################################################################################
# get date of this observation and convert to weekday (1-7 Monday-Sunday)
# you can also the star time of each swath to distinguish each swath
from datetime import date
sample_date = test.time_reference.split("T")[0]
sample_date =  date.fromisoformat(sample_date)
sample_weekday = sample_date.weekday() +1

# get year,month,day separately
sample_date = test.time_reference.split("T")[0]
year,month,day = (int(x) for x in sample_date.split('-'))
#####################################################################################################################
# plot with Basemap
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

data  = np.ma.masked_array(NO2_value,np.isnan(NO2_value))
m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat = 90,llcrnrlon=-180, urcrnrlon = 180)
m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180, 180., 45.), labels=[0, 0, 0, 1])
vmin1=0.0
vmax1=0.008
m.pcolormesh(lon_value, lat_value, data, latlon=True, vmin=vmin1, vmax=np.nanmax(data)*vmax1,cmap='jet')
cb = m.colorbar()
cb.set_label('NO2 VCD')
plt.autoscale()
fig = plt.gcf()
plt.title('{0}\n{1}'.format("TROPOMI NO2", "example plot"))
plt.show()

# End
#####################################################################################################################
