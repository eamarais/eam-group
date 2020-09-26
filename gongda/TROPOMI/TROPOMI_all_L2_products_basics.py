#####################################################################################################################
#####################################################################################################################
# Here we use files within 2019-12 to explore all TROPOMI level 2 offline products
# Homepage: http://www.tropomi.eu/data-products/level-2-products

# NO2 - nitrogendioxide_tropospheric_column
# HCHO - formaldehyde_tropospheric_vertical_column
# SO2 - sulfurdioxide_total_vertical_column
# CO - carbonmonoxide_total_column
# O3 - ozone_total_vertical_column
    
import os
import glob
import numpy as np
from netCDF4 import Dataset
#####################################################################################################################
# import all files for each prodcut during 2019-12
# read the data from a sample file from each product
# during the same sampling period, the number of files are almost the same for all species

# NO2
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/TROPOMI_sample_L2_products/NO2")
NO2_files = sorted(glob.glob('S5P_OFFL_L2__NO2*.nc'))
NO2_sample = Dataset(NO2_files[0], "r", format="NETCDF4")

print("#"*50,"NO2",sep='\n')
print("Number of files:",len(NO2_files))
print("First file:",NO2_files[0])
print("Last file: ",NO2_files[-1])

# HCHO
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/TROPOMI_sample_L2_products/HCHO")
HCHO_files = sorted(glob.glob('S5P_OFFL_L2__HCHO*.nc'))
HCHO_sample = Dataset(HCHO_files[0], "r", format="NETCDF4")

print("#"*50,"HCHO",sep='\n')
print("Number of files:",len(HCHO_files))
print("First file:",HCHO_files[0])
print("Last file: ",HCHO_files[-1])

# SO2
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/TROPOMI_sample_L2_products/SO2")
SO2_files = sorted(glob.glob('S5P_OFFL_L2__SO2*.nc'))
SO2_sample = Dataset(SO2_files[0], "r", format="NETCDF4")

print("#"*50,"SO2",sep='\n')
print("Number of files:",len(SO2_files))
print("First file:",SO2_files[0])
print("Last file: ",SO2_files[-1])

# CO
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/TROPOMI_sample_L2_products/CO")
CO_files = sorted(glob.glob('S5P_OFFL_L2__CO*.nc'))
CO_sample = Dataset(CO_files[0], "r", format="NETCDF4")

print("#"*50,"CO",sep='\n')
print("Number of files:",len(CO_files))
print("First file:",CO_files[0])
print("Last file: ",CO_files[-1])

# O3
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/TROPOMI_sample_L2_products/O3")
O3_files = sorted(glob.glob('S5P_OFFL_L2__O3*.nc'))
O3_sample = Dataset(O3_files[0], "r", format="NETCDF4")

print("#"*50,"O3",sep='\n')
print("Number of files:",len(O3_files))
print("First file:",O3_files[0])
print("Last file: ",O3_files[-1])
#####################################################################################################################
# explore data fields in each product
# subgroups are always listed at the bottom unless there is no more
# can also use data user guide or NASA Panoply to find the paths to data fields

# 1> "root group" only describes the data
# 2> variables are stored within "PRODUCT" and its subgroups 
# 3> "METADATA" group stores retrieval details

# NO2
print(NO2_sample,"#"*100,sep='\n')
print(NO2_sample.groups['PRODUCT'],"#"*100,sep='\n')                                   
print(NO2_sample.groups['PRODUCT']['SUPPORT_DATA'])                     
print(NO2_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'],"#"*100,sep='\n') 
print(NO2_sample.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'],"#"*100,sep='\n') 
print(NO2_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'],"#"*100,sep='\n')

# HCHO
print(HCHO_sample,"#"*100,sep='\n')
print(HCHO_sample.groups['PRODUCT'],"#"*100,sep='\n')                                   
print(HCHO_sample.groups['PRODUCT']['SUPPORT_DATA'])                     
print(HCHO_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'],"#"*100,sep='\n') 
print(HCHO_sample.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'],"#"*100,sep='\n') 
print(HCHO_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'],"#"*100,sep='\n')

# SO2
print(SO2_sample,"#"*100,sep='\n')
print(SO2_sample.groups['PRODUCT'],"#"*100,sep='\n')                                   
print(SO2_sample.groups['PRODUCT']['SUPPORT_DATA'])                     
print(SO2_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'],"#"*100,sep='\n') 
print(SO2_sample.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'],"#"*100,sep='\n') 
print(SO2_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'],"#"*100,sep='\n')

# CO
print(CO_sample,"#"*100,sep='\n')
print(CO_sample.groups['PRODUCT'],"#"*100,sep='\n')                                   
print(CO_sample.groups['PRODUCT']['SUPPORT_DATA'])                     
print(CO_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'],"#"*100,sep='\n') 
print(CO_sample.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'],"#"*100,sep='\n') 
print(CO_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'],"#"*100,sep='\n')

# O3
print(O3_sample,"#"*100,sep='\n')
print(O3_sample.groups['PRODUCT'],"#"*100,sep='\n')                                   
print(O3_sample.groups['PRODUCT']['SUPPORT_DATA'])                     
print(O3_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'],"#"*100,sep='\n') 
print(O3_sample.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'],"#"*100,sep='\n') 
print(O3_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'],"#"*100,sep='\n')
#####################################################################################################################
# ALl TROPOMI L2 offline prodcuts are structured in almost the same way
# ['PRODUCT'] group: VCD (with fill value and unit conversion), precision, quality flag, lat and lon
# ['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'] group: lat_bounds and lon_bounds
# ['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'] group: wind fields (NO2 and CO products only)

# NO2 (also understand dimensions for variables (e.g. longitude_bounds(time, scanline, ground_pixel, corner))
NO2_VCD = NO2_sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column']
NO2_VCD_fill_value = NO2_VCD._FillValue
NO2_VCD_unit_convert = NO2_VCD.multiplication_factor_to_convert_to_molecules_percm2
NO2_VCD_precision = NO2_sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision']
NO2_VCD_qa_flag = NO2_sample.groups['PRODUCT'].variables['qa_value']
NO2_lat = NO2_sample.groups['PRODUCT'].variables['latitude']
NO2_lon = NO2_sample.groups['PRODUCT'].variables['longitude']
NO2_lat_bounds = NO2_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds']
NO2_lon_bounds = NO2_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds']
NO2_wind_east = NO2_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind']
NO2_wind_north = NO2_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind']
NO2_Y_axis = NO2_sample.groups['PRODUCT'].variables['scanline'] # along-track dimension index
NO2_X_axis = NO2_sample.groups['PRODUCT'].variables['ground_pixel'] # across-track dimension index

# HCHO
HCHO_VCD = HCHO_sample.groups['PRODUCT'].variables['formaldehyde_tropospheric_vertical_column']
HCHO_VCD_fill_value = HCHO_VCD._FillValue
HCHO_VCD_unit_convert = HCHO_VCD.multiplication_factor_to_convert_to_molecules_percm2
HCHO_VCD_precision = HCHO_sample.groups['PRODUCT'].variables['formaldehyde_tropospheric_vertical_column_precision']
HCHO_VCD_qa_flag = HCHO_sample.groups['PRODUCT'].variables['qa_value']
HCHO_lat = HCHO_sample.groups['PRODUCT'].variables['latitude']
HCHO_lon = HCHO_sample.groups['PRODUCT'].variables['longitude']
HCHO_lat_bounds = HCHO_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds']
HCHO_lon_bounds = HCHO_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds']
HCHO_Y_axis = HCHO_sample.groups['PRODUCT'].variables['scanline'] 
HCHO_X_axis = HCHO_sample.groups['PRODUCT'].variables['ground_pixel'] 

# SO2
SO2_VCD = SO2_sample.groups['PRODUCT'].variables['sulfurdioxide_total_vertical_column']
SO2_VCD_fill_value = SO2_VCD._FillValue
SO2_VCD_unit_convert = SO2_VCD.multiplication_factor_to_convert_to_molecules_percm2
SO2_VCD_precision = SO2_sample.groups['PRODUCT'].variables['sulfurdioxide_total_vertical_column_precision']
SO2_VCD_qa_flag = SO2_sample.groups['PRODUCT'].variables['qa_value']
SO2_lat = SO2_sample.groups['PRODUCT'].variables['latitude']
SO2_lon = SO2_sample.groups['PRODUCT'].variables['longitude']
SO2_lat_bounds = SO2_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds']
SO2_lon_bounds = SO2_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds']
SO2_Y_axis = SO2_sample.groups['PRODUCT'].variables['scanline'] 
SO2_X_axis = SO2_sample.groups['PRODUCT'].variables['ground_pixel'] 

# CO
CO_VCD = CO_sample.groups['PRODUCT'].variables['carbonmonoxide_total_column']
CO_VCD_fill_value = CO_VCD._FillValue
CO_VCD_unit_convert = CO_VCD.multiplication_factor_to_convert_to_molecules_percm2
CO_VCD_precision = CO_sample.groups['PRODUCT'].variables['carbonmonoxide_total_column_precision']
CO_VCD_qa_flag = CO_sample.groups['PRODUCT'].variables['qa_value']
CO_lat = CO_sample.groups['PRODUCT'].variables['latitude']
CO_lon = CO_sample.groups['PRODUCT'].variables['longitude']
CO_lat_bounds = CO_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds']
CO_lon_bounds = CO_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds']
CO_wind_east = CO_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind']
CO_wind_north = CO_sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind']
CO_Y_axis = CO_sample.groups['PRODUCT'].variables['scanline'] 
CO_X_axis = CO_sample.groups['PRODUCT'].variables['ground_pixel'] 

# O3
O3_VCD = O3_sample.groups['PRODUCT'].variables['ozone_total_vertical_column']
O3_VCD_fill_value = O3_VCD._FillValue
O3_VCD_unit_convert = O3_VCD.multiplication_factor_to_convert_to_molecules_percm2
O3_VCD_precision = O3_sample.groups['PRODUCT'].variables['ozone_total_vertical_column_precision']
O3_VCD_qa_flag = O3_sample.groups['PRODUCT'].variables['qa_value']
O3_lat = O3_sample.groups['PRODUCT'].variables['latitude']
O3_lon = O3_sample.groups['PRODUCT'].variables['longitude']
O3_lat_bounds = O3_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds']
O3_lon_bounds = O3_sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds']
O3_Y_axis = O3_sample.groups['PRODUCT'].variables['scanline'] 
O3_X_axis = O3_sample.groups['PRODUCT'].variables['ground_pixel'] 

# the "shape" attribute can help you understand the variables and sometimes debug
print(NO2_lon.shape)
print(NO2_lon_bounds.shape)
print(NO2_lon[0][:][1][0:2])
print(NO2_lon_bounds[0][:][1][0:2])
#####################################################################################################################
# Now extract numeric values only from each variable
# fill values for VCD and precesion are masked by default, use "np.array()" will remove the mask

# take NO2 as the example
NO2_lat_value  = np.array(NO2_sample.groups['PRODUCT'].variables['latitude'][0][:][:])
NO2_lon_value  = np.array(NO2_sample.groups['PRODUCT'].variables['longitude'][0][:][:])
NO2_Flag_value = np.array(NO2_sample.groups['PRODUCT'].variables['qa_value'][0][:][:])
NO2_Pre_value  = np.array(NO2_sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision'][0][:][:])
NO2_value  = np.array(NO2_sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:]) 

# replace fill values with "NA"
NO2_value[NO2_value==NO2_VCD_fill_value]=np.nan

# the NO2 value and precision value should always have the same unit
NO2_value = NO2_value*NO2_VCD_unit_convert
NO2_value = NO2_value/1e15
NO2_Pre_value = NO2_Pre_value*NO2_VCD_unit_convert
NO2_Pre_value = NO2_Pre_value/1e15

# check NO2 VCD data range
print("min NO2 VCD:",np.nanmin(NO2_value))
print("max NO VCD:",np.nanmax(NO2_value))
print("mean NO VCD:",np.nanmean(NO2_value))
print("median NO VCD:",np.nanmedian(NO2_value))
#####################################################################################################################
# get date of this observation and convert to weekday (1-7 Monday-Sunday)
from datetime import date
sample_date = NO2_sample.time_reference.split("T")[0]
sample_date =  date.fromisoformat(sample_date)
sample_weekday = sample_date.weekday() +1
print(sample_date)
print(sample_weekday)

# get year,month,day separately
sample_date = NO2_sample.time_reference.split("T")[0]
year,month,day = (int(x) for x in sample_date.split('-'))
print(year,month,day)

# you can also the star time of each swath to distinguish each swath (but currently not really used)
print(NO2_sample.time_coverage_start)
#####################################################################################################################
# plot a swath using Basemap
# use NO2 as the example

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
