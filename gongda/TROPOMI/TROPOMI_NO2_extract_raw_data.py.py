#####################################################################################################################
#####################################################################################################################
# Extract TROPOMI NO2 data over Africa and Europe

# process data month by month, to avoid exceeding the RAM
# test the codes using sample data (2019-09) in JupyterNotebook
# submit complete Python jobs using ".py" files

import os
import glob
import numpy  as np
import pandas as pd
from netCDF4 import Dataset
import xarray as xr

#####################################################################################################################
# select TROPOMI files over Africa and Europe from the global dataset based on the time window 
# use 2019-09 as the example study period

os.chdir("/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2")
TROPOMI_files = sorted(glob.glob('S5P_OFFL_L2__NO2____201909*T0[8-9]*_*_*_*_*_*.nc') + 
                       glob.glob('S5P_OFFL_L2__NO2____201909*T1[0-4]*_*_*_*_*_*.nc')) 

print("Number of files:",len(TROPOMI_files))
print("First file:",TROPOMI_files[0])
print("Last  file:",TROPOMI_files[-1])
#####################################################################################################################
# explore relevant data fields using a sample file
# if you can not find the paths to some data fields, read data user guide or use NASA Panoply

sample = Dataset(TROPOMI_files[0])

sample_lat  = sample.groups['PRODUCT'].variables['latitude']
sample_lon  = sample.groups['PRODUCT'].variables['longitude']
sample_NO2  = sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column']
sample_Pre  = sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision']
sample_Flag = sample.groups['PRODUCT'].variables['qa_value']
sample_wind_east  = sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind']
sample_wind_north = sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind']

# get the mask value and unit conversion factor of TROPOMI NO2
NO2_fill_value = sample_NO2._FillValue
NO2_unit_convert = sample_NO2.multiplication_factor_to_convert_to_molecules_percm2

#####################################################################################################################
# process all files over Africa and Europe during the month
# only get values of variables (refer to the sample above for info)

TROPOMI_data = [Dataset(file, "r", format="NETCDF4") for file in TROPOMI_files]
lat  = [data.groups['PRODUCT'].variables['latitude'][0][:][:]  for data in TROPOMI_data]
lon  = [data.groups['PRODUCT'].variables['longitude'][0][:][:] for data in TROPOMI_data]
NO2  = [np.array(data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:]) for data in TROPOMI_data]
Pre  = [data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision'][0][:][:] for data in TROPOMI_data]
Flag = [data.groups['PRODUCT'].variables['qa_value'][0][:][:]  for data in TROPOMI_data]
wind_east  = [data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind'][0][:][:] for data in TROPOMI_data]
wind_north = [data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind'][0][:][:] for data in TROPOMI_data]

# remove fill_values for NO2 VCD
for i in range(len(NO2)):
    NO2[i][NO2[i]==NO2_fill_value] = np.nan
    
# convert NO2 VCD unit to molecules_percm2
NO2 = [data*NO2_unit_convert for data in NO2]

# scale the NO2 values (so the NO2 unit will be: 1e15 molecules_percm2)
NO2 = [data/1e15 for data in NO2]
#####################################################################################################################
# "date" info (optional)
# "date" for understanding time series, monthly or seasonal averages
# "weekday" for understanding urban background weekly variations

date  = [data.time_reference.split("T")[0] for data in TROPOMI_data] 
year  = [int(data.split('-')[0]) for data in date]
month = [int(data.split('-')[1]) for data in date]
day   = [int(data.split('-')[2]) for data in date]

# get "weekday" for each "date"
import datetime

weekday = []

for i in range(len(date)):
       weekday.append(datetime.date(year[i], month[i], day[i]).weekday())
        
# assign the same "date" and "weekday" for all observations on the same day
date_full = []
for i in range(len(NO2)):
    date_full.append([date[i]]*len(NO2[i].ravel()))

weekday_full = []
for i in range(len(NO2)):
    weekday_full.append([weekday[i]]*len(NO2[i].ravel()))
#####################################################################################################################
# now the input data are already extracted from the raw TROPOMI files

# check the mememory used by TROPOMI netCDF files
import sys
print(sys.getsizeof(TROPOMI_files))
print(sys.getsizeof(TROPOMI_data))

# close these to release some RAM
del (TROPOMI_files,TROPOMI_data) 
#####################################################################################################################
# now re-arrange the input data and combine into a pandas data frame
# do this for each variable seperately to avoid crashing the Jupyter Notebook kernel 

def prepare_TROPOMI(data):
    """get each input variable into the desired format before combining all variables into a single pandas dataframe"""
    data = np.concatenate(data,axis=0)
    data = data.ravel()
    data = pd.Series(data)
    return data
    
lat = prepare_TROPOMI(lat)
lon = prepare_TROPOMI(lon)
NO2 = prepare_TROPOMI(NO2)
Pre = prepare_TROPOMI(Pre)
Flag = prepare_TROPOMI(Flag)
wind_east = prepare_TROPOMI(wind_east)
wind_north = prepare_TROPOMI(wind_north)
date_full = prepare_TROPOMI(date_full)
weekday_full = prepare_TROPOMI(weekday_full)
#####################################################################################################################
# this dataframe contains all obervations within the time window (AF + EU) during this month

NO2_data_all = pd.concat([lon,lat,NO2,Pre,Flag,wind_east,wind_north,date_full,weekday_full], axis=1)
NO2_data_all.columns = ['lon','lat','NO2','Pre','Flag','wind_east','wind_north','date','weekday'] 

print('total observations:',len(NO2_data_all))
NO2_data_all.head()
#####################################################################################################################
# since the varaibles are already combined in the dataframe
# delete the input variables to release RAM
del (lon,lat,NO2,Pre,Flag,wind_east,wind_north,date_full,weekday_full)
#####################################################################################################################
# save the outputs to avoid reading raw TROPOMI files again
# this makes it easier to analyse data from multiple months together
# and further processing will always change (different regrid resolutions, cities, wind conditions or quality flags)
# so saving raw observations at this step can provide a 'restart' file
# and this 'restart' file can also be used by any other projects

# reduce the volumes of the data before saving outputs
# 1> zoom to the study region (AF or EU)
# 2> apply a loose data quality control (Flag >= 0.5)
#####################################################################################################################
# save TROPOMI NO2 in EU (be consistent with GEOS-Chem nested EU)

# zoom to EU domain
NO2_EU = NO2_data_all.loc[(NO2_data_all['lat'] >= (32.75 - 0.25/2)) &
                          (NO2_data_all['lat'] <= (61.25 + 0.25/2))  &
                          (NO2_data_all['lon'] >= (-15 - 0.3125/2)) &
                          (NO2_data_all['lon'] <= (40 + 0.3125/2))]

# apply a loose quality flag
NO2_EU = NO2_EU.loc[NO2_EU['Flag'] >= 0.5]

# reset the indices, otherwise later you can not slice
NO2_EU = NO2_EU.reset_index(drop=True) 

# check the selected data
print('total observations in EU:',len(NO2_EU))
NO2_EU.head()

# output EU data to netCDF

# convert pandas data frame to xarray data format
NO2_EU_201909 = xr.Dataset.from_dataframe(NO2_EU)

# edit attributes of variables
NO2_EU_201909.attrs={'Data summary':'TROPOMI NO2 in EU with quality flag >= 0.5 (2019-09)'}
NO2_EU_201909['NO2'].attrs = {'unit':'1e15 molecules_percm2'}
NO2_EU_201909['lon'].attrs = {'unit':'Degrees_east'}
NO2_EU_201909['lat'].attrs={'unit':'Degrees_north'}
NO2_EU_201909['wind_east'].attrs={'unit':'m s-1'}
NO2_EU_201909['wind_north'].attrs={'unit':'m s-1'}

# check the out structure
print(NO2_EU_201909)

# save the output
os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_extracted')
NO2_EU_201909.to_netcdf('NO2_EU_201909.nc')

# delte NO2_EU to release RAM
del NO2_EU_201909
#####################################################################################################################
# save TROPOMI NO2 in Africa (lat [38,-35] lon[-18,52])

# zoom to Africa domain
NO2_AF = NO2_data_all.loc[(NO2_data_all['lat'] >= -37) &
                          (NO2_data_all['lat'] <= 39)  &
                          (NO2_data_all['lon'] >= -22) &
                          (NO2_data_all['lon'] <= 53)]

NO2_AF = NO2_AF.loc[NO2_AF['Flag'] >= 0.5]
NO2_AF = NO2_AF.reset_index(drop=True) 
print("number of rows under all wind conditions:",len(NO2_AF))
print(NO2_AF.head())

# output AF data to netCDF
NO2_AF_201909 = xr.Dataset.from_dataframe(NO2_AF)

NO2_AF_201909.attrs={'Data summary':'TROPOMI NO2 in AF with quality flag >= 0.5 (2019-09)'}
NO2_AF_201909['NO2'].attrs = {'unit':'1e15 molecules_percm2'}
NO2_AF_201909['lon'].attrs = {'unit':'Degrees_east'}
NO2_AF_201909['lat'].attrs={'unit':'Degrees_north'}
NO2_AF_201909['wind_east'].attrs={'unit':'m s-1'}
NO2_AF_201909['wind_north'].attrs={'unit':'m s-1'}

print(NO2_AF_201909)

# save the output
os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_extracted')
NO2_AF_201909.to_netcdf('NO2_AF_201909.nc')

# delte NO2_AF to release RAM
del NO2_AF_201909
#####################################################################################################################
