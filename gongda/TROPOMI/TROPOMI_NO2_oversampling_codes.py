#####################################################################################################################
#####################################################################################################################
# oversampling TROPOMI NO2 data over Africa and Europe

# Steps
# 1> extract data from raw TROPOMI NO2 files
# 2> save to text file in format to be read in for RegridPixel fortran routine

import os
import glob
import numpy  as np
import pandas as pd
from netCDF4 import Dataset

#####################################################################################################################
# select TROPOMI files over Africa and Europe from the global dataset based on the time window 
# use 2019-09 as the sample study period

os.chdir("/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_0_RAW")
TROPOMI_files = sorted(glob.glob('S5P_OFFL_L2__NO2____201909*T0[8-9]*_*_*_*_*_*.nc') + 
                       glob.glob('S5P_OFFL_L2__NO2____201909*T1[0-4]*_*_*_*_*_*.nc')) 

print("Number of files:",len(TROPOMI_files))
print("First file:",TROPOMI_files[0])
print("Last  file:",TROPOMI_files[-1])

#####################################################################################################################
# open a single file
sample = Dataset(TROPOMI_files[0])

# get NO2 VCD, mask value and unit conversion factor
sample_NO2  = sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column']
NO2_fill_value = sample_NO2._FillValue
NO2_unit_convert = sample_NO2.multiplication_factor_to_convert_to_molecules_percm2

# other relevant data fileds
# use data user guide or NASA Panoply to find the paths to data fields
# understand the dimensions of each variable
# for example: longitude_bounds(time, scanline, ground_pixel, corner)
sample_lat  = sample.groups['PRODUCT'].variables['latitude']
sample_lon  = sample.groups['PRODUCT'].variables['longitude']
sample_lon_bounds = sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds']
sample_lat_bounds = sample.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds']
sample_Pre  = sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision']
sample_Flag = sample.groups['PRODUCT'].variables['qa_value']
sample_wind_east  = sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind']
sample_wind_north = sample.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind']

# use examples to understand lon_bounds and lat_bounds
print(sample_lon.shape)
print(sample_lon_bounds.shape)
print(sample_lon[0][:][1][0:2])
print(sample_lon_bounds[0][:][1][0:2])

# get the date of observation and convert to weekday
from datetime import date

sample_date = sample.time_reference.split("T")[0]
sample_date =  date.fromisoformat(sample_date)
sample_weekday = sample_date.weekday() +1

print(sample_date)
print(sample_weekday)

#####################################################################################################################
# create an "object" to extract and hold all relevant variables from raw TROPOMI NO2 L2 swath files

# I'd suggest the following:
# Use def __init__ to define the arrays of interest (i.e., those arrays you will use in the code).
# Then include something like def read_tropomi as a separate function to read in the data (so moving what you have currently in __init__ to read_tropomi
# The next function can be process_tropomi where you apply quality flags, scale the data (unit conversion), and any other relevant processing.
# Then the next fuction can be write_to_txt_file to write the data to the text file to be read in for the oversampling routine.

# You'll cycle through all these functions in __main__. See https://github.com/eamarais/eam-group/blob/main/eloise/uptrop/tropomi_ut_no2.py for an example of this. The __main__ bit is near the
# end of the file.

class TROPOMI_NO2:
    def __init__(self,TROPOMI_file):
        
        # Read data:
        TROPOMI_data = Dataset(TROPOMI_file, "r", format="NETCDF4")
        
        # Extract NO2 VCD
        self.NO2 = np.array(TROPOMI_data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:])
        
        # Extract other relevant variables
        self.lat = TROPOMI_data.groups['PRODUCT'].variables['latitude'][0][:][:]
        self.lon = TROPOMI_data.groups['PRODUCT'].variables['longitude'][0][:][:]
        self.lat_bounds = TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds'][0][:][:][:]
        self.lon_bounds = TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds'][0][:][:][:]
        self.Pre = TROPOMI_data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision'][0][:][:]
        self.Flag = TROPOMI_data.groups['PRODUCT'].variables['qa_value'][0][:][:]
        self.wind_east  = TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind'][0][:][:]
        self.wind_north = TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind'][0][:][:]
        
        # Extract the date of observation and convert to weekday
        self.date = TROPOMI_data.time_reference.split("T")[0]
        self.date = date.fromisoformat(self.date)
        self.weekday = self.date.weekday() +1
        
        # Number of data ponits
        self.nvals = len(self.lon)*len(self.lon[0])
        
#####################################################################################################################
# now use the function to extract and hold all raw data
NO2_201909 = [TROPOMI_NO2(data) for data in TROPOMI_files]
        
#####################################################################################################################
# remove fill values and the data with qa_flag < 0.75

for i in range(len(NO2_201909)):
    for j in range(len(NO2_201909[i].NO2)):
        for k in range(len(NO2_201909[i].NO2[0])):
            if ((NO2_201909[i].NO2[j,k] == NO2_fill_value) or (NO2_201909[i].Flag[j,k] < 0.75)):
                NO2_201909[i].NO2[j,k] = np.nan
        
#####################################################################################################################
# problems to solve:
# 1> Lines 100-106 are very slow, even much slower than my previous approach (combine all variables into a single pandas dataframe, and then subset it)
# 2> Even if lines 100-106 are done, I still do not know how to subset the data
#    The codes may be like: subset_data = [data for data in NO2_201909 if data.NO2 != np.nan], but "data.NO2" is a data matrix, so the codes won't work
# 3> Till now, I have extracted all the information needed, but I do not know how to save to text file in format to be read in for RegridPixel fortran routine
#    I think the corresponding codes started from line 116 in your IASI script, but I do not understand how to output     
        
