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
# get NO2 mask value and unit conversion factor using a sample file
sample_file = Dataset(TROPOMI_files[0])
sample_NO2  = sample_file.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column']
NO2_fill_value = sample_NO2._FillValue
NO2_unit_convert = sample_NO2.multiplication_factor_to_convert_to_molecules_percm2
#####################################################################################################################
# define an "object" to extract and hold all relevant variables from raw TROPOMI NO2 L2 swath files

# I'd suggest the following:
# Use def __init__ to define the arrays of interest (i.e., those arrays you will use in the code).
# Then include something like def read_tropomi as a separate function to read in the data (so moving what you have currently in __init__ to read_tropomi
# The next function can be process_tropomi where you apply quality flags, scale the data (unit conversion), and any other relevant processing.
# Then the next fuction can be write_to_txt_file to write the data to the text file to be read in for the oversampling routine.

# You'll cycle through all these functions in __main__. See https://github.com/eamarais/eam-group/blob/main/eloise/uptrop/tropomi_ut_no2.py for an example of this. The __main__ bit is near the
# end of the file.
#####################################################################################################################


# I updated my codes in two different ways.
# Here I only included four variables here to shorten the codes.
#####################################################################################################################
# Method 1 
# This method can succesfully extract the raw values, filter fill values and apply the quality flag. But it is still not in the way as you suggested. 
# I tried to read data in a function like "read_tropomi", but it seems that I have to provide things like "self.NO2","self.lat" in the initialization section.

class TROPOMI_NO2:
    def __init__(self,TROPOMI_file):
        
        # Open the file
        self.TROPOMI_data = Dataset(TROPOMI_file, "r", format="NETCDF4")   
        
        # Provide fill value and unit conversion value
        NO2_fill_value = sample_NO2._FillValue
        NO2_unit_convert = sample_NO2.multiplication_factor_to_convert_to_molecules_percm2

        # Extract relevant variables from raw TROPOMI NO2 files        
        self.NO2 = np.array(self.TROPOMI_data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:])
        self.lat = self.TROPOMI_data.groups['PRODUCT'].variables['latitude'][0][:][:]
        self.lon = self.TROPOMI_data.groups['PRODUCT'].variables['longitude'][0][:][:]
        self.Flag = self.TROPOMI_data.groups['PRODUCT'].variables['qa_value'][0][:][:]
        
        # Initialize the 1D arrays to hold the processed information
        self.NO2 = self.convert_to_1D_array(self.NO2)
        self.lat = self.convert_to_1D_array(self.lat)
        self.lon = self.convert_to_1D_array(self.lon)
        self.Flag = self.convert_to_1D_array(self.Flag)  
                      
        # Remove fill values
        self.NO2[self.NO2 == NO2_fill_value] = np.nan
        
        # Conver the unit
        self.NO2 = self.NO2*NO2_unit_convert/1e15
        
        # apply quality flag
        for i in range(len(self.Flag)):
            if (self.Flag[i] < 0.75):
                self.NO2[i] = np.nan
            
    def convert_to_1D_array(self,data):
        
        # convert "self.xxx" from 2D matrix to 1D array
        data = np.concatenate(data,axis=0)
        data = data.ravel()
        data = pd.Series(data)
        return data      
#####################################################################################################################
# Method 2
# Here I created two functions:
# The "convert_to_1D_array(self,data)" worked successfully
# The "process_TROPOMI(self)" did not work. But there was no error, and print statement indicated that those codes were read.

class TROPOMI_NO2:
    def __init__(self,TROPOMI_file):
        
        # Open the file
        self.TROPOMI_data = Dataset(TROPOMI_file, "r", format="NETCDF4")     
        
        # Provide fill value and unit conversion value
        NO2_fill_value = sample_NO2._FillValue
        NO2_unit_convert = sample_NO2.multiplication_factor_to_convert_to_molecules_percm2

        # Extract relevant variables from raw TROPOMI NO2 files        
        self.NO2 = np.array(self.TROPOMI_data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:])
        self.lat = self.TROPOMI_data.groups['PRODUCT'].variables['latitude'][0][:][:]
        self.lon = self.TROPOMI_data.groups['PRODUCT'].variables['longitude'][0][:][:]
        self.Flag = self.TROPOMI_data.groups['PRODUCT'].variables['qa_value'][0][:][:]

        # Initialize the 1D arrays to hold the information extracted from the 2D matrices
        self.NO2 = self.convert_to_1D_array(self.NO2)
        self.lat = self.convert_to_1D_array(self.lat)
        self.lon = self.convert_to_1D_array(self.lon)
        
    def convert_to_1D_array(self,data):
        # convert "self.xxx" from 2D matrix to 1D array
        data = np.concatenate(data,axis=0)
        data = data.ravel()
        data = pd.Series(data)
        return  data       

    def process_TROPOMI(self):
        print("The codes are read at least until 'process_TROPOMI'")
        # Remove fill values
        self.NO2[self.NO2 == NO2_fill_value] = np.nan
        
        # Conver the unit
        self.NO2 = self.NO2*NO2_unit_convert/1e15
        
# I did a quick test to check each step (outputs are shown following the "#")
test = TROPOMI_NO2(TROPOMI_files[0])

test.NO2.head()
# 0    9.969210e+36
# 1    9.969210e+36
# 2    9.969210e+36
# 3    9.969210e+36
# 4    9.969210e+36
# dtype: float32

test.Flag.head()
# 0    0.0
# 1    0.0
# 2    0.0
# 3    0.0
# 4    0.0
# dtype: float32
  
test.NO2.shape
# (1877400,)

test.process_TROPOMI()
# The codes are read at least until 'process_TROPOMI' 

# This means the first function "convert_to_1D_array" worked. But the "process_TROPOMI(self)" did not.
# So I did not continue with quality flag.
#####################################################################################################################

# About saving to txt, I tried the codes below.
# The error message is: "TypeError: write_to_txt_file() missing 2 required positional arguments: 'fId' and 'cnt'"
# I do not understand how you provided "fid,cnt" in your IASI codes
# Also, how should I define the string?

def write_to_txt_file(self, fId, cnt): 
        
        amf=0.0
        
        # Loop over data, skip where "NO2 = np.nan":
        for w in range(len(self.NO2)):
            if (self.NO2[w] != np.nan):
              
                cnt += 1
                self.cnt += 1

                # Define string of data to print to file:
                tstr="{:8}".format(cnt)+("{:15.6f}"*13).format(self.lat[w],(* self.lat_bounds[w]),self.lon[w],(* self.lon_bounds[w]),amf)+("{:15.6E}"*2).format(self.NO2[w])

                # Print to file:
                fId.write(tstr)
                fId.write("\n")            
#####################################################################################################################   
# End
