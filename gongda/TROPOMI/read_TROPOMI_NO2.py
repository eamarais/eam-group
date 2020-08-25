#####################################################################################################################
#####################################################################################################################
# New codes to process TROPOMI NO2 data

import os
import glob
import numpy  as np
import pandas as pd
from netCDF4 import Dataset
#####################################################################################################################
# load raw TROPOMI files
os.chdir("/rds/projects/s/shiz-shi-aphh/TROPOMI_NO2_0_RAW")

# list all the years & months (yyyymm) during the sampling period
def iterate_months(start_ym, end_ym):
    for ym in range(int(start_ym), int(end_ym) + 1):
        if ym % 100 > 12 or ym % 100 == 0:
            continue
        yield str(ym)

yyyymm = list(iterate_months('201908', '202007'))
print(yyyymm)

# import files from each month separately
# select TROPOMI files over Africa and Europe from the global dataset based on the time window 
for i in range(len(yyyymm)):
    load_files_str = "TROPOMI_files_"+yyyymm[i]+" = sorted(glob.glob('S5P_OFFL_L2__NO2____'+yyyymm[i]+'*T0[8-9]*_*_*_*_*_*.nc') + \
                                                           glob.glob('S5P_OFFL_L2__NO2____'+yyyymm[i]+'*T1[0-4]*_*_*_*_*_*.nc'))"
    exec(load_files_str)
#####################################################################################################################
# get NO2 fill value and unit conversion value using a sample file
sample = Dataset(TROPOMI_files_201908[0])
sample_NO2  = sample.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column']
NO2_fill_value = sample_NO2._FillValue
NO2_unit_convert = sample_NO2.multiplication_factor_to_convert_to_molecules_percm2
#####################################################################################################################
#####################################################################################################################
# create a "class" of "objects" to extract, process and store data fields from raw TROPOMI files

class TROPOMI_NO2:
    def __init__(self,TROPOMI_file):
        
        # open the file and get variables
        self.TROPOMI_data = Dataset(TROPOMI_file, "r", format="NETCDF4")   
            
    def extract_raw_TROPOMI(self):    
        '''Extract relevant variables from raw TROPOMI NO2 files '''       
        self.NO2 = np.array(self.TROPOMI_data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:])
        self.lat = np.array(self.TROPOMI_data.groups['PRODUCT'].variables['latitude'][0][:][:])
        self.lat_bounds = np.array(self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds'][0][:][:][:])
        self.lat_cor1 = self.lat_bounds[:,:,0] # this method can read 4 corners, but what do these 4 corners mean? Is this order correct for the oversampling?
        self.lat_cor2 = self.lat_bounds[:,:,1]
        self.lat_cor3 = self.lat_bounds[:,:,2]
        self.lat_cor4 = self.lat_bounds[:,:,3]
        self.lon = np.array(self.TROPOMI_data.groups['PRODUCT'].variables['longitude'][0][:][:])
        self.lon_bounds = np.array(self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds'][0][:][:][:])
        self.lon_cor1 = self.lon_bounds[:,:,0]
        self.lon_cor2 = self.lon_bounds[:,:,1]
        self.lon_cor3 = self.lon_bounds[:,:,2]
        self.lon_cor4 = self.lon_bounds[:,:,3]       
        self.flag = self.TROPOMI_data.groups['PRODUCT'].variables['qa_value'][0][:][:]
        self.pre  = self.TROPOMI_data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision'][0][:][:]
        self.sza  = self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['solar_zenith_angle'][0][:][:]
        self.cld  = self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['cloud_fraction_crb_nitrogendioxide_window'][0][:][:]       
        self.wind_east  = self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind'][0][:][:]
        self.wind_north = self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind'][0][:][:]

        # now process the data fields which have been stored by the "TROPOMI_NO2" object          
    def process_raw_TROPOMI(self):
        '''some processing steps: removing the fill values, applying quality flags, converting units etc.'''
        self.NO2 = np.where(self.flag >= 0.75,self.NO2,np.nan) # np.where(condition,x,y): if condition == True, yield x, otherwise y
        self.NO2 = self.NO2*NO2_unit_conversion  # convert unit to "molecules_percm2"
        self.pre = self.pre*NO2_unit_conversion   
        self.NO2 = self.NO2/1e15 # convert to "1e15 molecules_percm2"
        self.pre = self.pre/1e15
        
    def convert_2D_to_1D(self):
        '''convert the variables to 1D, so you can skip the '''
        self.NO2 = self.NO2.ravel()
        self.lat = self.lat.ravel()
        self.lat_cor1 = self.lat_cor1.ravel()
        self.lat_cor2 = self.lat_cor2.ravel()
        self.lat_cor3 = self.lat_cor3.ravel()
        self.lat_cor4 = self.lat_cor4.ravel()
        self.lon = self.lon.ravel()
        self.lon_cor1 = self.lon_cor1.ravel()
        self.lon_cor2 = self.lon_cor2.ravel()
        self.lon_cor3 = self.lon_cor3.ravel()
        self.lon_cor4 = self.lon_cor4.ravel()
        self.pre = self.pre.ravel()
        self.sza = self.sza.ravel()
        self.cld = self.cld.ravel()
        
    def write_to_txt_file(self): 
        amf=0.0
        self.tstr = []
        # Loop over data, skip where "NO2 = np.nan":
        for w in range(len(self.lat)):   
            if (self.NO2[w] != np.nan):
                # Define string of data to print to file:
                # Eloise example: tstr=("{:15.6f}"*13).format(self.lat[w],(* self.lat_bounds[w,:]),self.lon[w],(* self.lon_bounds[w,:]),self.sza[w],self.cld[w],amf)+("{:15.6E}"*2).format(self.NO2[w],self.pre[w])
                self.tstr.append(("{:15.6f}"*13).format(self.lat[w],(* self.lat_cor1[w],self.lat_cor2[w],self.lat_cor3[w],self.lat_cor4[w]),self.lon[w],(* self.lat_cor1[w],self.lat_cor2[w],self.lat_cor3[w],self.lat_cor4[w]),self.sza[w],self.cld[w],amf)+("{:15.6E}"*2).format(self.NO2[w],self.pre[w]))
            else:
                self.tstr.append("NA")
#####################################################################################################################
# test your class  objects using a raw TROPOMI file
test = TROPOMI_NO2(TROPOMI_files_201908[7])

# the later functions are relying on what has been provided in advances
# you must call these functions in orders so that they can work
test.extract_raw_TROPOMI() 
test.process_raw_TROPOMI()
test.convert_2D_to_1D()
test.write_to_txt_file() # Error: 'numpy.float32' object is not iterable

# now wrap your functions so you can repeat all the functions together for all files
def read_TROPOMI(test):
    '''given a raw TROPOMI file, convert it to a "TROPOMI_NO2" object and get the results'''
    test = TROPOMI_NO2(test)
    test.extract_raw_TROPOMI() 
    test.process_raw_TROPOMI()
    test.convert_2D_to_1D()
    test.write_to_txt_file()
    return test

# use the main routine that can cycle through the function above in the right order
# To use the code from main, enter at the terminal: python TROPOMI_NO2_oversampling_codes.py.

if __name__ == "__main__":
    for i in range(len(yyyymm)):
    read_files_str = "TROPOMI_NO2_"+yyyymm[i]+" = read_TROPOMI(TROPOMI_files_"+yyyymm[i]+")"
    exec(read_files_str)

 # End
#####################################################################################################################   
