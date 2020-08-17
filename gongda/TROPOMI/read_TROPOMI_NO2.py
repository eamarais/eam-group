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
# create a "class" of "objects" to process and store variables from raw TROPOMI files
# given a raw TROPOMI file, output a pandas data frame containing variables of interest only

class TROPOMI_NO2:
    def __init__(self,TROPOMI_file):
        
        # first load the file and get variables
        self.TROPOMI_data = Dataset(TROPOMI_file, "r", format="NETCDF4")   
            
    def extract_raw_TROPOMI(self):    
        '''Extract relevant variables from raw TROPOMI NO2 files '''       
        self.NO2 = np.array(self.TROPOMI_data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:])
        self.lat = self.TROPOMI_data.groups['PRODUCT'].variables['latitude'][0][:][:]
        self.lat_bounds = self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds'][0][:][:][:]
        self.lon = self.TROPOMI_data.groups['PRODUCT'].variables['longitude'][0][:][:]
        self.lon_bounds = self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds'][0][:][:][:]
        self.Flag = self.TROPOMI_data.groups['PRODUCT'].variables['qa_value'][0][:][:]
        self.Pre = self.TROPOMI_data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision'][0][:][:]
        self.wind_east  = self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['eastward_wind'][0][:][:]
        self.wind_north = self.TROPOMI_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['northward_wind'][0][:][:]

        # now process the data fields (of this "TROPOMI_NO2" object)            
    def process_raw_TROPOMI(self):
        '''processing steps including: removing the fill values, applying quality flags, converting units etc.'''
        self.NO2[self.NO2 == NO2_fill_value] = np.nan # remove Fill values
        for i in range(self.NO2.shape[0]):            # apply quality flag = 0.75
            for j in range(self.NO2.shape[1]):
                if (self.Flag[i,j] < 0.75):
                    self.NO2[i,j] = np.nan
        self.NO2 = self.NO2*NO2_unit_convert/1e15     # convert the unit

    def convert_to_1D(self,data):
        '''convert TROPOMI NO2 attribute ("self.xxx") from a 2D matrix to 1D, so that you can combine them later'''
        data = np.concatenate(data,axis=0)
        data = data.ravel()
        data = pd.Series(data)
        return  data                
   
    # but this pandas df is not neccesary, if you know what to output, you can just perform on the 1D self.xxx directly
    def save_as_TROPOMI_dataframe(self): 
        '''save the variables to a single pandas dataframe'''
        NO2 = self.convert_to_1D(self.NO2)
        lon = self.convert_to_1D(self.lon)
        lat = self.convert_to_1D(self.lat)
        Flag = self.convert_to_1D(self.Flag)
        Pre = self.convert_to_1D(self.Pre)
        self.TROPOMI_results =  pd.concat([lon,lat,NO2,Flag,Pre], axis=1)
        self.TROPOMI_results.columns = ['lon','lat','NO2','Flag','Precision'] 
        self.TROPOMI_results = self.TROPOMI_results.dropna(subset = ["NO2"], inplace=False) # remove the row if "NO2" column is "NaN"
 #####################################################################################################################
# test your class  objects using a raw TROPOMI file
test = TROPOMI_NO2(TROPOMI_files_201908[7])

# the later functions are relying on what has been provided in advances
# you must call these functions in orders so that they can work
test.extract_raw_TROPOMI() 
test.process_raw_TROPOMI()
test.save_as_TROPOMI_dataframe() 

# show the final outputs after running all functions needed
test.TROPOMI_results
#####################################################################################################################
# now wrap your functions so that they can be mapped in multiprocessing
def read_TROPOMI(test):
    '''given a raw TROPOMI file, convert it to a "TROPOMI_NO2" object and get the results'''
    test = TROPOMI_NO2(test)
    test.extract_raw_TROPOMI() 
    test.process_raw_TROPOMI()
    test.save_as_TROPOMI_dataframe() 
    return test.TROPOMI_results
    
# set up multiprocessing using 16 CPUs
from multiprocessing import  Pool
pool = Pool(processes=16)

# parallel processing for files within a single month: TROPOMI_NO2_201908 = pool.map(read_TROPOMI,TROPOMI_files_201908)
# loop over all months
for i in range(len(yyyymm)):
    read_files_str = "TROPOMI_NO2_"+yyyymm[i]+" =pool.map(read_TROPOMI,TROPOMI_files_"+yyyymm[i]+")"
    exec(read_files_str)

# End
#####################################################################################################################
