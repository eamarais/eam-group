#####################################################################################################################
#####################################################################################################################
# Read TROPOMI NO2 L2 product and prepare for oversampling algorithms

# 1> select raw TROPOMI NO2 L2 files
# 2> extract data fileds of interest
# 3> perfrom processings (e.g. applying quality flag)
# 4> output strings for direct use in Fortran oversampling algorithms

import os
import glob
import numpy  as np
import pandas as pd
from netCDF4 import Dataset

import math
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
        '''convert the variables to 1D, so you can skip the unwanted observations'''
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
        
    def create_output_file(self):      
        """create an outputfile to store the results""" 
        # Get the date of observations
        self.date = self.TROPOMI_data.time_reference.split("T")[0]            
        # Define output direcotry and name outfile with date
        self.outfile = '/rds/projects/s/shiz-shi-aphh/TROPOMI_oversampling_input_'+self.date
        self.fId = open(self.outfile, "w+") 
        # Initialize line counts in the file
        self.line_count = 0
        
    def write_to_output_file(self): 
        """combine variables into a string and write strings to a single txt file"""
        amf=0.0
        # Loop over data, skip invalid NO2:
        for w in range(self.NO2.shape[0]):   
            if (math.isnan(self.NO2[w]) == False):
                tstr="{:8}".format(self.line_count)+("{:15.6f}"*13).format(self.lat[w],self.lat_cor1[w],self.lat_cor2[w],self.lat_cor3[w],self.lat_cor4[w],self.lon[w],self.lon_cor1[w],self.lon_cor2[w],self.lon_cor3[w],self.lon_cor4[w],self.sza[w],self.cld[w],amf)+("{:15.6E}"*2).format(self.NO2[w],self.pre[w])
                self.fId.write(tstr) # write tstr to the file 
                self.fId.write("\n") # progresses to the next line
                self.line_count += 1 # increment the line number
  

# test your class  objects using a raw TROPOMI file
test = TROPOMI_NO2(TROPOMI_files_201908[7])

# the later functions are relying on what has been provided in advances
# you must call these functions in orders so that they can work
test.extract_raw_TROPOMI() 
test.process_raw_TROPOMI()
test.convert_2D_to_1D()
test.create_output_file()

# now wrap your functions so you can loop over all files together
def prepare_TROPOMI_for_oversampling(raw_TROPOMI_file):
    '''given a raw TROPOMI file, convert it to a "TROPOMI_NO2" object and save results to text file'''
    test = TROPOMI_NO2(raw_TROPOMI_file)
    test.extract_raw_TROPOMI() 
    test.process_raw_TROPOMI()
    test.convert_2D_to_1D()
    test.create_output_file()
    test.write_to_txt_file()
    
       
# later I will edit the  codes to get TROPOMI files and if the function above is OK, I will start processeing multiple files.
test.write_to_txt_file()
#####################################################################################################################
        
