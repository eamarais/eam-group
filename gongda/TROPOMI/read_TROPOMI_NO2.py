#####################################################################################################################
#####################################################################################################################
# Read TROPOMI NO2 L2 product and prepare for oversampling algorithms
# This script allows the user to select sampling periods, domain and quality flag without editing the codes (unless you need a custom domain).
# Add a line like this when submitting Python job:
# "python TROPOMI_NO2_prepare_oversampling.py --qa_flag 0.5 --start_date 20190801 --end_date 20190802 --domain CH"

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
import argparse
#####################################################################################################################
# provide the NO2 VCD fill value and unit conversion factor obtained from a sample file
# values are the same for Precision
NO2_fill_value = 9.96921e+36
NO2_unit_conversion = 6.02214e+19

# initialize some parameters for later use, their values will change according to your arguments
start_date = int()  # first sampling date
end_date = int()    # last sampling date
timewindow1 = str() # initialize the string for timewindow
timewindow2 = str() # initialize the string for timewindow

qa_threshold = int(0)  # quality flag threshold 
lat_min = int(0)    # lat_min for the requested domain
lat_max = int(0)    # lat_max for the requested domain
lon_min = int(0)    # lon_min for the requested domain
lon_max = int(0)    # lon_max for the requested domain
#####################################################################################################################
# create a function to select files between two dates

def select_TROPOMI_files_between(start_date,end_date):
    '''Select files within the start date ('yyyymmdd') and end date ('yyyymmdd')'''
    # list all the dates between the start and the end
    from datetime import date, timedelta
    start_date = date(int(start_date[0:4]),int(start_date[4:6]),int(start_date[6:8]))
    end_date   = date(int(end_date[0:4]),int(end_date[4:6]),int(end_date[6:8]))
    delta = end_date - start_date 
    sampling_dates = []
    for i in range(delta.days + 1):
        sampling_dates.append((start_date + timedelta(days=i)).strftime('%Y%m%d'))
    
    # check the sampling dates
    print("#############################################","All sampling dates:",*sampling_dates,sep="\n")
    
    # get files on each day
    TROPOMI_files_each_day = []
    for i in range(len(sampling_dates)):
        TROPOMI_files_each_day.append(sorted(glob.glob('S5P_OFFL_L2__NO2____'+sampling_dates[i]+'*T'+timewindow1+'*_*_*_*_*_*.nc') + \
                                             glob.glob('S5P_OFFL_L2__NO2____'+sampling_dates[i]+'*T'+timewindow2+'*_*_*_*_*_*.nc')))
    
    # combine the files on each day into a single list
    import itertools
    
    # check the selected files
    print("Number of files:",len((list(itertools.chain.from_iterable(TROPOMI_files_each_day)))))
    print("Files processed:",*list(itertools.chain.from_iterable(TROPOMI_files_each_day)),sep = "\n") 
    
    # now return the selected files
    return list(itertools.chain.from_iterable(TROPOMI_files_each_day))
#####################################################################################################################
# create a "class" of "object" to extract, process and store data fields from the selected TROPOMI files

class TROPOMI_NO2:
    def __init__(self,TROPOMI_files_list):
        '''Process a list of raw TROPOMI NO2 files. This is a flexible selection of files to accommodate each specific task. 
           It can be a list of files between any dates, on certain dates, or grouped by weekdays/weekends.'''
        
        # open the raw files from this list 
        self.TROPOMI_data = [Dataset(file, "r", format="NETCDF4") for file in  TROPOMI_files_list]
            
    def extract_raw_TROPOMI(self):    
        '''Extract relevant variables from raw TROPOMI NO2 files '''       
        self.NO2 = [np.array(data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0][:][:]) for data in self.TROPOMI_data]
        self.lat = [np.array(data.groups['PRODUCT'].variables['latitude'][0][:][:]) for data in self.TROPOMI_data]
        self.lat_bounds = [np.array(data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['latitude_bounds'][0][:][:][:]) for data in self.TROPOMI_data]
        self.lat_cor1 = [data[:,:,0] for data in self.lat_bounds] 
        self.lat_cor2 = [data[:,:,1] for data in self.lat_bounds] 
        self.lat_cor3 = [data[:,:,2] for data in self.lat_bounds] 
        self.lat_cor4 = [data[:,:,3] for data in self.lat_bounds] 
        self.lon = [np.array(data.groups['PRODUCT'].variables['longitude'][0][:][:]) for data in self.TROPOMI_data]
        self.lon_bounds = [np.array(data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].variables['longitude_bounds'][0][:][:][:]) for data in self.TROPOMI_data]
        self.lon_cor1 = [data[:,:,0] for data in self.lon_bounds]
        self.lon_cor2 = [data[:,:,1] for data in self.lon_bounds]
        self.lon_cor3 = [data[:,:,2] for data in self.lon_bounds]
        self.lon_cor4 = [data[:,:,3] for data in self.lon_bounds]       
        self.flag = [np.array(data.groups['PRODUCT'].variables['qa_value'][0][:][:]) for data in self.TROPOMI_data]
        self.pre  = [data.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column_precision'][0][:][:] for data in self.TROPOMI_data]
        self.sza  = [data.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['solar_zenith_angle'][0][:][:] for data in self.TROPOMI_data]
        self.cld  = [data.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['cloud_fraction_crb_nitrogendioxide_window'][0][:][:] for data in self.TROPOMI_data]
            
        # now process the data fields which have been stored by the "TROPOMI_NO2" object          
    def process_raw_TROPOMI(self):
        '''For data from each file, do some processing: remove the fill values, apply quality flags, convert units etc.'''
        for i in range(len(self.TROPOMI_data)):
            self.NO2[i] = np.where(self.flag[i] >= 0.5,self.NO2[i],np.nan) # np.where(condition,x,y): if condition == True, yield x, otherwise y
            self.NO2[i] = self.NO2[i]*NO2_unit_conversion  # convert unit to "molecules_percm2"
            self.pre[i] = self.pre[i]*NO2_unit_conversion   
            self.NO2[i] = self.NO2[i]/1e15 # convert to "1e15 molecules_percm2"
            self.pre[i] = self.pre[i]/1e15
        
    def convert_2D_to_1D(self):
        '''For data from each file, convert the variables to 1D, so you can skip the unwanted observations later'''
        self.NO2 = [data.ravel() for data in self.NO2]
        self.lat = [data.ravel() for data in self.lat]
        self.lat_cor1 = [data.ravel() for data in self.lat_cor1]
        self.lat_cor2 = [data.ravel() for data in self.lat_cor2]
        self.lat_cor3 = [data.ravel() for data in self.lat_cor3]
        self.lat_cor4 = [data.ravel() for data in self.lat_cor4]
        self.lon = [data.ravel() for data in self.lon]
        self.lon_cor1 = [data.ravel() for data in self.lon_cor1]
        self.lon_cor2 = [data.ravel() for data in self.lon_cor2]
        self.lon_cor3 = [data.ravel() for data in self.lon_cor3]
        self.lon_cor4 = [data.ravel() for data in self.lon_cor4]
        self.pre = [data.ravel() for data in self.pre]
        self.sza = [data.ravel() for data in self.sza]
        self.cld = [data.ravel() for data in self.cld]
        
    def create_output_file(self):      
        """create a single outputfile to store the results from all selected files""" 
        # use a unique suffix to dinguish the output file from each task
        self.suffix = args.domain+'_NO2_'+'_'+args.start_date+'_'+args.end_date           
        # set output direcotry and name the output file
        self.outfile = '/rds/projects/s/shiz-shi-aphh/TROPOMI_NO2_2_Oversampling/Input/TROPOMI_oversampling_input_'+self.suffix
        self.fId = open(self.outfile, "w+") 
        # initialize line counts in the file
        self.line_count = 0
        
    def write_to_output_file(self): 
        """combine variables into a string and write strings to a single txt file"""
        amf=0.0
        # Loop over data, skip invalid NO2 from each data:
        for i in range(len(self.NO2)):
            for j in range(self.NO2[i].shape[0]):   
                if ((math.isnan(self.NO2[i][j]) == False) &
                    (self.lat[i][j] >= lat_min) &
                    (self.lat[i][j] <= lat_max) &
                    (self.lon[i][j] >= lon_min) &
                    (self.lon[i][j] <= lon_max)):
                    tstr="{:8}".format(self.line_count)+("{:15.6f}"*13).format(self.lat[i][j],self.lat_cor1[i][j],self.lat_cor2[i][j],self.lat_cor3[i][j],self.lat_cor4[i][j],self.lon[i][j],self.lon_cor1[i][j],self.lon_cor2[i][j],self.lon_cor3[i][j],self.lon_cor4[i][j],self.sza[i][j],self.cld[i][j],amf)+("{:15.6E}"*2).format(self.NO2[i][j],self.pre[i][j])
                    self.fId.write(tstr) # write tstr to the file 
                    self.fId.write("\n") # progresses to the next line
                    self.line_count += 1 # increment the line number
        self.fId.close # close the txt file at end of loop
#####################################################################################################################
# now wrap your processing steps so you can loop over all files together
def prepare_TROPOMI_for_oversampling(raw_TROPOMI_file_list):
    '''given a list of raw TROPOMI files, convert the list to a "TROPOMI_NO2" object,
       and save results from all files to a single text file.
    '''
    test = TROPOMI_NO2(raw_TROPOMI_file_list)
    test.extract_raw_TROPOMI() 
    test.process_raw_TROPOMI()
    test.convert_2D_to_1D()
    test.create_output_file()
    test.write_to_output_file()
#####################################################################################################################
# pass the arguments to specify all requests for this job

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--qa_flag", default='0.75', help="Can be 0.75,0.5")
    parser.add_argument("--start_date", default='20190801', help="Can be any date string in 'yyyymmdd' format")
    parser.add_argument("--end_date", default='20190801', help="Can be any date string in 'yyyymmdd' format")
    parser.add_argument("--domain", default='AF', help="Can be AF,EU,CH,IND,NA or a custom domain")
    args = parser.parse_args()
    
    # pass arguments for the sampling period
    start_date = args.start_date
    end_date = args.end_date
    
    # pass arguments for quality flag
    qa_threshold = args.qa_flag

    # pass arguments for domain (data time window + domain lat/lon boundaries)
    if args.domain == "AF":         # African domain
        lat_min = -40
        lat_max =  48
        lon_min = -20
        lon_max =  50
        timewindow1 = '0[8-9]'
        timewindow2 = '1[0-4]'
    elif args.domain == "EU":        # EU nested GEOS-Chem domain
        lat_min =  32.75 - 0.25/2
        lat_max =  61.25 + 0.25/2
        lon_min = -15.0  - 0.3125/2
        lon_max =  40.0  + 0.3125/2
        timewindow1 = '0[8-9]'
        timewindow2 = '1[0-4]'
    elif args.domain == "CH":        # CH nested GEOS-Chem domain
        lat_min = 15.0 - 0.25/2
        lat_max = 55.0 + 0.25/2
        lon_min = 70.0 - 0.3152/2
        lon_max = 140.0 + 0.3125/2
        timewindow1 = '0[3-7]'
        timewindow2 = 'None'
    elif args.domain == "IND":       # Indian domain
        lat_min = 7
        lat_max = 35
        lon_min = 67
        lon_max = 90
        timewindow1 = '0[6-8]'
        timewindow2 = 'None'
    elif args.domain == "NA":        # NA nested GEOS-Chem domain
        lat_min =  9.75 - 0.25/2
        lat_max =  60.0 + 0.25/2
        lon_min = -130.0 - 0.3125/2
        lon_max = -60.0 + 0.3125/2
        timewindow1 = '1[5-9]'
        timewindow2 = '2[0-2]'
    else:
        print("Invalid domain: avaiable domains are AF,EU,CH,IND,NA. Or you can add your custom domain in the script.")
     
   
    # after passing all arguments, the job is already specified (period,quality flag and domain) 
    
    # now move to directory where the raw files are
    os.chdir("/rds/projects/s/shiz-shi-aphh/TROPOMI_NO2_0_RAW")
    
    # after passing all arguments, the job is already specified (period,quality flag and domain) 
    
    # now move to directory where the raw files are
    os.chdir("/rds/projects/s/shiz-shi-aphh/TROPOMI_NO2_0_RAW")
    
    # now start the job
    # import files during the sampling period
    TROPOMI_files_list = select_TROPOMI_files_between(start_date,end_date)
    
    # process each and save outputs from all files to a single txt file
    prepare_TROPOMI_for_oversampling(TROPOMI_files_list)
    
    # close the txt file that was created for this list
    
    # check if the job is what you expected
    # the print statements will appear upon job completion
    print("#############################################")
    print("Here we check the arguments used in this job:")
    print("Sampling period:",start_date,"-",end_date)
    print("Quality flag threshold:",qa_threshold)
    print("Selected domain:",args.domain)
    print("Domain boundaries (lat_min,lat_max,lon_min,lon_max):",lat_min,lat_max,lon_min,lon_max)
    print("TROPOMI files time window:",timewindow1,"+",timewindow2)

    # End of job
    print("The job is done.")
    print("#############################################")
