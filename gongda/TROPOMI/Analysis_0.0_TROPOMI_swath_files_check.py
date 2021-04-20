#####################################################################################################################
#####################################################################################################################
# There are 14-15 swaths per day from TROPOMI, and the files are huge in size
# before downloading or processing the data, first find which files are actually over your study domains

import os
import glob
import numpy as np
import pandas as pd

#####################################################################################################################
# Based on previous analysis for the shipping project
# we have defined three study domains (lon_min,lat_min,lon_max,lat_max)

# EU + Mid East: -20,5,60,65
# Southeast Asia: 60,-10,120,20
# America: -130,-10,-30,50

#####################################################################################################################
# Explore TROPOMI NO2 and TROPOMI HCHO products on NASA GES DISC homepage (https://disc.gsfc.nasa.gov/datasets?page=1)

# TROPOMI NO2
# Product name: (S5P_L2__NO2___)     Resolution: (7km x 3.5km)   Version: 1  Period: (From 2018-04-30 to 2019-08-07)
# Product name: (S5P_L2__NO2____HiR) Resolution: (5.5km x 3.5km) Version: 1  Period: (From 2019-08-06 till now)

'''                           Sentinel-5P TROPOMI V1.4 Nitrogen Dioxide Improvement
The TROPOMI NO2 product has been upgraded to version 1.4 beginning from the orbit 16213 dated on November 29th, 2020. 
The tropospheric vertical column amount may increase by up to 50% due to the overall decrease of the auxiliary cloud pressure input. 
Please be aware of the version 1.4 improvement compared to prior versions in time series analyses!
The TROPOMI NO2 product continuity and consistency will be improved when the fully reprocessed version 2 product becomes available, which is expected to be in fall 2021.
                             (Apart from this, check the differences between other 1.x versions)
'''

# TROPOMI HCHO
# Product name: (S5P_L2__HCHO__1)      Resolution: (7km x 3.5km)    Version: 1  Period: (From 2018-05-14 to 2019-08-07)
# Product name: (S5P_L2__HCHO___HiR_1) Resolution: (5.5km x 3.5km)  Version: 1  Period: (From 2019-08-06 to 2020-07-14)
# Product name: (S5P_L2__HCHO___HiR_2) Resolution: (5.5km x 3.5km)  Version: 2  Period: (From 2020-07-13 till now)

#####################################################################################################################
# Then use the NASA data subset tool to get the download links of all the TROPOMI files over the study domains since 2018
# read the links into Python

os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/Domains_swaths/')

EU_swaths_files = sorted(glob.glob("EU*"))
AS_swaths_files = sorted(glob.glob("AS*"))
US_swaths_files = sorted(glob.glob("US*"))

print("EU swath files links:",*EU_swaths_files,"#"*30, sep="\n")
print("AS swath files links:",*AS_swaths_files,"#"*30, sep="\n")
print("US swath files links:",*US_swaths_files,"#"*30, sep="\n")

#####################################################################################################################
# get swath starting hours over EU

# read the EU links
EU_links = [pd.read_csv(file,sep="/",skiprows=3) for file in EU_swaths_files]

# only keep the columns for product name and the swath file name
EU_links = [df.iloc[:,[5,8]] for df in EU_links]

# add a new column to hold information of the starting hour of the swath
for i in range(len(EU_links)):
    EU_links[i]['swath_start_hour'] = np.nan
    
# extract the hour info from each download link
for i in range(len(EU_links)):
    for j in range(EU_links[i].shape[0]):
        EU_links[i]['swath_start_hour'][j] = str(EU_links[i].iloc[j,1].split("T")[1])[0:2]
        
# get unique swath starting hours over EU
EU_swath_hours = pd.concat([EU_links[0]['swath_start_hour'],
                            EU_links[1]['swath_start_hour'],
                            EU_links[2]['swath_start_hour'],
                            EU_links[3]['swath_start_hour'],
                            EU_links[4]['swath_start_hour'],])

EU_swath_hours = EU_swath_hours.reset_index(drop=True)
EU_swath_hours = sorted(EU_swath_hours.unique())
EU_swath_hours

#####################################################################################################################
# get swath starting hours over AS

# read the AS links
AS_links = [pd.read_csv(file,sep="/",skiprows=3) for file in AS_swaths_files]

# only keep the columns for product name and the swath file name
AS_links = [df.iloc[:,[5,8]] for df in AS_links]

# add a new column to hold information of the starting hour of the swath
for i in range(len(AS_links)):
    AS_links[i]['swath_start_hour'] = np.nan
    
# extract the hour info from each download link
for i in range(len(AS_links)):
    for j in range(AS_links[i].shape[0]):
        AS_links[i]['swath_start_hour'][j] = str(AS_links[i].iloc[j,1].split("T")[1])[0:2]
        
# get unique swath starting hours over AS
AS_swath_hours = pd.concat([AS_links[0]['swath_start_hour'],
                            AS_links[1]['swath_start_hour'],
                            AS_links[2]['swath_start_hour'],
                            AS_links[3]['swath_start_hour'],
                            AS_links[4]['swath_start_hour'],])

AS_swath_hours = AS_swath_hours.reset_index(drop=True)
AS_swath_hours = sorted(AS_swath_hours.unique())
AS_swath_hours

#####################################################################################################################
# get swath starting hours over US

# read the US links
US_links = [pd.read_csv(file,sep="/",skiprows=3) for file in US_swaths_files]

# only keep the columns for product name and the swath file name
US_links = [df.iloc[:,[5,8]] for df in US_links]

# add a new column to hold information of the starting hour of the swath
for i in range(len(US_links)):
    US_links[i]['swath_start_hour'] = np.nan
    
# extract the hour info from each download link
for i in range(len(US_links)):
    for j in range(US_links[i].shape[0]):
        US_links[i]['swath_start_hour'][j] = str(US_links[i].iloc[j,1].split("T")[1])[0:2]
        
# get unique swath starting hours over US
US_swath_hours = pd.concat([US_links[0]['swath_start_hour'],
                            US_links[1]['swath_start_hour'],
                            US_links[2]['swath_start_hour'],
                            US_links[3]['swath_start_hour'],
                            US_links[4]['swath_start_hour'],])

US_swath_hours = US_swath_hours.reset_index(drop=True)
US_swath_hours = sorted(US_swath_hours.unique())
US_swath_hours

#####################################################################################################################
# print out the results for these domains
print("EU TROPOMI swath files staring hours:", EU_swath_hours,sep="\n")
print("AS TROPOMI swath files staring hours:", AS_swath_hours,sep="\n")
print("US TROPOMI swath files staring hours:", US_swath_hours,sep="\n")

#####################################################################################################################
# test the swaths covering over all these domains
all_domains = EU_swath_hours + AS_swath_hours + US_swath_hours
all_domains = np.array(all_domains)
np.unique(all_domains)

#####################################################################################################################
# So for the shipping project, we still need to download all the global data
# But we can specify the hour windows when processing the files for each domain

# Some files may be missing while downloading
# Some files can be faulty, they trend to have a weried file size and normally they crash the data processing jobs
# Here we use the codes to test if there are massive missing files and identify any potential faulty files
# Then you can compare the your downloaded files with the url list from NASA

# provide the information for target TPOPOMI files directory
species = 'NO2'
year = '2019'
month = '03'
average_file_size = 325

'''
average_file_size on NASA webpage

NO2: 325 MB before 2019-08-06
     437 MB after 2019-08-07
HCHO: 545 MB before 2019-08-06
      697 MB after 2019-08-07
'''

# move to the target files directory
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/0_TROPOMI_'+str(species)+'_raw/'+str(year)+'/'+str(species)+'_'+str(year)+str(month)+'/')
print(os.getcwd())

# first get all file names
all_files = sorted(glob.glob("*.nc"))
print("Total number of files:", len(all_files))

# retrieve file sizes from file names to identify any potential faulty files (with over 20% difference in file size)
def convert_unit(size_in_bytes,output_unit):
        """ Convert the size from bytes to KB, MB or GB"""
        if output_unit == 'KB':
            return size_in_bytes/1024
        elif output_unit == 'MB':
            return size_in_bytes/(1024**2)
        elif output_unit == 'GB':
            return size_in_bytes/(1024**3)
        else:
            return size_in_bytes
        
def get_file_size(file_name,output_unit):
    """ Get file size in units like KB, MB or GB"""
    import os
    import numpy as np
    size = os.path.getsize(file_name)
    size = convert_unit(size,output_unit)
    size = np.round(size,2)
    return size      

# print out the any potenial faulty files
print("check files with unexpected size:")

for i in range(len(all_files)):
    if abs(get_file_size(all_files[i],'MB')/average_file_size - 1) >= 0.2:
            print(all_files[i])
            print(get_file_size(all_files[i],'MB'),'MB')
    if i == len(all_files)-1:
            print("all files checked")
            
# but think about how to save out the printed results into a file so you can keep records of what has been removed
# and if you can combine the file names, you can use Python codes to remove them directly
# for now, use the manual approach to record and remove via Linux

# after excluding the unwanted files
# extract the date information of each remaining swath file
# so we can check if there are massive missing files
# but note there is minor difference in NO2 and HCHO file names

# NO2

# extract the date information of each remaining swath file
all_files = sorted(glob.glob("*.nc"))
datetime = [file_name.split("_")[8] for file_name in all_files]
date = [datetime_string[0:8] for datetime_string in datetime]

# finally, check if there are massive missing files
print("Total number of days:", len(np.unique(date)))
print("Total number of files:", len(all_files))
print("Number of files on each day:")
from collections import Counter
Counter(date)

# HCHO

# extract the date information of each remaining swath file
all_files = sorted(glob.glob("*.nc"))
datetime = [file_name.split("_")[7] for file_name in all_files]
date = [datetime_string[0:8] for datetime_string in datetime]

# finally, check if there are massive missing files
print("Total number of days:", len(np.unique(date)))
print("Total number of files:", len(all_files))
print("Number of files on each day:")
from collections import Counter
Counter(date)