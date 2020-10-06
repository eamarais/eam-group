##############################################################################
# Merge TROPOMI oversampling input data from multiple months/sampling periods
# This is to generate a single TROPOMI oversampling input file for a long sampling period which has massive files

# For example, following "TROPOMI_0_oversampling_preparation.py":
# Reading raw TROPOMI L2 files over Africa over year, and outputing the strings together will require over 800G RAM, which exceeds the capacity of the HPC
# So generate the TROPOMI oversampling input file from each month seperately, then merge all of them into a single one for the year

import os
import glob
import pandas as pd
##############################################################################
# move to the directory which holds the oversampling input data from each month
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/TROPOMI_oversampling_data/Oversampling_Input_Data')

# list the monthly input files of interest
TROPOMI_NO2_Africa_201808_201907_files = sorted(glob.glob(...))

# read all the monthly input files
TROPOMI_NO2_Africa_201808_201907 = [pd.read_csv(file,delim_whitespace=True, header=None) for file in TROPOMI_NO2_Africa_201808_201907_files] 

# merge the monthly data into a single data frame
TROPOMI_NO2_Africa_201808_201907 = pd.concat(TROPOMI_NO2_Africa_201808_201907)

# reset the index (row numbers)
TROPOMI_NO2_Africa_201808_201907 = TROPOMI_NO2_Africa_201808_201907.reset_index(drop=True)

# insert the row numbers into the output file
TROPOMI_NO2_Africa_201808_201907.iloc[:,0] = TROPOMI_NO2_Africa_201808_201907.index
##############################################################################
# read each row of this dataframe, and output in the format required by the Fortran oversampling codes
def write_to_output_file(input_data,outfile_name): 
        """read each row of the data, and output in the format required by the Fortran oversampling codes"""  
        # create an output file, set path and output file name
        input_data.outfile = '/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/TROPOMI_oversampling_data/Oversampling_Input_Data/'+outfile_name
        # open the file
        input_data.fId = open(input_data.outfile, "w+") 
        # initialize line counts in the file
        input_data.line_count = 0
        # loop over rows, output information to the strings in the way required by the Fortran oversampling codes
        for i in range(input_data.shape[0]):
                    tstr="{:8}".format(input_data.iloc[i,0])+("{:15.6f}"*13).format(input_data.iloc[i,1],input_data.iloc[i,2],input_data.iloc[i,3],input_data.iloc[i,4],input_data.iloc[i,5],input_data.iloc[i,6],input_data.iloc[i,7],input_data.iloc[i,8],input_data.iloc[i,9],input_data.iloc[i,10],input_data.iloc[i,11],input_data.iloc[i,12],input_data.iloc[i,13])+("{:15.6E}"*2).format(input_data.iloc[i,14],input_data.iloc[i,15])
                    input_data.fId.write(tstr) # write tstr to the file 
                    input_data.fId.write("\n") # progresses to the next line
                    input_data.line_count += 1 # increment the line number
        input_data.fId.close # close the file at end of loop
##############################################################################        
# use a small dataset (6 swaths files on 2018-08-01) to test the function 
# this data is also used when comparing the "standard error" vs "standard deviation"
test1 = pd.read_csv('test_SE',delim_whitespace=True, header=None)

# output the file in the format required by Fortran oversampling routine
write_to_output_file(test1,'test_output')

# read the output file that you've just created
test2 = pd.read_csv('test_output',delim_whitespace=True, header=None)

# check if the raw input file and the one you've created are identical
test1.equals(test2)
##############################################################################
# generate the TROPOMI oversampling input file for one year over Africa
write_to_output_file(TROPOMI_NO2_Africa_201808_201907,'TROPOMI_oversampling_input_AF_NO2_0.75_20180801_20190731')

# End
##############################################################################
