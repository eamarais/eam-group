#############################################################################################
# This script is to read netCDF files, and save the data fields into a single Excel spreadsheet or multiple csv files.
# The netCDF files are normally suggested, but sometimes the collaborators have no experience of handling netCDF files.
# Here I use the MEIC emission inventory netcdf files as the example.

import os
import pandas as pd
import xarray as xr
#############################################################################################
#############################################################################################
# simple example: MEIC emissions at 05x0666 (emissions from all sectors are merged: files already prepared for direct use in GEOS-Chem)
os.chdir('/rds/projects/2018/maraisea-glu-01/RDS/RDS_Data/BTH_project/Inventory/MEIC_05x0666/')

# open netcdf files
MEIC_OC = xr.open_dataset("MEIC_OC.05x0666.nc")
MEIC_BC = xr.open_dataset("MEIC_BC.05x0666.nc")
print(MEIC_OC,MEIC_BC,sep="\n###################")

# convert xarray data array to pandas dataframe
def xarray_to_pandas(data):
    data = data.to_dataframe()
    data.reset_index(inplace=True)
    return data

# why the for loop doesn't work here?
MEIC_OC_df = xarray_to_pandas(MEIC_OC)
MEIC_BC_df = xarray_to_pandas(MEIC_BC)

# extract values in 2017 and group by month
MEIC_OC_2017 = []
MEIC_BC_2017 = []

for i in range(12):
    MEIC_OC_2017.append(MEIC_OC_df[(pd.DatetimeIndex(MEIC_OC_df['time']).year == 2017) & 
                                   (pd.DatetimeIndex(MEIC_OC_df['time']).month == i +1)])
    MEIC_BC_2017.append(MEIC_BC_df[(pd.DatetimeIndex(MEIC_BC_df['time']).year == 2017) & 
                                   (pd.DatetimeIndex(MEIC_BC_df['time']).month == i +1)])
    
# think about a better/safer way to drop rows where all data fields are "NaN"
for i in range(len(MEIC_OC_2017)):
    MEIC_OC_2017[i] = MEIC_OC_2017[i][MEIC_OC_2017[i]['OC_agriculture'] >= 0]
    MEIC_BC_2017[i] = MEIC_BC_2017[i][MEIC_BC_2017[i]['BC_agriculture'] >= 0]
    
# reset index
MEIC_OC_2017 = [x.reset_index(drop=True) for x in MEIC_OC_2017]
MEIC_BC_2017 = [x.reset_index(drop=True) for x in MEIC_BC_2017]

# save results to a single xlsx file
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/BTH/geoschem')

yymm=list(range(201701,201713))

writer=pd.ExcelWriter(r"MEIC_OC_2017_05x0666.xlsx")
for i,data in enumerate(MEIC_OC_2017):
    data.to_excel(writer,sheet_name="{0}".format(yymm[i]))

writer.save()

writer=pd.ExcelWriter(r"MEIC_BC_2017_05x0666.xlsx")
for i,data in enumerate(MEIC_BC_2017):
    data.to_excel(writer,sheet_name="{0}".format(yymm[i]))

writer.save()

# save results to multiple csv files
for i in range(12):
    MEIC_OC_2017[i].to_csv("MEIC_OC_05x0666"+str(yymm[i])+".csv",index=False,sep=',')
    MEIC_BC_2017[i].to_csv("MEIC_BC_05x0666"+str(yymm[i])+".csv",index=False,sep=',')
#############################################################################################
#############################################################################################
# complicated example: MEIC emissions at 025x025 (emissions from all sectors are seperated: raw files)
os.chdir('/rds/projects/2018/maraisea-glu-01/RDS/RDS_Data/BTH_project/Inventory/MEIC_025x025')

# first, how to remove all the items defined previously in this job?
background = [file for file in globals().keys()]
del background

# import multiple files
import glob
import re

MEIC_OC_IND = glob.glob("*industry-OC.nc")
MEIC_OC_POW = glob.glob("*power-OC.nc")
MEIC_OC_TRA = glob.glob("*transportation-OC.nc")
MEIC_OC_RES = glob.glob("*residential-OC.nc")
MEIC_OC_AGR = glob.glob("*agriculture-OC.nc")

# group the items, so you can perform the same functions for all
# if I can clean the working space, I think I will be able to group the items which name start with "MEIC_"?
# or are there other ways to combine the files with similar names more efficiently? As there will be more data fields and items to be defined.

all_MEIC_OC = [MEIC_OC_IND,MEIC_OC_POW,MEIC_OC_TRA,MEIC_OC_RES,MEIC_OC_AGR]

# sort all files numerically
for file in all_MEIC_OC:
    file.sort(key=lambda var:[int(x) if x.isdigit() else x for x in re.findall(r'[^0-9]|[0-9]+', var)])

# check the sorted files
print('number of files:',len(MEIC_OC_IND),MEIC_OC_IND[0],MEIC_OC_IND[-1],sep="  ")

# To read all files together "all_MEIC_OC = [xr.open_dataset(file) for x in all_MEIC_OC for file in x]"
# but this is not a good option here, because the emission rate fields are named "z" in the raw files for all the species 
# and there is no info/attribute within the file to distinguish each species

# so for now, I extract emission from each sector seperately
MEIC_OC_IND = [xr.open_dataset(file) for file in MEIC_OC_IND]
MEIC_OC_POW = [xr.open_dataset(file) for file in MEIC_OC_POW]
MEIC_OC_TRA = [xr.open_dataset(file) for file in MEIC_OC_TRA]
MEIC_OC_RES = [xr.open_dataset(file) for file in MEIC_OC_RES]
MEIC_OC_AGR = [xr.open_dataset(file) for file in MEIC_OC_AGR]

# convert xarray data array to pandas dataframe 
MEIC_OC_IND_df = [xarray_to_pandas(data) for data in MEIC_OC_IND]
MEIC_OC_POW_df = [xarray_to_pandas(data) for data in MEIC_OC_POW]
MEIC_OC_TRA_df = [xarray_to_pandas(data) for data in MEIC_OC_TRA]
MEIC_OC_RES_df = [xarray_to_pandas(data) for data in MEIC_OC_RES]
MEIC_OC_AGR_df = [xarray_to_pandas(data) for data in MEIC_OC_AGR]

# check one example
print(MEIC_OC_IND_df[0].head())

# but why the loop below does not work?
# for x in all_MEIC:
#     x = [xr.open_dataset(file) for file in x]

# so the lat and lon are not provided in the raw file, we have to generate those on our own
lon = np.arange(70+0.25/2,150,0.25)
lat = np.arange(60-0.25/2,10,-.25)
print(len(lon)*len(lat))

def expand_grid(lon, lat):
    xG, yG = np.meshgrid(lon, lat) # create the actual grid
    xG = xG.flatten() # make the grid 1d
    yG = yG.flatten() # same
    return pd.DataFrame({'lon':xG, 'lat':yG}) # return a dataframe

grid_points = expand_grid(lon,lat)

# insert lat and lon to the dataframes
MEIC_OC_IND_df = [pd.concat([data,grid_points], axis=1, sort=False) for data in MEIC_OC_IND_df]
MEIC_OC_POW_df = [pd.concat([data,grid_points], axis=1, sort=False) for data in MEIC_OC_POW_df]
MEIC_OC_TRA_df = [pd.concat([data,grid_points], axis=1, sort=False) for data in MEIC_OC_TRA_df]
MEIC_OC_RES_df = [pd.concat([data,grid_points], axis=1, sort=False) for data in MEIC_OC_RES_df]
MEIC_OC_AGR_df = [pd.concat([data,grid_points], axis=1, sort=False) for data in MEIC_OC_AGR_df]

# assign corresponding names to emissions from each sector
for i in range(len(MEIC_OC_IND_df)):
    MEIC_OC_IND_df[i].rename(columns={'z':'OC_industry'},inplace=True)
    MEIC_OC_POW_df[i].rename(columns={'z':'OC_power'},inplace=True)
    MEIC_OC_TRA_df[i].rename(columns={'z':'OC_transportation'},inplace=True)
    MEIC_OC_RES_df[i].rename(columns={'z':'OC_residential'},inplace=True)
    MEIC_OC_AGR_df[i].rename(columns={'z':'OC_agriculture'},inplace=True)

# combine emissions from all sectors
MEIC_OC = []

for i in range(len(MEIC_OC_IND_df)):
    MEIC_OC.append(pd.concat([MEIC_OC_IND_df[i],MEIC_OC_POW_df[i],MEIC_OC_TRA_df[i],MEIC_OC_RES_df[i],MEIC_OC_AGR_df[i]],axis=1))
    
# remove dulicated columns
MEIC_OC = [df.loc[:,~df.columns.duplicated()] for df in MEIC_OC]

# re-arrange the data
MEIC_OC = [df[['lon','lat','OC_agriculture','OC_industry','OC_power','OC_residential','OC_transportation']] for df in MEIC_OC]

# remove the row if the fill value of "-9999" is deteced. But how to do this in a straightforward way?
MEIC_OC = [df.replace(-9999,np.nan) for df in MEIC_OC]
MEIC_OC = [df.dropna() for df in MEIC_OC]

# save results to a single xlsx
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/BTH/geoschem')

yymm=list(range(201701,201713))

writer=pd.ExcelWriter(r"MEIC_OC_2017_025x025.xlsx")
for i,data in enumerate(MEIC_OC):
    data.to_excel(writer,sheet_name="{0}".format(yymm[i]))

writer.save()

# save results to multiple csvs
for i in range(len(yymmdd)):
    MEIC_OC[i].to_csv("MEIC_OC_2017_025x025"+str(yymm[i])+".csv",index=False,sep=',')
    
# End
#############################################################################################
