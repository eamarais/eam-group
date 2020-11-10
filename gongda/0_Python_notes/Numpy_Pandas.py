#####################################################################
# Numpy

import numpy as np

# array/matrix/sequence
a = np.array([1,2,3,4,5])
b = np.zeros(3,dtype=float) 
c = np.zeros((3,5),dtype=int)
d = np.full((3,5),'NA')
e = np.arange(0,20,2)
f = np.linspace(0,1,5)

print(a,b,c,d,e,f,sep="\n")
print(c.ndim,c.shape,c.size,c.dtype,end='')

# dtype is important
x = np.array([1,2,3,4,5])
x[0] = 3.14
print(x)

y = np.array([1,2,3,4,5],dtype = float)
y[0] = 3.14
print(y)

# make an empty N-d array
test = np.zeros((3,5))
#####################################################################
# Pandas

import os
import glob
import pandas as pd

os.chdir("/../../../your working directory/")

# read csv files
surface_files = sorted(glob.glob("china_sites_201610*.csv"))
surface_data  = [pd.read_csv(file) for file in surface_files] 

# only select rows for "NO2" measurements at all sites
NO2  = [data[data['type'] == 'NO2'] for data in surface_data] 

# merge the list of data frames into a single dataframe by rows
NO2_total  = pd.concat(NO2) 
print(len(NO2_total))       
print(NO2_total.head())     
print(NO2_total.tail())    

# subset the dataframe using "df.iloc[rows,columns]"
NO2_test = NO2_total.iloc[0:3,0:6] 
print(NO2_test)

# merge two data frames by column
df1 = NO2_total.iloc[0:3,0:4]
df2 = NO2_total.iloc[0:3,0:6]
NO2_test  = pd.merge(df1,df2,left_on = 'date',right_on = 'date') 
print(NO2_test)

# select a column by column name
NO2_date_1 = NO2_total['date']
NO2_date_2 = NO2_total.date

# subset the dataframe by column names
Info  = ['date','hour','type']
Sites = ['1001A','1002A','1003A'] 
NO2_test = NO2_total[Info + Sites]
print(NO2_test)

# subset the dataframe by conditions
NO2_test_1 = NO2_total[(NO2_total['1001A'] < 80) & (NO2_total['1002A'] < 60)] 
NO2_test_2 = NO2_total[(NO2_total['1001A'] < 80) | (NO2_total['1002A'] < 60)] 

print(NO2_test_1)
print(NO2_test_2)

# delete columns by column names
NO2_test = NO2_total.drop(['date','hour','type'],1)
print(NO2_test.head())

# add a column
NO2_test = NO2_total.iloc[0:3,0:6] 
NO2_test['new_column'] = 'test'
print(NO2_test)

# convert YYYYMMDD to YYYY-MM-DD
NO2_date = pd.to_datetime(NO2_total.date, format = '%Y%m%d')
print(NO2_date.head())

# another way to handle datetime
df = pd.DataFrame()
df['year'] = [2019,2020,2021]
df['month'] = [1,2,3]
df['day'] = [28,29,30]
df = df.astype(str)
df['date'] = df['year'] +'-'+ df['month'] + '-' + df['day']
print(df)

# save the pandas dataframe to csv
NO2_total.to_csv('NO2_total.csv', index=False)

# save the pandas dataframe to netcdf
import xarray as xr

# convert pandas data frame to xarray data format
test = xr.Dataset.from_dataframe(NO2_total)

# see the data format is already changed
test

# add attributes of variables
test['NO2'].attrs = {'unit':'ug/m3'}
test['Site_Longitude'].attrs = {'unit':'Degrees_east'}
test['Site_Latitude'].attrs={'unit':'Degrees_north'}

# save xarray data to netcdf format
test.to_netcdf('China_NO2.nc')

# check if the output is what you expected
China_NO2_nc = xr.open_dataset('China_NO2.nc')

# some fun tips of Pandas: https://towardsdatascience.com/30-examples-to-master-pandas-f8a2da751fa4

# End
#####################################################################
