###############################################################################
###############################################################################
# Download data from Air Quality System (AQS) API provided by US EPA

# The US EPA provides historical data in "Pre-generated Data files", but the records are not always up-to-date
# This script uses a web scraping technique to request data of interest from AQS API

import os
import pandas as pd
import json
from urllib.request import urlopen
###############################################################################
# set up working directory
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/1_co_author_papers/COVID-ML/MODEL_INPUT_DATA/USA/OpenAQ_validation")
###############################################################################
# US EPA pages will help you to understand how the API works
# EPA AQ data home page： https://www.epa.gov/outdoor-air-quality-data
# AQS API: https://aqs.epa.gov/aqsweb/documents/data_api.html
# Map of sites: https://gispub.epa.gov/airnow/?mlayer=ozonepm&clayer=none&panel=0%203%3E%20TROPOMI%20codes%20-%20visualisation

# you need target state,county,site codes, species parameter codes to download the data
# go through the links below to understand how the network is structured

# list of states in the US
https://aqs.epa.gov/data/api/list/states?email=test@aqs.api&key=test
# list of counties in a sample state
https://aqs.epa.gov/data/api/list/countiesByState?email=test@aqs.api&key=test&state=06
# list of sites in a sample county
https://aqs.epa.gov/data/api/list/sitesByCounty?email=test@aqs.api&key=test&state=06&county=037
# request data at a sample site
https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param=42602&bdate=20190101&edate=20190131&state=06&county=037&site=2005
###############################################################################
# now build a function to download json data from the target API link
# mearsurements are in "Data" ("Header" returns the status of this request)

def get_data(download_link):
    """open the target download link, return the data or error message (if there is no data)"""
    response = urlopen(download_link).read().decode("utf-8")
    responseJson = json.loads(response)
    if (len(responseJson.get("Data")) == 0):
        return responseJson.get("Header")
    else:
        return responseJson.get("Data") 

# request data at the same sample site
# you need to input your registered email and key
test_link = "https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param=42602&bdate=20190101&edate=20190131&state=06&county=037&site=2005"
test_results = get_data(test_link)
test_results

# the function returns a list of dictionaries
print(type(test_results[0]))
print("dictionary keys:",test_results[0].keys(),sep="\n")
print("number of dictionaries:",len(test_results))

# check the first record
print(test_results[0])
###############################################################################
###############################################################################
# New York City seems to be missing from the historical data provided by EPA
# so here we use the script to do a systematic check for New York City

# New York state code: 36
# New York City county code: 061
# NO2: 42602
NY_sites = get_data("https://aqs.epa.gov/data/api/list/sitesByCounty?email=test@aqs.api&key=test&state=36&county=061")

# There seems to be a lot of missing sites ('value_represented': None)
NY_sites

# get the site code list in New York City
NY_sites_codes = []

for i in range(len(NY_sites)):
        NY_sites_codes.append(NY_sites[i]['code']) 
        
print(NY_sites_codes)

# generate corresponding links for each site code
NY_sites_links = []

for i in range(len(NY_sites_codes)):
    NY_sites_links.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param=42602&bdate=20190101&edate=20190131&state=36&county=061&site="+str(NY_sites_codes[i]))

# print out all the download links
for i in range(len(NY_sites_links)):
    print(NY_sites_links[i])
    
# get data from each link
NY_data = [get_data(link) for link in NY_sites_links]   

# print out the results
for i in range(len(NY_data)):
    print(NY_data[i])
    
# so there is no data for NO2 in New York City in 2019    
###############################################################################
# Now download NO2,SO2,CO,O3,PM2.5,PM10 from in Los Angeles
# first summarize the codes needed for iput variables (parameter,bdate,edate,state code,county code,site code)

# store information using regular expressions

import re

parameters = '''
             NO2: 42602
             SO2: 42401
             CO: 42101
             O3: 44201
             PM2.5 FRM/FEM Mass: 88101
             PM2.5 non FRM/FEM Mass: 88502
             PM10: 81102
             '''
date = '''
       begin date: 20190101
       end date: 20191231
       '''

state_code = '''
             California: 06
             '''

county_code = '''
              Los Angeles: 037
              '''

site_code = '''
             Compton: 1302
             Lancaster: 9033
             North Main: 1103
             LAX: 5005
             Glendora: 0016
             '''

# only keep codes on the right (tell Python to find "digits" after ": ")
parameters = re.findall(r'(?<=:\s)\d+',parameters)
date = re.findall(r'(?<=:\s)\d+',date)
state_code = re.findall(r'(?<=:\s)\d+',state_code)
county_code = re.findall(r'(?<=:\s)\d+',county_code)
site_code = re.findall(r'(?<=:\s)\d+',site_code)

# create download links to all species at each site
LA_Compton = []
LA_Lancaster = []
LA_North_Main = []
LA_LAX = []
LA_Glendora = []

for i in range(len(parameters)):
    LA_Compton.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate="+str(date[0])+"&edate="+str(date[1])+"&state="+str(state_code[0])+"&county="+str(county_code[0])+"&site="+str(site_code[0]))
    LA_Lancaster.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate="+str(date[0])+"&edate="+str(date[1])+"&state="+str(state_code[0])+"&county="+str(county_code[0])+"&site="+str(site_code[1]))
    LA_North_Main.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate="+str(date[0])+"&edate="+str(date[1])+"&state="+str(state_code[0])+"&county="+str(county_code[0])+"&site="+str(site_code[2]))
    LA_LAX.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate="+str(date[0])+"&edate="+str(date[1])+"&state="+str(state_code[0])+"&county="+str(county_code[0])+"&site="+str(site_code[3]))
    LA_Glendora.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate="+str(date[0])+"&edate="+str(date[1])+"&state="+str(state_code[0])+"&county="+str(county_code[0])+"&site="+str(site_code[4]))
    
# create the list of species in the same order (for later use when naming the output files)
species = ['NO2','SO2','CO','O3','PM2.5 FRM FEM Mass','PM2.5 non FRM FEM Mass','PM10'] 

# get all species at each site in 2019
LA_Compton_results = [get_data(link) for link in LA_Compton]
LA_Lancaster_results = [get_data(link) for link in LA_Lancaster]
LA_North_Main_results = [get_data(link) for link in LA_North_Main]
LA_LAX_results = [get_data(link) for link in LA_LAX]
LA_Glendora_results = [get_data(link) for link in LA_Glendora]

# if there is data, measurements during the samplng period are returned as a list of dictionaries
# now build a function to convert the list of dictionary to a pandas dataframes

def save_raw_EPA_results_to_df(raw_EPA_data_results):
    """For each request, the API returns measurements of one species at one site during one sampling period.
       Results are returned as a list of dictionaries. This function converts the results from each request to a single pandas dataframe.
    """
    test_data = [pd.DataFrame([data]) for data in raw_EPA_data_results]
    test_df = pd.concat(test_data,ignore_index=True)
    return test_df

# if there is no data, an error message is kept
# the length will be "1"

# now convert results to pandas if there is data
# 1> Comton site
LA_Compton_df = []

for i in range(len(parameters)):
    if (len(LA_Compton_results[i]) > 1):
        LA_Compton_df.append(save_raw_EPA_results_to_df(LA_Compton_results[i]))
    else:
        LA_Compton_df.append("There is no observation for "+str(species[i]))

# 2> Lancaster site
LA_Lancaster_df = []

for i in range(len(parameters)):
    if (len(LA_Lancaster_results[i]) > 1):
        LA_Lancaster_df.append(save_raw_EPA_results_to_df(LA_Lancaster_results[i]))
    else:
        LA_Lancaster_df.append("There is no observation for "+str(species[i]))

# 3> North Main street site
LA_North_Main_df = []

for i in range(len(parameters)):
    if (len(LA_North_Main_results[i]) > 1):
        LA_North_Main_df.append(save_raw_EPA_results_to_df(LA_North_Main_results[i]))
    else:
        LA_North_Main_df.append("There is no observation for "+str(species[i]))

# 4> LAX site
LA_LAX_df = []

for i in range(len(parameters)):
    if (len(LA_LAX_results[i]) > 1):
        LA_LAX_df.append(save_raw_EPA_results_to_df(LA_LAX_results[i]))
    else:
        LA_LAX_df.append("There is no observation for "+str(species[i]))

# 5> Glendora site
LA_Glendora_df = []

for i in range(len(parameters)):
    if (len(LA_Glendora_results[i]) > 1):
        LA_Glendora_df.append(save_raw_EPA_results_to_df(LA_Glendora_results[i]))
    else:
        LA_Glendora_df.append("There is no observation for "+str(species[i]))

###############################################################################
# output the results as csv files

for i in range(len(parameters)):
    if (len(LA_Compton_results[i]) > 0):
        LA_Compton_df[i].to_csv("LA_Compton_2019_"+str(species[i]+".csv"))
        
for i in range(len(parameters)):
    if (len(LA_Lancaster_results[i]) > 0):
        LA_Lancaster_df[i].to_csv("LA_Lancaster_2019_"+str(species[i]+".csv"))
 
for i in range(len(parameters)):
    if (len(LA_North_Main_results[i]) > 0):
        LA_North_Main_df[i].to_csv("LA_North_Main_2019_"+str(species[i]+".csv"))

for i in range(len(parameters)):
    if (len(LA_LAX_results[i]) > 0):
        LA_LAX_df[i].to_csv("LA_LAX_2019_"+str(species[i]+".csv"))
        
for i in range(len(parameters)):
    if (len(LA_Glendora_results[i]) > 0):
        LA_Glendora_df[i].to_csv("LA_Glendora_2019_"+str(species[i]+".csv"))
        
# For the next version, I aim to use arugments to make the codes more flexible and sharable. 
# Input variables like "parameter","state code" can be imported as arguments.
 ###############################################################################
