###############################################################################
###############################################################################
# Use web scraping to download data from Air Quality System (AQS) API provided by US EPA

# It is recommended to read instructions from US EPA to understand how this API works before using this script
# US EPA air quality home page: https://www.epa.gov/outdoor-air-quality-data
# AQS API: https://aqs.epa.gov/aqsweb/documents/data_api.html

import os
import pandas as pd
import json
from urllib.request import urlopen

# set up working directory
os.chdir("Y:\\Study\\Research_Data\\COVID-ML\\MODEL_INPUT_DATA\\USA\\OpenAQ_validation")

# build a function to download json data from the target API link
# keep "Data" only ("Header" descrips the status of this request)

def get_data(download_link):
    """open the input ipAddress and save data of interest"""
    response = urlopen(download_link).read().decode("utf-8")
    responseJson = json.loads(response)
    return responseJson.get("Data") 
    
# try a sample link
test_link = "https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param=42602&bdate=20190101&edate=20190131&state=06&county=037&site=1103"
test_results = get_data(test_link)

# the function returns a list of dictionaries
print(type(test_results))
print(type(test_results[0]))
print("dictionary keys:",test_results[0].keys(),sep="\n")
print("number of dictionaries:",len(test_results))

# check the first record
print(test_results[0])

# check the last record
print(test_results[-1])

# see the full results
test_results

# save this list of dictionaries to pandas dataframe
test_data = [pd.DataFrame([data]) for data in test_results]
test_df = pd.concat(test_data,ignore_index=True)
test_df

###############################################################################
# now use the function above to download massive data

# download NO2,SO2,CO,O3,PM2.5,PM10 from 5 sites in Los Angeles in 2019
"""
summary of the codes needed to requeset data:
1> email: your registered email address
2> key: the key you received after registration
3> parameters
   NO2: 42602
   SO2: 42401
   CO: 42101
   O3: 44201
   PM2.5 FRM/FEM Mass: 88101
   PM2.5 non FRM/FEM Mass: 88502
   PM10: 81102
4> bdate: begin date
5> edate: end date
6> state code: 06 for CA
7> county code: 037 for LA
8> site code:
   Compton: 1302
   Lancaster: 9033
   North Main: 1103 
   LAX: 5005
   Glendora: 0016
"""

# create download links to all species at each site in 2019
# mannually type 0016 for now
sites = [1302,9033,1103,5005]
parameters = [42602,42401,42101,44201,88101,88502,81102]
species = ['NO2','SO2','CO','O3','PM2.5 FRM FEM Mass','PM2.5 non FRM FEM Mass','PM10'] # this is used later

LA_Compton = []
LA_Lancaster = []
LA_North_Main = []
LA_LAX = []
LA_Glendora = []

for i in range(len(parameters)):
    LA_Compton.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate=20190101&edate=20191231&state=06&county=037&site="+str(sites[0]))
    LA_Lancaster.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate=20190101&edate=20191231&state=06&county=037&site="+str(sites[1]))
    LA_North_Main.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate=20190101&edate=20191231&state=06&county=037&site="+str(sites[2]))
    LA_LAX.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate=20190101&edate=20191231&state=06&county=037&site="+str(sites[3]))
    LA_Glendora.append("https://aqs.epa.gov/data/api/sampleData/bySite?email=gxl642@student.bham.ac.uk&key=carmelmallard48&param="+str(parameters[i])+"&bdate=20190101&edate=20191231&state=06&county=037&site=0016")
    
    # get all species at each site in 2019
# no idea why list compreshension does not work here: LA_Compton_results = [get_data(link) for link in LA_Compton]

LA_Compton_results = []
LA_Lancaster_results = []
LA_North_Main_results = []
LA_LAX_results = []
LA_Glendora_results = []

for i in range(len(parameters)):
    LA_Compton_results.append(get_data(LA_Compton[i]))
    LA_Lancaster_results.append(get_data(LA_Lancaster[i]))
    LA_North_Main_results.append(get_data(LA_North_Main[i]))
    LA_LAX_results.append(get_data(LA_LAX[i]))
    LA_Glendora_results.append(get_data(LA_Glendora[i]))
    
    # data of each species at each site are returned as a list of dictionaries
# build a function to combine the results to pandas dataframes

def save_raw_EPA_results_to_df(raw_EPA_data_results):
    """For each request, the API returns measurements of one species at one site during one sampling period.
       Results are returned as a list of dictionaries. This function converts the results from each request to a single pandas dataframe.
    """
    test_data = [pd.DataFrame([data]) for data in raw_EPA_data_results]
    test_df = pd.concat(test_data,ignore_index=True)
    return test_df
    
    # convert results to pandas for all species and sites
LA_Compton_df = []

for i in range(len(parameters)):
    if (len(LA_Compton_results[i]) > 0):
        LA_Compton_df.append(save_raw_EPA_results_to_df(LA_Compton_results[i]))
    else:
        LA_Compton_df.append("There is no observation for "+str(species[i]))

LA_Lancaster_df = []

for i in range(len(parameters)):
    if (len(LA_Lancaster_results[i]) > 0):
        LA_Lancaster_df.append(save_raw_EPA_results_to_df(LA_Lancaster_results[i]))
    else:
        LA_Lancaster_df.append("There is no observation for "+str(species[i]))

LA_North_Main_df = []

for i in range(len(parameters)):
    if (len(LA_North_Main_results[i]) > 0):
        LA_North_Main_df.append(save_raw_EPA_results_to_df(LA_North_Main_results[i]))
    else:
        LA_North_Main_df.append("There is no observation for "+str(species[i]))

LA_LAX_df = []

for i in range(len(parameters)):
    if (len(LA_LAX_results[i]) > 0):
        LA_LAX_df.append(save_raw_EPA_results_to_df(LA_LAX_results[i]))
    else:
        LA_LAX_df.append("There is no observation for "+str(species[i]))

LA_Glendora_df = []

for i in range(len(parameters)):
    if (len(LA_Glendora_results[i]) > 0):
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
    if (len(LA_LAX_results[i]) > 0):
        LA_LAX_df[i].to_csv("LA_LAX_2019_"+str(species[i]+".csv"))
###############################################################################
