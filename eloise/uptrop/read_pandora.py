#!/usr/bin/python

# Import relevant packages:
import glob
import sys
import os
import re

import numpy as np
import pandas as pd


def get_column_description_index(filename):
    """
    Returns a dictionary of {description:column index} for a pandora file
    """
    # See https://regex101.com/r/gAjFtL/1 for an in-depth explanation of this regex
    # Returns two groups; group 1 is the column number, group 2 is the description.
    searcher = re.compile(r"Column ([0-9]+): (.*)")
    with open(filename, 'r', encoding="Latin-1") as pandora_file:
        pandora_text = pandora_file.read()
    groups = searcher.findall(pandora_text)
    column_dict = {column_description: int(column_index) for column_index, column_description in groups}
    return column_dict


def get_start_of_data(filepath):
    """Gets the line number of the start of the data itself"""
    with open(filepath, encoding="Latin-1") as pandora_file:
        # Line numbers are 1-indexed....
        line_number = 1
        # Look for the dotted line twice
        while not pandora_file.readline().startswith("-------"):
            line_number += 1
        while not pandora_file.readline().startswith("-------"):
            line_number += 1
    # Increment one more time to get the index of the first line of actual data
    line_number += 1
    return line_number


def get_column_from_description(column_dict, description):
    """Returns the column index matching a substring of description"""
    index = [value for key, value in column_dict.items() if description in key][0]
    if index is None:
        return
    else:
        return index


def get_lat_lon(filename):
    """Returns a dictionary of lat, lon extracted from the pandora file """
    lat = None
    lon = None
    with open(filename, 'r', encoding="Latin-1") as pandora_file:
        while (lat == None or lon == None):
            current_line = pandora_file.readline()
            if current_line.startswith("Location latitude [deg]:"):
                lat = float(current_line.split(":")[1])
            elif current_line.startswith("Location longitude [deg]:"):
                lon = float(current_line.split(":")[1])
            elif current_line.startswith("--------"):
                # TODO: Maybe change for exception if this might happen
                print("Lat/lon not found in file {}".format(filename))
                return
    return {"lat": lat, "lon": lon}


def readpandora(filename,no2col):

    """
    Reads a Pandora file to a pandas dataframe
    """

    loc = get_lat_lon(filename)

    column_dict = get_column_description_index(filename)
    dateind = get_column_from_description(column_dict, 'UT date and time for center of m')
    jdayind = get_column_from_description(column_dict, 'Fractional days since 1-Jan-2000')
    # SZA:
    szaind = get_column_from_description(column_dict,  'Solar zenith angle for center of')
    # NO2 column:
    # (a) Total:
    if no2col == 'Tot':
        no2ind = get_column_from_description(column_dict,  'Nitrogen dioxide total vertical ')
    # (b) Troposphere:
    if no2col == 'Trop':
        no2ind = get_column_from_description(column_dict,  'Nitrogen dioxide tropospheric v')
    # NO2 column error:
    # (a) Total:
    if no2col == 'Tot':
        errind = get_column_from_description(column_dict,  'Uncertainty of nitrogen dioxide total ver')
    # (b) Troposphere:
    if no2col == 'Trop':
        errind = get_column_from_description(column_dict,  'Uncertainty of nitrogen dioxide troposph')
    # Data quality flag:
    qaflagind = get_column_from_description(column_dict,  'L2 data quality flag for nitrog')

    # Level 2 fit flag:
    # There are two entries (min and max) of this in the tropospheric
    # column file. The minimum is being used here.
    fitflagind = get_column_from_description(column_dict,'Level 2 Fit data quality flag')

    data_start = get_start_of_data(filename)

    names = ["ut_date", "jday", "sza", "no2", "no2err", "qaflag", "fitflag"]
    columns = [dateind, jdayind,szaind,no2ind,errind, qaflagind, fitflagind]
    columns = [column -1 for column in columns]  # Pandora columns are 1-indexed, Pandas are 0

    df = pd.read_csv(filename,
                     sep=" ",
                     skiprows=data_start,
                     usecols=columns,
                     names=names,
                     parse_dates=['ut_date']
                     )

    date_df = pd.DataFrame({
        "year": df.ut_date.dt.year,
        "month": df.ut_date.dt.month,
        "day": df.ut_date.dt.day,
        "hour_utc": df.ut_date.dt.hour,
        "minute": df.ut_date.dt.minute
    })

    df = pd.concat((date_df, df), axis=1)
    df = df.drop("ut_date", axis=1)

    # Output:
    return (loc, df)
