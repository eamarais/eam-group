import pandas as pd
import netCDF4
import argparse
import dateutil
import numpy as np
from os import path
import datetime as dt

def can_be_float(string):
    """Returns True if a string can be converted to a float,
    False if not"""
    try:
        float(string)
        return True
    except ValueError:
        return False


def extract_notes_from_lat(lat_val):
    """Extracts notes from a latitude column"""
    if can_be_float(lat_val):
        return None
    else:
        return lat_val


def build_datetime_column(df):
    """Creates a date column from year, month, day and time columns"""
    hour_column = (df['Time (UTC)']
                   .fillna(0)
                   .map(lambda dec_hour : pd.Timedelta(hours=np.floor(dec_hour))))
    minute_column = (df['Time (UTC)']
                     .fillna(0)
                     .map(lambda dec_hour : pd.Timedelta(minutes=(dec_hour%1)*60)))
    date_column = pd.to_datetime(df[['Day','Month','Year']])
    return date_column + hour_column + minute_column


def calc_lat_lon_count(df):
    lat_min = df.Latitude.min()
    lat_max = df.Latitude.max()
    lon_min = df.Longitude.min()
    lon_max = df.Longitude.max()
    precision = 0.001   # Hard-coding for now

    lat_count = (lat_max - lat_min)/precision
    lon_count = (lon_max - lon_min)/precision

    return lat_count, lon_count


def calc_time_axis(df):
    """Returns the number of minutes between midnight of the earliest date
    and midnight of the day after the latest date in the dataset"""
    time_min = df.Timestamp.min()
    time_max = df.Timestamp.max()
    
    diff_in_days =(time_max.date() - time_min.date())
    diff_in_minutes = (diff_in_days.days + 1)*24*60
    return diff_in_minutes


def make_new_corads_netcdf(df, path):
    ds = netCDF4.Dataset(path, 'w', clobber="True")
    
    ds.createDimension('n', len(df.index))
    ds.createDimension('str_len', 10)

    for column_name, data in df.iteritems():
        short_name = column_name.partition('/')[0].strip()
 
        if data.dtype == object:
            new_variable = ds.createVariable(short_name, 'S1', ('n', 'str_len'))
            new_variable._Encoding = 'ascii'
            new_variable[:] = data.array.astype('S10')
        else:
            new_variable = ds.createVariable(short_name, 'f8', ('n'))
            new_variable[:] = data.array
        
        new_variable.long_name = column_name
        if column_name.endswith('kg'):
            new_variable.units = "kg"
        if column_name in ["Latitude", "Longitude"]:
            new_variable.units = "Degrees"
        
    
    ds.close()


def clean_dataframe(df):
    df.rename(columns=lambda x: x.strip(), inplace=True)
    df['Time (UTC)'] = df['Time (UTC)'].where(df['Time (UTC)'].map(can_be_float)).astype(float)
    df['Timestamp'] = build_datetime_column(df)
    df['Notes'] = df['Latitude'].map(extract_notes_from_lat)
    df['Latitude'] = df['Latitude'].where(df['Latitude'].map(can_be_float)).astype(float)
    return df


def main(excel_path, out_dir):
    debris_df = pd.read_excel(excel_path, sheet_name='Re-entries')
    debris_df = clean_dataframe(debris_df)
    make_new_corads_netcdf(debris_df, path.join(out_dir, "reentries.nc4")) 
    
    

    launch_df = pd.read_excel(excel_path, sheet_name='Launch emissions')
    launch_total_df = launch_df.tail(1)
    launch_df.drop(launch_df.tail(1).index,inplace=True)
    launch_df = clean_dataframe(launch_df)
    make_new_corads_netcdf(launch_df, path.join(out_dir, "launches.nc4"))
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Converts your rocketry spreadsheet to netCDF")
    parser.add_argument("--excel_path", default = "Launch-Re-entry-emis-only-03-08.xlsx")
    parser.add_argument("--out_dir", default = '.')
    args = parser.parse_args()
    main(args.excel_path, args.out_dir)

