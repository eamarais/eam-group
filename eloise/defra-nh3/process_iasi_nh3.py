#!/usr/bin/python

# Import relevant packages:
import os
import sys
import glob
import argparse

import numpy as np
import math

import netCDF4 as nc4
from netCDF4 import Dataset

from geographiclib.geodesic import Geodesic

# Fill value (although not currently used)
FILL_VAL = 9.969209968386869e+36

class EdgeOfData(Exception):
    pass

class IasiData:
    def __init__(self, iasi_file):

        # Read data:
        fh = Dataset(iasi_file, mode='r')
        
        # Extract relevant information:
        self.time = fh.variables['time'][:]
        self.lat = fh.variables['latitude'][:]
        self.lon = fh.variables['longitude'][:]
        self.col = fh.variables['column'][:]
        tcol_unc = fh.variables['error'][:]
        self.cld = fh.variables['CLcov'][:]
        self.sza = fh.variables['angle'][:]
        nvals = len(self.lon)

        # Convert error from % to absolute:
        self.col_unc = tcol_unc * 1e-2 * self.col
        
        # Initialize lon and lat edges (should go elsewhere):
        self.lon_edges = np.zeros((nvals,4))
        self.lat_edges = np.zeros((nvals,4))
        self.cnt = 0

    def get_number_of_data_points(self):
        nvals = len(self.lon)
        return nvals

    def get_lat_lon_edges(self,i,nvals):

        # Calculate distance between neighbouring points to get
        # lat and long coordinates:
        # Crude assumption is that these are halfway between neighbouring pixels:

        # Skip first and last entry (can't calculate edges):
        # Could approximate these as being similar to nearest neighbour
        if i==0 or i==nvals-1: 
            self.lon_edges[i,:] = np.nan
            self.lat_edges[i,:] = np.nan

        else:

            # Initialize:
            lon_edges=np.zeros(4)
            lat_edges=np.zeros(4)

            # CONER 1:
            # Define the path from 1 to 2
            l = Geodesic.WGS84.InverseLine(self.lat[i],self.lon[i],self.lat[i+1],self.lon[i-1])
            # Compute the midpoint:
            m = l.Position(0.5 * l.s13)
            lon_edges[0]=m['lon2']
            lat_edges[0]=m['lat2']

            # CONER 2:
            # Define the path from 1 to 2
            l = Geodesic.WGS84.InverseLine(self.lat[i],self.lon[i],self.lat[i+1],self.lon[i+1])
            # Compute the midpoint:
            m = l.Position(0.5 * l.s13)
            lon_edges[1]=m['lon2']
            lat_edges[1]=m['lat2']

            # CONER 3:
            # Define the path from 1 to 2
            l = Geodesic.WGS84.InverseLine(self.lat[i],self.lon[i],self.lat[i-1],self.lon[i+1])
            # Compute the midpoint:
            m = l.Position(0.5 * l.s13)
            lon_edges[2]=m['lon2']
            lat_edges[2]=m['lat2']

            # CONER 4:
            # Define the path from 1 to 2
            l = Geodesic.WGS84.InverseLine(self.lat[i],self.lon[i],self.lat[i-1],self.lon[i-1])
            # Compute the midpoint:
            m = l.Position(0.5 * l.s13)
            lon_edges[3]=m['lon2']
            lat_edges[3]=m['lat2']

            # Error checks:
            zero_diff_lon = np.where( abs(np.subtract(self.lon[i],lon_edges)) < 1e-5 )[0]
            zero_diff_lat = np.where( abs(np.subtract(self.lat[i],lat_edges)) < 1e-5 )[0]
            if ( len(zero_diff_lon) !=0 or len(zero_diff_lat) !=0  ):
                self.lon_edges[i,:] = np.nan
                self.lat_edges[i,:] = np.nan
            elif ( len(np.unique(lon_edges)) !=4 or len(np.unique(lat_edges)) !=4 ):
                self.lon_edges[i,:] = np.nan
                self.lat_edges[i,:] = np.nan
            else:
                self.lon_edges[i,:] = lon_edges
                self.lat_edges[i,:] = lat_edges

    def write_to_txt_file(self, fId, cnt):

        # Save to text file in format to be read in for RegridPixel fortran routine
        #file_out = open(file, "w+")
        
        # Just set AMF to 0.0 (not needed, but read in to Fortran routine):
        amf=0.0

        # Loop over data:
        for w in range(len(self.lat)):

            # Check for and remove outliers (currently set at 5e17):
            # Could try 3e17 if this doesn't work.
            if abs(self.col[w]) > 5e17: continue
            if abs(self.col[w])==np.nan or abs(self.col[w])==np.inf:
                print('Column value is NAN',flush=True)
                sys.exit(1)
                continue
            if abs(self.col_unc[w])==np.nan or abs(self.col_unc[w])==np.inf:
                print('Column uncertainty is NAN',flush=True)
                sys.exit(1)
                continue
            if self.col_unc[w]==0.0:
                print('Uncertainty is zero',flush=True)
                sys.exit(1)
                continue

            # Check for and remove scenes with very large errors:
            # Doesn't really improve the data. Makes the IASI NH3
            # data over the UK worse.
            #if abs(self.col_unc[w] / self.col[w]) > 1.0: continue

            # Skip points with large difference in midpoint and edges:
            err_ind = np.where(abs(self.lon[w]-self.lon_edges[w,:])>1)[0]

            # Skip first and last:
            if not math.isnan(self.lat_edges[w,0]) and not len(err_ind)>0:

                cnt += 1
                self.cnt += 1

                # Define string of data to print to file:
                tstr="{:8}".format(cnt)+("{:15.6f}"*13).format(self.lat[w],(* self.lat_edges[w,:]),self.lon[w],(* self.lon_edges[w,:]),self.sza[w],self.cld[w],amf)+("{:15.6E}"*2).format(self.col[w],self.col_unc[w])

                # Print to file:
                fId.write(tstr)
                fId.write("\n")

def get_iasi_files(iasidir,stryr,strmon):

    # Loop over years:
    for y in range(len(stryr)):

        iasi_glob_string = os.path.join(iasidir,stryr[y],'nh3nn_v3R_'+stryr[y]+strmon+'*')
        tfiles = glob.glob(iasi_glob_string)
        tfiles = sorted(tfiles)
        if y == 0: 
            iasifiles=tfiles
        else:
            for f in tfiles:
                iasifiles.append(f)

    print('Found {} tropomi files for month {} '.format(len(iasifiles), strmon),flush=True)

    return iasifiles

# Eventually put file/directory def in dedicated section.

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    #parser.add_argument("--season", default='jja')
    parser.add_argument("--month", default='01')
    args = parser.parse_args()

    # Define months based on input:
    #season = args.season
    #if season=='jja': months=['06','07','08']
    #if season=='son': months=['09','10','11']
    #if season=='djf': months=['12','01','02']
    #if season=='mam': months=['03','03','05']

    # Define month:
    month = args.month 
    
    # Initialize count:
    count = 0

    # Define output filename:
    outfile = '/data/uptrop/Projects/DEFRA-NH3/Data/IASI_Oversampled/iasi_nh3_'+month+'_2008-2018'
    file_id = open(outfile, "w+")

    iasi_dir = '/data/uptrop/nobackup/iasi_nh3/ANNI_v3R_metopa/'

    years=['2008','2009','2010','2011','2012','2013',\
           '2014','2015','2016','2017','2018']

    iasi_files = get_iasi_files(iasi_dir, years, month)

    # Loop over files:
    for f in iasi_files:

        print('Processing: ', f,flush=True)

        iasi_data = IasiData(f)

        nvals = iasi_data.get_number_of_data_points()
        
        # Loop over points to get lon and lat edges:
        for i in range(nvals):
            iasi_data.get_lat_lon_edges(i,nvals)

        iasi_data.write_to_txt_file(file_id,count)

        # Increment:
        count += iasi_data.cnt

    # Close file:
    file_id.close
