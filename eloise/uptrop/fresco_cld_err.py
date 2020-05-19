#!/usr/bin/python

''' Read and regrid cloud information from the official TROPOMI CLOUD product 
    and from the FRESCO cloud product in the TROPOMI NO2 data product file.

    Data are regridded to a 2x2.5 degree global grid and saved as monthly 
    means in NetCDF format for scenes with TROPOMI FRESCO cloud fraction >
    0.7 and clouds top (centroid) pressure in the upper troposphere 
    (180-450 hPa).

    (1) Initial version created on 15 March 2020 by Eloise Marais
        (eloise.marais@le.ac.uk).

    '''

# Import relevant packages:
import glob
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from gamap_colormap import WhGrYlRd
import argparse
from os import path

# Turn off warnings:
# np.warnings.filterwarnings('ignore')

# TODO: Check with Eloise how constant these are going to be
# Define constants/maybe parameters?
FILL_VAL = 9.96921e+36
# DEFINE GRID:
# Define model grid:
MIN_LAT = -88.
MAX_LAT = 88.
MIN_LON = -178.
MAX_LON = 178.
DELTA_LAT = 2  # 4#0.5#4#0.5
DELTA_LON = 2.5  # 5#0.5#5#0.5


class FileMismatchException(Exception):
    """
    Raised when there are an unequal number of cloud and N02 files
    """
    pass


class LatLonException(Exception):
    """
    Raised when there is a LatLon mismatch
    """
    pass


class CloudVariableStore:
    """
    Class containing the running results of various datas for this project; also saving to netCDF and plotting.
    """
    def __init__(self, data_shape):
        """
        Creates an empty CloudVariableStore of shape data_shape
        :param data_shape: Tuple of (presumably) two elements for shape of cloud data
        """
        # TODO Find out from Eloise what these stand for for better variable names
        self.gknmi_cf = np.zeros(data_shape)
        self.gknmi_ct = np.zeros(data_shape)
        self.gknmi_cb = np.zeros(data_shape)
        self.gknmi_cnt = np.zeros(data_shape)
        self.gdlr_cf = np.zeros(data_shape)
        self.gdlr_ct = np.zeros(data_shape)
        self.gdlr_cb = np.zeros(data_shape)
        self.gdlr_od = np.zeros(data_shape)
        self.gdlr_cnt = np.zeros(data_shape)
        self.gdiff_cf = np.zeros(data_shape)
        self.gdiff_ct = np.zeros(data_shape)
        self.nobs_dlr = 0
        self.nobs_fresco = 0
        self.shape = data_shape

    def update_pixel(self, tropomi_data, trop_i, trop_j):
        """
        Updates the appropriate pixe of running cloud variables with the tropomi data at trop_i, trop_j.
        """
        # TODO: Make this vectorised
        # Skip where FRESCO cloud product is NAN:
        if np.isnan(tropomi_data.tffrc[trop_i, trop_j]):
            print("FRESCO cloud product is NaN, skipping...")
            return
        # Skip if there is also no valid DLR data due to
        # poor data quality or missing values:
        if np.isnan(tropomi_data.tdfrc[trop_i, trop_j]):
            print("No DLR data, skipping....")
            return

        # Find corresponding gridsquare
        # TODO: Find a better place to put out_lon
        p = np.argmin(abs(out_lon - tropomi_data.tdlons[trop_i, trop_j]))
        q = np.argmin(abs(out_lat - tropomi_data.tdlats[trop_i, trop_j]))

        # Error checks:
        # Check that lat and lon values are the same (except for
        # those at 180/-180 that can have different signs):
        if (tropomi_data.tdlons[trop_i, trop_j] != tropomi_data.tflons[trop_i, trop_j]) \
                and (abs(tropomi_data.tdlons[trop_i, trop_j]) != 180.0):
            print('Longitudes not the same')
            print(tropomi_data.tdlons[trop_i, trop_j], tropomi_data.tflons[trop_i, trop_j])
            print(p)
            print(out_lon[p])
            raise Exception("Longitudes not the same")
        if tropomi_data.tdlats[trop_i, trop_j] != tropomi_data.tflats[trop_i, trop_j]:
            print('Latitudes not the same')
            raise Exception("Latitudes not the same")

        # Add data to output arrays:
        # if np.isnan(tffrc[i, j]): continue
        self.gknmi_cf[p, q] += tropomi_data.tffrc[trop_i, trop_j]
        self.gknmi_ct[p, q] += np.divide(tropomi_data.tftop[trop_i, trop_j], 100)
        self.gknmi_cnt[p, q] += 1.0
        # if np.isnan(tdfrc[i, j]): continue
        self.gdlr_cf[p, q] += tropomi_data.tdfrc[trop_i, trop_j]
        self.gdlr_ct[p, q] += np.divide(tropomi_data.tdtop[trop_i, trop_j], 100)
        self.gdlr_cb[p, q] += np.divide(tropomi_data.tdbase[trop_i, trop_j], 100)
        self.gdlr_od[p, q] += tropomi_data.tdoptd[trop_i, trop_j]
        self.gdlr_cnt[p, q] += 1.0

    def update_nobs(self, tropomi_data):
        """
        Given a tropomi_data object, updates the number of observations.
        """
        nobs_dlr, nobs_fresco = tropomi_data.get_nobs()
        self.nobs_dlr += nobs_dlr
        self.nobs_fresco += nobs_fresco

    def calc_cloud_statistics(self):
        """
        For each pixel in the final model, calculate the cloud fraction, cloud top pressure, cloud albedo and (for
        DLR) cloud base pressure. Cloud pressure is given in Pa.
        """
        # Get means and differences (only means for cloud base height):
        # (1) Cloud fraction:
        self.gknmi_cf = np.where(self.gknmi_cnt == 0, np.nan, np.divide(self.gknmi_cf, self.gknmi_cnt))
        self.gdlr_cf = np.where(self.gdlr_cnt == 0, np.nan, np.divide(self.gdlr_cf, self.gdlr_cnt))
        self.gdiff_cf = np.subtract(self.gknmi_cf, self.gdlr_cf)
        # (2) Cloud top pressure:
        self.gknmi_ct = np.where(self.gknmi_cnt == 0, np.nan, np.divide(self.gknmi_ct, self.gknmi_cnt))
        self.gdlr_ct = np.where(self.gdlr_cnt == 0, np.nan, np.divide(self.gdlr_ct, self.gdlr_cnt))
        self.gdiff_ct = np.subtract(self.gknmi_ct, self.gdlr_ct)
        # (3) Cloud albedo/optical depth:
        self.gdlr_od = np.where(self.gdlr_cnt == 0, np.nan, np.divide(self.gdlr_od, self.gdlr_cnt))
        # (4) Cloud base pressure (DLR only):
        self.gdlr_cb = np.divide(self.gdlr_cb, self.gdlr_cnt)
        self.gdlr_cb[self.gdlr_cnt == 0.0] = np.nan

    def write_to_netcdf(self, out_dir):
        """Given an out_directory, writes totalled data to netcdf"""
        # Write data to file:
        outfile = os.path.join(out_dir, 'fresco-dlr-cloud-products-' + MMName + '-' + StrYY + '-2x25-v2.nc')
        ncfile = Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
        ncfile.createDimension('xdim', len(out_lon))
        ncfile.createDimension('ydim', len(out_lat))
        # Global attributes:
        ncfile.title = 'FRESCO and DLR TROPOMI monthly mean cloud properties for ' + \
                       StrMM + ' ' + StrYY
        ncfile.subtitle = 'Data written to file in ' + StrMM + ' ' + StrYY
        ncfile.anything = 'Verions used are v010107 for the DLR product and v010302 for the FRESCO product'
        # Longitudes:
        lon = ncfile.createVariable('lon', np.dtype('float'), ('xdim',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        lon[:] = out_lon
        # Latitudes:
        lat = ncfile.createVariable('lat', np.dtype('float'), ('ydim',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lat[:] = out_lat
        # DLR cloud fraction:
        data1 = ncfile.createVariable('dlr_cld_frac', np.dtype('float'), ('xdim', 'ydim'))
        data1.units = 'unitless'
        data1.long_name = 'DLR CAL cloud fraction'
        data1[:] = self.gdlr_cf
        # DLR cloud top pressure:
        data2 = ncfile.createVariable('dlr_cld_top_pres', np.dtype('float'), ('xdim', 'ydim'))
        data2.units = 'hPa'
        data2.long_name = 'DLR CAL cloud top pressure'
        data2[:] = self.gdlr_ct
        # DLR cloud base pressure:
        data3 = ncfile.createVariable('dlr_cld_base_pres', np.dtype('float'), ('xdim', 'ydim'))
        data3.units = 'hPa'
        data3.long_name = 'DLR CAL cloud base pressure'
        data3[:] = self.gdlr_cb
        # DLR cloud optical depth:
        data4 = ncfile.createVariable('dlr_cld_optd', np.dtype('float'), ('xdim', 'ydim'))
        data4.units = 'unitless'
        data4.long_name = 'DLR CAL cloud optical thickness'
        data4[:] = self.gdlr_od
        # FRESCO cloud fraction:
        data5 = ncfile.createVariable('knmi_cld_frac', np.dtype('float'), ('xdim', 'ydim'))
        data5.units = 'unitless'
        data5.long_name = 'FRESCO cloud fraction'
        data5[:] = self.gknmi_cf
        # FRESCO cloud top pressure:
        data6 = ncfile.createVariable('knmi_cld_top_pres', np.dtype('float'), ('xdim', 'ydim'))
        data6.units = 'hPa'
        data6.long_name = 'FRESCO cloud top pressure'
        data6[:] = self.gknmi_ct
        # FRESCO data points in each gridsquare:
        data7 = ncfile.createVariable('knmi_points', np.dtype('float'), ('xdim', 'ydim'))
        data7.units = 'unitless'
        data7.long_name = 'Number of FRESCO data points'
        data7[:] = self.gknmi_cnt
        # DLR data points in each gridsquare:
        data8 = ncfile.createVariable('dlr_points', np.dtype('float'), ('xdim', 'ydim'))
        data8.units = 'unitless'
        data8.long_name = 'Number of DLR data points'
        data8[:] = self.gdlr_cnt
        # close the file.
        ncfile.close()

    def plot_clouds_products(self, plot_dir):
        """For cloud fraction, cloud top pressure, plots the products for DLR, FRESCO and DLR-FRESCO"""
        # TODO: Get MMName and StrYY from filename instead of the global at the bottom of the script
        # PLOT THE DATA:
        m = Basemap(resolution='l', projection='merc', lat_0=0, lon_0=0,
                    llcrnrlon=MIN_LON, llcrnrlat=-70,
                    urcrnrlon=MAX_LON, urcrnrlat=70)
        xi, yi = m(X, Y)
        # (1) Cloud fraction:
        plt.figure(1)
        plt.subplot(3, 1, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.gknmi_cf), vmin=0.7, vmax=1.3, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('FRESCO Cloud fraction')
        plt.subplot(3, 1, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_cf), vmin=0.7, vmax=1.3, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR Cloud fraction')
        plt.subplot(3, 1, 3)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdiff_cf), vmin=-0.2, vmax=0.2, cmap='bwr')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Difference (FRESCO-DLR)')
        plt.savefig(path.join(plot_dir, 'fresco-vs-dlr-cloud-frac-' + MMName + '-' + StrYY + '-v2.ps'), \
                    format='ps')
        # (2) Cloud top pressure:
        plt.figure(2)
        plt.subplot(3, 1, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.gknmi_ct), vmin=150, vmax=500, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('FRESCO Cloud top pressure [hPa]')
        plt.subplot(3, 1, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_ct), vmin=150, vmax=500, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR Cloud top pressure [hPa]')
        plt.subplot(3, 1, 3)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdiff_ct), vmin=-30, vmax=30, cmap='bwr')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Difference (FRESCO-DLR)')
        plt.savefig(path.join(plot_dir, 'fresco-vs-dlr-cloud-top-press-' + MMName + '-' + StrYY + '-v1.ps'),
                    format='ps')
        # (3) Cloud optical depth/albedo (DLR only):
        plt.figure(3)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_od), vmin=0, vmax=200, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR cloud optical thickness')
        plt.savefig(path.join(plot_dir, 'dlr-cloud-optical-depth-' + MMName + '-' + StrYY + '-v1.ps'),
                    format='ps')
        # (4) Cloud base pressure (DLR only):
        plt.figure(4)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_cb), vmin=150, vmax=900, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR base pressure [hPa]')
        plt.savefig(path.join(plot_dir,'dlr-cloud-base-press-' + MMName + '-' + StrYY + '-v1.ps'),
                    format='ps')
        # (5) Number of points (both):
        plt.figure(4)
        plt.subplot(2, 1, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.gknmi_cnt), vmin=0, vmax=5000, cmap=WhGrYlRd)
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('FRESCO No. of obs')
        plt.subplot(2, 1, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_cnt), vmin=0, vmax=5000, cmap=WhGrYlRd)
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR No. of obs')
        plt.savefig(path.join(plot_dir, 'fresco-dlr-number-of-obs-' + MMName + '-' + StrYY + '-v1.ps'),
                    format='ps')
        plt.show()


class TropomiData:
    """
    Class for holding the data for an individual Tropomi file. Applies data filters on creation.
    """
    # TODO: The filter values in get_nobs and filter_t*file are the same _I think_. Put them in attributes, instead of code

    def __init__(self, td_file_path, tf_file_path):
        forb = tf_file_path[104:109]
        dorb = td_file_path[106:111]
        # Check orbit/swath is the same. If not, skip this iteration:
        if (forb != dorb):
            print("Orbit is not swath")
            return
        self.read_tdfile(td_file_path)
        self.read_tffile(tf_file_path)
        self.filter_tdfile()
        self.filter_tffile()
        self.check_parity(dorb)
        self.shape = self.tdlons.shape

    def read_tdfile(self, tdfile):
        """Read DLR cloud data:"""
        dlr_cloud_data = Dataset(tdfile, mode='r')
        # Extract data of interest:
        self.tdlons = dlr_cloud_data.groups['PRODUCT'].variables['longitude'][:][0, :, :]
        self.tdlats = dlr_cloud_data.groups['PRODUCT'].variables['latitude'][:][0, :, :]
        tfrc = dlr_cloud_data.groups['PRODUCT'].variables['cloud_fraction'][:]
        self.tdfrc = tfrc.data[0, :, :]
        ttop = dlr_cloud_data.groups['PRODUCT'].variables['cloud_top_pressure'][:]
        self.tdtop = ttop.data[0, :, :]
        tbase = dlr_cloud_data.groups['PRODUCT'].variables['cloud_base_pressure'][:]
        self.tdbase = tbase.data[0, :, :]
        tqval = dlr_cloud_data.groups['PRODUCT'].variables['qa_value'][:]
        self.tdqval = tqval.data[0, :, :]  # for more accurate CAL product
        toptd = dlr_cloud_data.groups['PRODUCT'].variables['cloud_optical_thickness'][:]
        self.tdoptd = toptd.data[0, :, :]
        tsnow = dlr_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['snow_ice_flag_nise'][:]
        self.tdsnow = tsnow.data[0, :, :]
        dlr_cloud_data.close()   # Note to self; does NetCDF have a context manager?

    def read_tffile(self, tffile):
        """Read in FRESCO cloud data:"""
        dlr_cloud_data = Dataset(tffile, mode='r')
        # Extract data of interest:
        self.tflons = dlr_cloud_data.groups['PRODUCT'].variables['longitude'][:][0, :, :]
        self.tflats = dlr_cloud_data.groups['PRODUCT'].variables['latitude'][:][0, :, :]
        tfrc = dlr_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['cloud_fraction_crb'][:]
        self.tffrc = tfrc.data[0, :, :]
        talb = dlr_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['cloud_albedo_crb'][:]
        self.tfalb = talb.data[0, :, :]
        ttop = dlr_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['cloud_pressure_crb'][:]
        self.tftop = ttop.data[0, :, :]
        tqval = dlr_cloud_data.groups['PRODUCT'].variables['qa_value'][:]
        self.tfqval = tqval.data[0, :, :]
        tsnow = dlr_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['snow_ice_flag'][:]
        self.tfsnow = tsnow.data[0, :, :]
        dlr_cloud_data.close()

    def filter_tdfile(self):
        # Convert all valid snow/ice free flag values (0,255) to 0.
        self.tdsnow = np.where(self.tdsnow == 255, 0, self.tdsnow)

        # Set missing/poor quality/irrelevant data to NAN:
        # Apply cloud fraction filter, but set it to 0.7 rather
        #     than 0.9 to account for variability around this threshold
        #     in the two cloud products.
        # inicnt=np.count_nonzero(~np.isnan(tdfrc))
        # print(inicnt)
        # tdfrc=np.where(tdfrc<0.7, np.nan, tdfrc)
        # Fill values (do for cloud fraction and cloud top pressure, as missing
        # values for cloud fraction may not be the same as missing values for
        # other data, as these are obtained with different algorithms):
        # necessarily the same as the fill values for the other data):
        self.tdfrc = np.where(self.tdfrc == FILL_VAL, np.nan, self.tdfrc)
        self.tdtop = np.where(self.tdtop == FILL_VAL, np.nan, self.tdtop)
        # Set cloud fraciton to nan for scenes with missing cloud top
        # pressure data:
        self.tdfrc = np.where(np.isnan(self.tdtop), np.nan, self.tdfrc)
        # Poor data quality:
        self.tdfrc = np.where(self.tdqval < 0.5, np.nan, self.tdfrc)
        # Snow/ice cover:
        self.tdfrc = np.where(self.tdsnow != 0, np.nan, self.tdfrc)
        # Apply filter to remaining data:
        self.tdtop = np.where(np.isnan(self.tdfrc), np.nan, self.tdtop)
        self.tdbase = np.where(np.isnan(self.tdfrc), np.nan, self.tdbase)
        self.tdoptd = np.where(np.isnan(self.tdfrc), np.nan, self.tdoptd)
        # Error check:
        if np.nanmax(self.tdtop) == FILL_VAL:
            raise Exception('Not all missing values converted to NANs')
        if np.nanmax(self.tdbase) == FILL_VAL:
            raise Exception('Not all missing values converted to NANs')
        if np.nanmax(self.tdoptd) == FILL_VAL:
            raise Exception('Not all missing values converted to NANs')

    def filter_tffile(self):
        # Convert all valid snow/ice free flag values (0,255) to 0.
        self.tfsnow = np.where(self.tfsnow == 255, 0, self.tfsnow)

        # Set missing/poor quality/irrelevant data to NAN:
        # Apply cloud fraction filter, but set it to 0.7 rather
        #     than 0.9 to account for variability around this threshold
        #     in the two cloud products.
        # inicnt=np.count_nonzero(~np.isnan(tffrc))
        # print(inicnt)
        self.tffrc = np.where(self.tffrc < 0.7, np.nan, self.tffrc)
        self.tffrc = np.where(self.tffrc == FILL_VAL, np.nan, self.tffrc)
        # QA Flags. Threshold of 0.45 suggested by Henk Eskes in email
        # exchange on 18 Jan 2020:
        self.tffrc = np.where(self.tfqval < 0.45, np.nan, self.tffrc)
        # Apply cloud top pressure filter to only consider clouds above
        #   500 hPa and below 150 hPa (more generous than the 450-200 hPa
        #   range to account for variability around this threshold in the
        #   two cloud products:
        self.tffrc = np.where(self.tftop > 450 * 1e2, np.nan, self.tffrc)
        self.tffrc = np.where(self.tftop < 180 * 1e2, np.nan, self.tffrc)
        # Snow/ice cover:
        self.tffrc = np.where(self.tfsnow != 0, np.nan, self.tffrc)
        # Apply filter to remaining data:
        self.tftop = np.where(np.isnan(self.tffrc), np.nan, self.tftop)

    def check_parity(self,dorb):
        # Skip files if the number of indices are not equal:
        if self.tdlons.shape != self.tflons.shape:
            print('Indices ne for ', dorb)
            return  # CONTINUE
            # m=min(md,mf)
            # n=min(nd,nf)

    def get_nobs(self):
        dlr_ind = np.where((self.tdqval < 0.5) & (self.tdfrc >= 0.7) & (self.tdtop >= 18000)
                          & (self.tdtop <= 45000) & (self.tdsnow == 0))[0]
        fr_ind = np.where((self.tfqval < 0.45) & (self.tffrc >= 0.7) & (self.tftop >= 18000)
                          & (self.tftop <= 45000) & (self.tfsnow == 0))[0]
        return len(dlr_ind), len(fr_ind)


def process_file(tdfile, tffile, running_total_container):
    """Processes a paired dtr and fresco product, and adds the pixels to the running total in the
    running_total_container"""
    # Track progress:
    print('===> Processing: ', tdfile)
    file_data_container = TropomiData(tdfile, tffile)

    # REGRID THE DATA:
    for i in range(file_data_container.shape[0]):
        for j in range(file_data_container.shape[1]):
            running_total_container.update_pixel(file_data_container, i, j)
    running_total_container.update_nobs(file_data_container)
    return


def get_files_for_month(sen_5_p_dir, month_index):
    """
    For a given month index (jan-may 2020 being 1-5, jun-dec 2019 being 6-12), returns every Tropomi and Fresco
    filepath. Also sets the globals StrMM, StrYY, MMName (used in the plotting method)(at least until I fix it)
    """
    # TODO Roll the string manufacturing into the CloudVariableStore class
    global StrMM, StrYY, MMName


    # Input parameter (to selet month of interest):
    # Define month of interest as string 2 characters in length:
    # TODO: Change this entire section into datetime
    StrMM = str(month_index)
    # Define string of year and first 3 letters of month name based on above entry:
    if StrMM == '01': StrYY, MMName = '2020', 'jan'
    if StrMM == '02': StrYY, MMName = '2020', 'feb'
    if StrMM == '03': StrYY, MMName = '2020', 'mar'
    if StrMM == '04': StrYY, MMName = '2020', 'apr'
    if StrMM == '05': StrYY, MMName = '2020', 'may'
    if StrMM == '06': StrYY, MMName = '2019', 'jun'
    if StrMM == '07': StrYY, MMName = '2019', 'jul'
    if StrMM == '08': StrYY, MMName = '2019', 'aug'
    if StrMM == '09': StrYY, MMName = '2019', 'sep'
    if StrMM == '10': StrYY, MMName = '2019', 'oct'
    if StrMM == '11': StrYY, MMName = '2019', 'nov'
    if StrMM == '12': StrYY, MMName = '2019', 'dec'
    # Define days:
    ndays = 31
    dd = list(range(1, ndays + 1))
    strdd = [''] * ndays
    cnt = 0
    for d in dd:
        strdd[cnt] = '0' + str(d) if d < 10 else str(d)
        cnt = cnt + 1

    # dir structure specific to HPC, filename from TROPOMI
    # Get DLR data file names:
    td_file_list = glob.glob(
        path.join(sen_5_p_dir,
                  'CLOUD_OFFL',
                  StrYY,
                  StrMM,
                  'S5P_OFFL_L2__CLOUD__' + StrYY + StrMM + '*'))
    td_file_list = sorted(td_file_list)
    # Get FRESCO file names:
    tf_file_list = glob.glob(
        path.join(sen_5_p_dir,
                  'NO2_OFFL',
                  StrYY,
                  StrMM ,
                  'S5P_OFFL_L2__NO2____' + StrYY + StrMM + '*'))
    tf_file_list = sorted(tf_file_list)
    # Check that number of files are equal. If not, exit the programme:
    if len(td_file_list) != len(tf_file_list):
        print('DLR files = ', len(td_file_list))
        print('FRESCO files = ', len(tf_file_list))
        print('unequal number of files')
        raise FileMismatchException("Unequal number of DLR and FRESCO files, check archive.")
    if len(td_file_list) == 0:
        raise FileNotFoundError("No files found in {}".format(sen_5_p_dir))

    return td_file_list, tf_file_list


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extracts and plots cloud data")
    parser.add_argument("--s5p_data_dir", default='/data/uptrop/nobackup/tropomi/Data/')
    parser.add_argument("--output_dir", default = '~/eos_library/cloud_test_output')
    parser.add_argument("--month", default="10")
    parser.add_argument("--plot_dir", default = "~/eos_library/cloud_test_plots")
    args = parser.parse_args()

    # TODO: Make this a member of CloudVariableStore
    out_lon=np.arange(MIN_LON, MAX_LON, DELTA_LON)
    out_lat=np.arange(MIN_LAT, MAX_LAT, DELTA_LAT)
    # Convert output lats and long to 2D:
    X, Y = np.meshgrid(out_lon, out_lat, indexing='ij')

    td_file_list, tf_file_list = get_files_for_month(args.s5p_data_dir, args.month)
    running_cloud_total = CloudVariableStore(X.shape)

    # Loop over files:
    for td_file, tf_file in zip(td_file_list, tf_file_list):
        print("Processing files {} and {}".format(tf_file, td_file))
        process_file(td_file, tf_file, running_cloud_total)

    running_cloud_total.calc_cloud_statistics()

    # Print number of observations to screen:
    print('No. of FRESCO obs for '+MMName+' = ',running_cloud_total.nobs_fresco)
    print('No. of DLR obs for '+MMName+' = ',running_cloud_total.nobs_dlr)
    print("Writing to NetCDF at {}".format(args.output_dir))
    running_cloud_total.write_to_netcdf(args.output_dir)
    print("Creating plots at {}".format(args.plot_dir))
    running_cloud_total.plot_clouds_products(args.plot_dir)
    print("Processing complete.")
