#!/usr/bin/python

'''Code to compare TROPOMI and Pandora column NO2 at high altitude
   sites to assess skill of TROPOMI at reproducing Pandora observations
   of free tropospheric NO2. 

   Code is set up to process Pandora total or tropospheric column NO2
   at the Mauna Loa, Izana, or Altzomoni sites.
   '''

# Import relevant packages:
import glob
import sys
import os
import netCDF4 as nc4
from netCDF4 import Dataset
import numpy as np
from bootstrap import rma

from read_pandora import readpandora
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from constants import DU_TO_MOLECULES_PER_CM2 as du2moleccm2

# Turn off warnings:
np.warnings.filterwarnings('ignore')

# Define final array of coincident data for each day at Pandora site:
# JR: This can be changed to produce an appendable array, as the size
#     will vary depending on which site is being processed.
nvals = 366
pan_ml = np.zeros(nvals)
s5p_ml = np.zeros(nvals)
s5p_ch = np.zeros(nvals)
s5p_cf = np.zeros(nvals)
pan_wgt = np.zeros(nvals)
s5p_wgt = np.zeros(nvals)
pan_cnt = np.zeros(nvals)
s5p_cnt = np.zeros(nvals)
daycnt = 0

# TODO: Change this to datetime. I say this in every single script, but never get around to it....
# Define strings of months:
StrMon = ['06', '07', '08', '09', '10', '11', '12', '01', '02', '03', '04', '05']
NMon = len(StrMon)
Month = [6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5]
StrYY = ['19', '19', '19', '19', '19', '19', '19', '20', '20', '20', '20', '20']
Year = [2019, 2019, 2019, 2019, 2019, 2019, 2019, 2020, 2020, 2020, 2020, 2020]
DayInMon = [30, 31, 31, 30, 31, 30, 31, 31, 29, 31, 30, 31]

# Arguments to input:
# Apply bias correction?
apply_bias_correction = False

# Define column comparing to:
no2col = 'Tot'  # Either Tot or Trop; default it Trop

# Define cloud product to be used:
cldprd = 'fresco'  # options are fresco, dlr-ocra; default is fresco

# Define Pandora measurement site to compare to TROPOMI:
SSite = 'mauna_loa'  # options are izana,mauna_loa,altzomoni; default is izana

# Define degree range to sample TROPOMI around Pandora site:
strdiffdeg = '02'  # options are: 0.3,0.2,0.1,0.05; default is 0.2 ("02")


class TropomiData:
    def __init__(self, filepath):
        # Read file:
        fh = Dataset(filepath, mode='r')
        # Extract data of interest (lon, lat, clouds, NO2 total column & error):
        glons = fh.groups['PRODUCT'].variables['longitude'][:]
        self.tlons = glons.data[0, :, :]
        glats = fh.groups['PRODUCT'].variables['latitude'][:]
        self.tlats = glats.data[0, :, :]
        self.xdim = len(self.tlats[:, 0])
        self.ydim = len(self.tlats[0, :])
        # Factor to convert from mol/m3 to molecules/cm2:
        self.no2sfac = fh.groups['PRODUCT']. \
            variables['nitrogendioxide_tropospheric' \
                      '_column'].multiplication_factor_to_convert_to_molecules_percm2
        # Get delta-time (along x index):
        gdtime = fh.groups['PRODUCT'].variables['delta_time'][:]
        self.tdtime = gdtime.data[0, :]
        # Get start (reference time):
        greftime = fh.groups['PRODUCT'].variables['time_utc'][:]
        self.treftime = greftime[0, :]
        # Extract UTC hours and minutes:
        gomi_dd = [x[8:10] for x in self.treftime]
        gomi_utc_hh = [x[11:13] for x in self.treftime]
        gomi_min = [x[14:16] for x in self.treftime]
        gomi_utc_hh = [int(i) for i in gomi_utc_hh]
        gomi_min = [int(i) for i in gomi_min]
        gomi_dd = [int(i) for i in gomi_dd]
        # Convert time from 1D to 2D:
        self.tomi_min = np.zeros((self.xdim, self.ydim))
        self.tomi_utc_hh = np.zeros((self.xdim, self.ydim))
        self.tomi_dd = np.zeros((self.xdim, self.ydim))
        for i in range(self.xdim):
            self.tomi_min[i, :] = gomi_min[i]
            self.tomi_utc_hh[i, :] = gomi_utc_hh[i]
            self.tomi_dd[i, :] = gomi_dd[i]
        # Get QA flag scale factor:
        self.qasfac = fh.groups['PRODUCT'].variables['qa_value'].scale_factor
        # QA value:
        self.qaval = fh.groups['PRODUCT'].variables['qa_value'][0, :, :]
        # NO2 fill/missing value:
        self.fillval = fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column']._FillValue
        # Total vertical column NO2 column:
        self.gtotno2 = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].variables['nitrogendioxide_total_column'][:]
        # Preserve in case use in future:
        # gtotno2=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
        #         variables['nitrogendioxide_summed_total_column'][:]
        self.ttotno2 = self.gtotno2.data[0, :, :]
        # Total slant column:
        gscdno2 = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].variables[
                      'nitrogendioxide_slant_column_density'][:]
        self.tscdno2 = gscdno2.data[0, :, :]
        # Precision of total slant column:
        gscdno2err = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'] \
                         .variables['nitrogendioxide_slant_column_density_''precision'][:]
        self.tscdno2err = gscdno2err.data[0, :, :]
        # Tropospheric vertical column :
        gtropno2 = fh.groups['PRODUCT'].variables['nitrogendioxide_' \
                                                  'tropospheric_column'][:]
        self.ttropno2 = gtropno2.data[0, :, :]
        # Summed column precision:
        # Preserve in case use in future:
        # ttotno2err=fh.groups['PRODUCT']['SUPPORT_DATA']\
        #            ['DETAILED_RESULTS'].\
        #            variables['nitrogendioxide_summed_total_column_'\
        #                      'precision'][0,:,:]
        # Tropospheric column:
        self.ttropno2err = fh.groups['PRODUCT'].variables['nitrogendioxide_' \
                                                     'tropospheric_column_' \
                                                     'precision'][0, :, :]
        # Total columnn:
        self.ttotno2err = fh.groups['PRODUCT']['SUPPORT_DATA'] \
                         ['DETAILED_RESULTS']. \
                         variables['nitrogendioxide_total_column_precision'] \
            [0, :, :]
        # Statospheric column:
        gstratno2 = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                        variables['nitrogendioxide_stratospheric_column'][:]
        self.tstratno2 = gstratno2.data[0, :, :]
        # Statospheric column error:
        self.stratno2err = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                          variables['nitrogendioxide_stratospheric_column_precision'][0, :, :]
        # Surface pressure:
        gsurfp = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                     variables['surface_pressure'][:]
        self.tsurfp = gsurfp.data[0, :, :]
        # Solar zenith angle (degrees):
        tsza = fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']. \
                   variables['solar_zenith_angle'][:]
        self.sza = tsza[0, :, :]
        # Viewing zenith angle (degrees):
        tvza = fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']. \
                   variables['viewing_zenith_angle'][:]
        self.vza = tvza[0, :, :]
        # Stratospheric AMF:
        gstratamf = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                        variables['air_mass_factor_stratosphere'][:]
        self.tstratamf = gstratamf.data[0, :, :]
        fh.close()

    def preprocess(self):
        # Calculate the geometric AMF:
        self.tamf_geo = np.add((np.reciprocal(np.cos(np.deg2rad(self.sza)))), \
                          (np.reciprocal(np.cos(np.deg2rad(self.vza)))))
        # Calculate the total column with a geometric AMF:
        # Step 1: calculate stratospheric SCD (not in data product):
        self.tscdstrat = np.multiply(self.tstratno2, self.tstratamf)
        # Step 2: calculate tropospheric NO2 SCD:
        self.ttropscd = np.subtract(self.tscdno2, self.tscdstrat)
        # Step 3: calculate tropospheric NO2 VCD:
        self.tgeotropvcd = np.divide(self.ttropscd, self.amf_geo)
        if (apply_bias_correction):
            # Correct for bias in the tropospheric column based on
            # assessment of TROPOMI with Pandora
            self.tgeotropvcd = np.divide(self.tgeotropvcd, 2.0)
            # Apply bias correction to stratosphere (underestimated by 10%):
            self.tstratno2 = np.multiply(self.tstratno2, 0.9)
        # Step 4: sum up stratospheric and tropospheric NO2 VCDs:
        self.tgeototvcd = np.add(self.tgeotropvcd, self.tstratno2)
        # Calculate total VCD column error by adding in quadrature
        # individual contributions:
        self.ttotvcd_geo_err = np.sqrt(np.add(np.square(self.stratno2err), \
                                         np.square(self.tscdno2err)))
        # Estimate the tropospheric NO2 error as the total error
        # weighted by the relative contribution of the troposphere
        # to the total column. This can be done as components that
        # contribute to the error are the same:
        self.ttropvcd_geo_err = np.multiply(self.ttotvcd_geo_err, \
                                       (np.divide(self.tgeotropvcd, self.tgeototvcd)))

    def apply_cloud_filter(self, no2col, cloud_product):

        # Select which NO2 data to use based on no2col selection:
        if (no2col == 'Tot'):
            tno2val = self.tgeototvcd
            tno2err = self.ttotvcd_geo_err
        elif (no2col == 'Trop'):
            tno2val = self.tgeotropvcd
            tno2err = self.ttropvcd_geo_err
            stratcol = self.tstratno2
            totcol = self.tgeototvcd
        else:
            raise BadNo2ColException
        # Account for files where mask is missing (only appears to be one):
        if len(self.gtotno2.mask.shape) == 0:
            tno2val = np.where(tno2val == self.fillval, np.nan, tno2val)
        else:
            tno2val[self.gtotno2.mask[0, :, :]] = float("nan")
        # Find relevant data only:
        # Filter out low quality retrieval scenes (0.45 suggested
        # by Henk Eskes at KNMI):
        tno2val = np.where(self.qaval < 0.45, np.nan, tno2val)

        # Also set scenes with snow/ice to nan. Not likely for the tropical
        # sites selected for this comparison, but included this here in
        # case of future comparisons that in midlatitudes or poles:
        tno2val = np.where(cloud_product.tsnow != 0, np.nan, tno2val)
        # Convert NO2 from mol/m3 to molec/cm2:
        self.tno2val = np.multiply(tno2val, self.no2sfac)
        self.tno2err = np.multiply(tno2err, self.no2sfac)
        # Trim to remove data where relevant NO2 data is not NAN:
        lons = self.tlons[~np.isnan(tno2val)]
        lats = self.tlats[~np.isnan(tno2val)]
        self.no2err = tno2err[~np.isnan(tno2val)]
        self.omi_utc_hh = self.tomi_utc_hh[~np.isnan(tno2val)]
        self.omi_min = self.tomi_min[~np.isnan(tno2val)]
        self.omi_dd = self.tomi_dd[~np.isnan(tno2val)]
        self.cldfrac = cloud_product.tcldfrac[~np.isnan(tno2val)]
        self.cldpres = cloud_product.tcldpres[~np.isnan(tno2val)]
        self.no2val = tno2val[~np.isnan(tno2val)]
        if (no2col == 'Trop'):
            self.stratcol = stratcol[~np.isnan(tno2val)]
            self.totcol = totcol[~np.isnan(tno2val)]
        # Combine hour and minute into xx.xx format:
        self.tomi_hhmm = self.omi_utc_hh + np.divide(self.omi_min, 60.)


class CloudData:
    def __init__(self, filepath, product_type):
        if product_type == "dlr-ocra":
            self.read_ocra_data(filepath)
        elif product_type == "fresco":
            self.read_fresco_data(filepath)

    def read_ocra_data(self,filepath):
        # Read data:
        fh = Dataset(filepath, mode='r')
        # TODO: Watch out for those string indexes. Change when format is understood.
        # Check that date is the same as the no2 file:
        strdate = filepath[-66:-51]
        # TODO: Move check elsewhere
        if strdate != tomi_files_on_day[-66:-51]:
            print('NO2 file, Cloud file: ' + strdate+ ", " + \
                  strdate, flush=True)
            print('EXITING: Files are not for the same date!', \
                  flush=True)
            sys.exit()
        # Get cloud fraction and cloud top pressure:
        gcldfrac = fh.groups['PRODUCT'].variables['cloud_fraction'][:]
        self.tcldfrac = gcldfrac.data[0, :, :]
        gcldpres = fh.groups['PRODUCT'].variables['cloud_top_pressure'][:]
        self.tcldpres = np.ma.getdata(gcldpres[0, :, :])  # extract data from masked array
        # QA value:
        self.cldqa = fh.groups['PRODUCT'].variables['qa_value'][0, :, :]
        # Snow/ice flag:
        self.gsnow = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                    variables['snow_ice_flag'][:]
        self.tsnow = self.gsnow.data[0, :, :]
        # Set poor quality cloud data to nan:
        self.tcldfrac = np.where(self.cldqa < 0.5, np.nan, self.tcldfrac)
        self.tcldpres = np.where(self.cldqa < 0.5, np.nan, self.tcldpres)
        # Set clouds over snow/ice scenes to nan:
        self.tcldfrac = np.where(self.tsnow != 0, np.nan, self.tcldfrac)
        self.tcldpres = np.where(self.tsnow != 0, np.nan, self.tcldpres)
        # TODO: Move check elsewhere
        # Check that data shapes are equal:
        if self.tcldfrac.shape != sza.shape:
            print('Cloud product and NO2 indices ne!', flush=True)
            print(tcldfrac.shape, sza.shape, flush=True)
            print('Skipping this swath', flush=True)
            fh.close()
            raise BadCloudShapeException
        # Close file:
        fh.close()

    def read_fresco_data(self,filepath):
        fh = Dataset(filepath, mode='r')
        # Cloud input data (cldfrac, cldalb, cldpres):
        gcldfrac = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                       variables['cloud_fraction_crb'][:]
        self.tcldfrac = gcldfrac.data[0, :, :]
        gcldpres = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                       variables['cloud_pressure_crb'][:]
        self.tcldpres = np.ma.getdata(gcldpres[0, :, :])  #
        # Snow/ice flag:
        gsnow = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                    variables['snow_ice_flag'][:]
        # Apparent scene pressure:
        gscenep = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                      variables['apparent_scene_pressure'][:]
        self.tscenep = gscenep.data[0, :, :]
        self.tsnow = gsnow.data[0, :, :]
        # Convert all valid snow/ice free flag values (252,255) to 0.
        # Ocean values:
        self.tsnow = np.where(self.tsnow == 255, 0, self.tsnow)
        # Coastline values (listed as potential "suspect" in the ATBD
        # document (page 67):
        self.tsnow = np.where(self.tsnow == 252, 0, self.tsnow)
        # Less then 1% snow/ice cover:
        self.tsnow = np.where(self.tsnow < 1, 0, self.tsnow)
        # Snow/ice misclassified as clouds:
        self.tsnow = np.where(((self.tsnow > 80) & (self.tsnow < 104)
                               & (self.tscenep > (0.98 * self.tsurfp))),
                                0, self.tsnow)
        # Set clouds over snow/ice scenes to nan:
        self.tcldfrac = np.where(self.tsnow != 0, np.nan, self.tcldfrac)
        self.tcldpres = np.where(self.tsnow != 0, np.nan, self.tcldpres)
        # close file:
        fh.close()


def plot_data():
    global days
    # Plot time series:
    plt.figure(1, figsize=(10, 5))
    x = np.arange(0, daycnt, 1)
    days = x
    plt.errorbar(x, pan_ml[0:daycnt] * 1e-14, yerr=pan_wgt[0:daycnt] * 1e-14, \
                 fmt='.k', color='black', capsize=5, capthick=2, \
                 ecolor='black', markersize=12, label='Pandora')
    plt.errorbar(x, s5p_ml[0:daycnt] * 1e-14, yerr=s5p_wgt[0:daycnt] * 1e-14, \
                 fmt='.k', color='blue', capsize=5, capthick=2, \
                 ecolor='blue', markeredgecolor='blue', \
                 markerfacecolor='blue', markersize=12, label='TROPOMI')
    plt.ylim(ymin, ymax)
    plt.xlabel('Days since 1 June 2019')
    plt.ylabel('$NO_2$ total VCD [$10^{14}$ molecules $cm^2$]')
    leg = plt.legend(loc='lower left', fontsize='large')
    leg.get_frame().set_linewidth(0.0)
    # plt.savefig('./Images/tropomi-'+SSite+'-pandora-no2-timeseries-v1-jun2019-apr2020.ps', \
    #            format='ps',transparent=True,bbox_inches='tight',dpi=100)
    # Plot scatterplot:
    tx = pan_ml[0:daycnt]
    ty = s5p_ml[0:daycnt]
    nas = np.logical_or(np.isnan(tx), np.isnan(ty))
    print('No. of coincident points = ', len(tx[~nas]))
    r = stats.pearsonr(tx[~nas], ty[~nas])
    print('Correlation: ', r[0])
    # Get mean difference:
    Diff = np.subtract(np.mean(ty[~nas]), np.mean(tx[~nas]))
    print('TROPOMI minus Pandora (10^14) = ', Diff * 1e-14)
    NMB = 100. * np.divide(Diff, np.mean(tx[~nas]))
    print('TROPOMI NMB (%) = ', NMB)
    # RMA regression:
    result = rma(tx[~nas] * 1e-14, ty[~nas] * 1e-14, len(tx[~nas]), 10000)
    print('Intercept (10^14): ', result[1])
    print('Slope: ', result[0])
    fig = plt.figure(2)
    plt.figure(2, figsize=(6, 5))
    ax = fig.add_subplot(1, 1, 1)
    plt.plot(1e-14 * tx, 1e-14 * ty, 'o', color='black')
    plt.xlim(0, 60)
    plt.ylim(0, 60)
    plt.xlabel('Pandora $NO_2$ total VCD [$10^{14}$ molecules $cm^2$]')
    plt.ylabel('TROPOMI $NO_2$ total VCD [$10^{14}$ molecules $cm^2$]')
    xvals = np.arange(0, 60, 2)
    yvals = result[1] + xvals * result[0]
    plt.plot(xvals, yvals, '-')
    add2plt = ("y = {a:.3f}x + {b:.3f}". \
               format(a=result[0], b=result[1]))
    plt.text(0.1, 0.9, add2plt, fontsize=10, \
             ha='left', va='center', transform=ax.transAxes)
    add2plt = ("r = {a:.3f}".format(a=r[0]))
    plt.text(0.1, 0.84, add2plt, fontsize=10, \
             ha='left', va='center', transform=ax.transAxes)
    # plt.savefig('./Images/tropomi-'+SSite+'-pandora-no2-scatterplot-v1-jun2019-apr2020.ps', \
    #            format='ps',transparent=True,bbox_inches='tight',dpi=100)
    plt.show()


def write_to_netcdf():
    global time
    # Save the data to NetCDF:
    ncout = Dataset('./Data/tropomi-pandora-comparison-' + SSite + '-' + cldprd + \
                    '-' + no2col + '-' + strdiffdeg + 'deg-bias-corr-v1.nc', mode='w', format='NETCDF4')
    # Set array sizes:
    TDim = daycnt
    ncout.createDimension('time', TDim)
    # create time axis
    time = ncout.createVariable('time', np.float32, ('time',))
    time.units = 'days since 2019-06-01'
    time.long_name = 'time in days since 2019-06-01'
    time[:] = days
    panno2 = ncout.createVariable('panno2', np.float32, ('time',))
    panno2.units = 'molecules/cm2'
    panno2.long_name = 'Pandora error-weighted daily mean total column NO2 coincident with TROPOMI overpass'
    panno2[:] = pan_ml[0:daycnt]
    panerr = ncout.createVariable('panerr', np.float32, ('time',))
    panerr.units = 'molecules/cm2'
    panerr.long_name = 'Pandora weighted error of daily mean total columns of NO2 coincident with TROPOMI overpass'
    panerr[:] = pan_wgt[0:daycnt]
    pancnt = ncout.createVariable('pancnt', np.float32, ('time',))
    pancnt.units = 'unitless'
    pancnt.long_name = 'Number of Pandora observations used to obtain weighted mean'
    pancnt[:] = pan_cnt[0:daycnt]
    satno2 = ncout.createVariable('satno2', np.float32, ('time',))
    satno2.units = 'molecules/cm2'
    satno2.long_name = 'S5P/TROPOMI NO2 OFFL error-weighted daily mean total column NO2 coincident with Pandora'
    satno2[:] = s5p_ml[0:daycnt]
    satcldh = ncout.createVariable('satcldh', np.float32, ('time',))
    satcldh.units = 'hPa'
    satcldh.long_name = 'S5P/TROPOMI mean cloud top pressure at Pandora site'
    satcldh[:] = s5p_ch[0:daycnt]
    satcldf = ncout.createVariable('satcldf', np.float32, ('time',))
    satcldf.units = 'hPa'
    satcldf.long_name = 'S5P/TROPOMI mean cloud fraction at Pandora site'
    satcldf[:] = s5p_cf[0:daycnt]
    saterr = ncout.createVariable('saterr', np.float32, ('time',))
    saterr.units = 'molecules/cm2'
    saterr.long_name = 'S5P/TROPOMI NO2 OFFL weighted error of daily mean total columns of NO2 coincident with the Pandora site'
    saterr[:] = s5p_wgt[0:daycnt]
    satcnt = ncout.createVariable('satcnt', np.float32, ('time',))
    satcnt.units = 'unitless'
    satcnt.long_name = 'Number of S5P/TROPOMI observations used to obtain weighted mean'
    satcnt[:] = s5p_cnt[0:daycnt]
    ncout.close()


def daily_weighted_means():
    # Get daily error-weighted means:
    pan_ml[0:daycnt] = pan_ml[0:daycnt] / pan_wgt[0:daycnt]
    pan_wgt[0:daycnt] = np.divide(1, np.sqrt(pan_wgt[0:daycnt]))
    s5p_ml[0:daycnt] = s5p_ml[0:daycnt] / s5p_wgt[0:daycnt]
    s5p_ch[0:daycnt] = s5p_ch[0:daycnt] / s5p_cnt[0:daycnt]
    s5p_cf[0:daycnt] = s5p_cf[0:daycnt] / s5p_cnt[0:daycnt]
    s5p_wgt[0:daycnt] = np.divide(1, np.sqrt(s5p_wgt[0:daycnt]))
    print('Min & max relative errors (Pandora): ', np.nanmin(np.divide(pan_wgt, pan_ml)), \
          np.nanmax(np.divide(pan_wgt, pan_ml)))
    print('Min & max relative errors (TROPOMI): ', np.nanmin(np.divide(s5p_wgt, s5p_ml)), \
          np.nanmax(np.divide(s5p_wgt, s5p_ml)))


def add_pandora_to_day():
    # Find relevant Pandora data for this year, month and day:
    # Pandora flag threshold selected is from https://www.atmos-meas-tech.net/13/205/2020/amt-13-205-2020.pdf
    panind = np.argwhere((panyy == Year[m]) & (panmon == Month[m]) \
                         & (pandd == (d + 1)) & (panno2 > -8e99) \
                         & (panqaflag <= 11) \
                         & (panqaflag != 2) & (pan_hhmm >= hhsite[n] - 0.5) \
                         & (pan_hhmm <= hhsite[n] + 0.5))
    # Proceed if there are Pandora data points:
    if len(panind) == 0:
        raise PanIndException
    # Create arrays of relevant data and convert from DU to molec/cm2:
    tno2 = np.multiply(panno2[panind], du2moleccm2)
    terr = np.multiply(panno2err[panind], du2moleccm2)
    tqa = panqaflag[panind]
    # Add Pandora total NO2 to final array:
    for w in range(len(panind)):
        pan_ml[daycnt] += np.divide(tno2[w], np.square(terr[w]))
        pan_wgt[daycnt] += np.divide(1.0, np.square(terr[w]))
        pan_cnt[daycnt] += len(panind)




def add_trop_data_to_total():
    global tomiind
    # Find coincident data for this file:
    difflon = abs(np.subtract(lons, panlon))
    difflat = abs(np.subtract(lats, panlat))
    # Use distanc (degrees) to find coincident data.
    # For Pandora 'Trop' data, only consider TROPOMI scenes where the
    # total column exceeds the stratospheric column:
    if (no2col == 'Tot'):
        tomiind = np.argwhere((difflon <= diffdeg) & (difflat <= diffdeg) \
                              & (no2val != np.nan) & (omi_dd == (d + 1)))
    if (no2col == 'Trop'):
        tomiind = np.argwhere((difflon <= diffdeg) & (difflat <= diffdeg) \
                              & (no2val != np.nan) & (omi_dd == (d + 1)) \
                              & (stratcol < totcol))
    # Skip if no data:
    if (len(tomiind) == 0):
        raise NoTropDataException
    # TODO: Here's the total-ing object
    # Add TROPOMI total NO2 to final array of daily means:
    s5p_ml[daycnt] = s5p_ml[daycnt] + \
                     sum(np.divide(no2val[tomiind], \
                                   np.square(no2err[tomiind])))
    s5p_wgt[daycnt] = s5p_wgt[daycnt] + \
                      sum(np.divide(1.0, \
                                    np.square(no2err[tomiind])))
    s5p_ch[daycnt] = s5p_ch[daycnt] + sum(cldpres[tomiind] * 1e-2)
    s5p_cf[daycnt] = s5p_cf[daycnt] + sum(cldfrac[tomiind])
    s5p_cnt[daycnt] = s5p_cnt[daycnt] + len(tomiind)





def get_nhrs():
    global hhsite
    # Get min and max TROPOMI UTC for this orbit:
    # Choose min and max time window of TROPOMI 0.2 degrees
    # around Pandora site:
    minhh = np.nanmin(omi_utc_hh[tomiind])
    maxhh = np.nanmax(omi_utc_hh[tomiind])
    mintime = np.nanmin(tomi_hhmm[tomiind])
    maxtime = np.nanmax(tomi_hhmm[tomiind])
    if (minhh == maxhh):
        hhsite = [mintime]
    else:
        hhsite = [mintime, maxtime]
    nhrs = len(hhsite)
    return nhrs


def get_tropomi_files_on_day(d):
    # Get string of day:
    StrDD = str(d + 1)
    if (d + 1) <= 9:
        StrDD = '0' + StrDD
    # Get string of files for this day:
    tomi_files_on_day = glob.glob(tomidir + 'NO2_OFFL/20' + StrYY[m] + '/' + StrMon + '/' + \
                                  'S5P_OFFL_L2__NO2____20' + StrYY[m] + StrMon + StrDD + '*')
    # Track progress:
    print('Processing day in month: ', StrDD)
    # Order the files:
    tomi_files_on_day = sorted(tomi_files_on_day)
    return tomi_files_on_day


def get_ocra_files_on_day():
    # Get string of day:
    StrDD = str(d + 1)
    if (d + 1) <= 9:
        StrDD = '0' + StrDD
    clddir = tomidir + 'CLOUD_OFFL/20' + StrYY[m] + '/'
    cldfile = glob.glob(clddir + StrMon + '/' + 'S5P_OFFL_L2__CLOUD__20' + \
                        StrYY[m] + StrMon + StrDD + '*')
    # Order the files:
    cldfile = sorted(cldfile)
    return cldfile


def read_pandora_data():
    global panlat, panlon, panyy, panmon, pandd, pan_hhmm, panno2, panno2err, panqaflag
    # Read Pdndora data from external function:
    p = readpandora(panfile[0])
    # Extract latitude and longitude:
    loc = p[0]
    panlat = loc['lat']
    panlon = loc['lon']
    # Extract data frame with relevant Pandora data:
    df = p[1]
    # Get variables names from column headers:
    varnames = df.columns.values
    # Rename Pandora data:
    panyy = df.year.values
    panmon = df.month.values
    pandd = df.day.values
    panhh_utc = df.hour_utc.values
    panmin = df.minute.values
    # Combine hour and minute into xx.xx format:
    pan_hhmm = panhh_utc + np.divide(panmin, 60.)
    # Change data at the date line (0-2 UTC) to (24-26 UTC) to aid sampling 30
    # minutes around the satellite overpass time at Mauna Loa. This won't
    # affect sampling over Izana, as it's at about 12 UTC.
    sind = np.argwhere((pan_hhmm >= 0.) & (pan_hhmm < 2.))
    pan_hhmm[sind] = pan_hhmm[sind] + 24.
    panjday = df.jday.values
    pansza = df.sza.values
    panno2 = df.no2.values
    panno2err = df.no2err.values
    panqaflag = df.qaflag.values
    panfitflag = df.fitflag.values
    # Create pseudo v1.8 data by decreasing Pandora column value and error by 90%.
    # Recommendation by Alexander Cede (email exchange) to account for lower
    # reference temperature at these sites that will be used in the future v1.8
    # retrieval rather than 254K used for sites that extend to the surface.
    # V1.8 data will be available in late 2020.
    panno2 = panno2 * 0.9
    panno2err = panno2err * 0.9
    # Get data length (i.e., length of each row):
    npanpnts = len(df)
    # Confirm processing correct site:
    print('Pandora Site: ', SSite)


if __name__ == "__main__":


    # Set degree range based on string entry.
    if ( strdiffdeg=='02' ):
        diffdeg=0.2
    if ( strdiffdeg=='03' ):
        diffdeg=0.3
    if ( strdiffdeg=='01' ):
        diffdeg=0.1
    if ( strdiffdeg=='005' ):
        diffdeg=0.05

    # Get Pandora site number:
    if ( SSite=='altzomoni' ):
        SiteNum='65'
        CSite='Altzomoni'
    if ( SSite=='izana' ):
        SiteNum='101'
        CSite='Izana'
    if ( SSite=='mauna_loa' ):
        SiteNum='59'
        CSite='MaunaLoaHI'

    # Conditions for choosing total or tropospheric column:
    if ( no2col=='Trop' ):
        fv='rnvh1p1-7'
        #maxval=3
        ymin=0
        ymax=25
    if ( no2col=='Tot'  ):
        #maxval=5
        fv='rnvs1p1-7'
        ymin=10
        ymax=50

    # Get Pandora filename (one file per site):
    pandir='/data/uptrop/nobackup/pandora/'+SSite+'/'
    panfile=glob.glob(pandir+'Pandora'+SiteNum+'s1_'+CSite+'_L2'+no2col+\
                      '_'+fv+'.txt')

    read_pandora_data()

    # Get OMI directory for this month:
    tomidir='/data/uptrop/nobackup/tropomi/Data/'

    # Loop over months:
    for m, StrMon in enumerate(StrMon):

        # Track progress:
        print('Processing month: ',StrMon)

        for daycnt, d in enumerate(range(DayInMon[m])):

            tomi_files_on_day = get_tropomi_files_on_day(d)

            # Get string of S5P TROPOMI cloud product file names:
            if cldprd=='dlr-ocra':
                cloud_files_on_day = get_ocra_files_on_day(d)

                # Check for inconsistent number of files:
                if len(orca_files_on_day) != len(tomi_files_on_day):
                    print('NO2 files = ', len(tomi_files_on_day), flush=True)
                    print('CLOUD files = ',len(orca_files_on_day),flush=True)
                    print('unequal number of files',flush=True)
                    raise UnequalFileException

            # Loop over files:
            for f, tomi_file_on_day, cloud_file_on_day \
                    in enumerate(zip(tomi_files_on_day, cloud_files_on_day)):
                trop_data = TropomiData(tomi_file_on_day)
                trop_data.preprocess()
                cloud_data = CloudData(cloud_file_on_day, cldprd)
                trop_data.apply_cloud_filter(no2col, cloud_data)
                # TODO: Got this far on Friday
                add_trop_data_to_total()
                nhrs = get_nhrs()
                # loop over TROPOMI hours at site:
                for n in range(nhrs):
                    add_pandora_to_day()

    daily_weighted_means()
    plot_data()
    write_to_netcdf()








