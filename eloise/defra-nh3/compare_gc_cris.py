#!/usr/bin/python

''' Python code to read and regrid CrIS NH3 and compare to GEOS-Chem
    Stream of consciosness code for now. Eventually break up into regrid
    CrIS and compare to GEOS-Chem'''

# import relevant packages:
import sys
import os
import glob
from os import path
import numpy as np
import netCDF4 as nc4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits import basemap    
import matplotlib as mpl
from scipy import stats

from sklearn.linear_model import TheilSenRegressor

from constants import AVOGADRO as na
from constants import MW_AIR as mmair
from constants import G as g
from bootstrap import rma

import cartopy
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

# Turn off warnings:
np.warnings.filterwarnings('ignore')

class CompareGCCris:
    '''Read in relevant CrIS NH3 data, filter for poor quality data.'''

    def __init__(self, cris_file, gc_file, out_lon, out_lat, grid_cris_nh3, grid_gc_nh3_og, grid_gc_nh3_ak, grid_gc_cld, grid_count):
        '''Read in CrIS data '''

        # Rename file names:
        self.cris_file = cris_file
        self.gc_file   = gc_file

        # Rename output lat and long:
        self.out_lon = out_lon
        self.out_lat = out_lat

        # Rename input data:
        self.grid_cris_nh3  = grid_cris_nh3
        self.grid_gc_nh3_og = grid_gc_nh3_og
        self.grid_gc_nh3_ak = grid_gc_nh3_ak
        self.grid_gc_cld    = grid_gc_cld
        self.grid_count     = grid_count

        # Read the CrIS NH3 data:
        self.read_cris_data()

        # Read the GEOS-Chem NH3 data:
        self.read_gc_data()

        # Filter the data:
        #self.filter_cris_data()

        # Aggregate data to target grid (for now, same as model):
        self.aggregate_data()

    def read_cris_data(self):
        '''Read and extract relevant variables'''

        # Read data variables:
        fh = Dataset(self.cris_file, mode='r')

        # Extract data variables of interest:
        self.tavgk = np.ma.getdata(fh.variables['avg_kernel'][:])
        self.dofs = np.ma.getdata(fh.variables['DOF'][:])

        # Day_night flag (day=1.0 and night=0.0) 
        self.tday_night = np.ma.getdata(fh.variables['Day_Night_Flag'][:])
        self.cris_lat = np.ma.getdata(fh.variables['Latitude'][:])
        self.cris_lon = np.ma.getdata(fh.variables['Longitude'][:])
        self.cris_pres = np.ma.getdata(fh.variables['pressure'][:])
        self.cris_alt  = np.ma.getdata(fh.variables['zm'][:])
        try:
            self.tqflag = np.ma.getdata(fh.variables['quality_flag'][:])
        except KeyError:
            self.tqflag = np.ma.getdata(fh.variables['Quality_Flag'][:])
        self.trvmr = np.ma.getdata(fh.variables['rvmr'][:])
        self.tsnr = np.ma.getdata(fh.variables['snr'][:])
        self.trvmr_pres = np.ma.getdata(fh.variables['rvmr_press'][:])
        self.cris_col = np.ma.getdata(fh.variables['tot_col'][:])
        self.cris_col_meas_err = np.ma.getdata(fh.variables['tot_col_meas_error'][:])
        self.cris_col_tot_err = np.ma.getdata(fh.variables['tot_col_total_error'][:])
        # A priori profile (ppmv)
        self.ta_priori = np.ma.getdata(fh.variables['xa'][:])
        # Cloud flag
        self.cld_flag = np.ma.getdata(fh.variables['Cloud_Flag'][:])
        
        # Temperature:
        # (1) Surface:
        self.tsurf  = np.ma.getdata(fh.variables['tsfc'][:])
        # (2) Atmospheric:
        self.tatmos = np.ma.getdata(fh.variables['tm'][:])

        fh.close()

        # Get the number of CrIS vertical layers:
        self.nlayers = len(self.tatmos[0,:])

    def read_gc_data(self):
        '''Read in relevant GEOS-Chem data'''

        fgc = Dataset(self.gc_file, mode='r')
        
        tlon,tlat,tgcnh3,tcldfr,tcldhgt,tadn,tbxhgt,tpedge,tdegk=\
          fgc.variables['LON'],fgc.variables['LAT'],\
          fgc.variables['IJ-AVG-S__NH3'],fgc.variables['TIME-SER__CF'],\
          fgc.variables['TIME-SER__CThgt'],fgc.variables['TIME-SER__AIRDEN'],\
          fgc.variables['BXHGHT-S__BXHEIGHT'],fgc.variables['PEDGE-S__PSURF'],\
          fgc.variables['DAO-3D-S__TMPU']
        self.gc_lon = tlon[:]
        self.gc_lat = tlat[:]
        self.gc_nh3 = tgcnh3[:]
        self.tcldfr=tcldfr[:]
        self.tcldhgt=tcldhgt[0,:,:]
        self.air_den = tadn[:]   # in molec/cm3
        self.box_hgt = tbxhgt[:]*1e2   # Converted from m to cm
        self.pres_edge = tpedge[:]
        #tdegk=tdegk[:]    # 3D temperature in deg K

        fgc.close()

    def read_uk_land_mask(self):
        '''Get UK land mask at GEOS-Chem EU nested resolution and domain'''

        mask_file = '/data/uptrop/Projects/covid-19/python/Data/uk-land-mask-025x03125.nc'

        fk = Dataset(mask_file, mode='r')

        tmask = fk.variables['uk-mask']
        uk_mask = tmask[:]
        tlon = fk.variables['lon']
        mask_lon = tlon[:]
        tlat = fk.variables['lat']
        mask_lat = tlat[:]

        uk_mask = np.transpose(uk_mask)
        #jind = np.where(mask_lat > 55)
        #print(jind)
        #print(mask_lat[jind])
        #uk_mask[:,jind] = 0

        fk.close()

        return uk_mask

    def aggregate_data(self):
        '''Isolate good quality daytime data'''

        # Daytime data only:
        self.cris_col = np.where( self.tday_night==1.0, self.cris_col, np.nan )
        # Good quality data only:
        # Using 3 and above 
        #self.cris_col = np.where( self.tqflag>=3, self.cris_col, np.nan )
        # Test effect of only considering data with DOFs >= 0.5:
        #self.cris_col = np.where( self.tsnr>=0.75, self.cris_col, np.nan )
        # Filter out cloudy scenes:
        #self.cris_col = np.where( self.cld_flag==0, self.cris_col, np.nan )
        # Test effect of only considering data with DOFs >= 0.5:
        self.cris_col = np.where( self.dofs>1., self.cris_col, np.nan )
        # Test the effect of only considering pixels with thermal 
        # contrast > -2K (as in Dammers et al., 2019):
        #self.cris_col = np.where( self.delta_t>0., self.cris_col, np.nan )
        # Calculate relative error (use measure error, could use total error)
        trel_err = np.divide( self.cris_col_meas_err, self.cris_col )
        #print(np.nanmin(trel_err),np.nanmax(trel_err))
        # Test effect of only considering data with relative error < 100%
        #self.cris_col = np.where( ((trel_err <= 0.5)&(~np.isnan(self.cris_col))), self.cris_col, np.nan )

        # Apply to all other data:
        for w in range(len(self.cris_col)):
            
            # Only consider valid data (daytime, good quality):
            if ~np.isnan(self.cris_col[w]):

                # Get thermal contrast:
                pmin_ind = np.argmax(self.cris_pres[w,:])
                #print(pmax_ind)
                pmax_ind = np.argmin(np.abs(self.cris_pres[w,:]-800.0))
                #print(self.cris_pres[w,pmax_ind],self.cris_pres[w,pmin_ind])
                #ndata = len(tatmos[:,0])
                #self.delta_t = np.zeros(ndata)
                #tdelta_t = self.tsurf[w] - np.mean(self.tatmos[w,pmin_ind:pmax_ind+1])
                tdelta_t = self.tsurf[w] - self.tatmos[w,pmin_ind]
                #tdelta_t = np.zeros(self.nlayers)
                #for l in range(nlayers):
                #    tdelta_t[:,l] = np.where(tatmos[:,l]!=-9.99500000e+02, np.subtract(tsurf,tatmos[:,l]), np.nan)
                #    ind=np.where((~np.isnan(tdelta_t[:,l])&(self.delta_t==0)))[0]
                #    if len(ind)>0:
                #        self.delta_t[ind] = tdelta_t[ind,l]

                # Skip if tdelta_t < -2 K:
                #if tdelta_t < 0.0: continue

                # Find GEOS-Chem grid square:
                p = np.argmin(abs(self.out_lon - self.cris_lon[w]))
                q = np.argmin(abs(self.out_lat - self.cris_lat[w]))
                self.p, self.q = p, q

                # Apply averaging kernel to model:
                cris_pres  = np.where( self.ta_priori[w,:]!=-999.5, self.cris_pres[w,:], np.nan)
                log_cris_avgk  = np.where( self.ta_priori[w,:]!=-999.5, self.tavgk[w,:,:], np.nan)
                #log_cris_avgk.transpose() # = np.transpose(log_cris_avgk)

                # Convert CrIS a priori from ppm to ppb:
                cris_prior = np.where( self.ta_priori[w,:]!=-999.5, self.ta_priori[w,:] * 1e3, np.nan)

                # Linearize the averaging kernel:
                cris_avgk = np.zeros((len(cris_prior),len(cris_prior)))
                #print(cris_prior.shape)
                #print(log_cris_avgk)
                for m in range(len(cris_prior)):
                    for n in range(len(cris_prior)):
                        cris_avgk[m,n] = min((cris_prior[m]/cris_prior[n]),3.0) * log_cris_avgk[m,n]

                #print(' === ')
                #print(cris_avgk)

# Linearization code from Hansen Cao at U. Boulder
#lin_avk0[idata,clevel,rlevel]=min(xa[bb[idata],rlevel]/xa[bb[idata],clevel],3.0) * avk[bb[idata],clevel,rlevel]

#lin_avk0: linearized averaging kernel
#avk: original log averaging kernel from CrIS data file
#xa: CrIS prior
#idata: pixel index of lin_avk0
#bb[idata]: pixel index of avk
#clevel, rlevel: column and row index of avk and lin_avk0
                #sys.exit()
                      
                self.apply_avgk(cris_pres, cris_avgk, cris_prior)

                # Sum up/aggregate the data of interest:
                #if ~np.isnan(self.gc_col_ak):  # and self.gc_col_ak>2e15:
                self.grid_cris_nh3[p][q]  += self.cris_col[w]
                self.grid_gc_nh3_og[p][q] += self.gc_col_og 
                self.grid_gc_nh3_ak[p][q] += self.gc_col_ak
                self.grid_count[p][q]     += 1.0
                if np.nansum( self.tcldfr[:,self.q,self.p] ) > 0.5:
                    self.grid_gc_cld[p][q] += 1.0

    def apply_avgk(self, cris_pressure, cris_avg_kern, cris_priori):
        '''Interpolate GEOS-Chem data to CrIS vertical grid and apply
           averaging kernel'''

        # Get GEOS-Chem pressure edges:
        tp_edge = self.pres_edge[:,self.q,self.p]

        # Get GEOS-Chemf NH3 for this grid:
        #self.gc_nh3[:,self.q,self.p] = self.gc_nh3[:,self.q,self.p] * 2. 
        tgc_nh3 = self.gc_nh3[0:len(tp_edge)-1,self.q,self.p]
        # Artificially increase NH3 in the bottom layer:

        # Calculate GEOS-Chem pressure mids:
        self.gc_pmid = np.zeros(len(tp_edge)-1)
        for l in range(len(tp_edge)-1):
            self.gc_pmid[l] = 0.5*(tp_edge[l]+tp_edge[l+1])

        # Interpolate GEOS-Chem onto CrIS vertical grid:
        # Need to flip arrays to be in ascending order:
        interp_gc_nh3 = np.flip(np.interp(np.flip(cris_pressure), np.flip(self.gc_pmid), np.flip(tgc_nh3)))

        # Get factor to convert data from ppbv to molec/cm2:
        dp = cris_pressure[:-1] - cris_pressure[1:]
        dp[0] = 0.
        ppb2col=(np.insert(dp,0,0)+np.append(dp,0))/2e9*na/g/mmair/1e2

        #print(interp_gc_nh3)

        # Apply averaging kernel to GEOS-Chem:
        #gc_with_avgk = np.exp(((np.identity(15) - cris_avg_kern)*np.log(cris_priori)) + cris_avg_kern*np.log(interp_gc_nh3))
        gc_with_avgk = (np.identity(15) - cris_avg_kern)*cris_priori + cris_avg_kern*interp_gc_nh3

        #print(len(cris_avg_kern)) 
        
        #gc_with_avgk = np.diagonal(gc_with_avgk)
        
        #print(gc_with_avgk)
        
        #print(len(ppb2col),len(interp_gc_nh3),len(gc_with_avgk))
        #sys.exit()

        # Get column value with and without averaging kernel:
        #bl_ind = np.where(cris_pressure<8.5E+2)[0]
        self.gc_col_og = np.nansum(interp_gc_nh3 * ppb2col)
        self.gc_col_ak = np.nansum(gc_with_avgk * ppb2col)

        # Remove any data obtained with very large difference between
        # GEOS-Chem with and without Avgk
        #if self.gc_col_ak/self.gc_col_og > 25:   #self.gc_col_og<=1e14: 
        #    #if self.gc_col_og > 1e15: 
        #    print(self.gc_col_og*1e-15,self.gc_col_ak/self.gc_col_og)
        #    self.gc_col_ak = np.nan

    def pres2alt(self,pressure):
        """Convert pressure provided in Pa to altitude in m.

        :param pressure: Pressure in Pa.

        :return: Height in m.
        """

        # Define additional constants:
        p_zero=101325    #Pa

        # Calculate height:
        height=(1 - ((pressure/p_zero)**0.190284)) * 145366.45 * 0.3048

        # Output:
        return (height)

    def get_average(self):
        '''Get mean over time period of interest'''
        self.grid_cris_nh3 = np.where( ((self.grid_count>0)&(uk_mask==1)), np.divide(self.grid_cris_nh3, self.grid_count), np.nan )
        self.grid_gc_nh3_og = np.where( ((self.grid_count>0)&(uk_mask==1)), np.divide(self.grid_gc_nh3_og, self.grid_count), np.nan )
        self.grid_gc_nh3_ak = np.where( ((self.grid_count>0)&(uk_mask==1)), np.divide(self.grid_gc_nh3_ak, self.grid_count), np.nan )
        self.grid_gc_cld = np.where( ((self.grid_count>0)&(uk_mask==1)), 100.*np.divide(self.grid_gc_cld, self.grid_count), np.nan )
        #self.grid_gc_cld = np.where( uk_mask==1, self.grid_gc_cld, np.nan )
        self.grid_count = np.where( uk_mask==1, self.grid_count, np.nan )

    def compare_gc_cris(self):
        '''Statistical comparison of the data (prints to screen)'''

        # Flatten data:
        self.cris1d  = self.grid_cris_nh3.flatten()
        self.gc1d    = self.grid_gc_nh3_ak.flatten()
        self.gc1d_og = self.grid_gc_nh3_og.flatten()
        self.cnt1d   = self.grid_count.flatten()

        # Get model standard deviation with and without the avgk:
        #print('Model std dev og (1e15): ',np.nanstd(self.gc1d_og)*1e-15)
        #print('Model std dev w/ avgk (1e15): ',np.nanstd(self.gc1d)*1e-15)

        # Find coincident data:
        stats_ind = np.where( (~np.isnan(self.cris1d))&(~np.isnan(self.gc1d))&(self.cnt1d>10) )[0]
                              
        # Correlation:
        self.r = stats.pearsonr(self.cris1d[stats_ind],self.gc1d[stats_ind])
        print('Correlation: ', self.r[0])

        # Regression:
        self.regres = rma(self.cris1d[stats_ind],self.gc1d[stats_ind],len(stats_ind),1000)
        print('RMA Slope: ', self.regres[0])
        print('RMA Intercept (10^15): ', self.regres[1]*1e-15)

        cris1d_reshape = self.cris1d[stats_ind].reshape(len(self.cris1d[stats_ind]),1)
        #print(tx_reshape.shape,tx[~nas].shape)
        #reg = TheilSenRegressor(random_state=0).fit(cris1d_reshape,(self.gc1d[stats_ind]))
        #print('TS Slope = ', reg.coef_)
        #print('TS Intercept (1e15) = ', 1e-15*reg.intercept_)

        # NMB:
        self.nmb = 100.0*(np.mean(self.gc1d[stats_ind])-np.mean(self.cris1d[stats_ind]))/np.mean(self.cris1d[stats_ind])
        print('Model mean NMB = ', self.nmb)
        print('Model mean [10^15]', 1e-15*np.mean(self.gc1d[stats_ind]))
        print('CrIS mean [10^15]', 1e-15*np.mean(self.cris1d[stats_ind]))
        print('Model (no avgk) mean [10^15]', 1e-15*np.mean(self.gc1d_og[stats_ind]))
        self.nmb = 100.0*(np.median(self.gc1d[stats_ind])-np.median(self.cris1d[stats_ind]))/np.median(self.cris1d[stats_ind])
        print('Model median normalized bias = ', self.nmb)

        # Data density
        print('Min no. of points = ', np.min(self.cnt1d[stats_ind]))
        print('Max no. of points = ', np.max(self.cnt1d[stats_ind]))

    def save_to_netcdf(self, month):
        '''Saves the seasonal_means to out_file as a netcdf4'''

        # Define filename:
        filedata = './Data/cris-gc-nh3-'+month+'-2013-2018-v1_6-v1.nc'

        ncout = Dataset(filedata, mode='w',format="NETCDF4")

        ncout.createDimension('lat', len(self.out_lat))
        ncout.createDimension('lon', len(self.out_lon))

        # create longitude axis:
        lon = ncout.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longiitude'
        lon[:] = self.out_lon

        # Create latitude axis:
        lat = ncout.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lat[:] = self.out_lat

        # Save data: 
        crisnh3 = ncout.createVariable('cris_nh3', np.float32, ('lon', 'lat'))
        crisnh3.units = 'molecules/cm2'
        crisnh3.long_name = 'CrIS multiyear (2013-2018) mean for '+mm+' at 0.25x0.3125'
        crisnh3[:] = self.grid_cris_nh3

        gcnh3og = ncout.createVariable('gcnh3_og', np.float32, ('lon', 'lat'))
        gcnh3og.units = 'molecules/cm2'
        gcnh3og.long_name = 'GEOS-Chem 2016 mean coincident with CrIS multiyear (2013-2018) mean for '+mm+' at 0.25x0.3125 (no averaging kernel)'
        gcnh3og[:] = self.grid_gc_nh3_og

        gcnh3ak = ncout.createVariable('gcnh3_ak', np.float32, ('lon', 'lat'))
        gcnh3ak.units = 'molecules/cm2'
        gcnh3ak.long_name = 'GEOS-Chem 2016 mean coincident with CrIS multiyear (2013-2018) mean for '+mm+' at 0.25x0.3125 (with averaging kernel)'
        gcnh3ak[:] = self.grid_gc_nh3_ak

        num_obs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
        num_obs.units = 'unitless'
        num_obs.long_name = 'Number of CrIS observations'
        num_obs[:] = self.grid_count

        ncout.close()

    def plot_data(self, month):
        '''Plot maps of CrIS NH3'''
        X, Y = np.meshgrid(self.out_lon, self.out_lat, indexing='ij')

        # Define plot parameters:
        nplots = 5
        max_val = [5.0, 1.5, 1.5, 75, 100]
        nbound = [25, 30, 30, 25, 25]
        vcd_unit = '[10$^{16}$ molecules cm$^{-2}$]'
        unit = [vcd_unit, vcd_unit, vcd_unit, 'unitless', 'unitless']
        plot_title = ['CrIS','GC without AK','GC with AK',
                      '% Cloud fraction > 0.5', 
                      'No. of CrIS observations']
        data_name = ['CrIS NH3','GC NH3 (og)','GC NH3 (avgk)','Cloud', 
                     'No. of obs']

        fig, ax = plt.subplots(1, 5, figsize=(22,4), subplot_kw=dict(projection=ccrs.PlateCarree()))

        #ind = np.where(np.isnan(self.grid_count))[0]

        # Plot the subplots:
        for i in range(nplots):

            # Define what data to plot:
            if i==0: plot_vals = np.squeeze(self.grid_cris_nh3*1e-16)
            if i==1: plot_vals = np.squeeze(self.grid_gc_nh3_og*1e-16)
            if i==2: plot_vals = np.squeeze(self.grid_gc_nh3_ak*1e-16)
            if i==3: plot_vals = np.squeeze(self.grid_gc_cld)
            if i==4: plot_vals = np.squeeze(self.grid_count)

            print('Std dev of ' + data_name[i] + ':', np.nanstd(plot_vals))
            print('Rel Std dev of ' + data_name[i] + ':', 
                  np.nanstd(plot_vals)/np.nanmean(plot_vals))

            if i==0: 
                tickval = [0,1,2,3,4,5]
                #tickval=tickval*2
            if ((i>0) & (i<=1)):
                tickval = [0,0.5,1,1.5]            
            if i>2:
                tickval = [0,20,40, 60 ,80 ,100]  #15,30,45,60,75]  
            if i==2:
                tickval = [0,0.5,1,1.5]  #[0,0.5,1,1.5,2,2.5]
        
            political = cartopy.feature.NaturalEarthFeature(
                category='cultural',
                name='admin_0_boundary_lines_land',
                scale='10m',
                facecolor='none')

            states_provinces = cartopy.feature.NaturalEarthFeature(
                category='cultural',
                name='admin_0_boundary_lines_map_units',
                scale='10m',
                facecolor='none')

            counties = cartopy.feature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='10m',
                facecolor='none')

            ax[i].add_feature(political, edgecolor='black')
            ax[i].add_feature(states_provinces, edgecolor='black')
            #ax[i].add_feature(counties, edgecolor='black')
            ax[i].coastlines(resolution='10m')

            ax[i].set_extent([-10, 3, 49.5, 60], crs=ccrs.PlateCarree())
            ax[i].set_title(plot_title[i])

            data_crs = ccrs.PlateCarree()
            bounds = np.linspace(0, max_val[i], nbound[i])
            #cmap = plt.cm.jet 

            c = ax[i].pcolormesh(X,Y, plot_vals, transform=data_crs,cmap='jet', 
                                 vmin=0, vmax=max_val[i])

            cb = fig.colorbar(c, ax=ax[i],label=unit[i], orientation='horizontal',
                              shrink=0.5,pad=0.01,boundaries=bounds,ticks=tickval )

            cb.ax.tick_params(labelsize=10, direction='in', length=6) 

        fileplot='./Images/cris-gc-nh3-'+month+'-2013-2018-v1.ps'
        #plt.savefig(fileplot, format='ps')
        plt.show()

if __name__ == "__main__":

    # Define output grid:
    LON_RES = 0.3125
    LAT_RES = 0.25
    OUT_LON = np.arange(-15.00, 40.00+LON_RES, LON_RES)
    OUT_LAT = np.arange( 32.75, 61.25+LAT_RES, LAT_RES)

    # Define output array dimensions:
    XDIM = len(OUT_LON)
    YDIM = len(OUT_LAT)

    # Define month of interest:
    MONTH = ["07"]#, "04", "05", "06", "07", "08", "09"]
    YEARS = ["2013", "2014", "2015", "2016", "2017", "2018"]
    NDAYS = 31

    # Define CrIS directory:
    #cris_dir = '/data/uptrop/nobackup/cris_nh3_ecc/UK/'  #'2013/05/'
    cris_dir = '/data/uptrop/nobackup/cris_nh3_ecc_v1_6/UK/'
    # Define GEOS-Chem directory:
    #gc_dir = '/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_scale_nh3_emis/nc_sat_files/'
    gc_dir = '/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/nc_sat_files/'

    # Loop over months:
    for mm in MONTH:

        # Initialize data arrays of interest for obtaining monthly
        # multiyear means:
        GRID_CRIS_NH3  = np.zeros((XDIM, YDIM))
        GRID_GC_NH3_OG = np.zeros((XDIM, YDIM))
        GRID_GC_NH3_AK = np.zeros((XDIM, YDIM))
        GRID_GC_CLD    = np.zeros((XDIM, YDIM))
        GRID_COUNTER   = np.zeros((XDIM, YDIM))

        # Track progress:
        print(' ==== >> Processing Month: ', mm)

        # Loop over years:
        for yy in YEARS:

            # Track progress:
            print(' ==== >> Processing Year: ', yy)

            # Loop over days:
            for dd in range(NDAYS):

                #if (yy=='2018')&(mm=='07'):
                #        if (dd+1)==26 or (dd+1)==27: continue

                # Define day:
                DAY = str(dd+1)
                if dd<10: DAY='0'+DAY

                # Consider filtering out 26-27 July 2018 (sustained high
                # temperatures that lead to 4-times higher IASI NH3).

                # Get CrIS files:
                cris_files = glob.glob( path.join( cris_dir, yy, mm, ('*'+mm+DAY+'.nc') ) )

                # Get GEOS-Chem file:
                gc_file = glob.glob( path.join( gc_dir, ('ts_12_15.EU.2016'+mm+DAY+'*') ) )

                # Loop over CrIS files:
                for f in cris_files:

                    print('Reading CrIS file: ', f)
                    print('Reading GC file: ', gc_file[0])

                    # Read relevant CrIS NH3 data:
                    cris_data = CompareGCCris(f, gc_file[0], OUT_LON, OUT_LAT, GRID_CRIS_NH3, GRID_GC_NH3_OG, GRID_GC_NH3_AK, GRID_GC_CLD, GRID_COUNTER)

        # Get UK land mask:
        uk_mask = cris_data.read_uk_land_mask()

        # Maybe replace with get_average:
        cris_data.get_average()    

        # Compare the data:
        cris_data.compare_gc_cris()

        # Save the data to NetCDF:
        #cris_data.save_to_netcdf(mm)

        # Plot the data:
        cris_data.plot_data(mm)

    # Read relevant GEOS-Chem NH3 data:
    #sys.exit()
