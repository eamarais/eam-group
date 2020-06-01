#!/usr/bin/python

'''
Description: uses jackknife resampling to estimate the reduced major axis
             regression slopes and intercepts and the errors of these.
             Relevant reference is: Hirsch and Gilroy, Water Res. Bull., 
             20(5), Oct 1984.

Input: x: 1D array of x values
       y: 1D array of y values
       n: number of elements of x array (equivalent to the number of elements 
          of the y array)
       ntrials: the number of trials of randomly selected data

Return: grad_arr, cept_arr: 1D arrays of slope and intercept estimates 
        the same size as the number of trials (ntrials). Mean of these is the
        slope and intercept value and standard deviation is the error on
        the slope and intercept.

Original code on which this is based: GAMAP package bootstrap.pro
'''

# Import relevant packages:
import sys
from scipy import stats
import numpy as np

def rma(x,y,n,ntrials):

    # Get correlation:
    r=stats.pearsonr(x,y)

    # Initialize:
    fac = 0.
    cnt = 0

    # Find fac based on the sign of the correlation coefficient:
    if ( r[0] >0.0 ): fac=1.0
    if ( r[0]==0.0 ): fac=0.0
    if ( r[0] <0.0 ): fac=-1.0

    if ( np.isnan(r[0]) ): 
        'R is NaN -- EXITING PROGRAMME'
        sys.exit()

    # Define output arrays:
    grad_arr=np.zeros(ntrials)
    cept_arr=np.zeros(ntrials)

    # Loop over trials:
    for w in range(ntrials):

        # Random selection of data indices:
        index=np.random.randint(n,size=n)
        # Avoid having the same index throughout the randomly selected array
        # leading to Sx=0 and slope=Inf by redoing the random number selection
        # if this is the case:
        if ( (len(np.where(index==index[0])[0]))==(len(index)) ):
            index=np.random.randint(n,size=n)
        #if ( w==0 ): print('index: ',index)
        #if ( w==5 ): print('index: ',index)

        # Define randomly selected x and y data:
        x_rdm=x[index]
        y_rdm=y[index]

        # Get shuffled x and y means:
        xbar=np.mean(x_rdm)
        ybar=np.mean(y_rdm)

        # Get the population standard deviation:
        Sx=np.sqrt(np.divide(np.sum(np.square(np.subtract(x_rdm,xbar))),\
                             float(n)))
        Sy=np.sqrt(np.divide(np.sum(np.square(np.subtract(y_rdm,ybar))),\
                             float(n)))

        if (Sy==0): continue

        # Get slope and intercept:
        grad=np.multiply(fac, np.divide(Sy,Sx))
        cept=np.subtract(ybar, np.multiply(grad,xbar))

        grad_arr[cnt]=grad
        cept_arr[cnt]=cept

        # Increment:
        cnt += 1

    # Get output values:
    slope=np.mean(np.trim_zeros(grad_arr))
    intercept=np.mean(np.trim_zeros(cept_arr))
    slope_err=np.std(np.trim_zeros(grad_arr))
    intercept_err=np.std(np.trim_zeros(cept_arr))

    # Return quantities are 1D arrays of gradient and intercept estimates:
    return(slope,intercept,slope_err,intercept_err)
