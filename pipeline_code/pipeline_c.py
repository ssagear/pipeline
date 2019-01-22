from __future__ import division

import bls, ktransit, math, pylab, os, batman, untrendy, csv

import numpy as np

from scipy.stats import gamma

from astropy.io import fits

import matplotlib.pyplot as plt
fig = plt.figure()
axes = plt.gca()

import pylab
from pylab import xlabel, ylabel, title, plot, show, ylim, xlim

from astropy.stats import sigma_clip

from ktransit import FitTransit

from scipy import signal


time = []
flux = []


file = 'hlsp_k2sff_k2_lightcurve_212826600-c06_kepler_v1_llc-default-aper.csv'

with open(r'/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/Cycle2_4_CSV/hlsp_k2sff_k2_lightcurve_212826600-c06_kepler_v1_llc-default-aper.csv') as f:
    reader = csv.reader(f)

    for row in reader:
        time.append(row[0])
        flux.append(row[1])

    for i in range(len(time)):
        time[i] = float(time[i])
        flux[i] = float(flux[i])

    #print(time)
    #print('\n')
    #print(flux)

#normalize flux values to 1.0
    flux_count = len(flux)
    sumflux = sum(flux)
    flux_avg = sumflux/flux_count
    
    flux = [i/flux_avg for i in flux]

#targets titled by EPIC names

    targetname = file[25:]
    targetname = targetname[:-35]
    #title('EPIC ' + targetname + ' Light Curve')
    print("TARGET: EPIC " + str(targetname))
        
    campaign_no = file[36:]
    campaign_no = campaign_no[:-31]
    print("CAMPAIGN " + str(campaign_no))
        
#normalize to 0.0
    flux = [i-1 for i in flux]

    flux = np.asarray(flux)

    #print(type(flux))


#SIGMA CLIPPING
    flux = sigma_clip(flux, sigma=3, iters=1)

#uncomment if extra time stamp
    #time.pop()

##########################
#IMPORT CDPP FROM FITS
##########################
    #f = fits.open('/Users/sheilasagear/OneDrive/K2_Research/Cycle2_FITS/CYCLE2FITS/hlsp_k2sff_k2_lightcurve_' + str(targetname) + '-c0' + str(campaign_no) + '_kepler_v1_llc.fits')

    #bestaper = f[1].data
    #besthead = f[1].header
        
    #QCDPP = besthead['QCDPP6']
    #print(QCDPP)
    #QCDPParr.append(QCDPP)

    
    #trend, ferr = untrendy.untrend(time, flux)

    #fluxDetrend = list()
                
    #for i in range(len(trend)):
    #    val = flux[i]-flux[i-1]
    #    fluxDetrend.append(val)

    u = [0.0]*len(time)
    v = [0.0]*len(time)

    u = np.array(u)
    v = np.array(v)

    #time, flux, u, v, number of freq bins (nf), min freq to test (fmin), freq spacing (df), number of bins (nb), min transit dur (qmi), max transit dur (qma)

    nf = 1000.0
    fmin = .02
    df = 0.001
    nbins = 700
    qmi = 0.001
    qma = 0.3
    
    results = bls.eebls(time, flux, u, v, nf, fmin, df, nbins, qmi, qma)
            
#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

    print('BLS period: ' + str(results[1]))

    

    SR_array = results[0]
    max_SR = max(SR_array)
    avg_SR = np.mean(SR_array)
    sd_SR = np.std(SR_array)

#normalize SR_array between 0 and 1
    SR_array = [i/max_SR for i in SR_array]

    #print(max_SR, avg_SR, sd_SR)

#np.arange(freq start, freq stop (must be calculated), df (freq step))
    #freq = np.arange(.2, 1.2, .001, dtype=None)

    freq = fmin + np.arange(nf) * df

    print(freq, SR_array)

    plt.clf()

    plt.plot(freq, SR_array)
    plt.title('EPIC ' + str(targetname))
    plt.show()
