
from __future__ import division

import bls, ktransit, math, pylab, os

from scipy.stats import gamma

import untrendy

from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

import csv

from matplotlib import pyplot
import pylab
from pylab import xlabel, ylabel, title, plot, show, ylim, xlim

from astropy.stats import sigma_clip

from ktransit import FitTransit

from scipy import signal

import PyAstronomy

from PyAstronomy import pyasl

tested_period = []
tested_rprs = []
time = []
flux = []


z = [0, 0, 0, 0]

z[0] = -3.22185
z[1] = 1.02655
z[2] = 1.6613
z[3] = .225603

p = np.poly1d(z)

file = 'hlsp_k2sff_k2_lightcurve_217976219-c07_kepler_v1_llc-default-aper.csv'

with open(r'/Users/sheilasagear/Dropbox/K2/K2_targets/Cycle2_4_CSV/hlsp_k2sff_k2_lightcurve_217976219-c07_kepler_v1_llc-default-aper.csv') as f:
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

    #print('flux = ' + str(type(flux)))

    flux = np.asarray(flux)
    time = np.asarray(time)

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
                
    #CDPP = besthead['QCDPP6']

    #print('QCDPP: ' + str(CDPP))



    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.scatter(time, flux, c='k', s=2)
    ax.set_xlabel('BJD')
    ax.set_ylabel('Flux')
    ax.set_title('EPIC ' + str(targetname))
    show()
    

    T0 = 3.0
    impact = 0.0
    rho = 1.5
    
    
    #period = np.random.uniform(low=1, high=26)
    period = 5.0
    tested_period.append(period)
        
    #rprs = np.random.uniform(low=.05, high=1.0)
    rprs = 0.3
    tested_rprs.append(rprs)

        
    M = ktransit.LCModel()
    M.add_star(
            rho=rho, # mean stellar density in cgs units
            ld1=0.2, # ld1--4 are limb darkening coefficients 
            ld2=0., # if only ld1 and ld2 are non-zero then a quadratic limb darkening law is used
            ld3=0.0, # if all four parameters are non-zero we use non-linear flavour limb darkening
            ld4=0.0, 
            dil=0.0, # a dilution factor: 0.0 -> transit not diluted, 0.5 -> transit 50% diluted
            zpt=0.0  # a photometric zeropoint, incase the normalisation was wonky
    )
    M.add_planet(
            T0=T0,     # a transit mid-time  
            period=period, # an orbital period in days
            impact=impact, # an impact parameter
            rprs=rprs,   # planet stellar radius ratio  
            ecosw=0.0,  # eccentricity vector
            esinw=0.0,
            occ=0.0)    # a secondary eclipse depth in ppm

    M.add_data(time=np.array(time[:]))      # integration time of each timestamp

    tmod = M.transitmodel # the out of transit data will be 0.0 unless you specify zpt

    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.plot(time, tmod, c='k')
    ax.set_xlabel('BJD')
    ax.set_ylabel('Injected flux')
    ax.set_title('Period: ' + str(period) + ' | RPRS: ' + str(rprs))
    show()




    merged_flux = tmod + flux

    fig.clf()
    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.scatter(time, merged_flux, c='k', s=2)
    ax.set_xlabel('BJD')
    ax.set_ylabel('Injected flux')
    ax.set_title('Merged Light Curve')
    show()
    

##########################
#BLS routine
##########################

    
    u = [0.0]*len(time)
    v = [0.0]*len(flux)

    u = np.array(u)
    v = np.array(v)

    nf = 1000.0
    fmin = .035
    df = 0.001
    nb = 300
    qmi = 0.001
    qma = 0.3

    results = bls.eebls(time, merged_flux, u, v, nf, fmin, df, nb, qmi, qma)
            
#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

    print('BLS period: ' + str(results[1]))


##########################
#Power Spectrum
##########################

    SR_array = results[0]
    max_SR = max(SR_array)
    avg_SR = np.mean(SR_array)
    sd_SR = np.std(SR_array)

#normalize SR_array between 0 and 1
    SR_array = [i/max_SR for i in SR_array]



##########################
#Calculate SDE
########################## 

    SDE = (max_SR-avg_SR)/sd_SR
        
    print('Signal Detection Efficiency: ' + str(SDE))

    if SDE >= 6:
        print('SDE above 6: Transit Detectable')
    else:
        print('SDE below 6: Transit Undetectable')

    #print(max_SR, avg_SR, sd_SR)

#np.arange(freq start, freq stop (must be calculated), df (freq step))
    #freq = np.arange(.2, 1.2, .001, dtype=None)

    freq = fmin + np.arange(nf) * df

    fig.clf()
    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.plot(freq, SR_array, c='k')
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Power')
    ax.set_title('Box Least-Squares: SDE ' + str(SDE))
    show()



#Power Spectrum with Period
    #SRper = [1/x for x in freq]
    

##########################
#Data Centering (from http://nbviewer.jupyter.org/github/ridlo/exoplanet_notebook/blob/master/bls_test01.ipynb)
##########################


    time1 = time[0]
    u = time-time1
    s = np.mean(merged_flux)
    v = merged_flux-s

    fig.clf()
    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.scatter(u, v, c='k', s=2)
    ax.set_xlabel('t - t0')
    ax.set_ylabel('flux')
    ax.set_title('Data centering')
    show()

##########################
#Folding and Binning (from http://nbviewer.jupyter.org/github/ridlo/exoplanet_notebook/blob/master/bls_test01.ipynb)
##########################
            
    f0 = 1.0/results[1]
    n = len(time)
    ibi = np.zeros(nb)
    y = np.zeros(nb)
    phase = np.linspace(0.0, 1.0, nb)
    
    for i in range(n):
        ph = u[i]*f0
        ph = ph-int(ph)
        j = int(nb*ph)
        ibi[j] = ibi[j] + 1.0
        y[j] = y[j] + v[i]

    fig.clf()
    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.scatter(phase, y/ibi, c='k', s=2)
    ax.set_xlabel('Phase')
    ax.set_ylabel('Binned FLux')
    ax.set_title('Period: ' + str(period) + ' bin: ' + str(nb))
    show()




##########################
#Detrending (untrendy)
##########################


    trend = untrendy.median(time, merged_flux)

    mergedfluxDetrend = []
    length = len(time)
        
    for x in range(length):
        mergedfluxDetrend.append(merged_flux[x]-trend[x])

    mergedfluxDetrend = np.asarray(mergedfluxDetrend)

#Detrended, not binned

    fig.clf()
    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.scatter(time, mergedfluxDetrend, c='k', s=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Detrended Flux')
    ax.set_title('Detrended light curve')
    show()


##########################
#Just Folding
##########################

    """

    start = time[0]
    end = time[len(time)-1]

    time_folded = [t/(results[1]) for t in time]
    time_folded = [i % 1 for i in time_folded]

    print(time_folded)

    """

    phases = PyAstronomy.pyasl.foldAt(time, results[1], getEpoch=False)

    sortPhase = np.argsort(phases)
    phases = phases[sortPhase]
    fluxFolded = mergedfluxDetrend[sortPhase]

    
#JUST PHASE FOLDED PLOT

    fig.clf()
    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.scatter(phases, fluxFolded, c='k', s=2)
    ax.set_xlabel('Phase')
    ax.set_ylabel('Binned FLux')
    ax.set_title('Phase Folded: ' + str(results[1]) + ' days')
    show()


###########################
#BLS Overlay
###########################

        
    high = results[3]*results[4]
    low = high - results[3]
    
    fit = np.zeros(nb) + high # H
    fit[results[5]:results[6]+1] = low # L

    depth = high-low


    fig.clf()
    fig, ax = plt.subplots(1, 1, figsize=[15,10])
    ax.scatter(phase, y/ibi, c='k', s=2)
    ax.plot(phase, fit, c='darkred')
    ax.set_xlabel('Phase')
    ax.set_ylabel('Binned FLux')
    ax.set_title('BLS period: ' + str(results[1]) + ' bin: ' + str(nb))
    show()

        

##########################
#Detrending Flux and Levenberg-Marquardt Fitting:
#if SDE>6 
##########################

#start values are correct values

    if True: #if SDE >= 6

        print('rprs guess ' + str(p(depth)))

        rprsguess = p(depth)

        fitT = FitTransit()
        fitT.add_guess_star(rho=rho)    
        fitT.add_guess_planet(
        period=results[1], impact=impact, 
        T0=T0, rprs=rprsguess)
        fitT.add_data(time=time, flux=mergedfluxDetrend)

        vary_star = ['rho']      # free stellar parameters
        vary_planet = (['period',       # free planetary parameters 
        'rprs'])                # free planet parameters are the same for every planet you model

        fitT.free_parameters(vary_star, vary_planet)
        fitT.do_fit()                   # run the fitting

        fitT.print_results()    

        bestFstellar = fitT.fitresultstellar.items()
        bestFrho = bestFstellar[0][1]#Best Fit Rho
                    
        bestFplanet = fitT.fitresultplanets.items()
        bestFperiod = bestFplanet[0][1]['period']#Best Fit Period
        bestFrprs = bestFplanet[0][1]['rprs']#Best Fit rprs


        #fig = ktransit.plot_results(time,mergedfluxDetrend,fitT.transitmodel)
        #fig.show()
        
        fig, ax = plt.subplots(1, 1, figsize=[11,5])
        ax.scatter(time, mergedfluxDetrend, color='k', s=2)
        ax.plot(time, fitT.transitmodel, color='darkred')
        #ax.set_ylim(-0.07,0.07)
        ax.set_xlabel('Time (BJD - 2450000)')
        ax.set_title('EPIC ' + str(targetname))
        bbox_props = dict(boxstyle="square,pad=0.3", facecolor='none', edgecolor='black')
        ax.annotate('Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), xy = (time[5], 0.07), bbox = bbox_props)

        show()
        

##########################
#END LOOP
##########################


##########################
#INITIAL LIGHT CURVE with trendline
##########################

q = untrendy.median(time, flux)
        
fig.clf()
fig, ax = plt.subplots(1, 1, figsize=[11,8])

ax.scatter(time, flux, color='k', s=2)
ax.plot(time, q, color='darkred')
ax.set_xlabel('Time')
ax.set_ylabel('Corrected Flux (Normalized)')
ax.set_title('EPIC ' + str(targetname))

show()



f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(10, 6))

bbox_props = dict(boxstyle="square,pad=0.4", facecolor='none', edgecolor='black')

ax1.scatter(time, flux, color='k', s=2)
ax1.plot(time, q, color='darkred')
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Normalized Flux')
ax1.set_title('EPIC ' + str(targetname))


ax2.plot(freq, SR_array, color='black')
ax2.set_ylabel('Power')
ax2.set_xlabel('Frequency')
ax2.set_title('Box Least-Squares')
ax2.set_ylim(0, 1.3)
ax2.text(0.8, .7, 'SDE: ' + str(SDE) + '\n' + 'BLS Period: ' + str(results[1]), bbox=bbox_props, fontsize=7)


ax3.scatter(phases, fluxFolded, color='k', s=2)
ax3.set_xlabel('Phase')
ax3.set_ylabel('Normalized Flux')
ax3.set_title('Phase Folded: Period ' + str(results[1]))


ax4.scatter(time, mergedfluxDetrend, color='k', s=2)
ax4.plot(time, fitT.transitmodel, color='darkred')
ax4.set_ylim(-0.25,0.25)
bbox_props = dict(boxstyle="square,pad=0.3", facecolor='none', edgecolor='black')
ax4.set_xlabel('Time (days)')
ax4.set_ylabel('Detrended Flux')
ax4.set_title('Levenberg-Marquardt')
ax4.text(time[1], 0.12, 'Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), bbox=bbox_props, fontsize=7)


plt.tight_layout()



show()



