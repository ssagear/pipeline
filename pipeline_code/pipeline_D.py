
from __future__ import division

import bls, ktransit, math, pylab, os, batman

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

import obspy
from obspy.signal.detrend import polynomial

import PyAstronomy

from PyAstronomy import pyasl




time = []
flux = []

file = 'hlsp_k2sff_k2_lightcurve_229227268-c06_kepler_v1_llc-default-aper.csv'

with open(r'/Users/sheilasagear/Dropbox/K2/K2_targets/Cycle2_4_CSV/hlsp_k2sff_k2_lightcurve_229227268-c06_kepler_v1_llc-default-aper.csv') as f:
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
    flux = sigma_clip(flux, sigma=4, iters=1)

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



##########################
#BLS routine
##########################

    
    u = [0.0]*len(time)
    v = [0.0]*len(flux)

    u = np.array(u)
    v = np.array(v)

    nf = 200.0
    fmin = 0.02 #up to 50 days
    df = 0.01 #down to .5 days
    nb = 200
    qmi = 0.001
    qma = 0.3

    results = bls.eebls(time, flux, u, v, nf, fmin, df, nb, qmi, qma)
            
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

    #plt.clf()

    #plt.plot(freq, SR_array)
    #plt.show()

    
#fitting????


#UNCOMMENT TO SAVE SIGNAL RESIDUE PLOT

    SRper = [1/x for x in freq]

    """
    fig, ax = plt.subplots(1, 1, figsize=[11,5])
    ax.scatter(SRper, SR_array, color='black')
    ax.set_ylabel('Power')
    ax.set_xlabel('Period (days)')
    ax.set_title('Box Least-Squares')
    #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5")
        
    #show()
    """
    
    
##########################
#Calculate SDE
########################## 

    SDE = (max_SR-avg_SR)/sd_SR
        
    print('Signal Detection Efficiency: ' + str(SDE))

    if SDE >= 6:
        print('SDE above 6: Transit Detectable')
    else:
        print('SDE below 6: Transit Undetectable')






#centering


##########################
#Folding and Binning
##########################

    #pylab.cla()
    
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
    
    #plt.scatter(phase, y/ibi, s=3)





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



    trend = untrendy.median(time, flux)

    fluxDetrend = []

    length = len(time)

    for i in range(length):
        fluxDetrend.append(flux[i]-trend[i])

    fluxDetrend = np.asarray(fluxDetrend)
    

    

    phases = PyAstronomy.pyasl.foldAt(time, results[1], getEpoch=False)

    sortPhase = np.argsort(phases)
    phases = phases[sortPhase]
    fluxFolded = fluxDetrend[sortPhase]

    #replace p with results[1]

        
    #pylab.cla()

#UNCOMMENT TO SAVE PHASE FOLDED PLOT

    #pylab.clf()
    #xlabel('Phase')
    #ylabel('Corrected Flux (normalized to 0)')
    #title('EPIC (name) Phase Folded Light Curve')
    #plt.plot(phases, fluxFolded)
    #pylab.show()

###########################
#BLS Overlay
###########################

        
    high = results[3]*results[4]
    low = high - results[3]
    
    fit = np.zeros(nb) + high # H
    fit[results[5]:results[6]+1] = low # L
        
    #plt.plot(phase, fit)
    #plt.xlabel(r"Phase")
    #plt.ylabel(r"Mean value of flux")
    #plt.title("SDE " + str(SDE) + "; BLS period " + str(results[1]))
    #plt.ylim(-.2, .2)

    """
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(11,5))
    
    ax.scatter(phases, fluxFolded, s=2, c='k')
    ax.plot(phase, fit, color='red')
    ax.set_xlabel("Phase")
    ax.set_ylabel("Mean value of flux")
    ax.set_title("SDE " + str(SDE) + "; BLS period " + str(results[1]))

    show()
    """
    

    #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/Pipeline/EPIC229227254-2sigclip/Period" + str(p) + "Radius" + str(r) + "/folded_pltPer" + str(p) + "Rad" + str(r) + 'BLSoverlay.png')

    #pylab.show()

    #ylabel('Corrected Flux (normalized to 0)')
    #title('EPIC ' + targetname + ' Light Curve')

    #print('Low = ' + str(low))
        

##########################
#Detrending Flux and Levenberg-Marquardt Fitting:
#if SDE>6 
##########################

#start values are correct values

    if True: #SDE >= 6:
            
        #polynomial detrending
        #fluxDetrend = polynomial(flux, order=5)

        #untrendy
        #fluxDetrend, ferr = untrendy.median(time, flux)
        #fluxDetrend = [x-1.0 for x in fluxDetrend
                
        #for i in range(len(trend)):
        #    val = flux[i]-flux[i-1]
        #    fluxDetrend.append(val)


        fitT = FitTransit()
        
        fitT.add_guess_star(rho=1.5)
        
        fitT.add_guess_planet(
                period=results[1], impact=0.0, 
                T0=3.0, rprs=.2)#need a guess rprs
        fitT.add_data(time=time, flux=fluxDetrend)
                    

        vary_star = ['rho']      # not sure how to avoid free stellar parameters? ideally would not vary star at all
        vary_planet = (['period', 'rprs', 'impact', 'T0'])
                
        fitT.free_parameters(vary_star, vary_planet)
        fitT.do_fit()

        #print(fitT.fitresultstellar.items())
        #print(fitT.fitresultplanets.items())

        bestFstellar = fitT.fitresultstellar.items()
        bestFrho = bestFstellar[0][1]#Best Fit Rho
                    
        bestFplanet = fitT.fitresultplanets.items()
        bestFperiod = bestFplanet[0][1]['period']#Best Fit Period
        bestFrprs = bestFplanet[0][1]['rprs']#Best Fit rprs
                    

        fitT.print_results()
                    
#save figure
        #fig = ktransit.plot_results(time,fluxDetrend,fitT.transitmodel)
        #fig.savefig('/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/ktransitfits_untrendy/' + str(targetname) + 'fitPer' + str(p) + 'Rprs' + str(r) + '.png')


        #if not os.path.exists("/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/EPIC" + str(targetname)):
        #    os.makedirs("/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/EPIC" + str(targetname))
        """        
        pylab.cla()
        fig, ax = plt.subplots(1, 1, figsize=[11,5])
        ax.scatter(time, fluxDetrend, color='k', s=2)
        ax.plot(time, fitT.transitmodel, color='mediumaquamarine')
        ax.set_ylim(-0.07,0.07)
        ax.set_xlabel('Time (BJD - 2450000)')
        ax.set_title('EPIC ' + str(targetname))
        bbox_props = dict(boxstyle="square,pad=0.3", facecolor='none', edgecolor='black')
        ax.annotate('Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), xy = (time[5], 0.04), bbox = bbox_props)
        fig.savefig('/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/EPIC' + str(targetname) + '/EPIC' + str(targetname) + 'KRANSITfit.png')

        #show()
        """

##########################
#END LOOP
##########################


##########################
#INITIAL LIGHT CURVE
##########################

"""
pylab.cla()

fig, ax = plt.subplots(1, 1, figsize=[11,8])

ax.scatter(time, flux, color='k', s=2)
ax.plot(time, q, color='mediumaquamarine')
ax.set_xlabel('Time')
ax.set_ylabel('Corrected Flux (Normalized)')
ax.set_title('EPIC ' + str(targetname))

#show()
"""


#color for outliers

q = untrendy.median(time, flux)


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(11, 5))
#fig.subplots_adjust(hspace=.01)

bbox_props = dict(boxstyle="square,pad=0.4", facecolor='none', edgecolor='black')

ax1.scatter(time, flux, color='k', s=2)
ax1.plot(time, q, color='mediumaquamarine')
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Normalized Flux')
ax1.set_title('EPIC ' + str(targetname))


ax2.plot(freq, SR_array, color='black')
ax2.set_ylabel('Power')
ax2.set_xlabel('Frequency')
ax2.set_title('Box Least-Squares')
ax2.set_ylim(0, 1.3)
ax2.text(0.2, 1.05, 'SDE: ' + str(SDE) + '\n' + 'BLS Period: ' + str(results[1]), bbox=bbox_props, fontsize=7)
#fitting?????


ax3.scatter(phases, fluxFolded, color='k', s=2)
ax3.set_xlabel('Phase')
ax3.set_ylabel('Normalized Flux')
ax3.set_title('Phase Folded: Period ' + str(results[1]))

"""
ax4.scatter(time, fluxDetrend, color='k', s=2)
ax4.plot(time, fitT.transitmodel, color='mediumaquamarine')
ax4.set_ylim(-0.25,0.25)

ax4.set_xlabel('Time (days)')
ax4.set_ylabel('Detrended Flux')
ax4.set_title('Levenberg-Marquardt')
#ax4.annotate('Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), xy = (time[1], 0.04), bbox = bbox_props)
ax4.text(time[1], 0.12, 'Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), bbox=bbox_props, fontsize=7)
"""


plt.tight_layout()
show()

#f.savefig('/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID//Cycle2_4_Pipeline_PNG/EPIC' + str(targetname) + 'pipeline_3sigclip.png')

