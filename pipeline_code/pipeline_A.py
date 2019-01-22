
from __future__ import division

import bls, ktransit, math, pylab, os, batman

from scipy.stats import gamma

import untrendy

from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
fig = plt.figure()
axes = plt.gca()



#from mpl_toolkits.mplot3d import Axes3D
#axe = Axes3D(plt.gcf())
#ax = fig.add_subplot(111, projection='3d')

from matplotlib import pyplot
import pylab
from pylab import xlabel, ylabel, title, plot, show, ylim, xlim

from astropy.stats import sigma_clip

#fige = pylab.figure()
#ax = Axes3D(fig)

from ktransit import FitTransit

from scipy import signal

import obspy
from obspy.signal.detrend import polynomial


##########################
#BEGIN LOOP
##########################

#LATER: loop through C6 and C7 light curves; open and inject range of transits into each 

file = 'hlsp_k2sff_k2_lightcurve_212820594-c06_kepler_v1_llc-default-aper.txt'

with open(r'/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/hlsp_k2sff_k2_lightcurve_212820594-c06_kepler_v1_llc-default-aper.txt') as f:
    data1 = f.read()

#converts txt file to string (and removes heading: 30 characters) 
#in order to isolate each (x,y) coordinate as an element in list 'data'
    datastr = str(data1)
    datastr = datastr[30:]
    data1 = datastr.split('\n')


#removes comma after each (x,y) coordinate;
#isolates x and y values as indicies of list 'data'
    index = -1
    while (index < len(data1)):
        tempstring = str(data1[index])
        data1[index] = tempstring.rstrip(',')
        data1[index] = tempstring.split(', ')
        index+=1
    data1.pop
        
    data2 = sum(data1, [])
        
    index = 0
    while (index < len(data2)):
        if index % 2 == 1:
            data2[index] = data2[index].rstrip(',')
        index+=1
    data2.pop()
           
#converts str data points to float
    data_final = [float(i) for i in data2]
    
#defines x and y values by index of 'data'
    time = data_final[0::2]
    flux = data_final[1::2]

#normalize flux values to 1.0
    flux_count = len(flux)
    sumflux = sum(flux)
    flux_avg = sumflux/flux_count
    
    flux = [i/flux_avg for i in flux]

#targets titled by EPIC names

    targetname = file[25:]
    targetname = targetname[:-35]
    title('EPIC ' + targetname + ' Light Curve')
    print("TARGET: EPIC " + str(targetname))
        
    campaign_no = file[37:]
    campaign_no = campaign_no[:-31]
    print("CAMPAIGN " + str(campaign_no))
        
#normalize to 0.0
    flux = [i-1 for i in flux]

    print('flux = ' + str(type(flux)))

    flux = np.asarray(flux)

    print(type(flux))


#SIGMA CLIPPING
    #flux = sigma_clip(flux, sigma=3, iters=1)

#uncomment if extra time stamp
    #time.pop()
    
##########################
#IMPORT CDPP FROM FITS
##########################
    f = fits.open('/Users/sheilasagear/OneDrive/K2_Research/Cycle2_FITS/CYCLE2FITS/hlsp_k2sff_k2_lightcurve_' + str(targetname) + '-c0' + str(campaign_no) + '_kepler_v1_llc.fits')

    bestaper = f[1].data
    besthead = f[1].header
                
    CDPP = besthead['QCDPP6']

    print(CDPP)



##########################
#BLS routine
##########################

    u = [0.0]*len(time)
    v = [0.0]*len(time)

    u = np.array(u)
    v = np.array(v)

    nbins = 200
        
    results = bls.eebls(time, flux, u, v, 1000.0, .3, .001, nbins, .001, .3)

#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

    print('Period from BLS is: ' + str(results[1]))
        
##########################
#SR plot
##########################

    SR_array = results[0]

    max_SR = max(SR_array)
    avg_SR = np.mean(SR_array)
    sd_SR = np.std(SR_array)

#normalize SR_array between 0 and 1
    SR_array = [i/max_SR for i in SR_array]

    #print(max_SR, avg_SR, sd_SR)

#np.arange(freq start, freq stop (must be calculated), df (freq step)
    freq = np.arange(.3, 1.3, .001, dtype=None)

    pylab.cla()

#UNCOMMENT TO SAVE SIGNAL RESIDUE PLOT

    
    xlabel('Frequency')
    ylabel('Signal Residue')    
    title('EPIC ' + str(targetname) + ' Merged Signal Residue')
    SR_freq_plot = plt.plot(freq, SR_array)
    #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5")
        
    #pylab.show(SR_freq_plot)

        
##########################
#Calculate SDE
########################## 

    SDE = (max_SR-avg_SR)/sd_SR
        
    print('Signal Detection Efficiency: ' + str(SDE))
    if SDE >= 6:
        print('SDE above 6: Transit Detectable')
        detect_BLS_per.append(results[1])
    else:
        print('SDE below 6: Transit Undetectable')
                

########
#data centering?
########



    

##########################
#Folding and Binning
##########################

    #pylab.cla()
            
    f0 = 1.0/results[1]
    n = len(time)
    ibi = np.zeros(nbins)
    y = np.zeros(nbins)
    phase = np.linspace(0.0, 1.0, nbins)
    
    for i in range(n):
        ph = u[i]*f0
        ph = ph-int(ph)
        j = int(nbins*ph)
        ibi[j] = ibi[j] + 1.0
        y[j] = y[j] + v[i]

    #binned and folded plot
    pylab.cla()
    plt.scatter(phase, y/ibi, s=3)
    plt.title("EPIC " + str(targetname) + " folded and binned LC")
    #pylab.show()

        

###########################
#BLS Overlay
###########################

        
    high = results[3]*results[4]
    low = high - results[3]
    
    fit = np.zeros(nbins) + high # H
    fit[results[5]:results[6]+1] = low # L
        
    plt.plot(phase, fit)
    plt.xlabel(r"Phase")
    plt.ylabel(r"Mean value of flux")
    plt.title("SDE " + str(SDE) + "; BLS period " + str(results[1]))
    plt.ylim(-.2, .2)

    #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/Pipeline/EPIC229227254-2sigclip/Period" + str(p) + "Radius" + str(r) + "/folded_pltPer" + str(p) + "Rad" + str(r) + 'BLSoverlay.png')

    #pylab.show()


##########################
#Detrending Flux and Levenberg-Marquardt Fitting:
#if SDE>6 
##########################

#start values are correct values

    if True: #SDE >= 6:
            
        #polynomial detrending
        fluxDetrend = polynomial(flux, order=5)

        #untrendy
        #fluxDetrend, ferr = untrendy.untrend(time, merged_flux)
            

        fitT = FitTransit()
        fitT.add_guess_star(rho=1.5)    
        fitT.add_guess_planet(
                period=1.0, impact=0.0, 
                T0=3.0, rprs=.2)#need a guess rprs
        fitT.add_data(time=time, flux=fluxDetrend)
                    

        vary_star = ['rho']      # not sure how to avoid free stellar parameters? ideally would not vary star at all
        vary_planet = (['period', 'rprs'])
                
        fitT.free_parameters(vary_star, vary_planet)
        fitT.do_fit()

        print(fitT.fitresultstellar.items())
        print(fitT.fitresultplanets.items())

        bestFstellar = fitT.fitresultstellar.items()
        bestFrho = bestFstellar[0][1]#Best Fit Rho
                    
        bestFplanet = fitT.fitresultplanets.items()
        bestFperiod = bestFplanet[0][1]['period']#Best Fit Period
        bestFrprs = bestFplanet[0][1]['rprs']#Best Fit rprs
                    

        fitT.print_results()
                    
#save figure
        fig = ktransit.plot_results(time,fluxDetrend,fitT.transitmodel)
        #fig.savefig('/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/ktransitfits_untrendy/' + str(targetname) + 'fitPer' + str(p) + 'Rprs' + str(r) + '.png')

        show()

##########################
#END LOOP
##########################

            
######################################
#UNCOMMMENT TO SAVE DIFF, SDE, PER, RAD TEXT OUTPUT
#CHANGE PATH TO APPROPRIATE BLANK TEXT FILE
######################################


#np.savetxt('/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/SDE_text/data213244700.txt', np.column_stack((detect_diff, detect_SDE, detect_BLS_per, detect_per, detect_rad)), fmt='%.5f')


##########################
#INITIAL LIGHT CURVE
##########################

pylab.cla()

xlabel('Time')
ylabel('Corrected Flux (normalized to 0)')
title('EPIC ' + targetname + ' Light Curve')
        

plot (time,flux)
show()
