
import bls, ktransit, math, pylab, os, batman

import untrendy

from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
fig = plt.figure()
axes = plt.gca()

from mpl_toolkits.mplot3d import Axes3D
axe = Axes3D(plt.gcf())
ax = fig.add_subplot(111, projection='3d')

from matplotlib import pyplot
import pylab
from pylab import xlabel, ylabel, title, plot, show, ylim, xlim

from astropy.stats import sigma_clip

fige = pylab.figure()
ax = Axes3D(fig)

from ktransit import FitTransit

from scipy import signal

import obspy
from obspy.signal.detrend import polynomial

##########################
#INITIALIZATION
##########################

#Depending on step size and limits:
#detect_x and undetect_x are a maximum of 95*199 = 18905
#Radius ratio from .05 to 1 with a step size of .01
#Period from .1 to 20 days with a step size of .1

detect_count = 0
detect_SDE = []*18905
detect_per = []*18905
detect_rad = []*18905
detect_BLS_per = []*18905

undetect_count = 0
undetect_SDE = []*18905
undetect_per = []*18905
undetect_rad = []*18905

detect_diff = []*18905
undetect_diff = []*18905

diff = [[0 for x in range(96)] for y in range(200)]

rindex = -1
pindex = -1

maskindex = []


rad_diff_values = [[] for y in range(96)]
per_diff_values = [[] for y in range(200)]

###########################################
#These limits and step size can be changed!
rad_range = np.arange(.3, .4, .05)
per_range = np.arange(3, 4, .5)
###########################################


total_rad = len(rad_range)-1
total_per = len(per_range)-1
total = total_rad*total_per

recovered = 0

##########################
#BEGIN LOOP
##########################

#LATER: loop through C6 and C7 light curves; open and inject range of transits into each 

#count number of files being graphed
filecount = sum(1 for line in open('all_txt_files.txt'))

#opens list of file names, creates list of file names
with open('all_txt_files.txt') as f:
    content = f.readlines()

#converts each file name to string, removes \n from end, opens each object's txt file (loop)
#and stores (x,y) in list 'data'
tempindex = 0
while (tempindex < filecount):
    file = str(content[tempindex])
    file = file.strip('\n')
    file = file.strip('\r')
    with open(file) as f:
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
        x_list = data_final[0::2]
        y_list = data_final[1::2]

#normalizes flux values to 1.0 (at avg of flux values)
        y_count = len(y_list)
        sumy_list = sum(y_list)
        y_avg = sumy_list/y_count
        
        y_list = [i/y_avg for i in y_list]
        
        xlabel('Time')
        ylabel('Corrected Flux (normalized to 1.0)')

#targets titled by EPIC names
        targetname = file[25:]
        targetname = targetname[:-35]
        #title('EPIC ' + targetname + ' Light Curve')


    for r in rad_range:
        rindex += 1
        pindex = 0
        for p in per_range:
            pindex += 1

######################################
#CREATE DIRECTORY FOR SAVED PLOTS
######################################

            #if not os.path.exists("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/Pipeline/EPIC229227254-2sigclip/Period" + str(p) + "Radius" + str(r)):
                #os.makedirs("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/Pipeline/EPIC229227254-2sigclip/Period" + str(p) + "Radius" + str(r))        


##########################
#Import time, flux
##########################

#INSERT PATH OF LIGHT CURVE TEXT FILE FROM MIKULSKI ARCHIVE
            #with open(r'/Users/sheilasagear/OneDrive/K2_Research/Corrected Light Curve by EPIC ID (.txt)/hlsp_k2sff_k2_lightcurve_229227248-c06_kepler_v1_llc-default-aper.txt') as f:
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

            data_final = [float(i) for i in data2]
    
#defines x and y values by index of 'data'
            time = data_final[0::2]
            flux = data_final[1::2]
    
    
#normalize flux values to 1.0
            flux_count = len(flux)
            sumflux = sum(flux)
            flux_avg = sumflux/flux_count
    
            flux = [i/flux_avg for i in flux]
#normalize to 0.0
            flux = [i-1 for i in flux]


#SIGMA CLIPPING
            flux = sigma_clip(flux, sigma=3, iters=1)

#uncomment if extra time stamp
            time.pop()




##########################
#IMPORT CDPP FROM FITS
##########################

            f = fits.open('/Users/sheilasagear/Downloads/hlsp_k2sff_k2_lightcurve_229227253-c06_kepler_v1_llc.fits')

            bestaper = f[1].data
            besthead = f[1].header

            CDPP = besthead['QCDPP6']

            print(CDPP)


        
##########################
#Create ktransit Data: CODE FROM GITHUB
##########################

            num_time = len(time)

            M = ktransit.LCModel()
            M.add_star(
            rho=1.5, # mean stellar density in cgs units
            ld1=0.2, # ld1--4 are limb darkening coefficients 
            ld2=0.4, # if only ld1 and ld2 are non-zero then a quadratic limb darkening law is used
            ld3=0.0, # if all four parameters are non-zero we use non-linear flavour limb darkening
            ld4=0.0, 
            dil=0.0, # a dilution factor: 0.0 -> transit not diluted, 0.5 -> transit 50% diluted
            zpt=0.0  # a photometric zeropoint, incase the normalisation was wonky
    )
            M.add_planet(
            T0=3.0,     # a transit mid-time  
            period=p, # an orbital period in days
            impact=0.0, # an impact parameter
            rprs=r,   # planet stellar radius ratio  
            ecosw=0.0,  # eccentricity vector
            esinw=0.0,
            occ=0.0)    # a secondary eclipse depth in ppm
            

            M.add_data(time=np.array(time[:])),

            tmod = M.transitmodel# the out of transit data will be 0.0 unless you specify zpt

            #pylab.cla()

            #plot(time, tmod)
            #show()
        
######################################
#UNCOMMENT TO SAVE KTRANSIT LIGHT CURVE
######################################

            """
            xlabel('Time')
            ylabel('Corrected Flux (normalized to 0)')    
            title('KTransit Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
            graph = plt.plot(M.time,tmod)
            pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5/Period" + str(p) + "Radius" + str(r) + "/ktransitLCPer" + str(p) + "Rad" + str(r) + '.png')
            """

            #pylab.show(graph)


##########################
#Inject ktransit LC into K2 data
##########################


            merged_flux = tmod + flux
            #plus flux

            #pylab.cla()

######################################
#UNCOMMENT TO SAVE MERGED LIGHT CURVE
######################################

            """
            xlabel('Time')        
            ylabel('Merged Flux (normalized to 0)')
            title('EPIC 212820594 Merged Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
            merged_LC = plt.scatter(time, merged_flux)
            pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5/Period" + str(p) + "Radius" + str(r) + "/merged_fluxPer" + str(p) + "Rad" + str(r) + '.png')
            """

            #plt.show(merged_LC)


##########################
#BLS routine
##########################

            u = [0.0]*len(time)
            v = [0.0]*len(time)

            u = np.array(u)
            v = np.array(v)

            nbins = 200
        
            results = bls.eebls(time, merged_flux, u, v, 1000.0, .3, .001, nbins, .001, .3)

#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

            print('Planet radius/star radius:' + str(r))
            print('Period from BLS is: ' + str(results[1]))
            print('Set period is: ' + str(p))
        
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

            #pylab.cla()

#UNCOMMENT TO SAVE SIGNAL RESIDUE PLOT

            """
            xlabel('Frequency')
            ylabel('Signal Residue')    
            title('EPIC 212820594 Merged Signal Residue: Period' + str(p) + ' Radius Ratio ' + str(r))
            SR_freq_plot = plt.plot(freq, SR_array)
            pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5/Period" + str(p) + "Radius" + str(r) + "/SR_freq_plotPer" + str(p) + "Rad" + str(r) + '.png')
            """
        
            #pylab.show(SR_freq_plot)

        
##########################
#Calculate SDE
########################## 

            SDE = (max_SR-avg_SR)/sd_SR
        
            print('Signal Detection Efficiency: ' + str(SDE))
            if SDE >= 6:
                print('SDE above 6: Transit Detectable')
                detect_count += 1
                detect_SDE.append(SDE)
                detect_per.append(p)
                detect_rad.append(r)
                detect_BLS_per.append(results[1])
                
                detect_diff.append(results[1]-p)
            else:
                print('SDE below 6: Transit Undetectable')
                undetect_count += 1
                undetect_SDE.append(SDE)
                undetect_per.append(p)
                undetect_rad.append(r)
                undetect_diff.append(results[1]-p)
                

            diff[pindex][rindex] = (results[1]-p)
            print('Difference in period is '+ str(diff[pindex][rindex]))
            print('\n')
        
            rad_diff_values[rindex].append(diff[pindex][rindex])
        
            per_diff_values[pindex].append(diff[pindex][rindex])
            

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
            #plt.scatter(phase, y/ibi)
            #pylab.show()

        

###########################
#BLS Overlay
###########################

        
            high = results[3]*results[4]
            low = high - results[3]
        
            fit = np.zeros(nbins) + high # H
            fit[results[5]:results[6]+1] = low # L
        
            #plt.plot(phase, fit)
            #plt.xlabel(r"Phase")
            #plt.ylabel(r"Mean value of flux")
            #plt.title("SDE " + str(SDE) + "; BLS period " + str(results[1]))
            #plt.ylim(-.1, .1)

            #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/Pipeline/EPIC229227254-2sigclip/Period" + str(p) + "Radius" + str(r) + "/folded_pltPer" + str(p) + "Rad" + str(r) + '.png')

            #pylab.show()


##########################
#Detrending Flux and Levenberg-Marquardt Fitting:
#if SDE>6 
##########################

#start values are correct values

            if True: #SDE >= 6:
            
            #polynomial detrending
                #mergedfluxDetrend = polynomial(merged_flux, order=5)

            #untrendy
                mergedfluxDetrend, ferr = untrendy.untrend(time, merged_flux)
            

                fitT = FitTransit()
                fitT.add_guess_star(rho=1.5)    
                fitT.add_guess_planet(
                    period=p, impact=0.0, 
                    T0=3.0, rprs=r)#need a guess rprs
                fitT.add_data(time=time, flux=mergedfluxDetrend)


                vary_star = ['rho']      # not sure how to avoid free stellar parameters? ideally would not vary star at all
                vary_planet = (['period', 'rprs'])

                fitT.free_parameters(vary_star, vary_planet)
                fitT.do_fit()

                fitT.print_results()

            #saves in current dir
                #fig = ktransit.plot_results(time,mergedfluxDetrend,fitT.transitmodel)
                #fig.savefig('fitPer' + str(p) + 'Rprs' + str(r) + '.png')



##########################
#CALCULATE SNR
##########################

                SNR = (r**2)*(((time[-1]-time[0])/p)**(1/2))*(1/CDPP)
                print(SNR)

            
##########################
#This needs to be changed
#based on what we consider
#to be a recovered transit:
#if best-fit period & rprs
#falls within a certain error?
            #if True:#period
                #if True:#radius
                    #recovered += 1
            recovered += 1
##########################


            
        tempindex+=1

##########################
#END LOOP
##########################








frac_recovered = recovered/total

#print(frac_recovered)


            
######################################
#UNCOMMMENT TO SAVE DIFF, SDE, PER, RAD TEXT OUTPUT
#CHANGE PATH TO APPROPRIATE BLANK TEXT FILE
######################################


#np.savetxt('/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/SDE_text/data213244700.txt', np.column_stack((detect_diff, detect_SDE, detect_BLS_per, detect_per, detect_rad)), fmt='%.5f')


#deletes extras from arrays

for i in range(len(detect_SDE)):
    if detect_SDE[i] == 0:
        del detect_SDE[i]

for i in range(len(detect_per)):
    if detect_per[i] == 0:
        del detect_per[i]

for i in range(len(detect_rad)):
    if detect_rad[i] == 0:
        del detect_rad[i]

for i in range(len(undetect_SDE)):
    if undetect_SDE[i] == 0:
        del undetect_SDE[i]

for i in range(len(undetect_per)):
    if undetect_per[i] == 0:
        del undetect_per[i]

for i in range(len(undetect_rad)):
    if undetect_rad[i] == 0:
        del undetect_rad[i]

#print(detect_count, detect_SDE, detect_per, detect_rad)
#print(undetect_count, undetect_SDE, undetect_per, undetect_rad)

##########################
#INITIAL LIGHT CURVE
##########################
"""
pylab.cla()

xlabel('Time')
ylabel('Corrected Flux (normalized to 0)')
title('EPIC (name)')
        

plot (time,flux)
show()
"""

##########################
#2D SCATTER ALL LIGHT CURVES SDE>6
##########################

"""
pylab.cla()
plt.scatter(detect_per, detect_rad)
axes.set_ylim([0,max(rad_range)])
axes.set_xlim([0,max(per_range)])
plt.xlabel('Orbital Period (days)')
plt.ylabel('Planet Radius/Star Radius')
plt.title('Detectable Transits (SDE > 6)')
plt.show()
"""

##########################
#3D SCATTER: ERROR, PERIOD, RADIUS FOR SDE>6
##########################

#print(len(detect_per))
#print(len(detect_rad))
#print(len(detect_diff))

"""
pylab.cla()
ax.scatter(detect_per, detect_rad, detect_diff, zdir='z', s=20, c=detect_SDE)
ax.set_xlabel('Orbital Period')
ax.set_ylabel('Planet Radius/Star Radius')
ax.set_zlabel('BLS Period - ktransit Period')
ax.set_title('Detectable Transits')

show()
"""

##########################
#3D SCATTER: ERROR, PERIOD, RADIUS FOR RECOVERED TRANSITS
##########################

"""
pylab.cla()
ax.scatter(per_range, rad_range, percent_recovered, zdir='z', s=20, c='blue')
ax.set_xlabel('Orbital Period')
ax.set_ylabel('Radius Ratio')
ax.set_zlabel('Fraction of Transits Recovered')
ax.set_title('Recovered Transits?')

show()

"""

print('SNR = ' + str(SNR))
print('% recovered for LC = ' +str(frac_recovered))
        
pylab.cla()
plot(SNR, frac_recovered)
show()
