

import bls, ktransit, math, numpy, pylab, os, batman

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


"""
def initialize(rad_min, rad_max, rad_step, per_min, per_max, per_step):

    #Depending on step size and limits:
    #detect_x and undetect_x are a maximum of 95*199 = 18905
    #Radius ratio from .05 to 1 with a step size of .01
    #Period from .1 to 20 days with a step size of .1

    rad_min = rad_min
    rad_max = rad_max
    rad_step = rad_step
    per_min = per_min
    per_max = per_max
    per_step = per_step
    
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

    #Limits and step size can be changed
    rad_range = numpy.arange(rad_min, rad_max, rad_step)
    per_range = numpy.arange(per_min, per_max, per_step)

"""

    
def makedir(pathoutput, r, p):
    if not os.path.exists(pathoutput + "/Period" + str(p) + "Radius" + str(r)):
        os.mkdirs(pathoutput + "/Period" + str(p) + "Radius" + str(r))

    return pathoutput + "/Period" + str(p) + "Radius" + str(r)


def importLC(pathLC):

##########################
#Import time, flux
##########################

#INSERT PATH OF LIGHT CURVE TEXT FILE FROM MIKULSKI ARCHIVE
    with open(pathLC) as f:
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


    time = time[0:len(time)-1]
    
    print(len(time), len(flux))
    
    return time, flux


def createktransit(time, p, r):
    
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
    T0=1.0,     # a transit mid-time  
        period=p, # an orbital period in days
        impact=0.0, # an impact parameter
        rprs=r,   # planet stellar radius ratio  
        ecosw=0.0,  # eccentricity vector
        esinw=0.0,
        occ=0.0)    # a secondary eclipse depth in ppm
        
    M.add_data(time=numpy.array(time[:])),

    tmod = M.transitmodel# the out of transit data will be 0.0 unless you specify zpt
    
    return tmod


def savektransitLC(currentdir, time, ktflux, p, r):

    xlabel('Time')
    ylabel('Corrected Flux (normalized to 0)')    
    title('KTransit Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
    graph = plt.plot(time,ktflux)
    pylab.savefig(str(currentdir) + "/ktransitLCPer" + str(p) + "Rad" + str(r) + '.png')


def injecttransit(ktflux, flux):
    
    return ktflux + flux


def savemergedLC(currentdir, EPICID, time, merged_flux, p, r):

    pylab.cla()
    xlabel('Time')        
    ylabel('Merged Flux (normalized to 0)')
    title(str(EPICID) + ' Merged Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
    merged_LC = plt.scatter(time, merged_flux)
    pylab.savefig(str(currentdir) + "/merged_fluxPer" + str(p) + "Rad" + str(r) + '.png')


def saveSRplot(currentdir, EPICID, SR_array, p, r):
    max_SR = max(SR_array)
    avg_SR = numpy.mean(SR_array)
    sd_SR = numpy.std(SR_array)

#normalize SR_array between 0 and 1
    SR_array = [i/max_SR for i in SR_array]

#numpy.arange(freq start, freq stop (must be calculated), df (freq step)
    freq = numpy.arange(.3, 1.3, .001, dtype=None)

    pylab.cla()

    xlabel('Frequency')
    ylabel('Signal Residue')    
    title(str(EPICID) + ' Merged Signal Residue: Period' + str(p) + ' Radius Ratio ' + str(r))
    SR_freq_plot = plt.plot(freq, SR_array)
    pylab.savefig(str(currentdir) + "/SR_freq_plotPer" + str(p) + "Rad" + str(r) + '.png')

    return max_SR, avg_SR, sd_SR


def SDEcount(max_SR, avg_SR, sd_SR, p, r, detect_count, detect_SDE, detect_per, detect_rad, detect_BLS_per, detect_diff, undetect_count, undetect_SDE, undetect_per, undetect_rad, undetect_diff, results1):

    SDE = (max_SR-avg_SR)/sd_SR
        
    print('Signal Detection Efficiency: ' + str(SDE))
    if SDE >= 6:
        print('SDE above 6: Transit Detectable')
        detect_count += 1
        detect_SDE.append(SDE)
        detect_per.append(p)
        detect_rad.append(r)
        detect_BLS_per.append(results1)
        detect_diff.append(results1-p)
    else:
        print('SDE below 6: Transit Undetectable')
        undetect_count += 1
        undetect_SDE.append(SDE)
        undetect_per.append(p)
        undetect_rad.append(r)
        undetect_diff.append(results1-p)

    return SDE, detect_count, detect_SDE, detect_per, detect_rad, detect_BLS_per, detect_diff, undetect_count, undetect_SDE, undetect_per, undetect_rad, undetect_diff

###changeresultsX
def savefoldbin (currentdir, SDE, time, nbins, results1, results3, results4, results5, results6, u, v, p, r):
    pylab.cla()
        
    f0 = 1.0/results1
    n = len(time)
    ibi = numpy.zeros(nbins)
    y = numpy.zeros(nbins)
    phase = numpy.linspace(0.0, 1.0, nbins)

    for i in range(n):
        ph = u[i]*f0
        ph = ph-int(ph)
        j = int(nbins*ph)
        ibi[j] = ibi[j] + 1.0
        y[j] = y[j] + v[i]

    plt.scatter(phase, y/ibi)
    #pylab.show()

    high = results3*results4
    low = high - results3
    
    fit = numpy.zeros(nbins) + high # H
    fit[results5:results6+1] = low # L

    plt.plot(phase, fit)
    plt.xlabel(r"Phase")
    plt.ylabel(r"Mean value of flux")
    plt.title("SDE " + str(SDE) + "; BLS period " + str(results1))
    plt.ylim(-.1, .1)

    pylab.savefig(str(currentdir) + "/folded_pltPer" + str(p) + "Rad" + str(r) + '.png')



def savetext(textdir, detect_diff, detect_SDE, detect_BLS_per, detect_per, detect_rad):
    numpy.savetxt(str(textdir), numpy.column_stack((detect_diff, detect_SDE, detect_BLS_per, detect_per, detect_rad)), fmt='%.5f')




    
def showinitialLC(time, flux):
    pylab.cla()

    xlabel('Time')
    ylabel('Corrected Flux (normalized to 0)')
    title('EPIC (name)')
        
    plot (time,flux)
    show()




def show2Dscatter(detect_per, detect_rad, rad_max, per_max):
    pylab.cla()
    plt.scatter(detect_per, detect_rad)
    axes.set_ylim([0,rad_max])
    axes.set_xlim([0,per_max])
    plt.xlabel('Orbital Period (days)')
    plt.ylabel('Planet Radius/Star Radius')
    plt.title('Detectable Transits (SDE > 6)')
    plt.show()

    

def show3Dscatter(detect_per, detect_rad, detect_diff, detect_SDE):
    pylab.cla()
    ax.scatter(detect_per, detect_rad, detect_diff, zdir='z', s=20, c=detect_SDE)
    ax.set_xlabel('Orbital Period')
    ax.set_ylabel('Planet Radius/Star Radius')
    ax.set_zlabel('BLS Period - ktransit Period')
    ax.set_title('Detectable Transits')

    show()



    

def main(rad_min, rad_max, rad_step, per_min, per_max, per_step):

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
    
    #Limits and step size can be changed
    rad_range = numpy.arange(rad_min, rad_max, rad_step)
    per_range = numpy.arange(per_min, per_max, per_step)

    rindex = -1
    pindex = -1

    for r in rad_range:
        rindex += 1
        pindex = 0
        for p in per_range:
            pindex += 1

        
            currentdir = makedir("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/Pipeline/EPIC229227250-2sigclip", r, p)

            time, flux = importLC('/Users/sheilasagear/OneDrive/K2_Research/Corrected Light Curve by EPIC ID (.txt)/hlsp_k2sff_k2_lightcurve_229227254-c06_kepler_v1_llc-default-aper.txt')

            flux = sigma_clip(flux, sigma=2, iters=1)

            ktflux = createktransit(time, p, r)
            
            savektransitLC(currentdir, time, ktflux, p, r)

            merged_flux = injecttransit(ktflux, flux)

            savemergedLC(currentdir, "EPICX", time, merged_flux, p, r)

        ##########################
        #BLS routine
        ##########################

            u = [0.0]*len(time)
            v = [0.0]*len(time)

            u = numpy.array(u)
            v = numpy.array(v)

            nbins = 200
        
            results = bls.eebls(time, merged_flux, u, v, 1000.0, .3, .001, nbins, .001, .3)


#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

            print('Planet radius/star radius:' + str(r))
            print('Period from BLS is: ' + str(results[1]))
            print('Set period is: ' + str(p))
        

            max_SR, avg_SR, sd_SR = saveSRplot(currentdir, 'EPICX', results[0], p, r)

        
            SDE, detect_count, detect_SDE, detect_per, detect_rad, detect_BLS_per, detect_diff, undetect_count, undetect_SDE, undetect_per, undetect_rad, undetect_diff = SDEcount(max_SR, avg_SR, sd_SR, p, r, detect_count, detect_SDE, detect_per, detect_rad, detect_BLS_per, detect_diff, undetect_count, undetect_SDE, undetect_per, undetect_rad, undetect_diff, results[1])


            savefoldbin(currentdir, SDE, time, nbins, results[1], results[3], results[4], results[5], results[6], u, v, p, r)

        ##########################
        #END LOOP
        ##########################


    savetext('/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/SDE_text/data213244700.txt', detect_diff, detect_SDE, detect_BLS_per, detect_per, detect_rad)

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


    showinitialLC(time, flux)

    show2Dscatter(detect_per, detect_rad, rad_max, per_max)

    show3Dscatter(detect_per, detect_rad, detect_diff, detect_SDE)


main(.05, .5, .05, .1, 5, .1)
