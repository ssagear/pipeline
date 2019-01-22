from __future__ import division

import numpy as np

import bls, ktransit, pylab, os, untrendy, math, commands

from astropy.io import fits

import scipy.constants as sc

from astropy.stats import sigma_clip

import matplotlib.pyplot as plt

import matplotlib

import lightkurve as lk

import csv

########################################
#time and flux arrays must be np arrays
########################################

natural_val = np.linspace(0, 1.415, num=8, endpoint=True)
natural_val = natural_val.round(decimals=2)
per_segments = [10**x for x in natural_val]
per_segments = [i.round(decimals=2) for i in per_segments]
rprs_segments = [x / 20 for x in range(8)]


def read_line(path, line=1):
    return commands.getoutput('head -%s %s | tail -1' % (line, path))

def loguniform(low=float(1), high=float(26)):
    return 10**(np.random.uniform(np.log10(low), np.log10(high)))


def randomInc(n):
    randNums = np.random.uniform(low=0.5, high=1, size=n)
    incs = np.arccos(2*randNums - 1)
    return incs

def data_K2(fitspath):

    targetname = fitspath[43:]
    targetname = targetname[:-18]

    print('TARGET ' + str(targetname))

    campaign_no = 18

    try:
        tpf = lk.search_targetpixelfile(str(targetname), campaign=18).download()
    except IOError:
        return 0, 0, 0, 0

    lc = tpf.to_lightcurve().normalize().flatten()
    time = lc.time
    flux = lc.flux

    flux_count = len(flux)
    sumflux = sum(flux)
    flux_avg = sumflux/flux_count

    flux = [i/flux_avg for i in flux]

    #normalize to 0.0
    flux = [i-1 for i in flux]

    flux = np.asarray(flux)
    time = np.asarray(time)

    #SIGMA CLIPPING
    flux = sigma_clip(flux, sigma=3, iters=1)
    flux = flux.filled(fill_value=0)

    return time, flux, targetname, campaign_no

def CDPP(fitspath):

    f = fits.open(fitspath)
    QCDPP = f[1].header['QCDPP6']
    return QCDPP

def createT(time, T0):#time array, transit mid-time, impact parameter

    #random period and radius: uniformly weighted
    #you can inject a specific period if you want
    period = loguniform(low=float(1), high=float(26))
    print('inj period: ' + str(period))
    rprs = np.random.uniform(low=.01, high=.4)
    rprs = float(rprs)
    print('inj rprs: ' + str(rprs))

    Mc = 9e28#kilograms - mass of L dwarf
    Rs = 7.1492e7#meters - radius of L dwarf
    periodsec = 86400*period
    a1 = sc.G*Mc*periodsec**2
    a2 = (4*sc.pi)**2
    a = np.cbrt(a1/a2)#semimajor axis, meters

    i = randomInc(1)
    impact = (a*math.cos(i))/Rs
    #impact = 0.0

    total_per.append(period)
    total_rprs.append(rprs)

    print('inj impact: ' + str(impact))

    #this is Tom Barclay's ktransit package I use for injection and fitting (https://github.com/mrtommyb/ktransit)
    M = ktransit.LCModel()
    M.add_star(
    rho=1.5, # mean stellar density in cgs units
    ld1=0.2, # ld1--4 are limb darkening coefficients
    ld2=0.4, # assuming quadratic limb darkening
    ld3=0.0,
    ld4=0.0,
    dil=0.0, # a dilution factor: 0.0 -> transit not diluted, 0.5 -> transit 50% diluted
    zpt=0.0  # photometric zeropoint
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

    return tmod, period, rprs, impact

#tmod = the noiseless LC with the transits you want to inject (flux array)


def addT(flux, tmod):#flux array, transit flux array you want to inject
    merged_flux = flux + tmod
    return merged_flux

#dt=x (2 in trappist ntbk) - untrendy median parameter
def detrend(time, merged_flux):
    trend = untrendy.median(time, merged_flux, dt=0.6)
    mergedfluxD = np.zeros(len(time))
    for i in range(len(time)):
        mergedfluxD[i] = merged_flux[i]-trend[i]

    return mergedfluxD


#dan foreman-mackey's python BLS from Kovacs' fortran subroutine
#https://github.com/dfm/python-bls
def BLS(time, mergedfluxD):#time array, detrended flux arr (with transits injected)

    u = [0.0]*len(time)
    v = [0.0]*len(time)
    u = np.array(u)
    v = np.array(v)

    #time, flux, u, v, number of freq bins (nf), min freq to test (fmin), freq spacing (df), number of bins (nb), min transit dur (qmi), max transit dur (qma)

    nf = 1000.0
    fmin = .035
    df = 0.001
    nbins = 300
    qmi = 0.001
    qma = 0.3

    results = bls.eebls(time, mergedfluxD, u, v, nf, fmin, df, nbins, qmi, qma)

    #RESULTS:
    #power, best_period, best_power, depth, q, in1, in2
    #0      1            2           3      4  5    6

    SR_array = results[0]
    max_SR = max(SR_array)
    avg_SR = np.mean(SR_array)
    sd_SR = np.std(SR_array)
    #normaze SR_array between 0 and 1
    SR_array = [i/max_SR for i in SR_array]

    freq = fmin + np.arange(nf)*df

    #Signal Detection Efficiency
    SDE = (max_SR-avg_SR)/sd_SR

    #depth
    high = results[3]*results[4]
    low = high - results[3]
    fit = np.zeros(nbins) + high # H
    fit[results[5]:results[6]+1] = low # L
    depth = high - low

    return results, SR_array, freq, SDE, depth
    #plotting freq vs. SR_array gives you a power spectrum
    #depth is the BLS transit depth



def fitT(time, mergedfluxD, guessper, guessrprs, T0):#time arr, merged (transits added) flux arr detrended, guess period and rp/rs from BLS, transit mid-time

    fitT = ktransit.FitTransit()
    fitT.add_guess_star(rho=1.5)
    fitT.add_guess_planet(
    period=guessper, impact=0.0,
    T0=T0, rprs=guessrprs)
    fitT.add_data(time=time, flux=mergedfluxD)

    vary_star = []      # free stellar parameters
    vary_planet = (['period',       # free planetary parameters
        'rprs'])                # free planet parameters are the same for every planet you model

    fitT.free_parameters(vary_star, vary_planet)
    fitT.do_fit()                   # run the fitting

    bestFplanet = fitT.fitresultplanets.items()
    bestFperiod = bestFplanet[0][1]['period']#Best Fit Period
    bestFrprs = bestFplanet[0][1]['rprs']#Best Fit rprs

    #fitT.print_results()

    return bestFperiod, bestFrprs




def is_recovered(period, fitper, rprs, fitrprs):#true (injected) period, L-M fitting period, true rp/rs, L-M fitting rp/rs
    if abs(fitper-period) < .05*period and abs(fitrprs-rprs) < .05*rprs:
        recovered = True
        recoveredarr.append(True)
        recovered_period.append(period)
        recovered_rprs.append(rprs)
        for p_i in range(len(per_segments)-1):
            for r_i in range(len(rprs_segments)-1):
                if period > per_segments[p_i] and period < per_segments[p_i+1]:
                    if rprs > rprs_segments[r_i] and rprs < rprs_segments[r_i+1]:
                        r_list[r_i][p_i] += 1
    else:
        recovered = False
        recoveredarr.append(False)
        unrecovered_period.append(period)
        unrecovered_rprs.append(rprs)
        for p_i in range(len(per_segments)-1):
            for r_i in range(len(rprs_segments)-1):
                if period > per_segments[p_i] and period < per_segments[p_i+1]:
                    if rprs > rprs_segments[r_i] and rprs < rprs_segments[r_i+1]:
                        u_list[r_i][p_i] += 1

    return recovered


def create_params(path):

    path = path

    time, flux, name, cnum = data_K2(path)

    print('TARGET: ' + str(name))

    QCDPP = CDPP(path)


    T0 = time[0]+((time[-1]-time[0])/2)
    tmod, period, rprs, impact = createT(time, T0)

    return time, flux, name, cnum, QCDPP, tmod, period, rprs, impact, T0

def inject(time, flux, tmod, T0):
    merged_flux = addT(flux, tmod)

    mergedfluxD = detrend(time, merged_flux)

    results, power, freq, SDE, depth = BLS(time, mergedfluxD)

    guessper = results[1]
    #print('guess period = ' + str(guessper))
    guessrprs = np.sqrt(depth)
    #print('guess rprs = ' + str(guessrprs))

    fitper, fitrprs = fitT(time, mergedfluxD, guessper, guessrprs, T0)

    isrec = is_recovered(period, fitper, rprs, fitrprs)
    print(isrec)

    return(isrec)



def plot(recovered_period, recovered_rprs, total_per, total_rprs, targetname):

    a = np.arange(26)

    ybins = [x / 20 for x in range(8)]


    print('total per ' + str(total_per))
    print('total_rprs ' + str(total_rprs))
    counts, _, _ = np.histogram2d(recovered_period, recovered_rprs, bins=(per_segments, ybins))
    counts_tot, _, _ = np.histogram2d(total_per, total_rprs, bins=(per_segments, ybins))

    for i in range(len(counts.T)):
        for j in range(len(counts.T[i])):
            counts.T[i][j] = counts.T[i][j]/counts_tot.T[i][j]
            if np.isnan(counts.T[i][j]):
                counts.T[i][j] = 0
    print(counts.T)
    countsT = np.flip(counts.T, axis=0)

    with open("/Users/sheilasagear/Dropbox/K2/heatmap_inclination_csv/heatmaps_inprogress/inclination/heatmap" + str(targetname) + "_fullplanetsXX.csv", "w") as heatmap_csv:
        np.savetxt(heatmap_csv, countsT, fmt='%f', delimiter=',')

    matplotlib.rcParams['xtick.minor.size'] = 0
    matplotlib.rcParams['xtick.minor.width'] = 0

    fig, ax = plt.subplots()
    heatmap = ax.pcolormesh(per_segments, ybins, counts.T, cmap='Blues_r')
    ax.set_xscale('log')

    ax.xaxis.set_ticks(per_segments)
    ax.xaxis.set_ticklabels(per_segments)

    ax.yaxis.set_ticks(ybins)
    ax.yaxis.set_ticklabels(ybins)

    ax.set_title('Fraction of Transits Recovered (EPIC ' + str(targetname) + ')')
    ax.set_xlabel('Period (days)')
    ax.set_ylabel('Radius (stellar radii)')


    cbar = plt.colorbar(heatmap)
    heatmap.set_clim(0.0, 1.0)
    cbar.set_label('Fraction recovered', rotation=270, labelpad=13)


    plt.savefig('/Users/sheilasagear/Dropbox/ssagear_k2/plots/in_progress/heatmaps_inprogress/inclination/heatmap_' + str(targetname) + '_fullplanetsXX.png')


path_to_list = '/Users/sheilasagear/Dropbox/K2/K2_fits/C18-search.txt'

fitslist = open(path_to_list)

filecount = sum(1 for line in fitslist)

for f in range(filecount):
    fitspath = read_line(path_to_list, f+1)

    time, flux, targetname, campaign_no = data_K2(fitspath)
    if isinstance(time, (int, long)):
        continue

    fluxD = detrend(time, flux)

    #QCDPP = CDPP(fitspath)

    results, SR_array, freq, SDE, depth = BLS(time, fluxD)

    print('BLS PERIOD: ' + str(results[1]))


    with open('/Users/sheilasagear/Dropbox/ssagear_k2/planet_search/bls_power.csv', 'a') as power_csv:
        writer = csv.writer(power_csv, delimiter=',')
        writer.writerow([campaign_no, targetname, SDE])

    plt.cla()

    plt.plot(freq, SR_array, color='black')
    plt.xlabel('Frequency (1/days)')
    plt.ylabel('Power')
    plt.ylim(0.2, 1.25)
    plt.title('EPIC ' + str(targetname) + ' BLS period: ' + str(results[1]))
    plt.savefig('/Users/sheilasagear/Dropbox/ssagear_k2/plots/C18_planet_search/BLS_power_spectra/' + str(targetname) + '_bls.png')

    plt.cla()

    lc = lk.lightcurve.LightCurve(time, fluxD)
    lc = lc.fold(period=results[1])
    cdpp = lc.estimate_cdpp()
    t, f = lc.time, lc.flux
    plt.scatter(t, f, color='black', s=2)
    plt.xlabel('Phase')
    plt.ylabel('Flux')
    plt.title('EPIC ' + str(targetname) + ' Folded LC: ' + str(results[1]))
    plt.savefig('/Users/sheilasagear/Dropbox/ssagear_k2/plots/C18_planet_search/BLS_power_spectra/' + str(targetname) + '_phase_fold.png')






#h
