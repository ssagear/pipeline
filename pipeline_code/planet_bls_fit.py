import lightkurve as lk
import matplotlib.pyplot as plt
import bls
import numpy as np

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

    return results, SR_array, freq, SDE, depth, high, low, u, v

tpf = lk.search_targetpixelfile(211435416, campaign=18).download()

lc = tpf.to_lightcurve().normalize().flatten()
time = lc.time
flux = lc.flux

lc = lc.fold(period=1.757469244)
t, f = lc.time, lc.flux

results, SR_array, freq, SDE, depth, high, low, u, v = BLS(time, flux)

plt.scatter(t,f, s=2, color='k')
plt.ylim(0.4,1.6)
plt.title('EPIC 211435416')
plt.show()
