import ktransit
import matplotlib.pyplot as plt
import numpy as np

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
        period=1.0, # an orbital period in days
        impact=0.1, # an impact parameter
        rprs=0.1,   # planet stellar radius ratio  
        ecosw=0.0,  # eccentricity vector
        esinw=0.0,
        occ=0.0)    # a secondary eclipse depth in ppm

M.add_planet() # you can add as many planets as you like (up to 10)

M.add_data(
        time=np.arange(0,10,0.0188),                                 # timestamps to evaluate the model on
        itime=np.zeros_like(np.arange(0,10,0.0188))+0.0188 )      # integration time of each timestamp

tmod = M.transitmodel # the out of transit data will be 0.0 unless you specify zpt
plt.plot(M.time,tmod)





####################

from ktransit import FitTransit
    # uncertainty on the data

fitT = FitTransit()
fitT.add_guess_star(rho=7.0)    
fitT.add_guess_planet(
        period=365.25, impact=0.0, 
        T0=0.0, rprs=0.009155)
fitT.add_data(time=M.time, flux=tmod)

vary_star = ['rho', 'zpt']      # free stellar parameters
vary_planet = (['period',       # free planetary parameters
        'T0', 'impact', 
        'rprs'])                # free planet parameters are the same for every planet you model

fitT.free_parameters(vary_star, vary_planet)
fitT.do_fit()                   # run the fitting

fitT.print_results()  


fig = ktransit.plot_results(M.time,tmod,fitT.transitmodel)
fig.savefig('my_beautiful_fit.png') # or fig.show()
