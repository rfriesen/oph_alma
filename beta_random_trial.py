import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib
import matplotlib.pyplot as plt

def planck(freq,tdust):
    denom = np.exp(c.h.cgs * freq/(c.k_B.cgs*tdust))-1.
    function = 2. *c.h.cgs * freq**3./c.c.cgs**2./denom
    return function

def find_beta(freq1,freq2,flux1,flux2,tdust):
    flux_ratio = np.log(flux2/flux1)
    planck_ratio = np.log(planck(freq2,tdust).cgs/planck(freq1,tdust).cgs)
    freq_ratio = np.log(freq2/freq1)
    beta = (flux_ratio - planck_ratio)/freq_ratio
    return beta

def beta_random_trial(freqs,fluxes,sigmas,tdust,ntrials):
    beta_array = np.zeros(ntrials)
    for i in range(100):
        new_fluxes = np.zeros(len(fluxes))
        for j in range(len(fluxes)):
            new_fluxes[j] = np.random.normal(loc=fluxes[j].cgs.value,scale=sigmas[j].cgs.value)
        beta = find_beta(freqs[0],freqs[1],new_fluxes[0],new_fluxes[1], tdust)
        beta_array[i] = beta

    print np.mean(beta_array)
    print np.std(beta_array)

# Measured source fluxes
# Assume 10% calibration uncertainty
fluxes_SM1 = [0.02875*u.Jy,0.1289*u.Jy,0.317*u.Jy]
freqs_SM1  = [107.*u.GHz,221.*u.GHz,359.*u.GHz]
#sigma_SM1   = [0.98*u.mJy,0.8*u.mJy,1.*u.mJy]
sigma_SM1 = [0.1* x for x in fluxes_SM1]

# Issue: outside of primary beam. Not corrected in model. 
fluxes_SM1_2 = [0.00756*u.Jy,0.0416*u.Jy]
freqs_SM1_2  = [107.*u.GHz,221.*u.GHz]
#sigma_SM1_2   = [0.00086*u.Jy,0.0009*u.Jy]
sigma_SM1_2 = [0.1* x for x in fluxes_SM1_2]

fluxes_1623 = [0.031*u.Jy,0.1423*u.Jy]
freqs_1623  = [107.*u.GHz,221.*u.GHz]
#sigma_1623  = [1.2*u.mJy,1.0*u.mJy]
sigma_1623 = [0.1* x for x in fluxes_1623]

fluxes_16263_2422_2 = [4.3*u.mJy,13.6*u.mJy] # 3mm, 1mm
freqs_16263 = [107.*u.GHz,221.*u.GHz]
#sigma_16263_2422_2 = [1.4*u.mJy,0.08*u.mJy]
sigma_16263_2422_2 = [1.4*u.mJy,13.6*0.1*u.mJy]

td = 30.*u.K

ntrials = 100
print 'SM1'
beta_random_trial(freqs_SM1,fluxes_SM1,sigma_SM1,td,ntrials)
print 'N3-mm'
beta_random_trial(freqs_SM1_2,fluxes_SM1_2,sigma_SM1_2,td,ntrials)
print 'GSS 30 IRS3'
beta_random_trial(freqs_1623,fluxes_1623,sigma_1623,td,ntrials)
print 'GSS 30 IRS1'
beta_random_trial(freqs_16263,fluxes_16263_2422_2,sigma_16263_2422_2,td,ntrials)

