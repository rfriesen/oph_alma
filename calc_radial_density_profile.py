'''
Script to calculate radial density profile from continuum image
Given input R.A., Dec. 
Starting with SM1N ALMA data but try to generalize. 
May want to do elliptical profile later. 
Use standard dense core dust parameters to determine mass
'''
import numpy as np
from astropy.coordinates import Angle,SkyCoord
from astropy.coordinates import ICRS, Galactic
import astropy.units as u
import astropy.constants as c
from astropy.wcs import WCS
from astropy.io import fits
from regions import CircleSkyRegion
import aplpy
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

def binArray(data, r, binsize):
    '''
    Return radial profile binned in terms of [binsize] pixels
    Use data as weights to calculate. 
    '''
    new_r = r/binsize
    new_r_int = new_r.astype(int)
    tbin = np.bincount(new_r_int.ravel(),data.ravel())
    t2bin = np.bincount(new_r_int.ravel(),data.ravel()**2.)
    nr = np.bincount(new_r_int.ravel())
    radial_profile=tbin/nr
    radial_sdev=np.sqrt(t2bin/nr-radial_profile**2.)
    rbin = np.bincount(new_r_int.ravel(),new_r.ravel())
    r_array = rbin/nr
    return r_array, radial_profile, radial_sdev


def radial_profile(hdu, peak_loc, radial_bin):
    wcs = WCS(hdu.header)
    # Get pixel location for peak emission
    x_peak, y_peak = peak_loc.to_pixel(wcs)
    # Distance from each point to peak emission in pixels
    # ALMA hdus have extra axes
    data = np.squeeze(hdu.data)
    x_arr, y_arr = np.indices(data.shape)
    r = np.sqrt((x_arr - x_peak)**2. + (y_arr - y_peak)**2.)
    radius, r_profile, r_sdev = binArray(data, r, radial_bin)
    return radius, r_profile, r_sdev

def get_radial_distances(hdu,peak_loc):
    wcs = WCS(hdu.header)
    # Get pixel location for peak emission
    x_peak, y_peak = peak_loc.to_pixel(wcs)
    # Distance from each point to peak emission in pixels
    # ALMA hdus have extra axes
    data = np.squeeze(hdu.data)
    x_arr, y_arr = np.indices(data.shape)
    r = np.sqrt((x_arr - x_peak)**2. + (y_arr - y_peak)**2.)
    return r

def gaussian_profile(r,sigma,A):
    # Assume r, sigma both in arcsec
    gaussian = A* np.exp(-1.*(r-0)**2./(2.*sigma**2.))
    return gaussian

def gauss(x,*p):
    A, sigma = p
    return A*np.exp((-(x)**2.)/(2.*sigma**2.))

def convol_gauss(x,*p):
    A, sigma = p
    beam_reff = np.sqrt(hdu.header['BMAJ']*hdu.header['BMIN'])*u.deg/2.3548
    beam_gaus = np.exp(-1.*(x-0)**2./(2.*beam_reff.to(u.arcsec).value**2.))
    #print x.shape, beam_gaus.shape
    return np.convolve(A*np.exp((-(x)**2.)/(2.*sigma**2.)),beam_gaus,mode='same')

def power_law(x,*p):
    A, s = p
    beam_reff = np.sqrt(hdu.header['BMAJ']*hdu.header['BMIN'])*u.deg/2.3548
    beam_gaus = np.exp(-1.*(x-0)**2./(2.*beam_reff.to(u.arcsec).value**2.))
    return np.convolve(A*x**((-1.)*s),beam_gaus,mode='same')
    #return A*x**(-1.*s)

def planck(freq,tdust):
    denom = np.exp(c.h.cgs * freq/(c.k_B.cgs*tdust))-1.
    function = 2. *c.h.cgs * freq**3./c.c.cgs**2./denom
    return function

def calc_column(flux,freq,kap0,nu0,beta,tdust,omega_beam):
    kapnu = kap0*(freq/nu0)**beta
    bnut  = planck(freq,tdust)
    nh2 = flux/(omega_beam.value * bnut * kapnu * 2.8 * c.m_p.cgs)
    return nh2.cgs

def calc_interior_flux(data,r_array,r,omega_beam,pix_scale):
    indices = np.where(r_array < r)
    flux_density = np.sum(data[indices])
    #flux = flux_density / omega_beam.value * u.Jy
    flux = flux_density/(omega_beam.value/(pix_scale.value**2.)) * u.Jy
    return flux

def calc_dust_mass(flux,freq,kap0,nu0,beta,tdust,distance):
    kapnu = kap0*(freq/nu0)**beta
    mass = flux * distance**2./(kapnu * planck(freq,tdust))
    return mass

def calc_density_profile(r_bin_arcsec, mass_profile, distance):
    density_profile = np.zeros(len(r_bin_arcsec))
    for i in range(len(r_bin_arcsec)-1):
        inner_au = r_bin_arcsec[i].value*distance.value * u.au
        outer_au = r_bin_arcsec[i+1].value*distance.value * u.au
        volume = 4./3.*np.pi * ((outer_au)**3. - inner_au**3.)
        shell_density = (mass_profile[i+1]-mass_profile[i])*u.g/volume.cgs/(2.8*c.m_p)
        density_profile[i+1] = shell_density.cgs.value
    return density_profile    
    
def func_powerlaw(x,m,c,c0):
    return c0 + x**m * c 

def func_powerlaw_fix(x,m,c):
    return x**m*c

source = 'SM1N'
image_file = 'continuum/SM1N_continuum.pbcor.fits'
peak_loc = SkyCoord('16h26m27.173s -24d23m21.963s',ICRS)
distance = 137.3 * u.pc

hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
data = np.squeeze(hdu.data)
r_bin = 1.6 # pixels
r_bin_array = [x * r_bin for x in range(30)]

r_array = get_radial_distances(hdu,peak_loc)
# Convert radius array in pixels to arcseconds from header info
pix_scale = np.abs(hdu.header['CDELT1']*u.deg)
r_array_arcsec = r_array * pix_scale.to(u.arcsec)
r_bin_arcsec = r_bin_array * pix_scale.to(u.arcsec)

# Get beam info to convert to N(H2) from arcsec:
bmaj = hdu.header['BMAJ']*u.deg
bmin = hdu.header['BMIN']*u.deg
omega_beam = 2.*np.pi*(bmaj*bmin/(2.3548**2.))
beam_reff = np.sqrt(bmaj*bmin)/2.3548
# For mass calculation:
kap0 = 0.1*u.cm**2./u.g
nu0 = 1.e12*u.Hz
beta = 1.7
tdust = 15. *u.K
freq = 221.*u.GHz

# Get fluxes
flux_profile = np.zeros(len(r_bin_arcsec))
mass_profile = np.zeros(len(r_bin_arcsec))
mass_profile_cgs = np.zeros(len(r_bin_arcsec))
mass_profile_td = np.zeros(len(r_bin_arcsec))
mass_profile_td_cgs = np.zeros(len(r_bin_arcsec))
density_profile = np.zeros(len(r_bin_arcsec))
for i in range(len(r_bin_arcsec)-1):
    r = r_bin_arcsec[i]
    interior_flux = calc_interior_flux(data,r_array_arcsec,r,omega_beam,pix_scale)
    flux_profile[i] = interior_flux.to(u.Jy).value
    interior_mass = calc_dust_mass(interior_flux,freq,kap0,nu0,beta,15.*u.K,distance)
    if i in [0,1]:
        tdust = 30.*u.K
    else:
        tdust = 15.*u.K
    interior_mass_td = calc_dust_mass(interior_flux,freq,kap0,nu0,beta,tdust,distance)
    mass_profile[i] = interior_mass.to(u.Msun).value
    mass_profile_td[i] = interior_mass_td.to(u.Msun).value
    mass_profile_cgs[i] = interior_mass.cgs.value
    mass_profile_td_cgs[i] = interior_mass_td.cgs.value


density_profile = calc_density_profile(r_bin_arcsec,mass_profile_cgs,distance)
density_profile_td = calc_density_profile(r_bin_arcsec,mass_profile_td_cgs,distance)

#print flux_profile

# POWER LAW! Fit.
r_bin_cm = (r_bin_arcsec.value * distance.value*u.au).cgs
popt,pcov = curve_fit(func_powerlaw,r_bin_cm.value[1:len(r_bin_arcsec)-1],
                      density_profile[1:len(r_bin_arcsec)-1],p0=np.asarray([-1.4,1e29,0]))
sigma = np.sqrt([pcov[0,0],pcov[1,1],pcov[2,2]])
#print popt
#print sigma

popt2,pcov2 = curve_fit(func_powerlaw_fix,r_bin_cm.value[1:12],
                        density_profile[1:12],p0=np.asarray([-1.4,1e29]))
sigma2 = np.sqrt([pcov2[0,0],pcov2[1,1]])
#print popt2
#print sigma2

# Singular isothermal sphere
rho_sis = c.k_B * tdust/(2.33*c.m_p)/(2.*np.pi*c.G*r_bin_cm**2)
n_sis = rho_sis/(2.8*c.m_p)

xlim_arcsec = [0.1,4.]
xlim_au = [x*distance.value for x in xlim_arcsec]

fig = plt.figure(figsize=(4,3.))
ax = plt.gca() #plt.subplot(211)
#plt.plot(r_bin_arcsec.value, mass_profile,'o',markersize=3,label='',color='black')
#plt.plot(r_bin_arcsec.value,flux_profile,'o',markersize=3,label='',color='black')
plt.plot(r_bin_arcsec.value, density_profile,'o',markersize=3,label='Isothermal',color='black')
plt.plot(r_bin_arcsec.value,density_profile_td,'o',markersize=3,label='Enhanced $T_d$',
         color='firebrick',zorder=1)
plt.plot(r_bin_arcsec.value,func_powerlaw_fix(r_bin_cm.value,*popt2),'--',zorder=1,label='Power law')
plt.plot(r_bin_arcsec.value,n_sis.cgs.value,alpha=0.5,zorder=1,color='red',label='SIS')

plt.xscale('log')
plt.yscale('log')
ax.set_xlabel('R (arcsec)',size=11)
#ax.set_ylabel('enclosed M (M$_\odot$)')
ax.set_ylabel('$n$ (cm$^{-3}$)',size=11)
ax.yaxis.set_ticks_position('both')
#ax.set_ylabel('enclosed flux (Jy)')
ax.set_xlim(xlim_arcsec)
#ax.set_ylim(0,1e9)
#ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticklabels(['0.01','0.1','1'])
ax.legend(frameon=False,fontsize=9)

# Add second x-axis with AU scales
ax2 = ax.twiny()
ax2.set_xlim(xlim_au)
ax2.set_xscale('log')
ax2.set_xlabel('R (au)',size=11)
import matplotlib.ticker
ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax2.set_xticks([20,100,500])
'''
ax2 = plt.subplot(212)
plt.plot(r_bin_arcsec.value, mass_profile,'o',markersize=3,label='',color='black')
plt.xscale('log')
#plt.yscale('log')
ax2.set_xlabel('R (arcsec)')
ax2.set_ylabel('enclosed M (M$_\odot$)')
ax2.set_xlim(0.1,4.)
ax2.set_xticklabels(['0.01','0.1','1'])
'''

fig.savefig('figures/SM1N_density_profile.pdf',bbox_inches='tight')
plt.close()

