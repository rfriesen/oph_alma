##################################################
# Script to fit Gaussians to ALMA continuum data
# Ophiuchus targets 
# ALMA Cycle 2
# Redoing script that was lost in computer failure
# (sigh)
# Also calculating masses for components
##################################################
from astropy import units as u
from astropy import constants as c
import numpy as np

def planck(freq,tdust):
    denom = np.exp(c.h.cgs * freq/(c.k_B.cgs*tdust))-1.
    function = 2. *c.h.cgs * freq**3./c.c.cgs**2./denom
    return function

def calc_dust_mass(flux,freq,distance,kap0,nu0,beta,tdust):
    # Calculate mass from dust emission assuming all inputs have correct units
    kapnu = kap0*(freq/nu0)**beta
    mass = flux * distance**2./(kapnu * planck(freq,tdust))
    return mass

def calc_dust_density(mass,distance,deconmaj,deconmin):
    # Calculate density from mass, fit results
    amaj_pc = deconmaj * distance /206265. # arcsec to pc
    amin_pc = deconmin * distance /206265.
    volume = 4./3.*np.pi * amaj_pc**2.*amin_pc
    density = mass/volume/(2.8 * c.m_p.cgs)
    return density

def calc_ra_hms(ra_deg):
    ra = ra_deg + 360. # The output uses a negative convention ??!?
    if not np.isscalar(ra):
        ra_h = np.zeros(len(ra))
        ra_m = np.zeros(len(ra))
        ra_s = np.zeros(len(ra))
        for l in range(len(ra)):
            ra_h[l] = np.int(ra[l]/15.)
            ra_m[l] = np.int((ra[l]/15. - ra_h[l])*60.)
            ra_s[l] = ((ra[l]/15. - ra_h[l])*60. - ra_m[l])*60.
    else:
        ra_h = np.int(ra/15.)
        ra_m = np.int((ra/15. - ra_h)*60.)
        ra_s = ((ra/15. - ra_h)*60. - ra_m)*60.
    return ra_h, ra_m, ra_s

def calc_de_dms(de_deg):
    de = de_deg
    if not np.isscalar(de):
        de_d = np.zeros(len(de))
        de_m = np.zeros(len(de))
        de_s = np.zeros(len(de))
        for l in range(len(de)):
            if de[l] < 0:
                de_neg = -1
                de[l] = np.abs(de[l])
            else:
                de_neg = 1
            de_d[l] = np.int(de[l])
            de_m[l] = np.int((de[l]-de_d[l])*60.)
            de_s[l] = ((de[l]-de_d[l])*60.-de_m[l])*60.
            de_d[l] = de_d[l] * de_neg
    else:
        if de < 0:
            de_neg = -1
            de = np.abs(de)
        else:
            de_neg = 1        
        de_d = np.int(de)
        de_m = np.int((de-de_d)*60.)
        de_s = ((de-de_d)*60.-de_m)*60.
        de_d = de_d * de_neg  
    return de_d, de_m, de_s

# Set up 
# Question here - use high resolution for AN6 fit? 
sources = ['16263_2422','AN6','B2-A7','SM1','SM1N']
sources_latex = ['16263{\_}2422','AN6','B2-A7','SM1','SM1N']
file_ext = ['_selfcal_2','_taper','_taper','_selfcal_3','']
rms_list = [0.81e-4,0.6e-4,0.6e-4,1.1e-4,0.5e-4]
# For masses
kap0 = 0.1*u.cm**2./u.g
nu0 = 1.e12*u.Hz
beta = 1.7
tdust = 15. *u.K
# For Oph
distance = 137.3 * u.pc

for i in range(len(sources)):
    filename = sources[i] + '_continuum'+ file_ext[i] + '.image.fits'
    region_file = sources[i] + '_imfit_region.crtf'
    output_summary = sources[i] + '_imfit.summary'
    estimates = sources[i] + '_imfit_estimates.txt'
    output_model = sources[i] + '_imfit.model'
    imfit(imagename=filename,
          region=region_file,
          rms=rms_list[i],
          summary=output_summary,
          estimates=estimates)
          #model=model)

# Works yay!
# Now to redo output table for paper argh
# Need to calculate masses but can't use astropy in casa?
# Need to follow import instructions in CASA
import numpy as np
f = open('imfit_output_table.tex','w')
# Latex table preamble
f.write('\\floattable \n')
f.write('\\begin{deluxetable}{lcccccccccccc} \n')
f.write('\\tabletypesize{\\footnotesize} \n')
#f.write('\rotate \n')
f.write('\\tablecolumns{13} \n')
f.write('\\tablewidth{0pt} \n')
f.write('\\tablecaption{Gaussian fit results \label{tab:imfit}} \n')
f.write('\\tablehead{ \n')
f.write('\colhead{Source} & \colhead{Comp} & \colhead{R.A.} & \colhead{Decl.} & \colhead{$S_\\nu$} & \colhead{Peak} & \colhead{$\sigma_\mathrm{maj}$} & \colhead{$\sigma_\mathrm{min}$} & \colhead{P.A.} & \colhead{$\sigma_\mathrm{maj}$} & \colhead{$\sigma_\mathrm{min}$} & \colhead{P.A.} & \colhead{$M$} \\\ \n')
f.write('\colhead{} & \colhead{} & \colhead{J2000} & \colhead{J2000} & \colhead{Jy} & \colhead{Jy\ beam$^{-1}$} & \colhead{\\arcsec} & \colhead{\\arcsec} & \colhead{\degr} & \colhead{\\arcsec} & \colhead{\\arcsec} & \colhead{\degr} & \colhead{M$_\odot$} \n')
f.write('} \n')
f.write('\startdata \n')
# Now include data with nice formatting
for k in range(len(sources)):
    i,ierr,peak,peakerr,ra,dec,raerr,decerr,conmaj,conmin,conpa,conmajerr,conminerr,conpaerr,deconmaj,deconmin,deconpa,deconmajerr,deconminerr,deconpaerr,freq = np.loadtxt(sources[k]+'_imfit.summary',usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21),unpack=True,skiprows=2)
    # Calculate RA, Dec. more nicely
    ra_h, ra_m, ra_s = calc_ra_hms(ra)
    de_d, de_m, de_s = calc_de_dms(dec)
    # Calculate masses
    # Add units to frequency, flux
    i = i * u.Jy
    freq = freq * u.Hz
    mass = calc_dust_mass(i,freq,distance,kap0,nu0,beta,tdust)
    # Calculate density
    density = calc_dust_density(mass,distance,deconmaj,deconmin)
    print density.to(u.cm**(-3))
    # Write to a table
    # For targets with 2 sources need to loop here
    if np.isscalar(i.value):
        f.write(sources_latex[k]+' & \\nodata & '+'{:2d}'.format(ra_h)+':'+\
                '{:2d}'.format(ra_m)+\
                ':'+'{:04.1f}'.format(ra_s)+' & '+'{:2d}'.format(de_d)+':'+\
                '{:2d}'.format(de_m)+':'+'{:04.1f}'.format(de_s)+' & '+\
                '{:6.4f}'.format(i.value)+'('+\
                '{:2d}'.format(np.int(ierr*1.e4))+\
                ') & '+'{:7.5f}'.format(peak)+'('+\
                '{:2d}'.format(np.int(peakerr*1.e5))+\
                ') & '+'{:4.2f}'.format(conmaj)+'('+\
                '{:2d}'.format(np.int(conmajerr*100))+') & '+\
                '{:4.2f}'.format(conmin)+'('+\
                '{:2d}'.format(np.int(conminerr*100.))+') & '+\
                '{:3.1f}'.format(conpa)+'('+\
                '{:1d}'.format(np.int(conpaerr*10))+\
                ') & '+'{:4.2f}'.format(deconmaj)+'('+\
                '{:2d}'.format(np.int(deconmajerr*100))+') & '+\
                '{:4.2f}'.format(deconmin)+'('+\
                '{:2d}'.format(np.int(deconminerr*100.))+') & '+\
                '{:3.1f}'.format(deconpa)+'('+\
                '{:3.1f}'.format(deconpaerr)+\
                ') & '+'{:6.4f}'.format(mass.to(u.Msun).value)+'\\\ \n')
    else:
        for j in range(len(i.value)):
            f.write(sources_latex[k]+' & '+'{:1d}'.format(j)+' & '+\
                    '{:2d}'.format(np.int(ra_h[j]))+':'+\
                    '{:2d}'.format(np.int(ra_m[j]))+\
                    ':'+'{:04.1f}'.format(ra_s[j])+' & '+\
                    '{:2d}'.format(np.int(de_d[j]))+':'+\
                    '{:2d}'.format(np.int(de_m[j]))+':'+\
                    '{:04.1f}'.format(de_s[j])+\
                    ' & '+'{:6.4f}'.format(i[j].value)+'('+\
                    '{:2d}'.format(np.int(ierr[j]*1.e4))+\
                    ') & '+'{:7.5f}'.format(peak[j])+'('+\
                    '{:2d}'.format(np.int(peakerr[j]*1.e5))+\
                    ') & '+'{:4.2f}'.format(conmaj[j])+'('+\
                    '{:2d}'.format(np.int(conmajerr[j]*100))+') & '+\
                    '{:4.2f}'.format(conmin[j])+'('+\
                    '{:2d}'.format(np.int(conminerr[j]*100.))+') & '+\
                    '{:3.1f}'.format(conpa[j])+'('+\
                    '{:1d}'.format(np.int(conpaerr[j]*10))+\
                    ') & '+'{:4.2f}'.format(deconmaj[j])+'('+\
                    '{:2d}'.format(np.int(deconmajerr[j]*100))+') & '+\
                    '{:4.2f}'.format(deconmin[j])+'('+\
                    '{:2d}'.format(np.int(deconminerr[j]*100.))+') & '+\
                    '{:3.1f}'.format(deconpa[j])+'('+\
                    '{:3.1f}'.format(deconpaerr[j])+\
                    ') & '+'{:6.4f}'.format(mass[j].to(u.Msun).value)+'\\\ \n')

f.write('\enddata \n')
f.write('\end{deluxetable}')
f.close()

