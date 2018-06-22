# Script to run gaussian fitting in CASA on continuum images
# Should be similar to uvmodelfit results, but include pb correction
# 
# For SM1N, AN6, B2-A7, use calibrated continuum ms
# For SM1, 16263_2422, use self-cal continuum ms
# 
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 

def planck(freq,tdust):
    denom = np.exp(c.h.cgs * freq/(c.k_B.cgs*tdust))-1.
    function = 2. *c.h.cgs * freq**3./c.c.cgs**2./denom
    return function

def calc_dust_mass(flux,freq,distance,kap0,nu0,beta,tdust):
    # Calculate mass from dust emission assuming all inputs have correct units
    kapnu = kap0*(freq/nu0)**beta
    mass = flux * distance**2./(kapnu * planck(freq,tdust))
    return mass

def calc_dust_density(mass,distance,amaj,amin):
    # Calculate density from mass, fit results
    amaj_pc = amaj.to(u.arcsec).value * distance /206265. # arcsec to pc
    amin_pc = amin.to(u.arcsec).value * distance /206265.
    volume = 4./3.*np.pi * amin_pc**2.*amaj_pc
    density = mass/volume/(2.33 * c.m_p.cgs)
    return density

def read_model_info(fit):
    # Get fit parameters and uncertainties. 
    flux = fit['flux']['value'][0]
    flux_err = fit['flux']['error'][0]
    flux_unit = fit['flux']['unit']
    major_axis = fit['shape']['majoraxis']['value']
    major_axis_err = fit['shape']['majoraxiserror']['value']
    major_axis_unit = fit['shape']['majoraxis']['unit']
    minor_axis = fit['shape']['minoraxis']['value']
    minor_axis_err = fit['shape']['minoraxiserror']['value']
    minor_axis_unit = fit['shape']['minoraxis']['unit']
    ra = fit['shape']['direction']['m0']['value']
    ra_unit = fit['shape']['direction']['m0']['unit']
    de = fit['shape']['direction']['m1']['value']
    de_unit = fit['shape']['direction']['m1']['unit']
    pa = fit['shape']['positionangle']['value']
    pa_err = fit['shape']['positionangleerror']['value']
    pa_unit = fit['shape']['positionangle']['unit']
    # Convert to useful units
    # Print warnings if units are not what assumed here
    if flux_unit != 'Jy':
        print 'Warning: non-standard flux unit'
    flux = flux * u.Jy
    flux_err = flux_err * u.Jy
    # Major axis
    if major_axis_unit == 'arcmin':
        major_axis = major_axis * u.arcmin
        major_axis_err = major_axis_err * u.arcmin
    elif major_axis_unit == 'arcsec':
        major_axis = major_axis * u.arcsec
        major_axis_err = major_axis_err * u.arcsec
    else:
        print 'Warning: non-standard major axis unit'
    # Minor axis
    if minor_axis_unit == 'arcmin':
        minor_axis = minor_axis * u.arcmin
        minor_axis_err = minor_axis_err * u.arcmin
    elif minor_axis_unit == 'arcsec':
        minor_axis = minor_axis * u.arcsec
        minor_axis_err = minor_axis_err * u.arcsec
    else:
        print 'Warning: non-standard minor axis unit'
    # Just assume PA is in degrees
    pa = pa *u.deg
    pa_err = pa_err * u.deg
    # RA and dec
    if ra_unit == 'rad':
        ra = ra * u.rad
        ra_deg = 360.*u.deg + (ra.to(u.deg))
    else:
        print 'Warning: non-standard ra unit'
    if de_unit == 'rad':
        de = de * u.rad
        de_deg = de.to(u.deg)
    # Convert to SkyCoords
    comp_coords = SkyCoord(ra_deg,de_deg,unit='deg')
    # Return a bunch of stuff
    return comp_coords, flux, flux_err, major_axis, major_axis_err, minor_axis, minor_axis_err, pa, pa_err

# Create dictionary for results
gauss_fit_results={}

# First: SM1
# Two components
cont_image = 'continuum/SM1_continuum_selfcal_3.pbcor'
est_file = 'models/SM1_imfit_estimates.txt'
res_image = 'models/SM1_continuum_selfcal_3.pbcor.imfit_res'
model_image = 'models/SM1_continuum_selfcal_3.pbcor.imfit_model'
logfile = 'models/SM1_imfit_log.txt'
complist = 'models/SM1_imfit_complist.cl'
dooff = False # Fit zero level offset?

fit_two_comp = imfit(imagename=cont_image,
                     estimates=est_file,
                     logfile=logfile,
                     residual=res_image,
                     model=model_image,
                     complist=complist,
                     dooff=dooff,
                     overwrite=True)

# Want 'deconvolved' results
results1 = read_model_info(fit_two_comp['deconvolved']['component0'])
results2 = read_model_info(fit_two_comp['deconvolved']['component1'])

SM1_results = {'comps':2,
               'centre':[results1[0],results2[0]],
               'flux':[results1[1],results2[1]],
               'flux_err':[results1[2],results2[2]],
               'major':[results1[3],results2[3]],
               'major_err':[results1[4],results2[4]],
               'minor':[results1[5],results2[5]],
               'minor_err':[results1[6],results2[6]],
               'pa':[results1[7],results2[7]],
               'pa_err':[results1[8],results2[8]]}

gauss_fit_results['SM1']=SM1_results

# Next: 16263-2422
# Two components
cont_image = 'continuum/16263_2422_continuum_selfcal_2.pbcor'
est_file = 'models/16263_2422_imfit_estimates.txt'
res_image = 'models/16263_2422_continuum_selfcal_2.pbcor.imfit_res'
model_image = 'models/16263_2422_continuum_selfcal_2.pbcor.imfit_model'
logfile = 'models/16263_2422_imfit_log.txt'
complist = 'models/16263_2422_imfit_complist.cl'
dooff = False # Fit zero level offset?

fit_two_comp = imfit(imagename=cont_image,
                     estimates=est_file,
                     logfile=logfile,
                     residual=res_image,
                     model=model_image,
                     complist=complist,
                     dooff=dooff,
                     overwrite=True)

results1 = read_model_info(fit_two_comp['deconvolved']['component0'])
results2 = read_model_info(fit_two_comp['deconvolved']['component1'])

r16263_2422_results = {'comps':2,
               'centre':[results1[0],results2[0]],
               'flux':[results1[1],results2[1]],
               'flux_err':[results1[2],results2[2]],
               'major':[results1[3],results2[3]],
               'major_err':[results1[4],results2[4]],
               'minor':[results1[5],results2[5]],
               'minor_err':[results1[6],results2[6]],
               'pa':[results1[7],results2[7]],
               'pa_err':[results1[8],results2[8]]}

gauss_fit_results['16263_2422']=r16263_2422_results

# Next: AN6
# One component
cont_image = 'continuum/AN6_continuum.pbcor'
est_file = 'models/AN6_imfit_estimates.txt'
res_image = 'models/AN6_continuum.pbcor.imfit_res'
model_image = 'models/AN6_continuum.pbcor.imfit_model'
logfile = 'models/AN6_imfit_log.txt'
complist = 'models/AN6_imfit_complist.cl'
dooff = False # Fit zero level offset?

fit_one_comp = imfit(imagename=cont_image,
                     estimates=est_file,
                     logfile=logfile,
                     box='240,270,310,325',
                     residual=res_image,
                     model=model_image,
                     complist=complist,
                     dooff=dooff,
                     overwrite=True)

results1 = read_model_info(fit_one_comp['deconvolved']['component0'])

AN6_results = {'comps':1,
               'centre':[results1[0]],
               'flux':[results1[1]],
               'flux_err':[results1[2]],
               'major':[results1[3]],
               'major_err':[results1[4]],
               'minor':[results1[5]],
               'minor_err':[results1[6]],
               'pa':[results1[7]],
               'pa_err':[results1[8]]}

gauss_fit_results['AN6']=AN6_results

# Write out table

sources = ['16263_2422','SM1','AN6']
sources_latex = [['GSS 30 IRS3','GSS 30 IRS1'],['SM1','N3-mm'],['N6-mm']]
rms_list = [0.81e-4,0.6e-4,0.6e-4,1.1e-4,0.5e-4]
# For masses
freq = 221.*u.GHz
kap0 = 0.1*u.cm**2./u.g
nu0 = 1.e12*u.Hz
beta_all = 1.7
beta_SM1 = 0.2
beta_16263 = 0.4
beta_SED = [[0.2,0.0],[0.2,0.5],[1.7]]
tdust = 30. *u.K
# For Oph
distance = 137.3 * u.pc


# Set up output file
f = open('models/imfit_table_v2.tex','w')
# Latex table preamble
f.write('\\floattable \n')
f.write('\\begin{deluxetable}{lccccccccccc} \n')
f.write('\\tabletypesize{\\footnotesize} \n')
#f.write('\rotate \n')
f.write('\\tablecolumns{12} \n')
f.write('\\tablewidth{0pt} \n')
f.write('\\tablecaption{Gaussian fit results for compact sources \label{tab:uvmodelfit}} \n')
f.write('\\tablehead{ \n')
f.write('\colhead{Source} & \colhead{R.A.} & \colhead{Decl.} & \colhead{$S_\\nu$} & \colhead{$\sigma_\mathrm{maj}$} & \colhead{$\sigma_\mathrm{min}$} & \colhead{P.A.} & \colhead{$M$} & \colhead{$n$} & \colhead{$\beta$} & \colhead{$M$\\tablenotemark{2}} & \colhead{$n$\\tablenotemark{a}} \\\ \n')
f.write('\colhead{} & \colhead{J2000} & \colhead{J2000} & \colhead{Jy} & \colhead{\\arcsec} & \colhead{\\arcsec} & \colhead{\degr} & \colhead{M$_\odot$} & \colhead{cm$^{-3}$} & \colhead{} & \colhead{M$_\odot$} & \colhead{cm$^{-3}$} \n')
f.write('} \n')
f.write('\startdata \n')

for i in range(len(sources)):
    source = sources[i]
    fit_dict = gauss_fit_results[source]
    ncomp = fit_dict['comps']
    for comp in range(ncomp):
        beta = beta_SED[i][comp]
        source_name = sources_latex[i][comp]
        #if comp == 0:
        #    if source == 'SM1':
        #        beta = beta_SM1
        #    elif source == '16263_2422':
        #        beta = beta_16263
        #else:
        #    beta = beta_all
        mass = calc_dust_mass(fit_dict['flux'][comp],freq,distance,kap0,nu0,beta,tdust)
        density = calc_dust_density(mass,distance,fit_dict['major'][comp],fit_dict['minor'][comp])
        #masserr = mass * (ierr/i)
        mass2 = calc_dust_mass(fit_dict['flux'][comp],freq,distance,kap0,nu0,beta_all,tdust)
        density2 = calc_dust_density(mass2,distance,fit_dict['major'][comp],fit_dict['minor'][comp])
        # Write results to latex file:
        f.write(source_name+' & '+\
                fit_dict['centre'][comp].to_string(style='hmsdms',sep=':')+\
                ' & {:6.4f}'.format(fit_dict['flux'][comp].value)+\
                '('+'{:2d}'.format(np.int(fit_dict['flux_err'][comp].value*1.e4))+') & '+\
                #') & '+'{:7.5f}'.format(peak)+'('+\
                #'{:2d}'.format(np.int(peakerr*1.e5))+') & '+\
                ' {:5.3f}'.format(fit_dict['major'][comp].to(u.arcsec).value)+\
                '('+'{:3d}'.format(np.int(fit_dict['major_err'][comp].value*1000))+') & '+\
                ' {:5.3f}'.format(fit_dict['minor'][comp].to(u.arcsec).value)+\
                '('+'{:3d}'.format(np.int(fit_dict['minor_err'][comp].value*1000.))+') & '+\
                ' {:3.1f}'.format(fit_dict['pa'][comp].value)+\
                '('+'{:3.1f}'.format(fit_dict['pa_err'][comp].value)+') & '+\
                ' {:6.4f}'.format(mass2.to(u.Msun).value)+\
                #'('+'{:1d}'.format(np.int(np.round(masserr.to(u.Msun).value*1e4)))+') & '+\
                ' & {:.1e}'.format(density2.to(u.cm**(-3)).value)+\
                ' & {:.1}'.format(beta)+\
                ' & {:6.4f}'.format(mass.to(u.Msun).value)+\
                #'('+'{:1d}'.format(np.int(np.round(masserr.to(u.Msun).value*1e4)))+') & '+\
                ' & {:.1e}'.format(density.to(u.cm**(-3)).value)+' \\\ \n')


f.write('\enddata \n')
f.write('\\tablenotetext{a}{Masses determined with $\\beta = 1.7$.} \n')
f.write('\\tablenotetext{b}{Masses determined with $\\beta$ derived in \\ref{sec:masses}.} \n')
f.write('\end{deluxetable}')
f.close()


