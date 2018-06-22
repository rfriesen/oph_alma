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
from config import plottingDictionary
from matplotlib.ticker import AutoMinorLocator

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

def gaussian_profile(r,sigma,A):
    # Assume r, sigma both in arcsec
    gaussian = A* np.exp(-1.*(r-0)**2./(2.*sigma**2.))
    return gaussian

def convol_gauss(x,*p):
    A, sigma = p
    beam_reff = np.sqrt(hdu.header['BMAJ']*hdu.header['BMIN'])*u.deg/2.3548
    beam_gaus = np.exp(-1.*(x-0)**2./(2.*beam_reff.to(u.arcsec).value**2.))
    #print x.shape, beam_gaus.shape
    return np.convolve(A*np.exp((-(x)**2.)/(2.*sigma**2.)),beam_gaus,mode='same')

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

def radial_profile(hdu, peak_loc, radial_bin, max_r):
    wcs = WCS(hdu.header)
    # Get pixel location for peak emission
    x_peak, y_peak = peak_loc.to_pixel(wcs)
    # Distance from each point to peak emission in pixels
    # ALMA hdus have extra axes
    data = np.squeeze(hdu.data)
    x_arr, y_arr = np.indices(data.shape)
    r = np.sqrt((x_arr - x_peak)**2. + (y_arr - y_peak)**2.)

    radius, r_profile, r_sdev = binArray(data, r, radial_bin)
    return radius[0:max_r], r_profile[0:max_r], r_sdev[0:max_r]

source = 'SM1N'
image_file = 'continuum/SM1N_continuum.pbcor.fits'
image_file_nopb = 'continuum/SM1N_continuum.image.fits'
rad_profile = 'SM1N_radial_profile.txt'
peak_loc = SkyCoord('16h26m27.173s -24d23m21.963s',ICRS)
plot_param = plottingDictionary[source]

# Get beam info:
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
bmaj = hdu.header['BMAJ']*u.deg
bmin = hdu.header['BMIN']*u.deg
omega_beam = 2.*np.pi*(bmaj*bmin/(2.3548**2.))
beam_reff = np.sqrt(bmaj*bmin)/2.3548

# And plot
# Flux values for figure
# SM1 - no taper
vmin = plot_param['vmin']
vmax = plot_param['vmax']
sigma = plot_param['sigma']
clevs = [5.*x*sigma + 5.*sigma for x in range(15)]
neglevs = [-5.*x*sigma - 5.*sigma for x in range(3)]
neglevs = neglevs[::-1]     
l1 = 1
labsize=9
# Pointing centre 16:26:27.200000 -24.23.22.4
rc = plot_param['rc']
dc = plot_param['dc']
# Field of view ~ 27.4 arcsec
csize = 0.003804
# Image width, height (degrees)
iw = plot_param['im_width']
ih = plot_param['im_width']
distance = plot_param['distance']

# Figure
fig = plt.figure(figsize=(4,7.5))
f0 = aplpy.FITSFigure(image_file_nopb,figure=fig,subplot=[0.18,0.62,0.73,0.32])
f0.recenter(rc,dc,width=iw,height=ih)
f0.show_colorscale(cmap='gist_heat',vmin=vmin,vmax=vmax)
f0.show_beam(edgecolor='white',facecolor='none',corner='top left')
f0.ticks.set_color('white')
f0.ticks.set_length(5)
#f0.axis_labels.hide()
#f0.axis_labels.set_font(size=labsize, weight='normal', \
#                        stretch='normal', family='sans-serif', \
#                        style='normal', variant='normal')
f0.tick_labels.set_style('colons')
f0.tick_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f0.ticks.set_xspacing(1.*15./3600)
f0.tick_labels.set_xformat('hh:mm:ss.s')
f0.tick_labels.set_yformat('dd:mm:ss')
ang_sep = (plot_param['scalebar_size'].to(u.au)/distance).to(u.arcsec,equivalencies=u.dimensionless_angles())
f0.add_scalebar(ang_sep,label='{0:.0f} au'.format(plot_param['scalebar_size'].value),
                color='white',corner=plot_param['scalebar_pos'])
f0.add_colorbar(location='top',axis_label_text='Jy/beam',width=0.1,ticks=[0,0.0005,0.001,0.0015])
f0.colorbar.set_font(family='sans_serif',size=9)
f0.show_contour(image_file_nopb, levels=clevs, colors='white', linewidths=l1*0.8)
f0.show_contour(image_file_nopb, levels=neglevs, colors='white', linewidths=l1*0.8,linestyles='dashed')
f0.show_markers(peak_loc.ra.value,peak_loc.dec.value,
                marker='+',edgecolor='blue',facecolor='blue',s=12)
f0.add_label(0.5,0.88,plot_param['image_label'],relative=True,fontsize=9,color='white',weight='bold',horizontalalignment='center')

# Radial profile
# Read in values
r_bin = 2. # pixels
max_r = 120

radius, r_profile, r_sdev = radial_profile(hdu, peak_loc, r_bin, max_r)
# Convert radius array in pixels to arcseconds from header info
pix_scale = np.abs(hdu.header['CDELT1']*u.deg)
radius_arcsec = radius * pix_scale.to(u.arcsec)

# Was reading in non-pbcor image fit 
#radius_arcsec,r_profile,r_sdev = np.loadtxt(rad_profile,unpack=True,usecols=(0,1,2))
# Convolve with beam: need to mirror this to make something reasonable. 
negx = radius_arcsec[::-1]*-1.
mirror_x = np.concatenate((negx,radius_arcsec))
neg_prof = r_profile[::-1]*1.e3
mirror_prof = np.concatenate((neg_prof,r_profile*1.e3))
neg_r_sdev = r_sdev[::-1]*1.e3
mirror_r_sdev = np.concatenate((neg_r_sdev,r_sdev*1.e3))
# Fit
p0 = [1.,1.]
coeff, var_matrix = curve_fit(convol_gauss, mirror_x, mirror_prof, p0=p0, 
                              sigma=mirror_r_sdev)
print coeff
# Need uncertainties
fit_sigma=coeff[1]
fit_amp = coeff[0]
fit_damp = np.sqrt(var_matrix[0,0])
fit_dsigma = np.sqrt(var_matrix[1,1])

# Font size
axis_font = {'fontname':'Arial','size':'9'}

# Minor ticks
minorLocator = AutoMinorLocator(4)

# Plot radial profile
ax2 = fig.add_axes([0.26,0.06,0.57,0.22])
#plt.errorbar(radius_arcsec.value, r_profile*1.e3, yerr=r_sdev*1.e3, 
#             fmt='o',markersize=3,label='',ecolor='gold',color='black')
ax2.plot(radius_arcsec.value,r_profile*1.e3,color='black',zorder=4,marker='o',markersize=2,label='Data')
ax2.fill_between(radius_arcsec.value,(r_profile-r_sdev)*1.e3,(r_profile+r_sdev)*1.e3,
                 facecolor='gray',alpha=0.5,interpolate=True)
plt.xscale('log')
ax2.set_xlabel('R (arcsec)',**axis_font)
ax2.set_ylabel('mJy beam$^{-1}$',**axis_font)
ax2.set_xlim(0.03,10.)
ax2.set_ylim(-0.3,2.1)
ax2.set_xticklabels(['0.001','0.01','0.1','1','10'],**axis_font)
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
plt.yticks(size=9)
# Overplot beam
# Want more finely sampled than radial profile
r_max = 25.
n_bins = 1000
r_beam = np.array([x*25./n_bins for x in range(n_bins)]) # assume arcsec
beam = gaussian_profile(r_beam,beam_reff.to(u.arcsec).value,1.)
ax2.plot(r_beam,beam*np.max(r_profile*1.e3),'-r',zorder=2,label='Beam')
# Overplot gaussian fit
fit = convol_gauss(mirror_x,fit_amp,fit_sigma)
ax2.plot(mirror_x,fit,'-b',zorder=2,label='Fit')
ax2.plot([0,10],[0,0],'--',zorder=1,color='gray')
ax2.legend(frameon=False,fontsize=8)
ax2.text(0.02,0.2,r'$\sigma$ = {:4.2f}" $\pm$ {:4.2f}"'.format(fit_sigma,fit_dsigma),
         transform=ax2.transAxes,fontsize=labsize-1)


# Plot radial uv data
uv_fit_file = 'models/{0}_model_uv_amps.txt'.format(source)
uv_dist_m, uv_amp, uv_amp_err = np.loadtxt(uv_fit_file,unpack=True)
ax3 = fig.add_axes([0.26,0.34,0.57,0.22])
plt.errorbar(uv_dist_m, uv_amp, yerr = uv_amp_err, fmt='o', 
             markersize=3,label='data',color='black',ecolor='black')
#plt.legend(loc = 1, numpoints = 1, fontsize=8, frameon=False)
plt.plot([0,1.e3],[0,0],linestyle=':',color='gray',zorder=0)
plt.ylim(-8,30)
plt.xlim(0,550.)
plt.xlabel(r'uv Distance (m)',**axis_font)
plt.ylabel('Real Visibility (mJy)',**axis_font)
plt.xticks(size=9)
plt.yticks(size=9)
ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
ax3.yaxis.set_ticks_position('both')
ax3.xaxis.set_ticks_position('both')
fig.savefig('figures/SM1N_radial_profile_v2.pdf')
plt.close()
