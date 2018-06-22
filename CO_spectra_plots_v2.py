# ##########
# Create spectrum figure for ALMA Cycle 2 paper
# ##########
import aplpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.units as u
from astropy import constants as c
from astropy.io import fits
from config import plottingDictionary
from matplotlib.ticker import AutoMinorLocator
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
from astropy import wcs

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

def align_yaxis(ax1,v1,ax2,v2):
    """adjust ax2 limit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0,v1))
    _, y2 = ax2.transData.transform((0,v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0,0)) - inv.transform((0,y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy,maxy+dy)

contDir = 'continuum'
coDir   = 'co_21'

jcmt12co = 'jcmt/CO_tmb.fits'
jcmt13co = 'jcmt/13co_tmb.fits'
jcmtc18o = 'jcmt/c18o_tmb.fits'

sources = ['SM1','16263_2422','16263_2422_1']
source_names  = ['SM1','GSS 30 IRS3','GSS 30 IRS1']
source_coords = [SkyCoord('16h26m27.856s -24d23m59.57s',ICRS),SkyCoord('16h26m21.72s -24d22m50.939s',ICRS),
                 SkyCoord('16h26m21.359s','-24d23m04.865s',ICRS)]
coImage_SM1   = '{0}/SM1/SM1_CO21_0.16_image.pbcor.fits'.format(coDir)
contImage_SM1 = '{0}/SM1_continuum_selfcal_3.pbcor.fits'.format(contDir)
coImage_16   = '{0}/16263_2422/16263_2422_CO21_image.pbcor.fits'.format(coDir)
contImage_16 = '{0}/16263_2422_continuum_selfcal_2.pbcor.fits'.format(contDir)
coImages   = [coImage_SM1,coImage_16,coImage_16]
contImages = [contImage_SM1,contImage_16,contImage_16]
co_spec    = ['{0}/SM1/co_21_avg_over_src.txt'.format(coDir),
              '{0}/16263_2422/co_21_avg_over_src.txt'.format(coDir),
              '{0}/16263_2422/co_21_avg_over_src.txt'.format(coDir)]
co_chans = [125,71,71]
co_vels  = [5.0,3.4,3.4]
nh3_vels = [3.65,3.23,3.23]
co_vels_label = ['{0} km s$^{{-1}}$'.format(x) for x in co_vels]
co_vmin  = -0.1
co_vmax  = [0.2,0.2,0.2]
cont_levs = [5,50,500]
vlsr_lims = [-25.,28.]
vlsr_lim_SM1 = [-10,15]

rc_zoom = 246.59042
dc_zoom = -24.38067

imw = 0.002
imh = 0.002

labsize = 9
l1 = 0.8

fig = plt.figure(figsize=(7.5,7.5))

for i in range(len(sources)):
    plotParam = plottingDictionary[sources[i]]
    # Need extra dim in slices because ALMA data have four axes
    if i == 2:
        f0 = aplpy.FITSFigure(coImages[i],figure=fig,slices=[co_chans[i],0],
                              subplot=[0.13,0.13+0.3*i,0.36,0.35])
    else:
        f0 = aplpy.FITSFigure(coImages[i],figure=fig,slices=[co_chans[i],0],
                              subplot=[0.202,0.13+0.3*i,0.22,0.33])
    f0.show_grayscale(vmin=co_vmin,vmax=co_vmax[i],stretch='linear')
    if i == 0:
        f0.recenter(plotParam['rc'],plotParam['dc'],width=imw,height=imh)
    else:
        f0.recenter(rc_zoom,dc_zoom,width=imw,height=imh)
        f0.add_colorbar(location='top',axis_label_text=r'Jy beam$^{-1}$',width=0.1,
                        ticks=[-0.1,0,0.1,0.2],pad=0.07)
        f0.colorbar.set_font(family='sans_serif',size=labsize-1)
        f0.colorbar.set_axis_label_font(family='sans_serif',size=labsize-1)
    f0.ticks.set_color('black')
    f0.ticks.set_xspacing(0.5*15./3600)
    f0.ticks.set_length(3)
    f0.axis_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
    f0.tick_labels.set_style('colons')
    f0.tick_labels.set_font(size=labsize, weight='normal', \
                            stretch='normal', family='sans-serif', \
                            style='normal', variant='normal')
    f0.tick_labels.set_xformat('hh:mm:ss')
    f0.tick_labels.set_yformat('dd:mm:ss')
    if i == 1:
        f0.axis_labels.hide_x()
    f0.show_beam(edgecolor='white',facecolor='none',corner='top left')
    f0.show_contour(contImages[i],levels=[x*plotParam['sigma'] for x in cont_levs],
                    colors='white',linewidths=l1)
    f0.add_label(0.1,0.7,co_vels_label[i],relative=True,color='white',size=(labsize-1),
                 horizontalalignment='left')
    f0.show_scalebar(0.00020,label='100 au',color='white',corner='bottom left',
                     fontsize=labsize)
    ax2 = fig.add_axes([0.55,0.13+0.3*i,0.35,0.27])
    vel, flux = np.loadtxt(co_spec[i],unpack=True,skiprows=8)
    plt.xlim(vlsr_lims[0],vlsr_lims[1])
    plt.yticks([-0.2,0,0.2,0.4,0.6])
    xminor_locator = AutoMinorLocator(5)
    yminor_locator = AutoMinorLocator(2)
    ax2.xaxis.set_minor_locator(xminor_locator)
    ax2.yaxis.set_minor_locator(yminor_locator)
    if i == 1:
        #ax2.set_xticklabels([])
        plt.ylim(-0.2,0.4)
    else:
        plt.xlabel(r'$v_\mathrm{LSR}$ (km s$^{-1})$',fontsize=labsize)
        plt.ylim(-0.1,0.7)
    ax2.tick_params(axis='both',which='major',labelsize=labsize)
    plt.plot(vel, flux, color='black',zorder=2)
    plt.plot([-50,50],[0,0],color='gray',linestyle=':',zorder=1)
    plt.plot([co_vels[i],co_vels[i]],[-1,1],color='blue',linestyle='-',zorder=1,alpha=0.5)
    plt.plot([nh3_vels[i],nh3_vels[i]],[-1,1],color='red',linestyle='-',zorder=1,alpha=0.5)
    plt.ylabel(r'Jy beam$^{-1}$',fontsize=labsize)
    plt.text(0.05,0.85,source_names[i],fontsize=labsize,horizontalalignment='left',
             transform=ax2.transAxes,)

fig.savefig('figures/CO_spectra.pdf')
plt.close('all')

# Look at CO spectra at sources
cube12co, header12co = fits.getdata(jcmt12co,header=True)
wcs12co = wcs.WCS(jcmt12co)
# Create velocity axis (already in km/s)
nchan = header12co['NAXIS3']
v_inc = header12co['CDELT3']
v_ctr = header12co['CRVAL3']
v_pix = header12co['CRPIX3']
velo_12co = (np.arange(nchan)*v_inc + v_inc + (v_ctr - v_inc*v_pix))

cube13co, header13co = fits.getdata(jcmt13co,header=True)
wcs13co = wcs.WCS(jcmt13co)
# Create velocity axis (already in km/s)
nchan = header13co['NAXIS3']
v_inc = header13co['CDELT3']
v_ctr = header13co['CRVAL3']
v_pix = header13co['CRPIX3']
velo_13co = (np.arange(nchan)*v_inc + v_inc + (v_ctr - v_inc*v_pix))

cubec18o, headerc18o =fits.getdata(jcmtc18o, header=True)
wcsc18o = wcs.WCS(jcmtc18o)
# Create velocity axis
nchan = headerc18o['NAXIS3']
v_inc = headerc18o['CDELT3']
v_ctr = headerc18o['CRVAL3']
v_pix = headerc18o['CRPIX3']
velo_c18o  = (np.arange(nchan)*v_inc + v_inc + (v_ctr - v_inc*v_pix))

rest_freq = 230.538e9 # CO 2-1 for ALMA data

fig = plt.figure(figsize=(7.5,7.5))

for i in range(len(sources)):
    source = sources[i]
    coords = source_coords[i]
    #coords_pix_12co = coords.to_pixel(wcs12co)
    coords_pix_13co = coords.to_pixel(wcs13co)
    coords_pix_c18o = coords.to_pixel(wcsc18o)
    spec_13co = cube13co[:,np.int(coords_pix_13co[1]),np.int(coords_pix_13co[0])]
    spec_c18o = cubec18o[:,np.int(coords_pix_13co[1]),np.int(coords_pix_13co[0])]
    #spec_12co = cube12co[:,np.int(coords_pix_12co[1]),np.int(coords_pix_12co[0])]
    plotParam = plottingDictionary[sources[i]]
    # Get spectra at source peaks rather than over beam:
    alma_12co, alma_header = fits.getdata(coImages[i],header=True)
    wcsalma = wcs.WCS(coImages[i])
    coords_pix_alma = coords.to_pixel(wcsalma)
    # Velocity:
    nchan = alma_header['NAXIS3']
    v_inc = alma_header['CDELT3']
    v_ctr = alma_header['CRVAL3']
    v_pix = alma_header['CRPIX3']
    freq_alma  = (np.arange(nchan)*v_inc + v_inc + (v_ctr - v_inc*v_pix))
    velo_alma = np.flip(c.c.cgs.value*(rest_freq-freq_alma)/rest_freq,0)/1e5
    spec_alma = np.flip(alma_12co[0,:,np.int(coords_pix_alma[1]),np.int(coords_pix_alma[0])],0)
    # Need extra dim in slices because ALMA data have four axes
    if i == 2:
        f0 = aplpy.FITSFigure(coImages[i],figure=fig,slices=[co_chans[i],0],
                              subplot=[0.13,0.13+0.3*i,0.33,0.28])
    else:
        f0 = aplpy.FITSFigure(coImages[i],figure=fig,slices=[co_chans[i],0],
                              subplot=[0.16,0.13+0.3*i,0.26,0.27])
    #f0.show_grayscale(vmin=co_vmin,vmax=co_vmax[i],stretch='linear')
    f0.show_colorscale(vmin=co_vmin,vmax=co_vmax[i],stretch='linear',cmap='magma')
    #f0.recenter(plotParam['rc'],plotParam['dc'],width=imw,height=imh)
    f0.recenter(source_coords[i].ra.deg,source_coords[i].dec.deg,width=imw,height=imh)
    if i == 2:
        f0.add_colorbar(location='top',axis_label_text=r'Jy beam$^{-1}$',width=0.1,
                        ticks=[-0.1,0,0.1,0.2],pad=0.07)
        f0.colorbar.set_font(family='sans_serif',size=labsize-1)
        f0.colorbar.set_axis_label_font(family='sans_serif',size=labsize-1)
    f0.ticks.set_color('white')
    f0.ticks.set_length(5)
    f0.ticks.set_xspacing(0.2*15./3600)
    f0.axis_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
    f0.tick_labels.set_style('colons')
    f0.tick_labels.set_font(size=labsize, weight='normal', \
                            stretch='normal', family='sans-serif', \
                            style='normal', variant='normal')
    f0.tick_labels.set_xformat('hh:mm:ss')
    f0.tick_labels.set_yformat('dd:mm:ss')
    if i in [1,2]:
        f0.axis_labels.hide_x()
    f0.show_beam(edgecolor='white',facecolor='none',corner='top left')
    f0.show_contour(contImages[i],levels=[x*plotParam['sigma'] for x in cont_levs],
                    colors='white',linewidths=l1)
    if i == 2:
        f0.add_label(0.9,0.85,co_vels_label[i],relative=True,color='white',size=(labsize-1),
                     horizontalalignment='right')
        f0.show_scalebar(0.00020,label='100 au',color='white',corner='bottom right',
                         fontsize=labsize)
    else:
        f0.add_label(0.1,0.7,co_vels_label[i],relative=True,color='white',size=(labsize-1),
                     horizontalalignment='left')
        f0.show_scalebar(0.00020,label='100 au',color='white',corner='bottom left',
                         fontsize=labsize)
    ax2 = fig.add_axes([0.55,0.13+0.3*i,0.35,0.26])
    vel, flux = np.loadtxt(co_spec[i],unpack=True,skiprows=8)
    xminor_locator = AutoMinorLocator(5)
    yminor_locator = AutoMinorLocator(2)
    ax2.xaxis.set_minor_locator(xminor_locator)
    ax2.yaxis.set_minor_locator(yminor_locator)
    if i == 1:
        #ax2.set_xticklabels([])
        plt.ylim(-0.25,0.4)
        plt.yticks([-0.2,-0.1,0,0.1,0.2,0.3,0.4])
    elif i == 0:
        plt.xlabel(r'$v_\mathrm{LSR}$ (km s$^{-1})$',fontsize=labsize)
        plt.ylim(-0.15,0.85)
        plt.yticks([-0.2,0,0.2,0.4,0.6,0.8])
    else:
        plt.ylim(-0.15,2.2)
        plt.yticks([0,0.5,1.0,1.5,2.0])
    ax2.tick_params(axis='both',which='major',labelsize=labsize)
    ax2.plot(velo_alma, spec_alma, color='black',zorder=2,label='$^{12}$CO (ALMA)')
    ax2.plot([-50,50],[0,0],color='gray',linestyle=':',zorder=1)
    #ax2.plot([co_vels[i],co_vels[i]],[-1,1],color='blue',linestyle='-',zorder=1,alpha=0.5)
    ax2.plot([nh3_vels[i],nh3_vels[i]],[-2,3],color='red',linestyle='--',linewidth=1.,zorder=1,alpha=0.5)
    plt.text(0.05,0.85,source_names[i],fontsize=labsize,horizontalalignment='left',
             transform=ax2.transAxes,)
    ax2.set_ylabel(r'Jy beam$^{-1}$',fontsize=labsize)
    if i == 2:
        plt.legend(frameon=False,fontsize=8,loc=1,handlelength=0.5)
    # Need to add second y-axis to deal with JCMT Tmb data
    ax3 = ax2.twinx()
    ax3.plot(velo_13co,spec_13co,color='darkgreen',alpha=0.5,zorder=1,linewidth=1,label='$^{13}$CO (JCMT)')
    ax3.plot(velo_c18o,spec_c18o,color='darkblue',alpha=0.5,zorder=1,linewidth=1,label='C$^{18}$O (JCMT)')
    #ax3.plot(velo_12co,spec_12co,color='darkorange',alpha=0.5,zorder=1,linewidth=1,label='$^{12}$CO')
    ax3.set_ylabel('$T_{MB}$ (K)')
    if source == 'SM1':
        ax2.set_xlim(vlsr_lim_SM1[0],vlsr_lim_SM1[1])
        ax3.set_xlim(vlsr_lim_SM1[0],vlsr_lim_SM1[1])
    elif source == '16263_2422':
        ax2.set_xlim(vlsr_lims[0],vlsr_lims[1])
        ax3.set_xlim(vlsr_lims[0],vlsr_lims[1])
        ax3.set_ylim(-5,40)
        #plt.legend(frameon=False,fontsize=8,handlelength=0.5,loc=1,bbox_to_anchor=(0.99,0.9))
    else:
        ax2.set_xlim(vlsr_lims[0],vlsr_lims[1])
        ax3.set_xlim(vlsr_lims[0],vlsr_lims[1])
        ax3.set_ylim(-5,30)
        plt.legend(frameon=False,fontsize=8,handlelength=0.5,loc=1,bbox_to_anchor=(0.985,0.9))
    align_yaxis(ax2,0,ax3,0)
    

fig.savefig('figures/CO_spectra_w_jcmt.pdf',bbox_inches='tight',dpi=300)
plt.close('all')


# Also plot channel maps? 
