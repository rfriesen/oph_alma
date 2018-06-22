# Create figure showing location of Cycle 2 ALMA-observed cores
# Inset to zoom in on Oph A

import aplpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

def gaussian_profile(r,sigma,A):
    # Assume r, sigma both in arcsec
    gaussian = A* np.exp(-1.*(r-0)**2./(2.*sigma**2.))
    return gaussian

# Files
dataDir  = 'continuum/'
bimaFile = 'bima/opha.bima.n2h.mom0.fits'
jcmtFile = 'jcmt/Ophiuchus_L1688_20140506_s850_noco_KMP.fits'
sm1nFile  = dataDir + 'SM1N_continuum.image.fits'
an6File   = dataDir + 'AN6_continuum_taper.image.fits'
# List of compact sources
ra_src, de_src = np.genfromtxt('compact_sources.txt',usecols=(1,2),unpack=True,comments='#',dtype='str')

# Flux, pointing, etc. values for figure
# Oph - all
oph_ctr = SkyCoord('16h26m57s -24d24m0s',ICRS)
oph_ra = oph_ctr.ra.deg
oph_de = oph_ctr.dec.deg
oph_iw = 0.35
oph_ih = 0.35
oph_vmin = -0.005
oph_vmax = 0.5
# JCMT 
j_stretch = 'linear'
j_vmin = -3
j_vmax = 13
# Try centering JCMT image on SM1
j_ih = 0.045
j_iw = 0.055
j_sig = 0.025
j_cl = [j_sig * x + j_sig for x in range(4)]
# BIMA contours
b_cl = [2,4,6,8,10,12,14]

# SM1N 
sm1n_stretch = 'linear'
sm1n_vmin = -0.001
sm1n_vmax = 0.003
sm1n_rc = 246.613333
sm1n_dc = -24.389556
#sm1n_rc = 246.61354 # From Cycle 0
#sm1n_dc = -24.39017
sm1n_c1 = [0.0009,0.0018,0.0027]
sm1n_compact = SkyCoord('16h26m27.173s -24d23m21.963s',ICRS)

l1 = 0.6

# SM1
sm1_stretch = 'log'
sm1_vmin = -0.001
sm1_vmax = 0.3
sm1_vmid = -0.002
sm1_cl = [0.006,.1,.2]
sm1_rc = 246.616
sm1_dc = -24.3998

# N6
n6_rc = 246.63167
n6_dc = -24.41444
# Compact source coord
n6_compact = SkyCoord('16:26:31.735 -24:24:52.0402',ICRS,unit=(u.hourangle,u.deg))

#16263_2422
coords_GSS30 = SkyCoord('16:26:21.64 -24:22:52.1',ICRS,unit=(u.hourangle,u.deg))
rc_16263 = 246.590167
dc_16263 = -24.381139

#B2-A7
coords_b2a7 = SkyCoord('16:27:32.5 -24:26:57.0',ICRS,unit=(u.hourangle,u.deg))

# 16267-2417
coords_16267 = SkyCoord('16:26:43.57 -24:17:23.2',ICRS,unit=(u.hourangle,u.deg))

# Field of view ~ 27.4 arcsec
csize = 0.003804
# Image width, height (degrees)
iw = 0.0078
ih = 0.0078

# Good for both ALMA plots
#csize = 0.00243
#iw = 0.0048
#ih = 0.0052
labsize=9
msize = 15
mark  = '+'
mline = 0.1
mc = 'orange'

fig = plt.figure(figsize=(8,6))
# Big image of Oph
f1 = aplpy.FITSFigure(jcmtFile, figure=fig, subplot=[0.1,0.2,0.5,0.7])
f1.recenter(oph_ra,oph_de,width=oph_iw,height=oph_ih)
f1.show_grayscale(vmin=oph_vmin,vmax=oph_vmax,stretch='sqrt',vmid=j_vmin*1.1,invert=True)
f1.ticks.set_color('black')
f1.axis_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f1.tick_labels.set_style('colons')
f1.tick_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f1.tick_labels.set_xformat('hh:mm:ss')
f1.tick_labels.set_yformat('dd:mm')
f1.show_contour(jcmtFile,levels=j_cl, colors='gray', linewidths=l1, zorder=1)
f1.show_rectangles(sm1_rc-0.001, sm1_dc+0.00,j_iw*1.1,j_ih*1.5,linewidths=l1*2.,edgecolor='black',facecolor='None')
f1.show_scalebar(0.042,label='0.1 pc',color='black',corner='top left')
f1.show_markers(sm1n_rc,sm1n_dc,c=mc,marker=mark,s=msize,linewidths=mline)
f1.show_markers(sm1_rc,sm1_dc,c=mc,marker=mark,s=msize,linewidths=mline)
f1.show_markers(n6_rc,n6_dc,c=mc,marker=mark,s=msize,linewidths=mline)
f1.show_markers(rc_16263,dc_16263,c=mc,marker=mark,s=msize,linewidths=mline)
f1.show_markers(coords_b2a7.ra.deg,coords_b2a7.dec.deg,c=mc,marker=mark,s=msize,linewidths=mline)
f1.show_markers(coords_16267.ra.deg,coords_16267.dec.deg,c=mc,marker=mark,s=msize,linewidths=mline)


f0 = aplpy.FITSFigure(bimaFile, figure=fig, subplot=[0.7,0.35,0.34,0.55])
#f0.aplpy.FITSFigure(jcmtFile,figure=fig,subplt=[0.13,0.15,0.8,0.8])
f0.recenter(sm1_rc-0.001, sm1_dc+0.005, width=j_iw, height=j_ih)
f0.show_grayscale(vmin=j_vmin,vmax=j_vmax,invert=True)
f0.ticks.set_color('black')
f0.axis_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
#f0.axis_labels.hide_y()
f0.tick_labels.set_style('colons')
f0.tick_labels.set_font(size=labsize-1, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f0.tick_labels.set_xformat('hh:mm:ss')
f0.tick_labels.set_yformat('dd:mm:ss')
f0.ticks.set_length(5)
f0.show_beam(edgecolor='black',facecolor='none',corner='top left')
f0.show_contour(bimaFile,levels=b_cl, colors=['white','white','lightgray','lightgray','gray','darkgray'], linewidths=l1)
#f0.show_contour(jcmtFile,levels=j_cl, colors='gray', linewidths=l1, zorder=1)
f0.show_circles(sm1n_rc,sm1n_dc,csize,edgecolor='yellow',facecolor='none',zorder=3)
f0.show_circles(sm1_rc,sm1_dc,csize,edgecolor='yellow',facecolor='none',zorder=3)
f0.show_circles(n6_rc,n6_dc,csize,edgecolor='yellow',facecolor='none',zorder=3)
f0.show_circles(rc_16263,dc_16263,csize,edgecolor='yellow',facecolor='none',zorder=3)
#f0.show_rectangles(sm1n_rc,sm1n_dc,iw,ih,edgecolor='red',
#                   facecolor='none',zorder=3,linewidth=2)
#f0.show_rectangles(sm1_rc,sm1_dc,iw,ih,edgecolor='red',
#                   facecolor='none',zorder=3,linewidth=2)
#f0.show_rectangles(n6_rc,n6_dc,iw,ih,edgecolor='red',
#                   facecolor='none',zorder=3,linewidth=2)
f0.show_scalebar(0.0023,label='1000 AU',color='black',corner='bottom right')
f0.add_label(sm1_rc+0.008,sm1_dc,'SM1', color='black',family='sans-serif',
             size='small',weight='bold')
f0.add_label(sm1n_rc+0.009,sm1n_dc,'SM1N',color='black',family='sans-serif',
             size='small',weight='bold')
f0.add_label(n6_rc+0.0055,n6_dc+0.0035,'N6', color='black',family='sans-serif',
             size='small',weight='bold')
f0.add_label(rc_16263+0.008,dc_16263+0.004,'GSS 30', color='black',family='sans-serif',
             size='small',weight='bold')
#f0.show_markers(246.60999,-24.4085,marker='*',edgecolor='red',s=100,zorder=3)
#f0.show_markers(246.61579,-24.40061,marker='+',edgecolor='black',s=100,zorder=3)
#f0.show_markers(246.61617,-24.39975,marker='D',edgecolor='black',s=30,zorder=3)
for i in range(len(ra_src)):
    coords = '{0} {1}'.format(ra_src[i],de_src[i])
    c = SkyCoord(coords,ICRS,unit=(u.hourangle,u.deg))
    f0.show_markers(c.ra.deg,c.dec.deg,marker='o',edgecolor='darkorange',facecolor='gold',s=17,zorder=4)

f0.show_markers(n6_compact.ra.deg,n6_compact.dec.deg,marker='o',edgecolor='red',
                facecolor='khaki',s=20,zorder=4)
f0.show_markers(sm1n_compact.ra.deg,sm1n_compact.dec.deg,marker='D',edgecolor='red',
                facecolor='khaki',s=15,zorder=4)

fig.savefig('figures/Oph_core_locs.pdf',bbox_inches='tight',dpi=250)
plt.close('all')

