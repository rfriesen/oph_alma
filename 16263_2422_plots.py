# Quick image script for Cycle 2 ALMA project
# 16163_2422 figures
import aplpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import units as u

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

from config import plottingDictionary

source = '16263_2422'
plot_param = plottingDictionary[source]

JCMTcont = 'jcmt/Ophiuchus_L1688_20140506_s850_noco_KMP.fits'
image = 'continuum/{0}_continuum_selfcal_2.image.fits'.format(source)
COimageBlue = 'co_21/{0}/{0}_CO21_image_blue.m0.fits'.format(source)
COimageRed  = 'co_21/{0}/{0}_CO21_image_red.m0.fits'.format(source)
COimageHVBlue = 'co_21/{0}/{0}_CO21_image_blue_hv.m0.fits'.format(source)
COimageHVRed  = 'co_21/{0}/{0}_CO21_image_red_hv.m0.fits'.format(source)

labsize=11
l1 = 0.8

image_file = plot_param['image_name']
image_file = 'continuum/{0}'.format(image_file)
vmin = plot_param['vmin']
vmax = plot_param['vmax']
vmid = plot_param['vmid']
rc = plot_param['rc']
dc = plot_param['dc']
iw = plot_param['im_width']
ih = plot_param['im_height']
sigma = plot_param['sigma']
clevs = [5.*sigma,50.*sigma,500.*sigma]

# Contours for CO
co_sigma = 0.03
co_levs = [6.*x*co_sigma + 6.*co_sigma for x in range(25)]
# Contours for high-velocity CO
co_hv_sigma = 0.01
co_hv_levs = [5.*x*co_hv_sigma + 5.*co_hv_sigma for x in range(15)]

# Contours for JCMT
jcmtlev = [0.03,0.05,0.07,0.09]

# Distance
d = 137. 
sb_500au = 500/d * u.arcsec

fig = plt.figure(figsize=(4.,7.))

f0 = aplpy.FITSFigure(image_file, figure=fig, subplot=[0.12,0.58,0.85,0.4])
f0.recenter(rc,dc,width=iw,height=ih)
f0.show_grayscale(vmin=vmin,vmax=vmax,vmid=vmid,stretch='log',invert=True)
f0.show_beam(edgecolor='black',facecolor='none',corner='top left')
f0.ticks.set_color('black')
f0.axis_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f0.axis_labels.hide_x()
f0.tick_labels.set_style('colons')
f0.tick_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f0.ticks.set_xspacing(1.*15./3600)
f0.ticks.set_length(5)
f0.tick_labels.set_xformat('hh:mm:ss.s')
f0.tick_labels.set_yformat('dd:mm:ss')
f0.add_scalebar(sb_500au.to(u.deg).value,label='500 AU',
                color='black',corner='bottom left')
f0.add_colorbar(location='right',axis_label_text='Jy/beam',width=0.1,
                ticks=[0,0.01,0.03,0.06,0.09],pad=0.08)
f0.colorbar.set_font(family='sans_serif',size=labsize-2)
f0.show_contour(COimageBlue,levels=co_levs,colors='blue',linewidths=l1,zorder=4)
f0.show_contour(COimageRed,levels=co_levs,colors='red',linewidths=l1,zorder=3)
#f0.show_contour(JCMTcont,levels=jcmtlev,colors='grey',linewidths=l1,zorder=1)
f0.show_contour(image, levels=clevs, colors='black', linewidths=l1*0.8,zorder=5)
f0.add_label(0.53,0.57,'IRS3',relative=True,size=9,color='black',horizontalalignment='left')
f0.add_label(0.65,0.23,'IRS1',relative=True,size=9,color='black',horizontalalignment='left')
#f0.add_label(0.9,0.1,'a)',relative=True,fontsize=9,color='white')


# High velocity CO
rc_zoom = 246.59042
dc_zoom = -24.38067
f1 = aplpy.FITSFigure(image_file,figure=fig,subplot=[0.12,0.13,0.85,0.4])
f1.recenter(rc_zoom,dc_zoom,width=iw*0.3,height=ih*0.3)
f1.show_grayscale(vmin=vmin,vmax=vmax,vmid=vmid,stretch='log',invert=True)
f1.show_beam(edgecolor='black',facecolor='none',corner='top left')
f1.ticks.set_color('black')
f1.axis_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f1.tick_labels.set_style('colons')
f1.tick_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f1.ticks.set_xspacing(0.5*15./3600)
f1.ticks.set_length(3)
f1.tick_labels.set_xformat('hh:mm:ss.s')
f1.tick_labels.set_yformat('dd:mm:ss')
f1.add_scalebar(sb_500au.to(u.deg).value,label='500 AU',
                color='black',corner='top right',linewidth=1.5)
f1.add_colorbar(location='right',axis_label_text='Jy/beam',width=0.1,
                ticks=[0,0.01,0.03,0.06,0.09],pad=0.08)
f1.colorbar.set_font(family='sans_serif',size=labsize-2)
f1.show_contour(image, levels=clevs, colors='black', linewidths=l1,zorder=5)
f1.show_contour(COimageHVBlue,levels=co_hv_levs,colors='blue',linewidths=l1*2.,zorder=4)
f1.show_contour(COimageHVRed,levels=[x*1. for x in co_hv_levs],colors='red',linewidths=l1*2.,zorder=3)
f1.add_label(0.63,0.56,'IRS3',relative=True,size=11,color='black',horizontalalignment='left')

f1.save('figures/16263_2422_ALMA_CO_v2.pdf',dpi=200)
plt.close('all')

'''
rc_zoom = 246.59042
dc_zoom = -24.38067
f1 = aplpy.FITSFigure(image)
f1.recenter(rc_zoom,dc_zoom,width=iw*0.3,height=ih*0.3)
f1.show_grayscale(vmin=vmin,vmax=vmax,vmid=vmid,stretch='log')
f1.show_beam(edgecolor='yellow',facecolor='black',corner='top left',linewidth=1.5)
f1.ticks.set_color('white')
f1.axis_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f1.tick_labels.set_style('colons')
f1.tick_labels.set_font(size=labsize, weight='normal', \
                        stretch='normal', family='sans-serif', \
                        style='normal', variant='normal')
f1.ticks.set_xspacing(0.5*15./3600)
f1.tick_labels.set_xformat('hh:mm:ss.s')
f1.tick_labels.set_yformat('dd:mm:ss')
#f1.show_circles(rc,dc,csize,edgecolor='white',facecolor='none',zorder=2)
f1.add_scalebar(sb_500au.to(u.deg).value,label='500 AU',
                color='yellow',corner='top right',linewidth=1.5)
f1.add_colorbar(location='right',axis_label_text='Jy/beam',width=0.3)
f1.colorbar.set_font(family='sans_serif',size=labsize-1)
#f0.show_contour(sm1NImage, levels=c1n, colors='white', linewidths=l1)
f1.show_contour(image, levels=clevs, colors='white', linewidths=l1*1.5,zorder=1)
#f1.show_contour(JCMTcont,levels=jcmtlev,colors='grey',linewidths=l1,zorder=1)
f1.show_contour(COimageHVBlue,levels=co_hv_levs,colors='blue',linewidths=l1*2.)
f1.show_contour(COimageHVRed,levels=[x*1.7 for x in co_hv_levs],colors='red',linewidths=l1*2.)

f1.save('figures/16263_2422_ALMA_CO_hv.pdf',dpi=200)
plt.close('all')
'''
