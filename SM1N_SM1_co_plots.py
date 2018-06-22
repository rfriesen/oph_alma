# Quick image script for Cycle 2 ALMA project
# 16163_2422 figures
import aplpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy import units as u

from config import plottingDictionary

sources = ['SM1','SM1N']

JCMTcont = '/media/DATAPART/data/jcmt/gbls/oph/Ophiuchus_L1688_20140506_s850_noco_KMP.fits'

labsize=11
l1 = 0.8

# Contours for CO
co_sigma = 0.03
co_levs_r = [3.*x*co_sigma + 3.*co_sigma for x in range(15)]
co_levs_b = [5.*x*co_sigma + 5.*co_sigma for x in range(15)]
# Multiplier for CO contours different for sources. 
# Contours for JCMT
jcmtlev = [0.03,0.05,0.07,0.09]

# Distance
d = 137. 
sb_500au = 500/d * u.arcsec

fig = plt.figure(figsize=(3.8,7.))

for i in range(len(sources)):
    source = sources[i]
    plot_param = plottingDictionary[source]
    if source == 'SM1N':
        COimageBlue = 'co_21/{0}/{0}_CO21_image_blue.m0.fits'.format(source)
        COimageRed  = 'co_21/{0}/{0}_CO21_image_red.m0.fits'.format(source)
    elif source == 'SM1':
        COimageBlue = 'co_21/{0}/{0}_CO21_0.4_image_blue.m0.fits'.format(source)
        COimageRed  = 'co_21/{0}/{0}_CO21_0.4_image_red.m0.fits'.format(source)
    else:
        print('Not SM1 or SM1N')
    image_file = plot_param['image_name']
    image = 'continuum/{0}'.format(image_file)
    vmin = plot_param['vmin']
    vmax = plot_param['vmax']
    rc = plot_param['rc']
    dc = plot_param['dc']
    iw = 0.008 #plot_param['im_width']
    ih = 0.008 #plot_param['im_height']
    sigma = plot_param['sigma']
    clev_mult = plot_param['clev_mult']
    if source in ['SM1','16263_2422']:
        clevs = [50.*sigma,500.*sigma]
        beam_loc = 'bottom left'
        cz = 4
    else:
        clevs = [clev_mult*(x+1)*sigma for x in range(15)]
        beam_loc = plot_param['beam_loc']
        cz = 1
    neglevs = [-1.*clev_mult*x*sigma - clev_mult*sigma for x in range(3)]
    neglevs = neglevs[::-1]
    f0 = aplpy.FITSFigure(image, figure=fig, subplot=[0.14,0.13+0.45*i,0.85,0.4])
    f0.recenter(rc,dc,width=iw,height=ih)
    if plot_param['scale'] == 'linear':
        f0.show_grayscale(vmin=vmin,vmax=vmax,invert=True)
    else:
        vmid = plot_param['vmid']
        f0.show_grayscale(vmin=vmin,vmax=vmax,vmid=vmid,stretch=plot_param['scale'],
                          invert=True)
    f0.show_beam(edgecolor='black',facecolor='none',corner=plot_param['beam_loc'])
    f0.ticks.set_color('black')
    f0.axis_labels.set_font(size=labsize, weight='normal', \
                            stretch='normal', family='sans-serif', \
                            style='normal', variant='normal')
    if i !=0:
        f0.axis_labels.hide_x()
    f0.tick_labels.set_style('colons')
    f0.tick_labels.set_font(size=labsize, weight='normal', \
                            stretch='normal', family='sans-serif', \
                            style='normal', variant='normal')
    f0.ticks.set_xspacing(0.5*15./3600)
    f0.tick_labels.set_xformat('hh:mm:ss.s')
    f0.tick_labels.set_yformat('dd:mm:ss')
    f0.add_scalebar(0.00028,label='100 AU',color='black',corner=plot_param['cbar_loc'])
    f0.add_colorbar(location='right',axis_label_text='Jy/beam',width=0.1,
                    ticks=plot_param['cbar_ticks'])
    f0.colorbar.set_font(family='sans_serif',size=labsize-1)
    f0.show_contour(image, levels=clevs, colors='black', linewidths=l1,zorder=cz)
    #f0.show_contour(image, levels=neglevs, colors='black', linewidths=l1,
    #                linestyles='dashed')
    f0.add_label(0.78,0.93,source,relative=True,fontsize=9,color='black')
    f0.show_contour(COimageBlue,levels=co_levs_b,colors='blue',linewidths=l1*1.25)
    f0.show_contour(COimageRed,levels=co_levs_r,colors='red',linewidths=l1*1.25)

fig.savefig('figures/{0}_{1}_CO.pdf'.format(sources[0],sources[1]),dpi=200,bbox_inches='tight')
plt.close('all')

