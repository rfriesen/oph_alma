# Quick image script for Cycle 2 ALMA project
# AN6 figures
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
import numpy as np
import aplpy
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

from config import plottingDictionary

fig = plt.figure(figsize=(4.,7.))

sources = ['N6','B2-A7']
srcPltName = ['AN6t','B2A7']
labsize = 10
l1 = 1.

for i in range(len(sources)):
    source = sources[i]
    plotConfig = srcPltName[i]
    #if source == '16263_2422':
    #    plot_param = plottingDictionary['L16263_2422']
    #else:
    plot_param = plottingDictionary[plotConfig]
    image_file = 'continuum/{0}'.format(plot_param['image_name'])
    vmin = plot_param['vmin']
    vmax = plot_param['vmax']
    rc = plot_param['rc']
    dc = plot_param['dc']
    iw = plot_param['im_width']
    ih = plot_param['im_height']
    sigma = plot_param['sigma']
    clevs = [2.*x*sigma + 3.*sigma for x in range(10)]

    f0 = aplpy.FITSFigure(image_file, 
                          figure=fig, subplot=[0.14,0.15+0.4*i,0.85,0.36])
    f0.recenter(rc,dc,width=iw,height=ih)
    if plot_param['scale'] == 'linear':
        f0.show_colorscale(cmap='gist_heat',vmin=vmin,vmax=vmax,stretch='linear')
    else:
        f0.show_colorscale(cmap='gist_heat',vmin=vmin,vmax=vmax,
                           stretch=plot_param['scale'],vmid=plot_param['vmid'])
    f0.show_beam(edgecolor='white',facecolor='none',corner=plot_param['beam_loc'])
    f0.ticks.set_color('white')
    f0.axis_labels.set_font(size=labsize, weight='normal', \
                            stretch='normal', family='sans-serif', \
                            style='normal', variant='normal')
    f0.tick_labels.set_style('colons')
    f0.tick_labels.set_font(size=labsize, weight='normal', \
                            stretch='normal', family='sans-serif', \
                            style='normal', variant='normal')
    f0.ticks.set_xspacing(1.0*15./3600)
    f0.tick_labels.set_xformat('hh:mm:ss.s')
    f0.tick_labels.set_yformat('dd:mm:ss')
    if i != 0:
        f0.axis_labels.hide_x()
    f0.add_scalebar(0.00028,label='100 AU',color='white',corner='bottom left')
#    if source == 'AN6':
    ticks=[-1.e-4,0,1.e-4,2.e-4,3.e-4,4.e-4,5.e-4]
#    else:
#        ticks=[0,0.01,0.03,0.06,0.09]
    f0.add_colorbar(location='right',axis_label_text='Jy/beam',width=0.1,ticks=ticks)
    f0.colorbar.set_font(family='sans_serif',size=8)
    f0.show_contour(image_file, levels=clevs, colors='white', linewidths=l1*0.8)
    f0.show_contour(image_file,levels=[-1.*x for x in clevs[::-1]], linewidths=l1*0.8,linestyles='dashed',colors='white')
    f0.add_label(0.9,0.88,plot_param['image_label'],relative=True,fontsize=9,color='white',weight='bold',horizontalalignment='right')

fig.savefig('figures/taper_plots.pdf',bbox_inches='tight')
plt.close('all')
