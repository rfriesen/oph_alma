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

fig = plt.figure(figsize=(4,8))

sources = ['AN6','SM1','16263_2422']
labsize = 8
l1 = 1.
distance = 137.3*u.pc

SM1N_peak_loc = SkyCoord('16h26m27.173s -24d23m21.963s',ICRS)
SM1_xray = SkyCoord('16h26m27.88s -24d23m59.1s',ICRS)
N3mm_xray = SkyCoord('16h26m27.48s -24d24m18.1',ICRS)
SM1_cm   = SkyCoord('16h26m27.79s -24d24m02.2s',ICRS)
GSS30_IRS3_cm = SkyCoord('16h26m21.73s -24d22m51.4s',ICRS)

for i in range(len(sources)):
    source = sources[i]
    #if source == '16263_2422':
    #    plot_param = plottingDictionary['L16263_2422']
    #else:
    plot_param = plottingDictionary[source]
    image_file = plot_param['image_name']
    vmin = plot_param['vmin']
    vmax = plot_param['vmax']
    rc = plot_param['rc']
    dc = plot_param['dc']
    iw = plot_param['im_width']
    ih = plot_param['im_height']
    sigma = plot_param['sigma']
    beam_color='white'
    ang_sep = (plot_param['scalebar_size'].to(u.au)/distance).to(u.arcsec,equivalencies=u.dimensionless_angles())
    tick_spacing = 0.5
    if source in ['SM1','16263_2422']:
        clevs = [50.*sigma,500.*sigma]
        beam_color = 'black'
        tick_spacing = 1.
    elif source == 'SM1N':
        clevs = [5.*(x+1)*sigma for x in range(10)]
    else:
        clevs = [5.*sigma,7.*sigma,10.*sigma]

    f0 = aplpy.FITSFigure('continuum/{0}'.format(image_file), 
                          figure=fig, subplot=[0.14,0.19+0.28*i,0.8,0.24])
    f0.recenter(rc,dc,width=iw,height=ih)
    if plot_param['scale'] == 'linear':
        f0.show_colorscale(cmap='gist_heat',vmin=vmin,vmax=vmax,stretch='linear')
    else:
        f0.show_colorscale(cmap='gist_heat',vmin=vmin,vmax=vmax,
                           stretch=plot_param['scale'],vmid=plot_param['vmid'])
    f0.show_beam(edgecolor=beam_color,facecolor='none',corner='top left')
    if source == 'AN6':
        f0.ticks.set_color('white')
    else:
        f0.ticks.set_color('black')
    f0.ticks.set_length(5)
    f0.axis_labels.set_font(size=labsize, weight='normal', \
                            stretch='normal', family='sans-serif', \
                            style='normal', variant='normal')
    f0.tick_labels.set_style('colons')
    f0.tick_labels.set_font(size=labsize, weight='normal', \
                            stretch='normal', family='sans-serif', \
                            style='normal', variant='normal')
    f0.ticks.set_xspacing(tick_spacing*15./3600)
    f0.tick_labels.set_xformat('hh:mm:ss.s')
    f0.tick_labels.set_yformat('dd:mm:ss')
    if i != 0:
        f0.axis_labels.hide_x()
    if source == 'AN6':
        ticks=[-1.e-4,0,1.e-4,2.e-4,3.e-4,4.e-4,5.e-4]
    else:
        ticks=[0,0.01,0.03,0.06,0.09]
    f0.show_markers(SM1_xray.ra.value,SM1_xray.dec.value,marker='+',
                    edgecolor='blue',facecolor='mediumseagreen',s=12,zorder=4)
    f0.show_markers(N3mm_xray.ra.value,N3mm_xray.dec.value,marker='+',
                    edgecolor='blue',facecolor='mediumseagreen',s=12,zorder=4)
    f0.show_markers(SM1_cm.ra.value,SM1_cm.dec.value,marker='x',
                    edgecolor='blue',facecolor='dodgerblue',s=12,zorder=4)
    f0.show_markers(GSS30_IRS3_cm.ra.value,GSS30_IRS3_cm.dec.value,marker='x',
                    edgecolor='blue',facecolor='dodgerblue',s=12,zorder=4)
    f0.add_scalebar(ang_sep,label='{0:.0f} au'.format(plot_param['scalebar_size'].value),
                    color='white',corner=plot_param['scalebar_pos'])
    f0.add_colorbar(location='right',axis_label_text='Jy/beam',width=0.1,ticks=ticks,pad=0.07)
    f0.colorbar.set_font(family='sans_serif',size=7)
    f0.show_contour('continuum/{0}'.format(image_file), levels=clevs, 
                    colors='white', linewidths=l1*0.8)
    if source == '16263_2422':
        f0.add_label(0.5,0.88,'GSS 30',relative=True,fontsize=9,color='white',weight='bold',
                     horizontalalignment='center')
        f0.add_label(0.55,0.58,'IRS3',relative=True,size=8,color='white',horizontalalignment='left')
        f0.add_label(0.52,0.2,'IRS1',relative=True,size=8,color='white',horizontalalignment='right')
    elif source == 'SM1':
        f0.add_label(0.5,0.88,plot_param['image_label'],relative=True,fontsize=9,color='white',weight='bold',
                     horizontalalignment='center')
        f0.add_label(0.55,0.55,'SM1',relative=True,size=8,color='white',horizontalalignment='left')
        f0.add_label(0.56,0.08,'N3-mm',relative=True,size=8,color='white',horizontalalignment='right')
    else:
        f0.add_label(0.5,0.88,plot_param['image_label'],relative=True,fontsize=9,color='white',weight='bold',horizontalalignment='center')

fig.savefig('figures/compact_sources_v2.pdf',bbox_inches='tight')
plt.close('all')
