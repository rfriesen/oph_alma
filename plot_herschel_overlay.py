######################################################################
# Image script for Cycle 2 ALMA project
# 16163_2422 figures
######################################################################
import aplpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

JCMT850m = 'jcmt/Ophiuchus_L1688_20140506_s850_noco_KMP.fits'
Hers70m  = 'herschel/pabraham_3_70micron_highpass.fits'
nh3_dir  = 'nh3'
nh3_opha = '{0}/opha.23693.atca.uvredo.res.m0.tb.fits'.format(nh3_dir)
nh3_ophb = '{0}/b.avg.tb.2pix.m0.fits'.format(nh3_dir)
nh3_iso1 = '{0}/16263-2422.23693.m0.tb.fits'.format(nh3_dir)
nh3_iso2 = '{0}/16267-2417.23693.m0.tb.fits'.format(nh3_dir)

sources = ['AN6','SM1N','SM1','B2-A7','16263_2422','16267_2417']
ra_list = [246.631667,246.613333,246.616083,246.8854167,246.590167,246.681542]
de_list = [-24.414444,-24.389556,-24.399861,-24.4491667,-24.381139,-24.289778]
vmin = 0.
vmax_list = [4.,4.,4.,0.01,10.,0.04]
vint = [1,1,1,0.003,2,0.01]

# For all targets
# Field of view ~ 27.4 arcsec
csize = 0.003804
# Image width, height (degrees)
iw = 0.01
ih = 0.01

# Contours for JCMT (needs adjusting for different sources?)
jcmtlev = [0.03,0.05,0.07,0.09,0.15,0.2,0.3,0.4]
# Contours for ATCA
tbsig = 1.0
ophalev = [x*2. * tbsig + 3. * tbsig for x in range(10)]
tbsig = 0.6
ophblev = [x*2. * tbsig + 3. * tbsig for x in range(10)]

labsize=9
l1=1

fig = plt.figure(figsize=(7.5,5))

for i in range(6):
    if i < 3:
        ypos = 0.05
        xpos = i*0.325+0.09
        cbar_label_text = ''
    else:
        ypos = 0.53
        xpos = (i-3)*0.325 + 0.09
        cbar_label_text = 'Jy pixel$^{-1}$'
    subplot=[xpos,ypos,0.25,0.37]
    f0 = aplpy.FITSFigure(Hers70m,figure=fig,subplot=subplot)
    f0.recenter(ra_list[i],de_list[i],width=iw,height=ih)
    f0.show_colorscale(cmap='gist_heat',vmin=vmin,vmax=vmax_list[i])
    f0.ticks.set_color('white')
    #plt.tick_params(direction='in')
    f0.axis_labels.hide()
    f0.tick_labels.set_style('colons')
    f0.tick_labels.set_font(size=labsize, weight='normal', \
                                stretch='normal', family='sans-serif', \
                                style='normal', variant='normal')
    f0.ticks.set_xspacing(1*15./3600)
    f0.tick_labels.set_xformat('hh:mm:ss')
    f0.tick_labels.set_yformat('dd:mm:ss')
    colorbar_step=vint[i]
    ticks = np.arange(np.ceil((vmax_list[i]-vmin)/colorbar_step)+2)*colorbar_step + np.int(vmin)
    f0.add_colorbar(location='top',axis_label_text=cbar_label_text,width=0.1,ticks=ticks)
    f0.colorbar.set_axis_label_font(family='sans_serif',size=labsize+1)
    f0.colorbar.set_font(family='sans_serif',size=labsize)
    f0.show_circles(ra_list[i],de_list[i],csize,edgecolor='white',facecolor='none',zorder=3)
    f0.show_contour(JCMT850m,levels=jcmtlev,colors='grey',linewidths=l1,zorder=1)
    # NH3 contours
    if sources[i] in ['AN6','SM1']:
        f0.show_contour(nh3_opha,levels=ophalev,colors='violet',linewidths=l1,zorder=2)
    elif sources[i] in ['SM1N']:
        f0.show_contour(nh3_opha,levels=[x*1.5 for x in ophalev],colors='violet',linewidths=l1,zorder=2)
    elif sources[i] in ['B2-A7']:
        f0.show_contour(nh3_ophb,levels=ophblev,colors='violet',linewidths=l1,zorder=2)
    elif sources[i] in ['16263_2422']:
        f0.show_contour(nh3_iso1,levels=ophalev,colors='violet',linewidths=l1,zorder=2)
    elif sources[i] in ['16267_2417']:
        f0.show_contour(nh3_iso2,levels=[x*0.8 for x in ophblev],colors='violet',linewidths=l1,zorder=2)

    f0.add_label(0.95,0.93,sources[i],relative=True,size=8,color='white',weight='bold',
                 horizontalalignment='right')

fig.savefig('figures/alma_cycle2_herschel.pdf')
plt.close(fig)
plt.close('all')
