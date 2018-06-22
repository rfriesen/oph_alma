#############################################
# Imaging the Continuuum
field = '3'

# Summary of the parameters
sourceName='16263_2422'
imagermode='csclean'

# Set the ms and continuum image name.
#contvis = '../calibrated_final_cont.ms'         
# Edit to output self-cal ms and save it!! 
contvis = 'calibrated_final_cont_statwt_{0}.ms'.format(sourceName)

# 1. No taper, multiscale
# Parameters
cell='0.075arcsec'
imsize = [600,600]
weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'

#contimagename = 'calibrated_final_cont_image_'+sourceName
contimagename = '{0}_continuum'.format(sourceName)

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=contvis)
#delmod(vis=contvis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=contvis,
      imagename=contimagename,
      #field=field,
      mode='mfs',
      psfmode='clark',
      multiscale=[0,16,40], # Note this is in pixels 0, 2xbeam, 5xbeam
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      interactive = True,
      imagermode = imagermode)

# After 300 iterations:
# std dev = 7 x 10-5 Jy/beam
# beam 0.655688" x 0.545758", PA 82.0 deg
# Bright compact source 0.103 Jy/beam
# Fainter compact source 6.17 x 10-3 Jy/beam
# Hm.. these values match self-cal results below.. ?
# Yep, self-cal applied. 

# 2. Self-cal
os.system('rm -rf phase.cal')
gaincal(vis=contvis,
        caltable='phase.cal',
        #field=field,
        solint='40s',
        calmode='p',
        refant='DV17',
        gaintype='G')

plotcal(caltable='phase.cal',
        xaxis='time',yaxis='phase',
        subplot=331,iteration='antenna',
        plotrange=[0,0,-30,30],markersize=5,fontsize=10,
        figfile='selfcal_phase_scan.png')

applycal(vis=contvis,
         field=field,
         gaintable=['phase.cal'],
         interp='linear')

# Split out self-cal corrected data
split(vis=contvis,
      outputvis=contvis+'_selfcal.ms',
      datacolumn='corrected')

# Clean again
clean(vis=contvis+'_selfcal.ms',
      imagename=contimagename+'_selfcal',
      field=field,
      mode='mfs',
      psfmode='clark',
      multiscale=[0,16,40], # Note this is in pixels 0, 2xbeam, 5xbeam
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      interactive = True,
      imagermode = imagermode)

# After 300 iterations, rms = 9.5 x 10-5 Jy/beam
# Bright compact source 1.03 x 10-1 Jy/beam
# Fainter compact source 6.18 x 10-3 Jy/beam

# Try round 2
os.system('rm -rf phase_2.cal')
gaincal(vis=contvis+'_selfcal.ms',
        caltable='phase_2.cal',
        field=field,
        solint='60s',
        calmode='p',
        refant='DV17',
        gaintype='G')

plotcal(caltable='phase_2.cal',
        xaxis='time',yaxis='phase',
        subplot=331,iteration='antenna',timerange='24-May-2015/0:4:00~24-May-2015/0:6:30',
        plotrange=[0,0,-30,30],markersize=5,fontsize=10,
        figfile='selfcal_phase_scan.png')
# Doesn't look like can do much with this. Ignore. Try amp self-cal. 
os.system('rm -rf amp.cal')
gaincal(vis=contvis+'_selfcal.ms',
        caltable='amp.cal',
        field=field,
        solint='60s',
        calmode='ap',
        refant='DV17',
        gaintype='G',
        solnorm=True)

plotcal(caltable='amp.cal',
        xaxis='time',yaxis='amp',subplot=331,
        iteration='antenna',plotrange=[0,0,0,0],
        markersize=5,fontsize=10)
applycal(vis=contvis+'_selfcal.ms',
         field=field,
         gaintable=['amp.cal'],
         interp='linear')

split(vis=contvis+'_selfcal.ms',
      outputvis=contvis+'_selfcal_2.ms',
      datacolumn='corrected')

# Clean again
clean(vis=contvis+'_selfcal_2.ms',
      imagename=contimagename+'_selfcal_2',
      field=field,
      mode='mfs',
      psfmode='clark',
      multiscale=[0,16,40], # Note this is in pixels 0, 2xbeam, 5xbeam
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      interactive = True,
      imagermode = imagermode)

# After 300 iterations, rms = 8.1 x 10-5 Jy/beam
# Bright compact source 1.03 x 10-1 Jy/beam
# Fainter compact source 6.13 x 10-3 Jy/beam

########################################
# Continuum Subtraction for Line Imaging

plotms(vis='../calibrated_final.ms',
       field=field,spw='0',
       xaxis='channel',yaxis='amp',
       avgtime='9999')

# Windows are fine here

fitspw = '0:0~1820;2040~3839,4:0~1820;2040~3839' # line-free channel for fitting continuum
linespw = '0,4' # line spectral windows. You can subtract the continuum from multiple spectral line windows at once.

finalvis='../calibrated_final.ms'

uvcontsub(vis=finalvis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False) # This value should not be changed.

# NOTE: Imaging the continuum produced by uvcontsub with
# want_cont=True will lead to extremely poor continuum images because
# of bandwidth smearing effects. For imaging the continuum, you should
# always create a line-free continuum data set using the process
# outlined above.

linevis = finalvis+'.contsub' 
listobs(vis=linevis, listfile=linevis+'.listobs')

##############################################
# Image line emission 
### Summary of the parameters
imagermode='csclean'
cell='0.075arcsec'
imsize = [600,600]
weighting = 'briggs'
robust=2.0
niter=1000
threshold = '3mJy'
start='-15km/s'  # Good if ignoring VLA 1623 outflow. 
width='0.25km/s'
nchan = 160
outframe='lsrk'
veltype='radio'
restfreq='230.538GHz'
spw='0,1'
mscales=[0,16,40]
gain=0.5

# Set the ms and line image name.
# Does cleaning keep going beyond threshold????
lineName = 'CO21' # name of transition (see science goals in OT for name) 
lineimagename = sourceName+'_'+lineName+'_image' # name of line image

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=linevis)
#delmod(vis=linevis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(lineimagename + ext)

# SO MUCH EMISSION. Also absorption at brightest source. Bad continuum subtraction? But channels look fine. 
clean(vis=linevis,
      imagename=lineimagename, 
      mask=linemaskname,
      field=field,
      spw=spw,
      mode='velocity',
      start=start,
      width=width,
      nchan=nchan, 
      outframe=outframe, 
      veltype=veltype, 
      restfreq=restfreq, 
      niter=niter,  
      threshold=threshold, 
      interactive=True,
      cell=cell,
      imsize=imsize, 
      weighting=weighting, 
      robust=robust,
      multiscale=mscales,
      gain=gain,
      imagermode=imagermode)

# After 1400 iterations, achieved a RMS level of 3.5e-3 Jy/beam over 0.5 km/s
# Peak intensity: 8.2e-1 Jy/beam, S/N ~ 230
# Beam: 0.76"x0.62", 85 deg

immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='39~65',
          outfile=lineimagename+'_blue.m0')

immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='82~107',
          outfile=lineimagename+'_red.m0')
# Full range
os.system('rm -rf '+lineimagename+'.m1')
immoments(imagename=lineimagename+'.image',
          moments=1,
          chans='20~115',
          includepix=[50.e-3,20],
          outfile=lineimagename+'.m1')

exportfits(imagename=lineimagename+'_blue.m0',fitsimage=lineimagename+'_blue.m0.fits')
exportfits(imagename=lineimagename+'_red.m0',fitsimage=lineimagename+'_red.m0.fits')

# If you'd like to redo your clean, but don't want to make a new mask
# use the following commands to save your original mask. This is an
# optional step.
linemaskname = sourceName+'_'+lineName+'_save.mask'
## rmtables(linemaskname) # uncomment if you want to overwrite the mask.
os.system('cp -ir ' + lineimagename + '.mask ' + linemaskname)

##############################################
# Apply a primary beam correction

import glob

myimages = glob.glob("*.image")

rmtables('*.pbcor')
for image in myimages:
    impbcor(imagename=image, pbimage=image.replace('.image','.flux'), outfile = image.replace('.image','.pbcor'))

##############################################
# Export the images

import glob

myimages = glob.glob("*.image")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("*.pbcor")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("*.flux")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True) 

################################################
# Adding mid-velocity outflow images
sourceName='16263_2422'
lineName = 'CO21' # name of transition (see science goals in OT for name) 
lineimagename = sourceName+'_'+lineName+'_image' # name of line image
os.system('rm -rf '+lineimagename+'_blue_mid.m0')
immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='52~62',
          outfile=lineimagename+'_blue_mid.m0')

os.system('rm -rf '+lineimagename+'_red_mid.m0')
immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='80~91',
          outfile=lineimagename+'_red_mid.m0')

exportfits(imagename=lineimagename+'_blue_mid.m0',fitsimage=lineimagename+'_blue_mid.m0.fits')
exportfits(imagename=lineimagename+'_red_mid.m0',fitsimage=lineimagename+'_red_mid.m0.fits')

#
