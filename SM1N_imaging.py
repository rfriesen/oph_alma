#############################################
# Imaging the Continuuum
field = '6'

# Summary of the parameters
sourceName='SM1N'
imagermode='csclean'

# Set the ms and continuum image name.
# Want to use original dataset so can apply self-cal results to line spws
scims   = '../data/calibrated_final.ms'
cont_spw = '1,2,3,5,6,7'
contvis = '../calibrated_final_cont.ms'         

# 1. No taper, multiscale
# Parameters
cell='0.075arcsec'
imsize = [600,600]
weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'

contimagename = sourceName+'_continuum'

# If necessary, run the following commands to get rid of older clean
# data.

clearcal(vis=scims)
delmod(vis=scims)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=scims, #contvis,
      imagename=contimagename,
      field=field,
      spw=cont_spw,
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
      gain=0.5,
      imagermode = imagermode)

# After 300 iterations:
# std dev = 4.6 x 10-5 Jy/beam
# beam 0.65" x 0.54", PA 82.6 deg
# Extended emission
# Peak flux density 1.6 x 10-3 Jy/beam
# Source is extended and I think not bright enough for self-cal. 

# Looking for compact emission
compact_imagename = sourceName+'_compact_continuum'
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(compact_imagename+ext)

clean(vis=scims, #contvis,
      imagename=compact_imagename,
      field=field,
      spw=cont_spw,
      mode='mfs',
      psfmode='clark',
      #multiscale=[0,16,40], # Note this is in pixels 0, 2xbeam, 5xbeam
      uvrange='>110m',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = 0.5,
      uvtaper=True,
      outertaper=['1.5arcsec'],
      niter = niter, 
      threshold = threshold, 
      interactive = True,
      gain=0.5,
      imagermode = imagermode)
# No evidence for compact emission above noise levels. *Maybe* 2-sigma. 

########################################
# Continuum Subtraction for Line Imaging

plotms(vis=scims,
       field=field,spw='0',
       xaxis='channel',yaxis='amp',
       avgtime='9999')

# Windows look good

fitspw = '0:0~1820;2040~3839,4:0~1820;2040~3839' # line-free channel for fitting continuum
linespw = '0,4' # line spectral windows. You can subtract the continuum from multiple spectral line windows at once.

uvcontsub(vis=scims,
          field=field,
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

linevis = scims+'.contsub' 
listobs(vis=linevis, listfile=linevis+'.listobs')

##############################################
# Image line emission 
### Summary of the parameters
# May need to image larger region - some bright CO emission to north
imagermode='csclean'
cell='0.075arcsec'
imsize = [600,600]
weighting = 'briggs'
robust=0.5
niter=10000
threshold = '3mJy'
start='-15km/s'  # Good if ignoring VLA 1623 outflow. 
width='0.4km/s'
nchan = 100
outframe='lsrk'
veltype='radio'
restfreq='230.538GHz'
spw='0,1'
mscales=[0,16,40]
gain=0.5

# Set the ms and line image name.
# Does cleaning keep going beyond threshold????
lineName = 'CO21' # name of transition (see science goals in OT for name) 
lineimagename = sourceName+'_'+lineName+'_0.4_image' # name of line image

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=linevis)
#delmod(vis=linevis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(lineimagename + ext)

clean(vis=linevis,
      imagename=lineimagename, 
      mask='SM1N_CO21_0.4_keep.mask',#linemaskname,
      #field=field,
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

# After 1400 iterations, achieved a RMS level of 4.0e-3 Jy/beam over 0.25 km/s
# Beam: 0.76"x0.62", 86.8` deg

immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='20~44',#'40~68',
          outfile=lineimagename+'_blue.m0')
# rms 0.027
immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='50~63',#'79~96',
          outfile=lineimagename+'_red.m0')
# rms 0.008
os.system('rm -rf '+lineimagename+'.m1')
immoments(imagename=lineimagename+'.image',
          moments=1,
          chans='20~63',
          includepix=[30.e-3,10],
          outfile=lineimagename+'.m1')
# High v
immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='20~33',
          outfile=lineimagename+'_blue_hv.m0')
# rms 0.027
immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='57~63',
          outfile=lineimagename+'_red_hv.m0')

exportfits(imagename=contimagename+'.image',fitsimage=contimagename+'.image.fits')
exportfits(imagename=lineimagename+'_blue.m0',fitsimage=lineimagename+'_blue.m0.fits')
exportfits(imagename=lineimagename+'_red.m0',fitsimage=lineimagename+'_red.m0.fits')
# If you'd like to redo your clean, but don't want to make a new mask
# use the following commands to save your original mask. This is an
# optional step.
linemaskname = 'field_'+field+'_line.mask'
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

