#############################################
# Imaging the Continuuum
field = '7'

# Summary of the parameters
sourceName='AN6'
imagermode='csclean'

# Set the ms and continuum image name.
contvis = '../calibrated_final_cont.ms'         

# 1. No taper, multiscale
# Parameters
cell='0.075arcsec'
imsize = [600,600]
weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'

contimagename = 'calibrated_final_cont_image_'+sourceName

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=contvis)
#delmod(vis=contvis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=contvis,
      imagename=contimagename,
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

# After 200 iterations:
# std dev = 0.4 mJy/beam
# beam 0.6569" x 0.54628", PA 81.7 deg
# Peak S/N ~ 7 in compact source (2.9 mJy/beam), some extended emission

# 2. Taper, multiscale
cell='0.075arcsec'
imsize = [600,600]
weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+'_taper'+ext)

clean(vis=contvis,
      imagename=contimagename+'_taper',
      field=field,
      mode='mfs',
      psfmode='clark',
      multiscale=[0,16,40],
      uvtaper=True,
      outertaper=['1.5arcsec'],
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      interactive = True,
      imagermode = imagermode)

# After 500 iterations
# std dev = 0.65 mJy/beam
# beam 1.4196" x 1.3555" PA 111.3 deg
# Peak S/N ~ 10 in compact source 5.5 mJy/beam
# Extended emission (S/N ~ 3-5) and negative bowls

########################################
# Continuum Subtraction for Line Imaging

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

# Move contsub file to this folder
os.system('mv '+finalvis+'.contsub calibrated_final_'+sourceName+'.contsub')
linevis = 'calibrated_final_'+sourceName+'.contsub'
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
start='-10km/s'  # Good if ignoring VLA 1623 outflow. 
width='0.25km/s'
nchan = 120
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

clean(vis=linevis,
      imagename=lineimagename, 
      #mask=linemaskname,
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

os.system('rm -rf '+lineimagename+'_red.m0')
immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='60~70',
          outfile=lineimagename+'_red.m0')

os.system('rm -rf '+lineimagename+'.m1')
immoments(imagename=lineimagename+'.image',
          moments=1,
          chans='62~68',
          includepix=[20.e-3,10],
          outfile=lineimagename+'.m1')

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

myimages = glob.glob("*.m0")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("*.pbcor")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("*.flux")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True) 

