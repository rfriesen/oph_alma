#############################################
# Imaging the Continuuum
field = '5'

# Summary of the parameters
sourceName='SM1'
imagermode='csclean'

# Set the ms and continuum image name.
# Can use an averaged dataset but need all spws present so can apply 
# self-cal results to line spws
# Just do this once and use final_cont.ms for other sources
contvis = 'calibrated_final_cont_allspw.ms'         
finalvis = 'calibrated_final.ms'
contspw = ''
scims   = 'calibrated_final.ms'

# If you have complex line emission and no dedicated continuum
# windows, you will need to flag the line channels prior to averaging.
flagmanager(vis=finalvis,mode='save',
            versionname='before_cont_flags')

# Flag the "line channels"
flagchannels='0:1820~2040,4:1820~2040' # for all fields
flagdata(vis=finalvis,mode='manual', # field='3',
          spw=flagchannels,flagbackup=False)

rmtables(contvis)
os.system('rm -rf ' + contvis + '.flagversions')
split(vis=finalvis,
      spw=contspw,      
      outputvis=contvis,
      width=[128,128,128,128,128,128,128,128], 
      datacolumn='data') # the split ms file includes all targets

listobs(vis=contvis, listfile=contvis+'.listobs')

# If you flagged any line channels, restore the previous flags
flagmanager(vis=finalvis,mode='restore',
            versionname='before_cont_flags')

# Inspect continuum for any problems. Figure out channels to use for 
# continuum imaging. 
plotms(vis=contvis,xaxis='uvdist',yaxis='amp',coloraxis='spw')


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

clearcal(vis=contvis)
delmod(vis=contvis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

contspw = '0:0~12;18~29,1,2,3,4:0~12;18~29,5,6,7' # Fix

clean(vis=contvis,
      imagename=contimagename,
      field=field,
      spw=contspw,
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
# std dev = 2.5 x 10-4 Jy/beam
# beam 0.6565" x 0.5461", PA 82.1 deg
# Bright compact source 9.6 x 10-2 Jy/beam
# Fainter compact source 5.2 x 10-3 Jy/beam
# Need self-cal

# 2. Self-cal
os.system('rm -rf phase.cal')
gaincal(vis=contvis,
        caltable='phase.cal',
        field=field,
        solint='60s',
        combine='scan',
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
      outputvis=contvis+'_'+sourceName+'_selfcal.ms',
      datacolumn='corrected')

# Clean again
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+'_selfcal'+ext)

clean(vis=contvis+'_'+sourceName+'_selfcal.ms',
      imagename=contimagename+'_selfcal',
      field=field,
      spw=contspw,
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
      gain=0.2,
      imagermode = imagermode)

# After 500 iterations, rms = 1.1 x 10-4 Jy/beam
# Bright compact source 1.02 x 10-1 Jy/beam
# Fainter compact source 5.5 x 10-3 Jy/beam

# Try round 2
os.system('rm -rf phase_2.cal')
gaincal(vis=contvis+'_'+sourceName+'_selfcal.ms',
        caltable='phase_2.cal',
        field=field,
        solint='60s',
        combine='scan',
        calmode='p',
        refant='DV17',
        gaintype='G')

plotcal(caltable='phase_2.cal',
        xaxis='time',yaxis='phase',
        subplot=331,iteration='antenna',timerange='24-May-2015/0:4:00~24-May-2015/0:6:30',
        plotrange=[0,0,-30,30],markersize=5,fontsize=10,
        figfile='selfcal_phase_scan.png')

applycal(vis=contvis+'_'+sourceName+'_selfcal.ms',
         field=field,
         gaintable=['phase_2.cal'],
         interp='linear')

# Split out self-cal corrected data
# Try adding field parameter to speed things up? 
split(vis=contvis+'_'+sourceName+'_selfcal.ms',
      field=field,
      outputvis=contvis+'_'+sourceName+'_selfcal_2.ms',
      datacolumn='corrected')

# Clean again
clean(vis=contvis+'_'+sourceName+'_selfcal_2.ms',
      imagename=contimagename+'_selfcal_2',
      #field=field,
      spw=contspw,
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

# After 800 iterations, rms = 1.1 x 10-4 Jy/beam
# Bright compact source 1.02 x 10-1 Jy/beam
# Fainter compact source 5.5 x 10-3 Jy/beam
# No difference. 

#################################################################
os.system('rm -rf amp.cal')
gaincal(vis=contvis+'_'+sourceName+'_selfcal_2.ms',
        caltable='amp.cal',
        #field=field,
        solint='60s',
        combine='scan',
        calmode='ap',
        refant='DV17',
        gaintype='G',
        solnorm=True)

plotcal(caltable='amp.cal',
        xaxis='time',yaxis='amp',subplot=331,
        iteration='antenna',plotrange=[0,0,0,0],
        markersize=5,fontsize=10)
applycal(vis=contvis+'_'+sourceName+'_selfcal_2.ms',
         #field=field,
         gaintable=['amp.cal'],
         interp='linear')

split(vis=contvis+'_'+sourceName+'_selfcal_2.ms',
      outputvis=contvis+'_'+sourceName+'_selfcal_3.ms',
      datacolumn='corrected')

# Clean again
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+'_selfcal_3'+ext)

clean(vis=contvis+'_'+sourceName+'_selfcal_3.ms',
      imagename=contimagename+'_selfcal_3',
      #field=field,
      spw=contspw,
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
      gain=0.2,
      imagermode = imagermode)

# After 300 iterations with higher gain, rms = 8.6 x 10-5 Jy/beam
# Bright compact source 1.02 x 10-1 Jy/beam
# Fainter compact source 5.5 x 10-3 Jy/beam

########################################
# Continuum Subtraction for Line Imaging
# Apply self-cal tables
applycal(vis=finalvis,
         field=field,
         gaintable=['phase.cal','phase_2.cal','amp.cal'],
         interp='linear')

plotms(vis=finalvis,
       field=field,spw='0',
       xaxis='channel',yaxis='amp',
       avgtime='9999')

# May need to adjust windows, self-absorption at SM1 position

fitspw = '0:0~1500;2400~3839,4:0~1500;2400~3839' # line-free channel for fitting continuum
linespw = '0,4' # line spectral windows. You can subtract the continuum from multiple spectral line windows at once.

uvcontsub(vis=finalvis,
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

linevis = finalvis+'.contsub' 
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
width='0.16km/s' # 0.4?
nchan = 250
outframe='lsrk'
veltype='radio'
restfreq='230.538GHz'
spw='0,1'
mscales=[0,16,40]
gain=0.5

# Set the ms and line image name.
# Does cleaning keep going beyond threshold????
lineName = 'CO21' # name of transition (see science goals in OT for name) 
lineimagename = sourceName+'_'+lineName+'_0.16_image' # name of line image

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=linevis)
#delmod(vis=linevis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(lineimagename + ext)

clean(vis=linevis,
      imagename=lineimagename, 
      mask=linemaskname,
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

# After 1400 iterations, achieved a RMS level of 4.3e-3 Jy/beam over 0.25 km/s
# Beam: 0.84"x0.72", 89.1 deg

immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='81~109',#'33~44', #'54~69',
          outfile=lineimagename+'_blue.m0')

immoments(imagename=lineimagename+'.image',
          moments=0,
          chans='125~159',#'50~62', #'80~96',
          outfile=lineimagename+'_red.m0')

os.system('rm -rf '+lineimagename+'.m1')
immoments(imagename=lineimagename+'.image',
          moments=1,
          chans='81~159', #'54~96',
          includepix=[30.e-3,10],
          outfile=lineimagename+'.m1')

os.system('rm -rf '+lineimagename+'.m1')
immoments(imagename=lineimagename+'.image',
          moments=1,
          chans='34~65', #'54~96',
          includepix=[30.e-3,20],
          outfile=lineimagename+'.m1')

exportfits(imagename=contimagename+'_selfcal_3.image',fitsimage=contimagename+'_selfcal_3.image.fits')
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

