# Script to convolve 1mm data to match 3mm 
from astropy.io import fits
import astropy.units as u

# Beam parameters
file3mm = '3mm/OphAMosaic.image.fits'
hdu3mm = fits.open(file3mm)
major = hdu3mm[0].header['BMAJ']*u.deg
minor = hdu3mm[0].header['BMIN']*u.deg
pa = hdu3mm[0].header['BPA']*u.deg

# SM1N continuum image
sm1n_cont_image = 'continuum/SM1N_continuum.pbcor'
ia.open(sm1n_cont_image)
im2 = ia.convolve2d(outfile=sm1n_cont_image+'.convol',axes=[0,1],
	type='gauss',major='{0}arcsec'.format(major.to(u.arcsec).value),
	minor='{0}arcsec'.format(minor.to(u.arcsec).value),
	pa='{0}deg'.format(pa.value),overwrite=True)
im2.done()
ia.close()

# Boom!
# Write to fits
exportfits(imagename=sm1n_cont_image+'.convol',fitsimage=sm1n_cont_image+'.convol.fits')

# SM1:
sm1_cont_image = 'continuum/SM1_continuum_selfcal_3.pbcor'
ia.open(sm1_cont_image)
im2 = ia.convolve2d(outfile=sm1_cont_image+'.convol',axes=[0,1],
	type='gauss',major='{0}arcsec'.format(major.to(u.arcsec).value),
        minor='{0}arcsec'.format(minor.to(u.arcsec).value),
        pa='{0}deg'.format(pa.value),overwrite=True)
im2.done()
ia.close()

# Boom!
# Write to fits
exportfits(imagename=sm1_cont_image+'.convol',fitsimage=sm1_cont_image+'.convol.fits')

# Want to do the same thing for ALMA Cycle 0 359 GHz data
# SM1N:
sm1n_cont_image = 'cycle0data/SM1N_contall_sc.image'
sm1n_pbcor_image = 'cycle0data/SM1N_contall_sc.image.pbcor'
# First need to correct for the primary beam
impbcor(imagename=sm1n_cont_image,
        pbimage='cycle0data/SM1N_contall_sc.flux',
        outfile=sm1n_pbcor_image, overwrite=True)

# Convolve
ia.open(sm1n_pbcor_image)
im2 = ia.convolve2d(outfile=sm1n_pbcor_image+'.convol',axes=[0,1],
	type='gauss',major='{0}arcsec'.format(major.to(u.arcsec).value),
	minor='{0}arcsec'.format(minor.to(u.arcsec).value),
	pa='{0}deg'.format(pa.value),overwrite=True)
im2.done()
ia.close()

# Boom!
# Write to fits
os.system('rm -rf {0}.{1}'.format(sm1n_pbcor_image,'convol.fits'))
exportfits(imagename=sm1n_pbcor_image+'.convol',fitsimage=sm1n_pbcor_image+'.convol.fits')

# SM1:
sm1_cont_image = 'cycle0data/SM1_contall_apcal.image'
sm1_pbcor_image = 'cycle0data/SM1_contall_apcal.image.pbcor'
# First need to correct for the primary beam
impbcor(imagename=sm1_cont_image,
        pbimage='cycle0data/SM1_contall_apcal.flux',
        outfile=sm1_pbcor_image, overwrite=True)

# Convolve
ia.open(sm1_pbcor_image)
im2 = ia.convolve2d(outfile=sm1_pbcor_image+'.convol',axes=[0,1],
	type='gauss',major='{0}arcsec'.format(major.to(u.arcsec).value),
	minor='{0}arcsec'.format(minor.to(u.arcsec).value),
	pa='{0}deg'.format(pa.value),overwrite=True)
im2.done()
ia.close()

# Boom!
# Write to fits
os.system('rm -rf {0}.{1}'.format(sm1_pbcor_image,'convol.fits'))
exportfits(imagename=sm1_pbcor_image+'.convol',fitsimage=sm1_pbcor_image+'.convol.fits')


# Beam parameters for other images
# 16263-2422:
file3mm = '3mm/162622-24225_selfcal3.fits'
hdu3mm = fits.open(file3mm)
major = hdu3mm[0].header['BMAJ']*u.deg
minor = hdu3mm[0].header['BMIN']*u.deg
pa = hdu3mm[0].header['BPA']*u.deg

mm_cont_image = 'continuum/16263_2422_continuum_selfcal_2.pbcor'

ia.open(mm_cont_image)
im2 = ia.convolve2d(outfile=mm_cont_image+'.convol',axes=[0,1],
	type='gauss',major='{0}arcsec'.format(major.to(u.arcsec).value),
	minor='{0}arcsec'.format(minor.to(u.arcsec).value),
	pa='{0}deg'.format(pa.value),overwrite=True)
im2.done()
ia.close()

# Write to fits
exportfits(imagename=mm_cont_image+'.convol',fitsimage=mm_cont_image+'.convol.fits')

# 16267_2417
file3mm = '3mm/162644-24173.fits'
hdu3mm = fits.open(file3mm)
major = hdu3mm[0].header['BMAJ']*u.deg
minor = hdu3mm[0].header['BMIN']*u.deg
pa = hdu3mm[0].header['BPA']*u.deg

mm_cont_image = 'continuum/16267_2417_continuum_taper.pbcor'
ia.open(mm_cont_image)
im2 = ia.convolve2d(outfile=mm_cont_image+'.convol',axes=[0,1],
	type='gauss',major='{0}arcsec'.format(major.to(u.arcsec).value),
	minor='{0}arcsec'.format(minor.to(u.arcsec).value),
	pa='{0}deg'.format(pa.value),overwrite=True)
im2.done()
ia.close()

# Write to fits
exportfits(imagename=mm_cont_image+'.convol',fitsimage=mm_cont_image+'.convol.fits')


