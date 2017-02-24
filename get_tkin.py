##################################################
# Script to get Tkin from GBT GAS data for
# targets in L1688. 
# Also ATCA/VLA data where available. 
##################################################
import numpy as np
from astropy.coordinates import Angle,SkyCoord
from astropy.coordinates import ICRS,Galactic
import astropy.units as u
from astropy import wcs
from astropy.io import fits
from regions import CircleSkyRegion

def mean_over_beam(coords,par_file,epar_file):
        breff, wcs = get_beam_reff(par_file)
        par_hdu = fits.open(par_file)
        epar_hdu = fits.open(epar_file)
        sky_region = CircleSkyRegion(coords,breff)
        pix_region = sky_region.to_pixel(wcs)
        mask = pix_region.to_mask()
        par_masked = mask.cutout(par_hdu[0].data)
        epar_masked = mask.cutout(epar_hdu[0].data)
        weighted_mean, weighted_mean_err = calc_weighted_mean(par_masked,epar_masked)
        return weighted_mean, weighted_mean_err

def calc_weighted_mean(par_data,epar_data):
        par_data[par_data == 0] = np.nan
        epar_data[epar_data == 0] = np.nan
        numer = np.nansum(par_data/(epar_data**2.))
        denom = np.nansum(1./(epar_data**2.))
        weighted_mean = numer/denom
        weighted_mean_err = 1./denom
        return weighted_mean, weighted_mean_err

def get_beam_reff(file):
        hdulist = fits.open(file)
        w = wcs.WCS(hdulist[0].header)
        bmaj = hdulist[0].header['BMAJ'] * u.deg
        bmin = hdulist[0].header['BMIN'] * u.deg
        breff = np.sqrt(bmaj*bmin)
        return breff, w

source_ra_set   = ['16h26m21.64s','16h26m43.57s','16h26m27.86s','16h26m27.20s','16h26m31.60s','16h27m32.50s']
source_de_set   = ['-24d22m52.1s','-24d17m23.2s','-24d23m59.5s','-24d23m22.4s','-24d24m52.0s','-24d26m27.0s']
source_coords = [SkyCoord(ra,de,frame='icrs') for ra,de in zip(source_ra_set,source_de_set)]

sources = ['16263_2422','16267_2417','SM1','SM1N','AN6','B2-A7']
tkin_dir  = 'tk_data'
gbt_tk_file = '{0}/L1688_Tkin_DR1_rebase3_flag.fits'.format(tkin_dir)
gbt_etk_file = '{0}/L1688_eTkin_DR1_rebase3_flag.fits'.format(tkin_dir)
atca_tk_file = '{0}/atca_opha_tk.mask.fits'.format(tkin_dir)
atca_etk_file = '{0}/atca_opha_dtk.mask.fits'.format(tkin_dir)

# ALMA footprint is smaller than the beam. Just want to average over GBT beam 
# area for now. Weight by uncertainty in Tk fit. 
# Use ellipses for this

f = open('{0}/GAS_Tkin.txt'.format(tkin_dir),'w')

for i in range(len(sources)):
        source_coord_i = source_coords[i]
        gbt_weighted_mean, gbt_weighted_mean_err = mean_over_beam(source_coord_i,gbt_tk_file,gbt_etk_file)
        if sources[i] in ['SM1','SM1N','AN6']:
                atca_weighted_mean, atca_weighted_mean_err = mean_over_beam(source_coord_i,atca_tk_file,atca_etk_file)
        else:
                atca_weighted_mean, atca_weighted_mean_err = 0,0
        f.write('{0:11s}  {1:23s}  {2:4.1f}  {3:6.3f} {4:4.1f} {5:5.2f} \n'.format(sources[i],source_coord_i.to_string(style='hmsdms',sep=':'),gbt_weighted_mean,gbt_weighted_mean_err,atca_weighted_mean,atca_weighted_mean_err))

f.close()
