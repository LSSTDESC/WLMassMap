# This module computes convergence maps
from optparse import OptionParser
import os
import yaml
import numpy as np
from astropy.io import fits
import pyssht as ssht
import massmappy.cy_mass_mapping as mm
import massmappy.cy_healpy_mass_mapping as hp_mm
import healpy as hp

from kaiser_squires import flat_KS_map, healpix_KS_map

def convergence_map(config):
    """
    Computes convergence map with specified algorithm
    """

    filename = config['input_filename']

    # Extracts the shear map
    gmap = fits.getdata(filename, 0)
    c = config['algorithm']

    if c['name'] == 'flat_ks':
        kappa_e, kappa_b = flat_KS_map(gmap)

    elif c['name'] == 'healpix_ks':
        kappa_e, kappa_b = healpix_KS_map(gmap, lmax=c['lmax'])

    elif c['name'] == 'SKS_mw':
        gmap1 = gmap[0]
        gmap2 = gmap[1]

        gamma_mw = gmap1 + 1j*gmap2
        k_mw = mm.reduced_shear_to_kappa_mw(gamma_mw,L=c['L'],Method='MW',sigma=(np.pi/(60.*180.*2.355))*(c['fwhm_arcmins']/c['L']), Iterate=c['Iterate'],return_count=c['return_count'])
        kappa_e = np.real(k_mw)
        kappa_b = np.imag(k_mw)

    elif c['name'] == 'SKS_hp':
        gmap1 = gmap[0].astype('<f8')
        gmap2 = gmap[1].astype('<f8')

        kappa_e, kappa_b = hp_mm.reduced_shear_to_kappa_hp(gmap1, gmap2, L=c['L'], Nside=c['nside'], sigma=(np.pi/(60.*180.*2.355))*(c['fwhm_arcmins']/c['L']), Iterate=c['Iterate'])

        if c['Convert_to_MW_pixels'] is True:
            alm_E_hp    = hp.map2alm(kappa_e, lmax=c['L']-1)
            alm_B_hp    = hp.map2alm(kappa_b, lmax=c['L']-1)

            alm_E_hp_mw    = mm.lm_hp2lm(alm_E_hp, c['L'])
            alm_B_hp_mw    = mm.lm_hp2lm(alm_B_hp, c['L'])

            kappa_e = ssht.inverse(alm_E_hp_mw, c['L'], Reality=True)
            kappa_b = ssht.inverse(alm_B_hp_mw, c['L'], Reality=True)

    else:
        raise NotImplementedError

    phdu = fits.PrimaryHDU(kappa_e)
    exthdu = fits.ImageHDU(kappa_b)

    hdulist = fits.HDUList([phdu, exthdu])

    filename = config['output_filename']
    hdulist.writeto(filename)

if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    convergence_map(config['convergence_map'])
