# This module computes convergence maps
from optparse import OptionParser
import os
import yaml
import numpy as np
from astropy.io import fits

from .kaiser_squires import flat_KS_map, healpix_KS_map

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
