# This module computes convergence maps
from optparse import OptionParser
import os
import yaml
import numpy as np
from astropy.io import fits

from kaiser_squires import flat_KS_map

def convergence_map(config):
    """
    Computes convergence map with specified algorithm
    """

    filename = os.path.join(config['input_directory'], config['input_filename'])

    # Extracts the shear map
    gmap = fits.getdata(filename, 0)

    if config['algorithm']['name'] == 'flat_ks':
        kappa_e, kappa_b = flat_KS_map(gmap)
    else:
        raise NotImplementedError

    phdu = fits.PrimaryHDU(kappa_e)
    exthdu = fits.ImageHDU(kappa_b)

    hdulist = fits.HDUList([phdu, exthdu])

    filename = os.path.join(config['output_directory'], config['output_filename'])
    hdulist.writeto(filename)

if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    convergence_map(config['convergence_map'])
