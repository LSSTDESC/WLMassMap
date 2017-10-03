# This module contains the code to create an estimated shear map from a shape
# catalog
import numpy as np
from optparse import OptionParser
import os
import yaml
from numpy.linalg import pinv
from projection import project_flat
from astropy.table import Table
from astropy.io import fits

def add_metacal_shear(catalog, delta_gamma=0.01):
    """
    Computes the responsivity for metacalibration measurements
    TODO: Add selection effects
    """
    R1 = (catalog['mcal_g_1p'] - catalog['mcal_g_1m']) / (2 * delta_gamma)
    R2 = (catalog['mcal_g_2p'] - catalog['mcal_g_2m']) / (2 * delta_gamma)
    R = np.stack([R1, R2],axis=1)

    # Averages the responsivity matrix over entire sample
    R = np.mean(R, axis=0)

    # Inverts the responsivity matrix
    Rinv = pinv(R)

    # Computes the estimated shear
    catalog['g'] = np.dot(Rinv, np.array(catalog['mcal_g']).T).T
    return catalog

def bin_shear_map(catalog, nx, ny):
    """
    Computes the shear map
    """
    # Bin the shear catalog
    g1map = np.bincount(catalog['pixel_index'],
                        weights=catalog['g'][:,0],
                        minlength=nx*ny).reshape((nx, ny))
    g2map = np.bincount(catalog['pixel_index'],
                        weights=catalog['g'][:,1],
                        minlength=nx*ny).reshape((nx, ny))
    Nmap  = np.bincount(catalog['pixel_index'],
                        minlength=nx*ny).reshape((nx, ny))

    # Normalize by number of galaxies
    nz_ind = Nmap > 0
    g1map[nz_ind] /= Nmap[nz_ind]
    g2map[nz_ind] /= Nmap[nz_ind]

    gmap = np.stack([g1map,g2map], axis=0)

    return gmap, Nmap

def shear_map(config):
    """
    Builds a shear map from a given shape catalog

    Parameters
    ----------
        config: dictionary
            Configuration dictionary read from yaml config file
    """
    filename = os.path.join(config['input_directory'], config['input_filename'])
    catalog = Table.read(filename)

    # Applies a first set of cuts
    sel = (catalog['flags']==0) & (catalog['mcal_s2n_r'] > 10) & (catalog['mcal_T']/catalog['psfrec_T'] > 0.5)
    catalog = catalog[sel]

    # Computes calibrated shear
    catalog = add_metacal_shear(catalog)

    # Applies projection to the catalog
    if config['projection']['type'] in ['Gnomonic']:
        c = config['projection']
        catalog, grid_ra, grid_dec = project_flat(catalog, c['nx'], c['ny'],
                    c['pixel_size'], c['center_ra'], c['center_dec'], c['type'])

        # Bins the projected catalog
        gmap, nmap = bin_shear_map(catalog, nx=c['nx'], ny=c['ny'])
    else:
        raise NotImplementedError

    # Saves the resulting map
    # In the case of a spherical map, only saves the shear map nmap
    # For projected map also saves the ra,dec of each pixels
    phdu = fits.PrimaryHDU(gmap)
    nhdu = fits.ImageHDU(nmap)
    rahdu = fits.ImageHDU(grid_ra)
    dechdu = fits.ImageHDU(grid_dec)

    hdulist = fits.HDUList([phdu, nhdu, rahdu, dechdu])

    filename = os.path.join(config['output_directory'], config['output_filename'])
    hdulist.writeto(filename)


if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    shear_map(config['shear_map'])
