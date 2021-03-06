# This module contains the code to create an estimated shear map from a shape
# catalog
import numpy as np
from optparse import OptionParser
import os
import yaml
from numpy.linalg import pinv
from .projection import project_flat, project_healpix
from astropy.table import Table
from astropy.io import fits
import healpy as hp

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

def bin_shear_map(catalog, nx=None, ny=None, npix=None):
    """
    Computes the shear map by binning the catalog according to pixel_index.
    Either nx,ny or npix must be provided.

    Parameters
    ----------
    catalog: table
        Input shape catalog with pixel_index column

    nx,ny: int, optional
        Number of pixels of a 2d flat map

    npix: int, optional
        Number of pixels of a spherical map (or other 1D pixelating scheme)

    Returns
    -------
    gmap: ndarray
        Shear map

    nmap: ndarray
        Number of galaxies per pixels
    """
    assert (npix is not None) or ((nx is not None) and (ny is not None))

    # Bin the shear catalog
    if npix is None:
        npix = nx*ny

    g1map = np.bincount(catalog['pixel_index'],
                        weights=catalog['g'][:,0],
                        minlength=npix)
    g2map = np.bincount(catalog['pixel_index'],
                        weights=catalog['g'][:,1],
                        minlength=npix)
    Nmap  = np.bincount(catalog['pixel_index'], minlength=npix)

    # Rebining as a 2D map if requested
    if nx is not None:
        g1map = g1map.reshape((nx, ny))
        g2map = g2map.reshape((nx, ny))
        Nmap = Nmap.reshape((nx, ny))

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
    filename = config['input_filename']
    catalog = Table.read(filename)

    # Computes calibrated shear
    catalog = add_metacal_shear(catalog)

    # Extracts projection configuration
    c = config['projection']
    projection2d = False

    # Applies projection to catalog
    if c['type'] in ['gnomonic']: # Any 2D flat projection
        projection2d =True

        catalog, grid_ra, grid_dec = project_flat(catalog, c['nx'], c['ny'],
                    c['pixel_size'], c['center_ra'], c['center_dec'], c['type'])

        # Bins the projected catalog
        gmap, nmap = bin_shear_map(catalog, nx=c['nx'], ny=c['ny'])

    elif c['type'] == 'healpix': # Any spherical projection
        catalog = project_healpix(catalog, nside=c['nside'])

        # Bins the projected catalog
        gmap, nmap = bin_shear_map(catalog, npix=hp.nside2npix(c['nside']))
    else:
        raise NotImplementedError

    # Saves the resulting map
    # In the case of a spherical map, only saves the shear map nmap
    # For projected map also saves the ra,dec of each pixels
    phdu = fits.PrimaryHDU(gmap)
    nhdu = fits.ImageHDU(nmap)

    # In the case of 2D map, we are also saving the coordinate grid
    if projection2d:
        rahdu = fits.ImageHDU(grid_ra)
        dechdu = fits.ImageHDU(grid_dec)
        hdulist = fits.HDUList([phdu, nhdu, rahdu, dechdu])
    else:
        hdulist = fits.HDUList([phdu, nhdu])


    filename = config['output_filename']
    hdulist.writeto(filename)


if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    shear_map(config['shear_map'])
