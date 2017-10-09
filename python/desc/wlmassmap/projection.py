# This module handles the projection of a catalog on a specific grid
import numpy as np
from scipy.stats import binned_statistic_2d
from projection_utils import radec2xy, xy2radec, eq2ang
import healpy as hp

def project_healpix(catalog, nside, hp_type='RING'):
    """
    Adds a HEALpix pixel index to all galaxies in the catalog

    Parameters
    ----------
    catalog: table
        Input shape catalog

    nside: int
        HEALpix nside parameter

    hp_type: string
        HEALpix pixel order ('RING', 'NESTED')

    Returns
    -------
    catalog: table
        Output shape catalog with pixel index column
    """
    theta, phi = eq2ang(catalog['ra'], catalog['dec'])
    catalog['pixel_index'] = hp.ang2pix(nside, theta, phi,
                                        nest=(hp_type=='NESTED'))
    return catalog

def project_flat(catalog, nx, ny, pixel_size, center_ra, center_dec, projection='gnomonic'):
    """
    Adds a pixel index for a Gnomonic projected map. Pixels are indexed
    starting from 0 according to ind = y * nx + x

    Parameters
    ----------

    catalog: table
        Input shape catalog

    nx: int
        Number of pixels along the x axis

    ny: int
        Number of pixels along the y axis

    pixel_size: float
        Size of the pixels [arcmin]

    center_ra: float
        RA coordinate of projection origin [degrees]

    center_dec: float
        DEC coordinate of projection origin [degrees]

    projection: string
        Type of 2D projection in ['Gnomonic'] (default:'Gnomonic')

    Returns
    -------
    catalog: table
        Output shape catalog with pixel index column

    grid_ra: 2d array
        Ra coordinates of pixels

    grid_dec: 2d array
        Dec coordinates of pixels
    """
    # Convert pixel_size to deg for consistency
    pixel_size = pixel_size / 60.

    # Create coordinate grid for the map
    edges_x = np.linspace(-nx/2*pixel_size, nx/2*pixel_size, nx + 1)
    edges_y = np.linspace(-ny/2*pixel_size, ny/2*pixel_size, ny + 1)

    # coordinate grid for center of pixels
    grid_x, grid_y = np.meshgrid(0.5*(edges_x[1:] +edges_x[:-1]),
                                 0.5*(edges_y[1:] +edges_y[:-1]))

    # Computes projected coordinates on 2D plane
    if projection == 'gnomonic':
        x,y = radec2xy(center_ra, center_dec, catalog['ra'], catalog['dec'])
        grid_ra, grid_dec = xy2radec(center_ra, center_dec,
                                     grid_x.flatten(), grid_y.flatten())
        grid_ra  = grid_ra.reshape(grid_x.shape)
        grid_dec = grid_dec.reshape(grid_y.shape)
    else:
        raise NotImplementedError

    # Restricts to the galaxies that fall within the patch
    sel = ((x >= edges_x[0]) & (x < edges_x[-1]) &
           (y >= edges_y[0]) & (y < edges_y[-1]))
    catalog = catalog[sel]

    # Use this function mostly for its 2d digitize equivalent
    s, xe, ye ,binnumber = binned_statistic_2d(x[sel], y[sel],
                                               np.ones_like(x[sel]),
                                               statistic='count',
                                               bins=(edges_x, edges_y),
                                               expand_binnumbers=True)
    binnumber -= 1
    pixel_index =  (binnumber[1,:] * nx + binnumber[0, :])
    catalog['pixel_index'] = pixel_index

    return catalog, grid_ra, grid_dec
