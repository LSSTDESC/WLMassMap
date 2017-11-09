# This module handles the projection of a catalog on a specific grid
import numpy as np
from scipy.stats import binned_statistic_2d
from projection_utils import radec2xy, xy2radec, eq2ang
import healpy as hp
import pyssht as ssht

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


def project_ssht_mw(catalog, L_mw, Method='MW'):
    g1 = np.zeros(ssht.sample_shape(L_mw,Method='MW'))
    g2 = np.zeros(ssht.sample_shape(L_mw,Method='MW'))
    shear_g1 = np.zeros(ssht.sample_shape(L_mw,Method='MW'))
    shear_g2 = np.zeros(ssht.sample_shape(L_mw,Method='MW'))
    Ngal = np.zeros(ssht.sample_shape(L_mw,Method='MW'))

    theta,phi=ssht.ra_dec_to_theta_phi(catalog['ra'],catalog['dec'],Degrees=True)
    for i in range(theta.size):
        pix_i = ssht.theta_to_index(theta[i],L_mw,Method='MW')
        pix_j = ssht.phi_to_index(phi[i],L_mw,Method='MW')
        shear_g1[pix_i,pix_j] += catalog['g'][i][0]
        shear_g2[pix_i,pix_j] += catalog['g'][i][1]
        Ngal[pix_i,pix_j] += 1

    #if there are greater than n_min galaxies in the box. weight the shear response.
    #if not, replace the estimate for that pixel to nan.
    for pix_i in range(shear_g1.shape[0]):
        for pix_j in range(shear_g1.shape[1]):
            if Ngal[pix_i][pix_j] != 0.0:
                g1[pix_i][pix_j] = shear_g1[pix_i][pix_j] / Ngal[pix_i][pix_j]
                g2[pix_i][pix_j] = shear_g2[pix_i][pix_j] / Ngal[pix_i][pix_j]

            else:
                g1[pix_i][pix_j] = np.nan
                g2[pix_i][pix_j] = np.nan

    gmap = np.stack([g1,g2], axis=0)

    return gmap, Ngal 

def project_ssht_hp(catalog,Nside):
    
    theta, phi = ssht.ra_dec_to_theta_phi(catalog['ra'],catalog['dec'],Degrees=True)
    Npix = hp.nside2npix(Nside)
    pixnum = hp.ang2pix(Nside,theta,phi)   # for healpix specifically
    e1map_hp = np.zeros(Npix)
    e2map_hp = np.zeros(Npix)
    Ngal_hp = np.zeros(Npix)
    #mask_hp = np.full((Npix,),np.nan)

    for i in range(pixnum.size):
        pix=pixnum[i]
        e1map_hp[pix] += catalog['g'][i][0]
        e2map_hp[pix] += catalog['g'][i][1]
        Ngal_hp[pix] += 1.0
        #mask_hp[pix] = 1.0

    for i in range(Npix):
        if Ngal_hp[i]!=0.0:
            e1map_hp[i] = e1map_hp[i] / Ngal_hp[i]
            e2map_hp[i] = e2map_hp[i] / Ngal_hp[i]

    else:
        e1map_hp[i]=hp.UNSEEN
        e2map_hp[i]=hp.UNSEEN

    e1map_hp[e1map_hp!=hp.UNSEEN] = e1map_hp[e1map_hp!=hp.UNSEEN]*(-1)
    e2map_hp[e2map_hp!=hp.UNSEEN] = e2map_hp[e2map_hp!=hp.UNSEEN]*(-1)

    gmap = np.stack([e1map_hp,e2map_hp], axis=0)

    return gmap, Ngal_hp



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
