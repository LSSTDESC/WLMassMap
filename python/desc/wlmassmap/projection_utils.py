import numpy as np

def eq2ang(ra, dec):
    """
    convert equatorial ra,dec in degrees to angular theta, phi in radians
    parameters

    ----------
    ra: scalar or array
        Right ascension in degrees
    dec: scalar or array
        Declination in degrees
    returns
    -------
    theta,phi: tuple
        theta = pi/2-dec*D2R # in [0,pi]
        phi   = ra*D2R       # in [0,2*pi]
    """
    dec = dec*np.pi/180
    ra  = ra*np.pi/180
    theta = np.pi/2 - dec
    phi  = ra
    return theta, phi

def radec2xy(ra0, dec0, ra, dec, radians=False):
    """
    Gnomonic projection of (ra, dec) coordinates to tangent plane about
    the point (ra0, dec0).

    Parameters
    ----------
    ra0 : float
        RA coordinate of projection origin [degrees]
    dec0 : float
        DEC coordinate of projection origin [degrees]
    ra : array_like
        RA value(s) to project [degrees]
    dec : array_like
        DEC value(s) to project [degrees]
    radians : bool, optional
        If True, return radian values with origin at (0,0) instead of RADEC

    Returns
    -------
    x : array_like
        Projected RA value(s) [degrees]
    y : array_like
        Projected DEC value(s) [degrees]
    """
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    # Convert the input coordinates to radians
    x0 = np.deg2rad(ra0)
    y0 = np.deg2rad(dec0)
    alpha = np.deg2rad(ra)
    delta = np.deg2rad(dec)

    # Compute projected values
    denom = (np.cos(y0) * np.cos(delta) * np.cos(alpha - x0) +
             np.sin(y0) * np.sin(delta))
    x = np.cos(delta) * np.sin(alpha - x0) / denom
    y = ((np.cos(y0) * np.sin(delta) -
          np.sin(y0) * np.cos(delta) * np.cos(alpha - x0)) / denom)

    if radians:
        return x, y

    return np.rad2deg(x), np.rad2deg(y)


def xy2radec(center_ra, center_dec, x, y):
    """
    Inverse projection from tangent plane (x, y) coordinates to the sphere
    about the point (center_ra, center_dec)

    Parameters
    ----------
    center_ra : float
        RA coordinate of projection center [degrees]
    center_dec : float
        DEC coordinate of projection center [degrees]
    x : float array_like
        RA value(s) to deproject [degrees]
    y : float array_like
        DEC value(s) to deproject [degrees]

    Returns
    -------
    ra : float array_like
        Deprojected RA value(s) [degrees]
    dec : float array_like
        Deprojected DEC value(s) [degrees]
    """
    # Convert projection center to radians
    x0 = np.deg2rad(center_ra)
    y0 = np.deg2rad(center_dec)

    x = np.atleast_1d(np.deg2rad(x))
    y = np.atleast_1d(np.deg2rad(y))

    # Compute deprojected coordinates
    z = np.sqrt(x * x + y * y)
    c = np.arctan(z)

    # Prevent division by zero
    factor = np.ones(len(z))
    inds = (z != 0)
    factor[inds] = y[inds] / z[inds]

    delta = np.arcsin(np.cos(c) * np.sin(y0) + factor * np.cos(y0) * np.sin(c))
    denom = z * np.cos(y0) * np.cos(c) - y * np.sin(y0) * np.sin(c)
    alpha = x0 + np.arctan2(x * np.sin(c), denom)

    return np.rad2deg(alpha), np.rad2deg(delta)
