from astropy.table import Table
import healpy as hp
import numpy as np


class subfield():
    """
    Class that represents one contiguous subfield of the survey
    """

    def __init__(self, subfield_name, randoms_filename=None, dset=None, nside=4096):
        """
        Initialise all the subfield information from the HSC randoms

        Note: can be improved to use the output of the stack when available
        """
        self.name = subfield_name

        if randoms_filename is not None:
            self.nside = nside
            randoms = Table.read(randoms_filename)
            self._init_from_randoms(randoms, nside=nside)
        elif dset is not None:
            self._init_from_hdf5(dset)

    def _init_from_randoms(self, randoms, nside):
        """
        Initialise the subfield from the randoms

        Parameters
        ----------

        nside: int
            Resolution of the mask
        """
        self.ra_range = [randoms['ra'].min(), randoms['ra'].max()]
        self.dec_range = [randoms['dec'].min(), randoms['dec'].max()]

        self.center_ra = 0.5*(self.ra_range[0] + self.ra_range[1])
        self.center_dec = 0.5*(self.dec_range[0] + self.dec_range[1])

        # Code generating the mask from the randoms
        x = randoms['ra']
        y = randoms['dec']
        radec = np.transpose([x, y])
        theta = np.deg2rad(90.0 - radec[:, 1])
        phi = np.deg2rad(radec[:, 0])
        gal_hppix = hp.ang2pix(nside, theta=theta, phi=phi, nest=False)
        npix = hp.nside2npix(nside)
        self.mask = np.bincount(gal_hppix, weights=None, minlength=npix)

    def _init_from_hdf5(self, dset):
        """
        Initialise the subfield from a saved hdf5 dataset
        """
        self.mask = dset[:]
        self.nside = dset.attrs['nside']
        self.ra_range = [dset.attrs['ra_min'], dset.attrs['ra_max']]
        self.dec_range = [dset.attrs['dec_min'], dset.attrs['dec_max']]
        self.center_ra = dset.attrs['center_ra']
        self.center_dec = dset.attrs['center_dec']

    def export_mask(self, filename):
        """
        Exports the mask as a fits file
        """
        hp.write_map('testmap4096.fits', countmap)

    def write(self, hdf_grp):
        """
        Writes all the data related to the subfield into an HDF5 group,
        """

        # Saves the mask
        dset = hdf_grp.create_dataset(self.name,
                                      self.mask.shape,
                                      dtype='i')
        dset[...] = self.mask
        # Adds some metadata
        dset.attrs['nside'] = self.nside
        dset.attrs['ra_min'] = self.ra_range[0]
        dset.attrs['ra_max'] = self.ra_range[1]
        dset.attrs['dec_min'] = self.dec_range[0]
        dset.attrs['dec_max'] = self.dec_range[1]
        dset.attrs['center_ra'] = self.center_ra
        dset.attrs['center_dec'] = self.center_dec
