from astropy.table import Table
from astropy.io import fits
from GCR import BaseGalaxyCatalog, register_reader
import os
import numpy as np
from astropy.cosmology import FlatLambdaCDM

__all__ = ['MiraTitanCatalog']

def _get_fits_data(fits_file):
    return fits_file[1].data

class MiraTitanCatalog(BaseGalaxyCatalog):
    """
    Mira Titan catalog reader for the DESCQA generic catalog reader
    """

    def _subclass_init(self,
            catalog_main_dir=os.curdir,
            cosmo_h=0.704,
            cosmo_Omega_M0=0.272,
            filename_template='MT_LSST_tomo{}.fits',
            nbins=10, **kwargs):

        self._quantity_modifiers ={
            'ra': 'ra_arcmin',
            'dec': 'dec_arcmin',
            'ra_true': 'ra_arcmin',
            'dec_true': 'dec_arcmin',
            'redshift': 'z_spec',
            'shear_1': 'shear1',
            'shear_2': 'shear2',
            'convergence': (lambda x: np.zeros_like(x), 'ra_arcmin'),
            'galaxy_id': (lambda x: np.zeros_like(x), 'ra_arcmin'),
        }

        self._nbins = nbins
        self._filename_template = filename_template
        self._catalog_main_dir = catalog_main_dir
        self._tomo_list = list(range(1,self._nbins+1))
        print("WARNING: Mira-Titan Catalog does not provide galaxy ids, magnification, or convergence, these fields will be set to 0")

    def _generate_native_quantity_list(self):
        with fits.open(os.path.join(self._catalog_main_dir, self._filename_template.format(1))) as d:
            native_quantities = d[1].data.dtype.names
        return native_quantities

    def _iter_native_dataset(self, pre_filters=None):
        for i in self._tomo_list:
            if pre_filters and not all(f[0](*([i]*(len(f)-1))) for f in pre_filters):
                continue

            fp = Table.read(os.path.join(self._catalog_main_dir, self._filename_template.format(i)))
            yield i, fp

    @staticmethod
    def _fetch_native_quantity(dataset, native_quantity):
        fid, data = dataset
        return data[native_quantity]


# Register reader
register_reader(MiraTitanCatalog)
