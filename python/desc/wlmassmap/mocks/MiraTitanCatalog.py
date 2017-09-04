from astropy.table import Table
from GCR import BaseGalaxyCatalog, register_reader
import os
from astropy.io import fits
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
            'ra_true': 'ra_arcmin',
            'dec_true': 'dec_arcmin',
            'redshift_true': 'z_spec',
            'shear_1': 'shear1',
            'shear_2': 'shear2'
        }

        self._nbins = nbins
        self._filename_template = filename_template
        self._catalog_main_dir = catalog_main_dir
        self._tomo_list = list(range(self._nbins))


    def _generate_native_quantity_list(self):
        native_quantities = {'original_healpixel'}
        for _, dataset in self._iter_native_dataset():
            for k, v in dataset.items():
                fields = _get_fits_data(v).dtype.fields
                for name, (dt, size) in _get_fits_data(v).dtype.fields.items():
                    if dt.shape:
                        for i in range(dt.shape[0]):
                            native_quantities.add((k, name, i))
                    else:
                        native_quantities.add((k, name))
            break
        return native_quantities

    def _iter_native_dataset(self, pre_filters=None):
        for i in self._tomo_list:
            if pre_filters and not all(f[0](*([i]*(len(f)-1))) for f in pre_filters):
                continue

            fp = fits.open(os.path.join(self._catalog_main_dir, self._filename_template.format(i)))
            yield i, fp
            fp.close()

    @staticmethod
    def _fetch_native_quantity(dataset, native_quantity):
        fid, fits_data = dataset
        data =  _get_fits_data(fits_data[native_quantity[0]])[native_quantity[1]]
        return data


# Register reader
register_reader(MiraTitanCatalog)

