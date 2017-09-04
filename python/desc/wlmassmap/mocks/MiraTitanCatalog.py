from astropy.table import Table
from GCR import BaseGalaxyCatalog, register_reader

__all__ = ['MiraTitanCatalog']

class MiraTitanCatalog(BaseGalaxyCatalog):
    """
    Mira Titan catalog reader for the DESCQA generic catalog reader
    """

    def _subclass_init(self, filename, **kwargs):

        self._quantity_modifiers ={
            'ra_true': 'ra_arcmin',
            'dec_true': 'dec_arcmin',
            'redshift_true': 'z_spec',
            'shear_1': 'shear1',
            'shear_2': 'shear2'
        }

    @staticmethod
    def _fetch_native_quantity(dataset, native_quantity):
        return dataset[native_quantity].value

# Register reader 
register_reader(MiraTitanCatalog)

