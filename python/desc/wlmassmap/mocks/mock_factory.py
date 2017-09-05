from GCR import load_catalog
import os
import yaml
from optparse import OptionParser
from .MiraTitanCatalog import MiraTitanCatalog

default_quantities = ['ra_true', 'dec_true', 'shear_1', 'shear_2',
                      'ellipticity_1_obs', 'ellipticity_2_obs',
                      'redshift', 'redshift_photometric']

def create_mock(config):
    """
    Extracts data from simulation and exports

    Parameters
    ----------
        config: dictionary
            Configuration dictionary read from yaml file
    """
    # Load a catalog using the GCR and store whether to use shear
    # or reduced shear
    catalog = load_catalog(config['input']['catalog'])

    filters = []

    # Adds an 'observed shear' column, either shear or reduced shear
    if config['reduced_shear']:
        catalog.add_quantity_modifier('shear_1_obs', (lambda x, y: x/(1+y), 'shear_1', 'convergence'))
        catalog.add_quantity_modifier('shear_2_obs', (lambda x, y: x/(1+y), 'shear_2', 'convergence'))
    else:
        catalog.add_quantity_modifier('shear_1_obs', 'shear_1')
        catalog.add_quantity_modifier('shear_2_obs', 'shear_2')

    # Adds shape noise
    if 'shape_noise' in config:
        if config['shape_noise']['type'] == 'Gaussian':
            sigma = config['shape_noise']['sigma']
            catalog.add_quantity_modifier('ellipticity_1_obs', (lambda x: x + sigma*randn(len(x)), 'shear_1_obs'))
            catalog.add_quantity_modifier('ellipticity_2_obs', (lambda x: x + sigma*randn(len(x)), 'shear_2_obs'))
        else:
            raise NotImplementedError

    if 'photoz' in config:
        if config['photoz']['type'] == 'Gaussian':
            sigma = config['photoz']['sigma']
            catalog.add_quantity_modifier('redshift_photometric', (lambda x: x + sigma * (1+x) *randn(len(x)), 'redshift'))
        else:
            raise NotImplementedError

    if 'mask' in config:
        if config['mask']['type'] == 'patch':
            ra_min, ra_max = config['mask']['ra_range']
            dec_min, dec_max = config['mask']['dec_range']
            filters.append((lambda ra: (ra > ra_min) & (ra < ra_max), 'ra_true'))
            filters.append((lambda dec: (dec > dec_min) & (dec < dec_max), 'dec_true'))
        else:
            raise NotImplementedError

    # Exports the catalog in an HDF5 file
    filename = os.path.join(config['output_directory'], config['output_filename'])
    catalog.get_quantities(default_quantities, return_hdf5=filename)


if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    create_mock(config)
