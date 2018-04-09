# This module extracts from the simulation source a given footprint will all
# fields required for downstream analysis
from GCRCatalogs import load_catalog
import os
import yaml
from optparse import OptionParser
import astropy.table as table

required_quantities = ['galaxy_id', 'ra', 'dec',
                       'ra_true', 'dec_true',
                       'shear_1', 'shear_2',
                       'convergence', 'redshift']

def extract_footprint(config):
    """
    Extracts data from simulation and exports the master catalog as an HDF5 table

    Parameters
    ----------
        config: dictionary
            Configuration dictionary read from yaml config file
    """

    # Load a catalog using the GCR a
    catalog = load_catalog(config['catalog'])

    filters = []

    # Applies some filters to the input simulation
    if 'footprint' in config:
        if config['footprint']['type'] == 'patch':
            ra_min, ra_max = config['footprint']['ra_range']
            dec_min, dec_max = config['footprint']['dec_range']

            # if simulation has magnification, use apparent (ra,dec)
            # otherwise fallback to true ra,dec
            if catalog.has_quantities(['ra','dec']):
                filters.append((lambda ra: (ra > ra_min) & (ra < ra_max), 'ra'))
                filters.append((lambda dec: (dec > dec_min) & (dec < dec_max), 'dec'))
            else:
                filters.append((lambda ra: (ra > ra_min) & (ra < ra_max), 'ra_true'))
                filters.append((lambda dec: (dec > dec_min) & (dec < dec_max), 'dec_true'))
        else:
            raise NotImplementedError

    # Exports the catalog in an HDF5 file
    filename = config['output_filename']
    quantitites = catalog.get_quantities(required_quantities, filters=filters)
    outdict = table.Table(catalog.get_quantities(required_quantities, filters=filters))
    outdict.write(filename,path="WLMassMap_data",format='hdf5')

if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    extract_footprint(config['extract_footprint'])
