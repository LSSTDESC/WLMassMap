import os
import yaml
from astropy.table import Table
from optparse import OptionParser
import h5py
import numpy as np
from numpy.random import randn

def metacal_format(catalog, e1, e2, g1, g2, R=np.diag([1,1]), delta_gamma=0.01):
    """
    Populate the catalog with metacal fields

    Parameters
    ----------
    catalog: Table
        Catalog to populate with format specific fields

    e1, e2: array_like
        Intrinsic ellipticity of the galaxies

    g1, g2: array_like
        Shear or reduced shear at the position of the galaxies

    R: (2,2) array
        Shear responsivity

    delta_gamma: float
        Shearing strength when measuring the responsivity
    """

    catalog['mcal_flags'] = 0

    # Model the measured shear as e = e|g=0 + R . g
    catalog['mcal_g'] = np.stack([e1, e2], axis=1) + np.dot(R, np.stack([g1, g2], axis=1).T).T

    catalog['mcal_g_1p'] = np.stack([e1, e2], axis=1) + np.dot(R, np.stack([g1 + delta_gamma, g2], axis=1).T).T
    catalog['mcal_g_1m'] = np.stack([e1, e2], axis=1) + np.dot(R, np.stack([g1 - delta_gamma, g2], axis=1).T).T

    catalog['mcal_g_2p'] = np.stack([e1, e2], axis=1) + np.dot(R, np.stack([g1, g2 + delta_gamma], axis=1).T).T
    catalog['mcal_g_2m'] = np.stack([e1, e2], axis=1) + np.dot(R, np.stack([g1, g2 - delta_gamma], axis=1).T).T


def mock_observation(config):
    """
    Create a mock shape catalog matching the format of a given shape measurement
    pipeline.

    Parameters
    ----------
        config: dictionary
            Configuration dictionary read from yaml config file
    """

    # Open the input ground_truth catalog
    filename = os.path.join(config['input_directory'], config['input_filename'])
    cat_gt = h5py.File(filename, 'r')

    # Initialize catalog with mandatory fields
    catalog = Table([cat_gt['galaxy_id'][...], cat_gt['ra'][...], cat_gt['dec'][...]],
                    names=['id', 'ra', 'dec'])

    # Extract the shear that will be used to form the mock shape measurement
    if config['reduced_shear']:
        factor = 1. / (1 + cat_gt['convergence'][...])
    else:
        factor = 1.
    g1 = factor * cat_gt['shear_1'][...]
    g2 = factor * cat_gt['shear_2'][...]

    # Computes some intrinsic shapes for the galaxies
    e1 = np.zeros_like(g1)
    e2 = np.zeros_like(g2)
    if 'shape_noise' in config:
        if config['shape_noise']['type'] == 'Gaussian':
            e1 = config['shape_noise']['sigma'] * randn(len(g1))
            e2 = config['shape_noise']['sigma'] * randn(len(g2))
        else
            raise NotImplementedError

    # Adds format specific catalog fields
    if config['format']['type'] == 'metacal':
        metacal_format(catalog, g1, g2, e1, e2, **config['format'])
    else:
        raise NotImplementedError

    # Exports the catalog in an HDF5 file
    filename = os.path.join(config['output_directory'], config['output_filename'])
    catalog.write(filename, overwrite=True)


if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    mock_observation(config['mock_observation'])
