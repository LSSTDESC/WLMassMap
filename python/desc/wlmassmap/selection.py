# This scripts applies a series of cuts to the input catalog
from optparse import OptionParser
import os
import yaml

import numpy as np

from astropy.table import Table
from astropy.io import fits

def selection(config):
    """
    Applies a series of cuts to the input catalog
    """

    filename = config['input_filename']
    shape = Table.read(filename)
    mask = np.ones(len(shape)).astype('bool')

    for cut in config['cuts']:
        mask &= eval(cut)

    shape = shape[mask]

    # Exports the catalog in
    filename = config['output_filename']
    shape.write(filename, overwrite=True)

if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    selection(config['selection'])
