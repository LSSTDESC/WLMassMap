# This scripts applies a series of cuts to the input catalog
from optparse import OptionParser
import os
import yaml

import numpy as np
import time
from astropy.table import Table
from astropy.io import fits
import fitsio

def selection(config):
    """
    Applies a series of cuts to the input catalog
    """
    print 'read data...'

    filename = config['input_filename']
    shape = fitsio.read(filename)
    mask = np.ones(len(shape)).astype('bool')

    print 'make cut...'
    end1 = time.time()
    print end1-start

    for cut in config['cuts']:
        mask &= eval(cut)

    shape = shape[mask]

    # Exports the catalog in
    filename = config['output_filename']
    shape.write(filename, overwrite=True)

if __name__ == "__main__":

    start = time.time()

    parser = OptionParser()
    (options, args) = parser.parse_args()

    with open(args[0]) as f:
        config = yaml.load(f.read())

    selection(config['selection'])

    end = time.time()
    print end-start
