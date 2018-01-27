# Examples of WLMassMap usage

This folder contains notebooks and example config files for the different components of the pipeline

## Producing a mock shape catalog

Make sure the python code of the package is in your path, for instance with:
```
$ export PYTHONPATH=$PYTHONPATH:[path to this repo]/python
```

First step, extract a *ground truth* catalog from the desired simulation:
```
$ python -m desc.wlmassmap.mocks.extract_footprint config.yaml
```
this will export a subsample of the original simulation, with only the relevant fields, stored as an HDF5 file.

Then, apply the `mock_observation` module to generate an *observed* catalog, with the same field names as the real shape measurement pipeline, and without the truth information (eg. convergence, redshift,etc..):
```
$ python -m desc.wlmassmap.mocks.mock_observation config.yaml
```

Check the [configuration file](config.yaml) to see the kind of options that can be set.

## Binning a catalog into a shear map and compute a convergence map

The shear_map.py module contains the code to apply the shear calibration,
project and bin the catalog:
```
$ python -m desc.wlmassmap.mocks.shear_map config.yaml
```

The next step is to compute a convergence map:
```
$ python -m desc.wlmassmap.mocks.convergence_map config.yaml
```

Check the [configuration file](config_map.yaml) to see the kind of options that can be set.
