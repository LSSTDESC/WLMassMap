# The WLMassMap package

**Goal:** Generate lensing mass map(s) from [Hyper Suprime-Cam (HSC)](http://www.naoj.org/Projects/HSC/) publicly released image data.

## Project activities

1. Use LSST Data Management (DM) to process HSC images to galaxy shape catalogs.
2. Use shape catalogs that we compute, and/or those available elsewhere, to produce projected mass maps from weak lensing shear.
3. Compare and validate mass mapping algorithms in conjunction with lensing mocks from DESC.

While we expect exciting mass map results from this project, a primary aim is to build key components of DESC weak lensing expertise and team that will be needed for LSST initial data analysis. 

See the [Issues](https://github.com/LSSTDESC/WLMassMap/issues) and [Milestones](https://github.com/LSSTDESC/WLMassMap/milestones) pages for further details on project activities.

## Code set-up and testing

Checking out the code and submodules:
```
 $ git clone --recursive git@github.com:EiffL/WLMassMap.git
```

From bash
```
$ source <WLMassMap install directory>/setup/setup.sh
$ nosetests <WLMassMap install directory>
```



## Mocks

See the doc/mocks directory for info on available mocks.

## HSC Data

Erin put an example HSC MEDS file here:

http://www.cosmo.bnl.gov/www/esheldon/hsc-meds-examples/

and it has been copied here on NERSC:

/project/projectdirs/lsst/desc-wl/wl-massmap/hsc/test-file/test-segrad5.fits.fz


to read it, install fitsio if you don't have it already

pip install fitsio
Also install the meds reader. The meds library isn't in pip, so git clone it and install.
https://github.com/esheldon/meds

python setup.py install


## People

* Jim Bosch
* Domique Boutigny
* Alan Heavens
* Rachel Mandelbaum
* Josh Meyers
* Javier Sanchez
* Michael Schneider
* Erin Sheldon
* Anze Slosar
* Joe Zuntz

## Demo


## License, etc.

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/DarkEnergyScienceCollaboration/WLMassMap/issues/new).
