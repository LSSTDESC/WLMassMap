# The WLMassMap package

**Goal:** Generate lensing mass map(s) from [Hyper Suprime-Cam (HSC)](http://www.naoj.org/Projects/HSC/) publicly released image data.

## Project activities

1. Use LSST Data Management (DM) to process HSC images to galaxy shape catalogs.
2. Use shape catalogs that we compute, and/or those available elsewhere, to produce projected mass maps from weak lensing shear.
3. Compare and validate mass mapping algorithms in conjunction with lensing mocks from DESC.

While we expect exciting mass map results from this project, a primary aim is to build key components of DESC weak lensing expertise and team that will be needed for LSST initial data analysis. 

See the [Issues](https://github.com/LSSTDESC/WLMassMap/issues) and [Milestones](https://github.com/LSSTDESC/WLMassMap/milestones) pages for further details on project activities.


## Mocks

See the doc/mocks directory for info on available mocks.


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

## Code set-up and testing
From bash
```
$ source <WLMassMap install directory>/setup/setup.sh
$ nosetests <WLMassMap install directory>
```

## Demo


## License, etc.

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/DarkEnergyScienceCollaboration/WLMassMap/issues/new).
