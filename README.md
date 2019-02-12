# The WLMassMap package

See the [Issues](https://github.com/LSSTDESC/WLMassMap/issues) and [Milestones](https://github.com/LSSTDESC/WLMassMap/milestones) pages for further details on project activities.

## Code set-up and testing

Checking out the code:
```
 $ git clone git@github.com:LSSTDESC/WLMassMap.git
```

Or, install code and dependencies with pip:
```
$ pip install --user git+git://github.com/LSSTDESC/WLMassMap.git
```

See [examples](https://github.com/LSSTDESC/WLMassMap/tree/master/examples) for instructions on how to use the code


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
