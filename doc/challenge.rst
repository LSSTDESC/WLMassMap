=================================================
WL mass mapping challenge
=================================================

This document describes our plans for a DESC mass mapping 'challenge' for 
different teams to compare algorithms on a common set of input shear catalogs.

We have a variety of groups within DESC WL that have developed algorithms for 
both 2D and 3D lensing mass reconstruction. We expect that different algorithms
will show varying features and limitations. By coordinating common analyses in a
challenge we aim to identify the features and limitations of each algorithm and
to develop additional pipeline capabilities for future DESC analyses.

We will use simulated mocks as the initial input data set so we may compare with 
known truths for the cosmology and mass structures. In a series of staged 
analyses we would like to increase the complexity of the mock catalogs including
some systematics and noise properties that would be blinded to the teams 
competing in the competition (although the cosmology will not be blinded here).

We will compare the results from each team using a standard set of metrics 
that are listed below. 

Purposes
=========
- Development of new pipelines and shear statistics
- Correlation with systematics
- Identification of screened/unscreened regions
- Identification of voids
- Preparation for application to HSC and LSST data sets

Challenge data sets
===================

Mira Titan Mocks octant at LSST density
---------------------------------------

The mocks, made by Joachim Harnois-Deraps, are available in the project space on NERSC here:
/project/projectdirs/lsst/desc-wl/wl-massmap/mira-titan

Consider variations with added information and complexity:

1. Shape noise and SNR noise added with secret values
2. Redshifts with errors


Metrics
========
M1 - RMS error in a selection of random points (suitably smoothed - healpix pixels)
M2 - Peak/Trough/Void counts (as a function of smoothing scales?)
M3 - RMS error on points that are masked out
M4 - Some kind of peak ellipticity metric ??
M5 - Power spectrum
J1 - "Jet" award for best presentation/colour scheme.


Challenge variations
====================

Challenge 1
-----------
 No mask
 2D reconstruction

Challenge 2
-----------
LSST mask
Tomographic reconstruction (truth photo-z's?)

Challenge 3
-----------
No mask
3D reconstruction

Challenge 4
------------
LSST mask
3D reconstruction
