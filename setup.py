#!/usr/bin/env python
"""
Weak lensing mass-mapping pipeline for LSST DESC
Copyright (c) 2017 LSST DESC
http://opensource.org/licenses/MIT
"""
from setuptools import setup

setup(
    name='wlmassmap',
    version='0.0.1',
    description='Weak lensing mass-mapping pipeline for LSST DESC',
    url='https://github.com/LSSTDESC/WLMassMap',
    maintainer='Francois Lanusse',
    maintainer_email='francois.lanusse@gmail.com',
    license='MIT',
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ],
    package_dir={'': 'python'},
    packages=['desc.wlmassmap', 'desc.wlmassmap.mocks'],
    install_requires=['numpy','h5py','astropy','GCR']
)
