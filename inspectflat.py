#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:07:26 2019

@author: renbin
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

dir = '/Users/renbin/ASICAP/'
hdu1=fits.open(dir+'haloflat1.fits')
hdu2=fits.open(dir+'haloflat2.fits')

haloflat1=hdu1[0].data
haloflat2=hdu2[0].data


