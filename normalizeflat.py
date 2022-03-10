#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:25:12 2019

@author: renbin
"""

import numpy as np
from astropy.io import fits

dir = '/Users/renbin/ASICAP/'
hdu=fits.open(dir+'ceilingflat_final2.fits')

flat=hdu[0].data

medflat=np.median(flat)
normflat = flat/medflat

hdu.close()

hdu = fits.PrimaryHDU(normflat)
hdul = fits.HDUList([hdu])
dir = '/Users/renbin/ASICAP/flat/'
hdu.writeto(dir+'normflat_ceiling.fits',overwrite=True)

