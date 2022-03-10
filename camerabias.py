#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 7, 2020

@author: renbin
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

dir = '/Users/renbin/astro/amase/labdata/qhy600m/bias/'


filelist=glob.glob(dir+"*.fit*")

ct=0
for file in filelist:
    with fits.open(file) as hdu:
        data = hdu[0].data
        if ct==0 :
            tot = data
        else:
            tot = np.dstack((tot,data))
        print("Finished reading and stacking "+file)
        ct = ct+1
medbias = np.median(tot,axis=2)
stdbias = np.std(tot,axis=2)*np.sqrt(ct)/np.sqrt(ct-1)

hdu = fits.PrimaryHDU(medbias)
hdu2 = fits.ImageHDU(stdbias)
hdul = fits.HDUList([hdu,hdu2])
hdul.writeto(dir+'avebias.fits',overwrite=True)
#hdu.writeto(dir+'biasstd.fits',overwrite=True)
