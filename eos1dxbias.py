#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 21:09:15 2019

@author: renbin
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

dir = '/Users/renbin/astro/specdragon/eos1dx/bias/fits/'


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

hdu = fits.PrimaryHDU(medbias)
hdul = fits.HDUList([hdu])
hdu.writeto(dir+'avebias.fits',overwrite=True)
