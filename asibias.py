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

#dir = '/Users/renbin/astro/amase/labdata/ASICAP/Apr19/'
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/May25/'

#dir = '/Users/renbin/astro/amase/labdata/ASICAP/Jun16/CapObj/2019-06-16_19_57_26Z/'
dir = '/Users/renbin/astro/amase/labdata/qhy163m/2019-12-21/'

names = np.empty(40,dtype=(np.unicode_,40))
imgtypes = np.empty(40,dtype=(np.unicode_,16))
f = open(dir+'list_bias','r')
ct=0
for line in f:
      line = line.strip()
      columns = line.split()
      if len(columns) == 2: 
          names[ct] = columns[0]
          imgtypes[ct] = columns[1]
          ct = ct+1
f.close()

bias=np.where(imgtypes == 'bias')

ct=0
for file in names[bias]:
    with fits.open(dir+file) as hdu:
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
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/bias/'
dir = '/Users/renbin/astro/amase/labdata/qhy163m/bias/'
hdu.writeto(dir+'avebias_2019-12-21.fits',overwrite=True)
