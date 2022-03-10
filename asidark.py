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

#dir = '/Users/renbin/ASICAP/Apr19/'
#dir = '/Users/renbin/ASICAP/May25/'

dir = '/Users/renbin/ASICAP/Jun16/CapObj/'

biasfile = '/Users/renbin/ASICAP/bias/avebias_1600cool.fits'

names = np.empty(40,dtype=(np.unicode_,60))
imgtypes = np.empty(40,dtype=(np.unicode_,16))
f = open(dir+'list_dark300','r')
ct=0
for line in f:
      line = line.strip()
      columns = line.split()
      if len(columns) == 2: 
          names[ct] = columns[0]
          imgtypes[ct] = columns[1]
          ct = ct+1
f.close()

dark=np.where(imgtypes == 'dark')

ct=0
for file in names[dark]:
    with fits.open(dir+file) as hdu:
        data = hdu[0].data
        if ct==0 :
            tot = data
        else:
            tot = np.dstack((tot,data))
        print("Finished reading and stacking "+file)
        ct = ct+1
meddark = np.median(tot,axis=2)

biashdu = fits.open(biasfile)
bias = biashdu[0].data

meddark = meddark-bias

hdu = fits.PrimaryHDU(meddark)
hdul = fits.HDUList([hdu])
dir = '/Users/renbin/ASICAP/'
hdu.writeto(dir+'dark300_1600cool.fits',overwrite=True)
