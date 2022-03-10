#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 14:51:58 2019

@author: renbin
"""

import numpy as np
from astropy.io import fits

dir = '/Users/renbin/ASICAP/'
names = np.empty(40,dtype=(np.unicode_,40))
imgtypes = np.empty(40,dtype=(np.unicode_,16))
f = open(dir+'listApr13latenight','r')
ct=0
for line in f:
      line = line.strip()
      columns = line.split()
      if len(columns) == 2: 
          names[ct] = columns[0]
          imgtypes[ct] = columns[1]
          ct = ct+1
f.close()

flat1=np.where(imgtypes == 'ceilingflat')
flat2=np.where(imgtypes == 'haloflat')

biashdu = fits.open(dir+'avebias.fits')
medbias = biashdu[0].data

ct=0
for file in names[flat1]:
    hdu=fits.open(dir+file)
    print('read in file:',file)
    data = hdu[0].data-medbias
    if ct ==0 :
        tot=data
    else:
        tot = np.dstack((tot,data))
    ct+=1
    hdu.close()
medflat1 = np.median(tot,axis=2)

hdu=fits.open(dir+file)
hdu[0].data = medflat1
hdu.writeto(dir+'ceilingflat_final2.fits',overwrite=True)
hdu.close()

ct=0
for file in names[flat2]:
    hdu=fits.open(dir+file)
    print('read in file:',file)
    data = hdu[0].data-medbias
    if ct ==0 :
        tot=data
    else:
        tot=np.dstack((tot,data))
    ct+=1
    hdu.close()
#medflat2 = np.median(tot,axis=2)
#hdu=fits.open(dir+file)
#hdu[0].data = medflat2
#hdu.writeto(dir+'haloflat_final.fits',overwrite=True)
#hdu.close()