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
f = open(dir+'listtoday','r')
ct=0
for line in f:
      line = line.strip()
      columns = line.split()
      names[ct] = columns[0]
      imgtypes[ct] = columns[1]
      ct = ct+1
f.close()

bias=np.where(imgtypes == 'bias')
flat1=np.where(imgtypes == 'flat1')
flat3=np.where(imgtypes == 'flat3')

ct=0
for file in names[bias]:
    hdu=fits.open(dir+file)
    data = hdu[0].data
    if ct==0 :
        tot = data
    else:
        tot = np.dstack((tot,data))
    ct = ct+1
    hdu.close()    
medbias = np.median(tot,axis=2)

hdu=fits.open(dir+file)
hdu[0].data = medbias
hdu.writeto(dir+'avebias.fits',overwrite=True)
hdu.close()

ct=0
for file in names[flat1]:
    hdu=fits.open(dir+file)
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
hdu.writeto(dir+'aveflat1.fits',overwrite=True)
hdu.close()

ct=0
for file in names[flat3]:
    hdu=fits.open(dir+file)
    data = hdu[0].data-medbias
    if ct ==0 :
        tot=data
    else:
        tot=np.dstack((tot,data))
    ct+=1
    hdu.close()
medflat3 = np.median(tot,axis=2)
hdu=fits.open(dir+file)
hdu[0].data = medflat3
hdu.writeto(dir+'aveflat3.fits',overwrite=True)
hdu.close()