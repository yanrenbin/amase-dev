#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 14:51:58 2019

@author: renbin
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

dir = '/Users/renbin/ASICAP/'
#names = np.empty(40,dtype=(np.unicode_,40))
#imgtypes = np.empty(40,dtype=(np.unicode_,16))
#f = open(dir+'listarcs','r')
#ct=0
#for line in f:
#      line = line.strip()
#      columns = line.split()
#      if len(columns) == 2: 
#          names[ct] = columns[0]
#          imgtypes[ct] = columns[1]
#          ct = ct+1
#f.close()

filelist =['arc1','arc2','arc3','arc4','arc5']
for i in range(0,len(filelist)): filelist[i] = filelist[i]+'med.fits'

for file in filelist:
    with fits.open(dir+file) as hdu:
        data = hdu[0].data
        spatial= np.sum(data,axis=1)
        rows = np.where(spatial > np.max(spatial)/2.)[0]
        spec = np.sum(data[rows[30:50],:],axis=0)

#        if ct ==0 :
#            tot=data
#        else:
#            tot = np.dstack((tot,data))
#        ct+=1
#        hdu.close()
#    meddata = np.median(tot,axis=2)
#    
#    hdu=fits.open(dir+file)
#    hdu[0].data = meddata
#    hdu.writeto(dir+flavor+'med.fits',overwrite=True)
#    hdu.close()