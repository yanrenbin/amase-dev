import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt

dir = '/Users/renbin/ASICAP/Apr16/'
names = np.empty(268,dtype=(np.unicode_,40))
imgtypes = np.empty(268,dtype=(np.unicode_,16))
lens = np.empty(268,dtype=(np.unicode_,10))
f = open(dir+'list_apr16','r')
ct=0
for line in f:
      line = line.strip()
      columns = line.split()
      if len(columns) == 3:
          names[ct] = columns[0]
          imgtypes[ct] = columns[1]
          lens[ct] = columns[2]
          ct = ct+1
f.close()

names = names[0:ct]
imgtypes = imgtypes[0:ct]
lens = lens[0:ct]
combtypes = np.array([x+'_'+y for x,y in zip(imgtypes,lens)])

#uniqtypes = np.unique(imgtypes)
#uniqlens = np.unique(lens)
uniqcomb = np.unique(combtypes)

for type in uniqcomb:
    ind = np.where(combtypes == type)[0]
    print('For ',type, len(ind))
    file=names[ind[-1]]
    hdu=fits.open(dir+file)
    exptime = hdu[0].header['EXPOINUS']
    hdu.close()
    hdu=fits.open(dir+'median_'+type+'.fits')
    hdr = hdu[0].header
    hdr.set('EXPOINUS',exptime)
    hdu[0].header = hdr
    hdu.writeto(dir+'median_'+type+'.fits',overwrite=True)
    hdu.close()
    

