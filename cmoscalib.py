import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

dir = '/Users/renbin/astro/mase/GFX50R/'
filenumber = np.arange(79,106)
filename = filenumber.astype('U60')
totarr = np.zeros([len(filename),500,500],dtype='int32')

for i in np.arange(len(filenumber)):
    filename[i] = dir+'DSCF'+filename[i].zfill(4)+'.fit'
    with fits.open(filename[i]) as hdu:
#        image = hdu[0].data
#        totarr[i,:,:] = hdu[0].data[2:502,0:500]
#	totarr[i,:,:] = hdu[0].data[2:6210,0:8280]
        totarr[i,:,:] = hdu[0].data[3000:3500,4000:4500]+32768

medflat = np.median(totarr,axis=0)
stdflat = np.std(totarr,axis=0)

biasnumber = np.arange(106,117) # use just 11 frames here for quick look. There are actually 29 frames (0106-0134)
biasname = biasnumber.astype('U60')

totbias = np.zeros([len(biasnumber),6208,8280],dtype='int32')

for i in np.arange(len(biasnumber)):
    biasname[i] = dir+'DSCF'+biasname[i].zfill(4)+'.fit'
    with fits.open(biasname[i]) as hdu:
#        image=hdu[0].data
        totbias[i,:,:] = hdu[0].data[2:6210,0:8280]
#        totbias[i,:,:] = hdu[0].data[3002:3502,4000:4500]+32768

medbias = np.median(totbias,axis=0)
hdu = fits.PrimaryHDU(medbias)
hdul = fits.HDUList([hdu])
biasdir = dir+'bias/'
hdu.writeto(biasdir+'avebias.fits',overwrite=True)

stdbias = np.std(totbias,axis=0)

gain = (medflat-medbias)/(stdbias**2+stdflat**2)
readnoise = gain*stdbias

ysize,xsize = totarr[0,:,:].shape
xsize = int(xsize)
ysize = int(ysize)
xodd = np.array(range(int(xsize/2)))*2+1
xeven = np.array(range(int(xsize/2)))*2
yodd = np.array(range(int(ysize/2)))*2+1
yeven = np.array(range(int(ysize/2)))*2
flag = np.zeros([ysize,xsize],dtype=int)
# red is 0, green is 1, blue is 2.
for y in yeven:
    flag[y,xeven] = 1
    flag[y,xodd] = 2
for y in yodd:
    flag[y,xodd] = 1
    flag[y,xeven] = 0


