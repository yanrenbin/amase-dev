import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import pdb

#biasdir = '/Users/renbin/astro/amase/labdata/ASICAP/bias/'
biasdir = '/Users/renbin/astro/amase/labdata/qhy163m/bias/'
flatdir = '/Users/renbin/astro/amase/labdata/ASICAP/flat/'

hdu = fits.open(biasdir+'avebias_2019-12-21.fits')
bias = hdu[0].data
hdu.close()
medbias = np.median(bias)

#hdu = fits.open(flatdir+'normflat_ceiling.fits')
#normflat = hdu[0].data
#hdu.close()
normflat = bias*0+1



#dir = '/Users/renbin/astro/amase/labdata/ASICAP/Apr16/'  # Sigma105
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/May25/'  # Nikon85
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/Jun5/'  # Canon85
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/Jun8/' # Fuji 110
dir = '/Users/renbin/astro/amase/labdata/qhy163m/2019-12-20/' # Nikkor58
names = np.empty(268,dtype=(np.unicode_,40))
imgtypes = np.empty(268,dtype=(np.unicode_,16))
lens = np.empty(268,dtype=(np.unicode_,10))
f = open(dir+'list_2019-12-20','r')

#hdu = fits.open(dir+'median_flat_fuji110.fits')
#flatfuji = hdu[0].data
#hdu.close()

#normfuji = np.percentile(flatfuji,98)
#flatfuji = flatfuji/normfuji

#hdu = fits.open(dir+'median_flat2_quartz.fits')
#flatquartz = hdu[0].data
#hdu.close()

#normquartz = np.percentile(flatquartz,98)
#flatquartz = flatquartz/normquartz

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
contonlyid = np.flatnonzero(np.core.defchararray.find(uniqcomb,'wave_')==0)

for type in uniqcomb[contonlyid]:
    ind = np.where(combtypes == type)[0]
    print('For ',type, len(ind))
    totarr = np.zeros(np.append(len(ind),bias.shape))
    ifile = 0
    exptime = np.zeros(len(ind))
    for file in names[ind]:
        hdu=fits.open(dir+file)
#        exptime[ifile] = hdu[0].header['EXPOINUS']
        exptime[ifile] = hdu[0].header['EXPTIME']

        print('  read in file:',file)
        med = np.median(hdu[0].data)
        kk = np.where((hdu[0].data > med-80) & (hdu[0].data < med+80))
        med = np.median((hdu[0].data)[kk])
        print('  median:',med)
        data = hdu[0].data-med - (bias-medbias)
        if med >4000:
            data = data/1.024
    #    if lens[ind[0]] == 'fuji110':
    #       print(' using fuji flat')
    #       totarr[ifile,:,:] = data/flatfuji
    #    if lens[ind[0]] == 'quartz':
    #       print(' using quartz flat')
    #       totarr[ifile,:,:] = data/flatquartz
        totarr[ifile,:,:] = data/normflat
        ifile+=1
        hdu.close()
    meddata = np.median(totarr,axis=0)
    if np.max(exptime)-np.min(exptime) > 0.001:
        print('  Inconsisten Exposure Time')
        pdb.set_trace()
#    spatial = np.sum(meddata,axis=1)
#    spectral = np.sum(meddata,axis=0)

#    linethresh = np.median(spectral)
#    t=np.where(spectral < linethresh)[0]
#    spatial_noline = np.sum(meddata[:,t],axis=1)
    hdu = fits.PrimaryHDU(meddata)
    hdr = hdu.header
#    hdr.set('EXPOINUS',exptime[0])
    hdr.set('EXPTIME',exptime[0])
    hdu.header = hdr
    hdul = fits.HDUList([hdu])
    hdu.writeto(dir+'median_'+type+'.fits',overwrite=True)
