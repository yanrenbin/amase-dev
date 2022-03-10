import numpy as np
from astropy.io import fits
import os
from astropy.table import Table
import matplotlib.pyplot as plt
import pdb
from astropy.convolution import convolve, Box1DKernel

#camera = 'sigma105'
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/Apr16/'
#hdu = fits.open(dir+'throughput.fits')
#res = hdu[1].data

#camera = 'canon85'
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/Jun1/'
#hdu = fits.open(dir+'throughput.fits')
#res = hdu[1].data
#hdu.close()
#dir2 = '/Users/renbin/astro/amase/labdata/ASICAP/Jun2/'
#hdu = fits.open(dir2+'throughput.fits')
#res2 = hdu[1].data
#res = np.concatenate((res[4:14],res2))

#camera = 'nikon85'
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/May25/'
#hdu = fits.open(dir+'throughput.fits')
#res = hdu[1].data
#hdu.close()

#camera = 'fuji110'
#dir = '/Users/renbin/astro/amase/labdata/ASICAP/Jun8/'
#hdu = fits.open(dir+'throughput.fits')
#res = hdu[1].data
#hdu.close()
camera = 'nikkor58'
dir = '/Users/renbin/astro/amase/labdata/qhy163m/2019-12-20/'
hdu = fits.open(dir+'throughput.fits')
res = hdu[1].data
hdu.close()

grpnames = res['grpname'].astype('U')
if camera == 'nikon85': 
    grplist = ['grp1','grp2','grp2f']
if camera == 'canon85':
    grplist = ['grp1','grp2','grp2f','grp3','grp4','grp5','grp6']
if camera == 'sigam105':
    grplist = ['grp1','grp2','grp3','grp4.f','grp5.f']
if camera == 'fuji110':
    grplist = ['grp1','grp2','grp3','grp3f','grp4']
if camera == 'nikkor58':
    grplist = ['grp1','grp2','grp3']

wave=np.array(range(3500,11500)).astype(float)
flux_quartz = np.zeros([len(grplist),8000])
flux_sigma = np.zeros([len(grplist),8000])

for i in range(0,len(grplist)):
    grpname = grplist[i]
    ind_quartz = np.where(grpnames == grpname+'_quartz')[0][0]
    ind_sigma = np.where(grpnames == grpname+'_'+camera)[0][0]
    ss = np.argsort(res['wave'][ind_quartz])
    flux_quartz[i,:] = np.interp(wave, res['wave'][ind_quartz][ss],
        res['flux'][ind_quartz][ss], left=0,right=0)
    ss = np.argsort(res['wave'][ind_sigma])
    flux_sigma[i,:] = np.interp(wave,res['wave'][ind_sigma][ss],
        res['flux'][ind_sigma][ss], left=0,right=0)

totsmratio = np.zeros(8000)
for i in range(0,len(grplist)):
    gd = np.where((flux_quartz[i,:]>0.01*np.max(flux_quartz[i,:])) & (flux_sigma[i,:]!=0))[0]
    ratio = flux_sigma[i,gd]/(flux_quartz[i,gd]+(flux_quartz[i,gd]==0).astype(int))*(0.965**4)
    smratio = convolve(ratio,Box1DKernel(101))
    if i==0: part =np.where((wave[gd] >= 3790) & (wave[gd]<5210))[0]
    if i==1: part =np.where((wave[gd] >= 5210) & (wave[gd]<6080))[0]
    if i==2: part =np.where((wave[gd] >= 6080) & (wave[gd]<7320))[0]
    totsmratio[gd[part]]=smratio[part]

plt.plot(wave,totsmratio,'-',lw=2)
plt.xlim([3790,7300])
plt.grid(True)
#plt.plot([3500,8000],[0.9,0.9],':')
#plt.plot([3500,8000],[0.5,0.5],':')
plt.xlabel('Wavelength (A)')
plt.ylabel('Transmission')
if camera == 'canon85':
    plt.title('Canon 85mm f/1.2L USM II')
if camera == 'sigma105':
    plt.title('Sigma 105mm f/1.4 ART')
if camera == 'nikon85':
    plt.title('Nikon 85mm f/1.4 G')
if camera == 'fuji110':
    plt.title('Fuji 110mm f/2')
if camera == 'nikkor58':
    plt.title('Nikkor 58mm f/0.95 S Noct')
plt.ylim([0,1])
plt.show()

