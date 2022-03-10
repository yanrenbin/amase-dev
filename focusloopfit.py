import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from astropy.io import ascii

def gauss2d(xdata_tuple,a,x0,y0,sigma,offset):
    (x,y) = xdata_tuple
    xo = float(x0)
    yo = float(y0)
    dist2 = (x-x0)**2+(y-y0)**2
    z = a*np.exp(-dist2/(2*sigma*sigma))+offset
    return z.ravel()

def gauss2d_zero(xdata_tuple,a,x0,y0,sigma):
    (x,y) = xdata_tuple
    xo = float(x0)
    yo = float(y0)
    dist2 = (x-x0)**2+(y-y0)**2
    z = a*np.exp(-dist2/(2*sigma*sigma))
    return z.ravel()

dir = '/Users/renbin/astro/amase/labdata/qhy163m/fiber50/2020-08-27/'
prefix = '11_37_52_fiber50_'

table = ascii.read(dir+'focusloop.txt')
table = table[4:12]
colla = table['Colla']

nframe=5
gain=2.37
data = np.zeros([len(table),120,160],dtype='float')
fig,axs = plt.subplots(len(table),5)
for i in range(len(table)):
    fileno = np.arange(nframe)+table['Framestart'][i]
    filegrp = np.array([dir+prefix+'00'+str(a)+'.fits' for a in fileno])
    for thisfile in filegrp:
        with fits.open(thisfile) as hdu:
            data[i]=data[i]+hdu[0].data
    data[i] = data[i]/16.*gain
    maxind = np.argmax(data[i])
    sz = data[i].shape
    xpos = maxind%sz[1]
    ypos = int(maxind/sz[1])

    hwidth=10

    x0 = np.arange(sz[1])
    y0 = np.arange(sz[0])
    xg,yg = np.meshgrid(x0,y0)

    other=np.where((xg < xpos-hwidth) | (xg > xpos+hwidth) | (yg < ypos-hwidth) | (yg > ypos+hwidth))

    background = np.median(data[i][other])

    subtracted = data[i]-background
    noise = np.sqrt(data[i]+2.2**2*3)
    cutout = subtracted[ypos-hwidth:ypos+hwidth+1,xpos-hwidth:xpos+hwidth+1]
    noisecutout = noise[ypos-hwidth:ypos+hwidth+1,xpos-hwidth:xpos+hwidth+1]

    neff = np.sum(cutout)**2/np.sum(cutout**2)
    x0 = hwidth
    y0 = hwidth
    sigma = np.sqrt(neff/4/math.pi)
    amp = np.max(cutout)

    xc = np.arange(hwidth*2+1)
    yc = np.arange(hwidth*2+1)
    xx,yy = np.meshgrid(xc,yc)
#    guess = (amp,x0,y0,sigma,0.0)
    guess=(amp,x0,y0,sigma)

#    popt,pcov = curve_fit(gauss2d, (xx,yy),cutout.ravel(),p0=guess)
    popt,pcov = curve_fit(gauss2d_zero, (xx,yy), cutout.ravel(), sigma=noisecutout.ravel(), p0=guess)
    yfit = gauss2d_zero((xx,yy),*popt)
    yfit = yfit.reshape(2*hwidth+1,2*hwidth+1)

    print(f'{i:2d}, {colla[i]:5.2f}, {np.max(cutout-yfit):9.1f},{np.min(cutout-yfit):9.1f},{np.sum(cutout):10.1f},{background:9.1f},{popt[0]:9.1f},{popt[3]:6.3f}')
    axs[i,0].plot(cutout[hwidth,:])
    axs[i,0].plot(yfit[hwidth,:])
    axs[i,0].plot(cutout[hwidth,:]-yfit[hwidth,:])
    axs[i,0].plot(xc,np.zeros(2*hwidth+1),'--')
    axs[i,0].set_xlabel('X')

    axs[i,1].plot(cutout[:,hwidth])
    axs[i,1].plot(yfit[:,hwidth])
    axs[i,1].plot(cutout[:,hwidth]-yfit[:,hwidth])
    axs[i,1].plot(yc,np.zeros(2*hwidth+1),'--')
    axs[i,1].set_xlabel('Y')

    pos2=axs[i,2].imshow(cutout-yfit,vmin=-1800,vmax=2700)
    fig.colorbar(pos2,ax=axs[i,2])
    pos3=axs[i,3].imshow(cutout,vmax=22000,vmin=0)
    fig.colorbar(pos3,ax=axs[i,3])
    pos4=axs[i,4].imshow(yfit,vmax=22000,vmin=0)
    fig.colorbar(pos4,ax=axs[i,4])

