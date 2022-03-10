import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

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

dir = '/Users/renbin/astro/amase/labdata/qhy163m/fiber50/'
prefix = '8_23_2020T17_23_56_fiber50_000'

nframe=5
fileno1 = np.arange(nframe)+38 # 6meter
fileno2 = np.arange(nframe)+33 # 5meter
fileno3 = np.arange(nframe)+28 # 4meter
fileno4 = np.arange(nframe)+43 # 2.7meter
filegrp1 = np.array([dir+prefix+str(a)+'.fits' for a in fileno1])
filegrp2 = np.array([dir+prefix+str(a)+'.fits' for a in fileno2])
filegrp3 = np.array([dir+prefix+str(a)+'.fits' for a in fileno3])
filegrp4 = np.array([dir+prefix+str(a)+'.fits' for a in fileno4])
files = np.vstack([filegrp1,filegrp2,filegrp3,filegrp4])

gain=2.37
fshape = files.shape

fig,axs = plt.subplots(fshape[0],5)

data = np.zeros([fshape[0],120,160],dtype='float')
for i in range(fshape[0]):
    for j in range(nframe):
        with fits.open(files[i,j]) as hdu:
            data[i] = data[i]+hdu[0].data
    data[i]=data[i]/16.*gain
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

    print(i, np.max(cutout-yfit),np.min(cutout-yfit),np.sum(cutout))
    print(i, background, *popt)
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

