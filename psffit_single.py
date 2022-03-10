import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import sys
import pdb

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

def doublegauss(xdata_tuple,a0,a1,x0,y0,sigma0,sigma1):
    (x,y) = xdata_tuple
    xo = float(x0)
    yo = float(y0)
    dist2 = (x-x0)**2+(y-y0)**2
    z = a0*np.exp(-dist2/(2*sigma0*sigma0))+a1*np.exp(-dist2/(2*sigma1*sigma1))
    return z.ravel()

def gauss2d_asym(xdata_tuple,a0,x0,y0,sigmax,sigmay):
    (x,y) = xdata_tuple
    xo = float(x0)
    yo = float(y0)
    z = a0*np.exp(-(x-x0)**2/(2*sigmax*sigmax)-(y-y0)**2/(2*sigmay*sigmay))
    return z.ravel()


def psffit_single(startframe, endframe, dir, prefix, visual=0):
    #if len(sys.argv)!=3:
    #   sys.exit(2)
    #print(sys.argv[1])
    #print(sys.argv[2])

    #dir = '/Users/renbin/astro/amase/labdata/qhy163m/fiber50/2020-08-28/'
    ##prefix = '10_04_15_fiber50_'
    #prefix = '12_04_18_fiber50_'

    #startframe = sys.argv[1]
    #endframe = sys.argv[2]
    nframe = int(endframe)-int(startframe)+1
    fileno = np.arange(nframe)+int(startframe)
    files = np.array([dir+prefix+str(a).zfill(5)+'.fits' for a in fileno])

    gain=2.37

    data = np.zeros([120,160],dtype='float')
    for j in range(nframe):
        with fits.open(files[j]) as hdu:
             data = data+hdu[0].data
    data=data/16.*gain

#    coeff=np.array([4.74076695e8,6.38476617e10])
#    data=(data**(1/0.35069)-coeff[1])/coeff[0]
    maxind = np.argmax(data)
    sz = data.shape
    xpos = maxind%sz[1]
    ypos = int(maxind/sz[1])

    hwidth=10

    x0 = np.arange(sz[1])
    y0 = np.arange(sz[0])
    xg,yg = np.meshgrid(x0,y0)

    other=np.where((xg < xpos-hwidth) | (xg > xpos+hwidth) | (yg < ypos-hwidth) | (yg > ypos+hwidth))

    background = np.median(data[other])

    subtracted = data-background
    noise = np.sqrt(data+2.2**2*3)
    cutout = subtracted[ypos-hwidth:ypos+hwidth+1,xpos-hwidth:xpos+hwidth+1]
    noisecutout = noise[ypos-hwidth:ypos+hwidth+1,xpos-hwidth:xpos+hwidth+1]

    neff = np.sum(cutout)**2/np.sum(cutout**2)
    print(neff)
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
    dof = cutout.size-len(guess)
    chi2 = np.sum((cutout-yfit)**2/(noisecutout**2))/dof
    hwhm = popt[3]*2.35/2
    popt[1]=popt[1]+xpos
    popt[2]=popt[2]+ypos
    x_up = np.interp(amp/2.,cutout[hwidth,int(hwidth-popt[3]*2.5):hwidth],xc[int(hwidth-popt[3]*2.5):hwidth])
#    pdb.set_trace()
    x_down = np.interp(amp/2.,np.flip(cutout[hwidth,hwidth:int(hwidth+popt[3]*2.5+1)]),np.flip(xc[hwidth:int(hwidth+popt[3]*2.5+1)]))
    y_up = np.interp(amp/2.,cutout[int(hwidth-popt[3]*2.5):hwidth,hwidth],yc[int(hwidth-popt[3]*2.5):hwidth])
    y_down = np.interp(amp/2.,np.flip(cutout[hwidth:int(hwidth+popt[3]*2.5+1),hwidth]),np.flip(yc[hwidth:int(hwidth+popt[3]*2.5+1)]))

    xfwhm = x_down-x_up
    yfwhm = y_down-y_up
    fwhm = (xfwhm+yfwhm)/2.

    guess2 = (amp, 0.05*amp,x0,y0,sigma,3*sigma)
    popt2,pcov2 = curve_fit(doublegauss, (xx,yy), cutout.ravel(),sigma=noisecutout.ravel(),p0=guess2)
    yfit2 = doublegauss((xx,yy),*popt2)
    yfit2 = yfit2.reshape(2*hwidth+1,2*hwidth+1)
    dof2 = cutout.size-len(guess2)
    doublechi2 = np.sum((cutout-yfit2)**2/(noisecutout**2))/dof2

    xtemp = np.arange(300)/100.*popt2[4]
    ytemp = popt2[0]*np.exp(-xtemp*xtemp/(2*popt2[4]*popt2[4]))+popt2[1]*np.exp(-xtemp*xtemp/(2*popt2[5]*popt2[5]))
    hwhm2 = np.interp((popt2[0]+popt2[1])/2.,ytemp[::-1],xtemp[::-1])


    guess3 = (amp, x0,y0,sigma, sigma)
    popt3,pcov3 = curve_fit(gauss2d_asym, (xx,yy), cutout.ravel(),sigma=noisecutout.ravel(),p0=guess3) 
    yfit3 = gauss2d_asym((xx,yy),*popt3)
    yfit3 = yfit3.reshape(2*hwidth+1,2*hwidth+1)
    dof3 = cutout.size-len(guess3)
    asymchi2 = np.sum((cutout-yfit3)**2/(noisecutout**2))/dof3

    if visual==1:
        fig,axs = plt.subplots(2,2)
        l1, =axs[0,0].plot(cutout[hwidth,:])
        l2, =axs[0,0].plot(yfit[hwidth,:])
        l3, =axs[0,0].plot(cutout[hwidth,:]-yfit[hwidth,:])
        axs[0,0].plot(xc,np.zeros(2*hwidth+1),'--')
        axs[0,0].set_xlabel('X')
        axs[0,0].set_title('Center row')
        axs[0,0].legend((l1,l2,l3),('data','model','residual'),loc='upper right')
        axs[0,1].plot(cutout[:,hwidth])
        axs[0,1].plot(yfit2[:,hwidth])
        axs[0,1].plot(cutout[:,hwidth]-yfit2[:,hwidth])
        axs[0,1].plot(yc,np.zeros(2*hwidth+1),'--')
        axs[0,1].set_xlabel('Y')
        axs[0,1].set_title('Center column')
        axs[0,1].annotate(str("G-FWHM: %4.2f pix" % (2*hwhm)),xy=(0.55,0.9),xycoords='axes fraction')
        axs[0,1].annotate(str("FWHM: %4.2f pix" % (fwhm)),xy=(0.55,0.8),xycoords='axes fraction')

        pos3=axs[1,0].imshow(cutout/popt[0],vmax=1,vmin=0)
        fig.colorbar(pos3,ax=axs[1,0])
        axs[1,0].set_xlabel('X')
        axs[1,0].set_ylabel('Y')
        axs[1,0].set_title('Data')
        pos2=axs[1,1].imshow((cutout-yfit)/popt[0],vmin=-0.1,vmax=0.1)
        fig.colorbar(pos2,ax=axs[1,1])
        axs[1,1].set_xlabel('X')
        axs[1,1].set_ylabel('Y')
        axs[1,1].set_title('Residual/Gaussian Amplitude')
#        pos4=axs[4].imshow(yfit2,vmax=22000,vmin=0)
#        fig.colorbar(pos4,ax=axs[4])
#        fig.show()
        rawdata=data
        return popt,chi2,hwhm,fig,rawdata
    else: return popt,chi2,hwhm


#print(np.max(cutout-yfit2),np.min(cutout-yfit2),np.sum(cutout))
#print(background, *popt)
#fig,axs = plt.subplots(1,5)
#axs[0].plot(cutout[hwidth,:])
#axs[0].plot(yfit2[hwidth,:])
#axs[0].plot(cutout[hwidth,:]-yfit2[hwidth,:])
#axs[0].plot(xc,np.zeros(2*hwidth+1),'--')
#axs[0].set_xlabel('X')
#
#axs[1].plot(cutout[:,hwidth])
#axs[1].plot(yfit2[:,hwidth])
#axs[1].plot(cutout[:,hwidth]-yfit2[:,hwidth])
#axs[1].plot(yc,np.zeros(2*hwidth+1),'--')
#axs[1].set_xlabel('Y')
#
#pos2=axs[2].imshow(cutout-yfit2,vmin=-1800,vmax=2700)
#fig.colorbar(pos2,ax=axs[2])
#pos3=axs[3].imshow(cutout,vmax=22000,vmin=0)
#fig.colorbar(pos3,ax=axs[3])
#pos4=axs[4].imshow(yfit2,vmax=22000,vmin=0)
#fig.colorbar(pos4,ax=axs[4])

