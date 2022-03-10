from astropy.io import ascii
import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import sys
from psffit_single import psffit_single
from astropy.io import ascii


dir = '/Users/renbin/astro/amase/labdata/qhy163m/fiber50/2020-08-28/'
prefix = '12_04_18_fiber50_'


table=ascii.read(dir+'../focusloop.txt')

uu,returnind,rever,cts=np.unique(table['pan'],return_index=True,return_inverse=True,return_counts=True)

n_pan = len(uu)

for i in np.arange(4)+1:
    indx=np.where(table['pan']==uu[i])[0]
    for j in indx:
        popt,chi2,fwhm,fig,data=psffit_single(table[j]['start'],table[j]['end'],dir=dir,prefix=prefix,visual=1)
        dx = (4656/2.-(table[j]['pan']+popt[1]))*3.8/1000.
        dy = (3522/2.-(table[j]['tilt']+popt[2]))*3.8/1000.
        df = table[j]['focus']*2.5
        fig.suptitle("pan: {0:5.3f} mm, tilt: {1:5.3f} mm, focus: {2:6.3f} mm".format(dx,dy,df))
        fig.set_figheight(10)
        fig.set_figwidth(10)
        filename="sweep_pan"+str(table[j]['pan']).zfill(4)+"_foc{0:4.2f}.png".format(table[j]['focus'])
        fig.savefig(filename)
        plt.close(fig)
