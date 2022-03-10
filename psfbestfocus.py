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


table=ascii.read(dir+'../tilt_shimmed_full.txt')

for i in range(len(table)):
    popt,chi2,fwhm,fig,data=psffit_single(table[i]['start'],table[i]['finish'],dir=dir,prefix=prefix,visual=1)
#    fig.suptitle('pan:'+str(table[i]['pan']+popt[1])+' tilt:'+str(table[i]['tilt']+popt[2])+' focus:'+str(table[i]['micrometer']), fontsize=14)
    dx = (4656/2.-(table[i]['pan']+popt[1]))*3.8/1000.
    dy = (3522/2.-(table[i]['tilt']+popt[2]))*3.8/1000.
    df = table[i]['micrometer']*2.5 
#    fig.suptitle(str("pan: %4i, tilt: %4i, focus: %4.2f" % (table[i]['pan']+popt[1], table[i]['tilt']+popt[2], table[i]['micrometer'])),fontsize=14)
    fig.suptitle("pan: {0:5.3f} mm, tilt: {1:5.3f} mm, focus: {2:6.3f} mm".format(dx,dy,df))
    fig.set_figheight(10)
    fig.set_figwidth(10)
    fig.savefig('bestfocus_pos'+str(i).zfill(2)+'.png')
    plt.close(fig)
#    if i==0: 
#       allfig = fig 
#    else:
#       allfig = np.hstack((allfig,fig))


