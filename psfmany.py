import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import sys
from psffit_single import psffit_single, doublegauss
from astropy.io import ascii


dir = '/Users/renbin/astro/amase/labdata/qhy163m/fiber50/2020-08-28/'
prefix = '12_04_18_fiber50_'

table =ascii.read(dir+'../focusloop.txt')

panarr = np.array([0,1048,2248,3448, 4496])
hwhmarr = np.zeros(len(table),dtype='float')
optarr = np.zeros([len(table),4],dtype='float')
chi2arr = np.zeros(len(table),dtype='float')

for i in range(5):
     indx = np.where(table['pan'] == panarr[i])[0]
     for j in indx:
         popt,chi,hwhm = psffit_single(table[j]['start'],table[j]['end'],dir=dir,prefix=prefix)
         optarr[j] = popt
         hwhmarr[j] = hwhm
         chi2arr[j] = chi
     res = np.polyfit(table[indx]['focus'],hwhmarr[indx],2)
     xfocus = (np.max(table[indx]['focus'])-np.min(table[indx]['focus']))*np.arange(100.)/100.+np.min(table[indx]['focus'])

     plt.plot(table[indx]['focus'],hwhmarr[indx],'x')
     plt.plot(xfocus,res[0]*xfocus**2+res[1]*xfocus+res[2])
#     res2 = np.polyfit(table[indx]['focus'],optarr[indx,4],2)
#     plt.plot(table[indx]['focus'],optarr[indx,4],'o')
#     plt.plot(xfocus,res2[0]*xfocus**2+res2[1]*xfocus+res2[2])
     plt.title('pan ='+str(panarr[i]))
     plt.show()



     
   
