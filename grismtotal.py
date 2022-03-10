import os
import math
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from astropy.io import ascii

n_4900=1.4629
n_6600=1.4563

rcwafile = '/Users/renbin/astro/amase/grating/vph353_107.8.txt'
rcwafile = '/Users/renbin/astro/amase/grating/vph317_102.4_n133.txt'
rcwafile = '/Users/renbin/astro/amase/grating/vph353_107.8_n140.txt'
rcwa = ascii.read(rcwafile)
eff_s = rcwa['col2']
eff_p = rcwa['col3']

refwave = 660.4*u.nm
wavearr = rcwa['col1']*refwave

line1=plt.plot(wavearr, eff_s,label='TE')
line2=plt.plot(wavearr, eff_p,label='TM')
line3=plt.plot(wavearr, (eff_s+eff_p)/2,label='Ave')
plt.legend()
plt.xlabel('Wavelength (A)')
plt.ylabel('Efficiency')
plt.title('VPH 35.27 deg, dlogn=0.07143, d=107.8')  
plt.show()


prismangle = 7.806*u.deg #5.171*u.deg
littrow = 50.14*u.deg
n1 = n_6600
beta1= np.arcsin(np.sin(littrow)/n1)

betaarr = np.arcsin(2*np.sin(beta1)/refwave*wavearr-np.sin(beta1))
alphaprimearr = betaarr+prismangle
betaprimearr = np.arcsin(n1*np.sin(alphaprimearr))

redangles=np.interp([651, 660., 677]*u.nm,wavearr, betaprimearr.to(u.deg))

plt.plot(wavearr,betaprimearr.to(u.deg))

refwave = 490.3*u.nm
wavearr = rcwa['col1']*refwave

prismangle = 7.741*u.deg #5.126*u.deg
littrow = 50.14*u.deg
n1 = n_4900
beta1= np.arcsin(np.sin(littrow)/n1)

betaarr = np.arcsin(2*np.sin(beta1)/refwave*wavearr-np.sin(beta1))
alphaprimearr = betaarr+prismangle
betaprimearr = np.arcsin(n1*np.sin(alphaprimearr))

blueangles=np.interp([482, 490., 505]*u.nm,wavearr, betaprimearr.to(u.deg))

plt.plot(wavearr,betaprimearr.to(u.deg))
plt.show()

line1=plt.plot(wavearr, eff_s,label='TE')
line2=plt.plot(wavearr, eff_p,label='TM')
line3=plt.plot(wavearr, (eff_s+eff_p)/2,label='Ave')
plt.legend()
plt.xlabel('Wavelength (A)')
plt.ylabel('Efficiency')
plt.title('VPH 35.27 deg, dlogn=0.07143, d=107.8')  
plt.show()

wavearr1 = rcwa['col1']*660.4*u.nm
wavearr2 = rcwa['col1']*490.3*u.nm
wavearr = (np.arange(301)+400.)*u.nm
eff1 = np.interp(wavearr,wavearr1,(eff_s+eff_p)/2.)
eff2 = np.interp(wavearr,wavearr2,(eff_s+eff_p)/2.)
line1=plt.plot(wavearr,(1-eff2)*eff1*0.92,label='Red grating')
line2=plt.plot(wavearr,(1-eff1)*eff2*0.92,label='Blue grating')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Ave Efficiency')
plt.legend()
plt.show()

