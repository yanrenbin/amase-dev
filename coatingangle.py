import os
import math
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import ascii

#Fused silica refractive index
idxtab = ascii.read('/Users/renbin/astro/amase/glass/RefractiveIndexINFO.csv')
#n_4900=1.4629
#n_6600=1.4563
bluewave = 500*u.nm #490.3*u.nm
redwave = 660.4*u.nm
resolution = 15000
littrow = 50.68*u.deg #incidence angle in air for a plane-parallel grating
fcoll = 200*u.mm
dfiber = 0.05*u.mm

Rbase = 2*fcoll/dfiber*np.tan(littrow)

n1 = np.interp(bluewave.value/1000.,idxtab['wave'],idxtab['n'])

beta0= littrow
beta1= np.arcsin(np.sin(beta0)/n1)  # incidence angle in the substrate with refrative index of n1
#gamma_arr = np.arange(95)/10.*u.deg # array of prism angles
gamma = 7.4*u.deg

betaprime = np.arcsin(n1*np.sin(beta1+gamma))

magnification = n1*np.tan(beta1)/np.tan(beta0)*np.cos(beta1+gamma)/np.cos(betaprime)

#angle=np.interp(resolution/Rbase, magnification, gamma_arr)
#bestbp = np.interp(resolution/Rbase, magnification,betaprime).to(u.deg)

lineden_blue = 2*np.sin(beta0)/bluewave
angulardispersion = (np.cos(beta1+gamma)/np.cos(betaprime)/np.cos(beta1)*u.rad).to(u.deg)*lineden_blue

print(" Wavelength 4903A:")
print("    Prism angle:",gamma)
print("    alphaprime:",(beta1+gamma).to(u.deg))
print("    betaprime: ",betaprime)
print("    overall grating tilt: ",betaprime-gamma)
print("    angular dispersion: ",angulardispersion)
print("    Resolution for 7.4 deg: ",magnification*Rbase)

wavearr = (np.arange(11)*5+470)*u.nm
refidxarr = np.interp(wavearr/u.nm/1000.,idxtab['wave'],idxtab['n'])
beta1arr = np.arcsin(lineden_blue*wavearr/refidxarr-np.sin(beta1))
alphaprimearr = beta1arr+gamma
betaprimearr = np.arcsin(n1*np.sin(alphaprimearr))

plt.plot(wavearr,betaprimearr.to(u.deg))

bluefile = open('blue_coating.tex', 'w')
for i in range(len(wavearr)):
    bluefile.write("{0:d} & {1:6.2f} & {2:6.2f}\\\\\n".format(int(wavearr[i].value),alphaprimearr[i].to(u.deg).value,betaprimearr[i].to(u.deg).value))
bluefile.close()

n1 = np.interp(redwave.value/1000.,idxtab['wave'],idxtab['n'])

beta0= littrow
beta1= np.arcsin(np.sin(beta0)/n1)

betaprime = np.arcsin(n1*np.sin(beta1+gamma))

magnification = n1*np.tan(beta1)/np.tan(beta0)*np.cos(beta1+gamma)/np.cos(betaprime)

#angle=np.interp(resolution/Rbase, magnification, gamma_arr)
#bestbp = np.interp(resolution/Rbase, magnification,betaprime).to(u.deg)

lineden_red = 2*np.sin(beta0)/redwave
angulardispersion = (np.cos(beta1+gamma)/np.cos(betaprime)/np.cos(beta1)*u.rad).to(u.deg)*lineden_red

print(" Wavelength 6600A:")
print("    Prism angle:",gamma)
print("    alphaprime:",(beta1+gamma).to(u.deg))
print("    betaprime: ",betaprime)
print("    overall grating tilt: ",betaprime-gamma)
print("    angular dispersion: ",angulardispersion)
print("    Resolution for 7.4 deg: ",magnification*Rbase)

wavearr = (np.arange(15)*5+620)*u.nm
refidxarr = np.interp(wavearr/u.nm/1000.,idxtab['wave'],idxtab['n'])
beta1arr = np.arcsin(lineden_red*wavearr/refidxarr-np.sin(beta1))
alphaprimearr = beta1arr+gamma
betaprimearr = np.arcsin(refidxarr*np.sin(alphaprimearr))
plt.plot(wavearr,betaprimearr.to(u.deg))

redfile = open('red_coating.tex', 'w')
for i in range(len(wavearr)):
#    print(f'{wavearr[i]/u.nm:3d},{alphaprimearr[i].to(u.deg)/u.deg:6.2f},{betaprimearr[i].to(u.deg)/u.deg:6.2f}')
    redfile.write('{0:d} & {1:6.2f} & {2:6.2f}\\\\ \n'.format(int(wavearr[i].value),alphaprimearr[i].to(u.deg).value,betaprimearr[i].to(u.deg).value))
redfile.close()
