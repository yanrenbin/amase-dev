import os
import math
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

#Fused silica refractive index
n_4900=1.4629
n_6600=1.4563
bluewave = 500.0*u.nm
redwave = 670*u.nm
resolution = 15000
littrow = 50.69*u.deg #incidence angle in air for a plane-parallel grating
fcoll = 200*u.mm
dfiber = 0.05*u.mm
n_air = 1.00027717

Rbase = 2*fcoll/dfiber*np.tan(littrow)

n1 = n_4900

beta0= littrow
beta1= np.arcsin(np.sin(beta0)/n1)  # incidence angle in the substrate with refrative index of n1
gamma_arr = np.arange(95)/10.*u.deg # array of prism angles

betaprime = np.arcsin(n1*np.sin(beta1+gamma_arr))

magnification = n1*np.tan(beta1)/np.tan(beta0)*np.cos(beta1+gamma_arr)/np.cos(betaprime)

angle=np.interp(resolution/Rbase, magnification, gamma_arr)
bestbp = np.interp(resolution/Rbase, magnification,betaprime).to(u.deg)

lineden_blue = 2*np.sin(beta0)/bluewave
angulardispersion = (np.cos(beta1+angle)/np.cos(bestbp)/np.cos(beta1)*u.rad).to(u.deg)*lineden_blue

print(" Wavelength %4i A:" % (bluewave.value*10))
print("    Prism angle:",angle)
print("    alphaprime:",(beta1+angle).to(u.deg))
print("    betaprime: ",bestbp)
print("    overall grating tilt: ",bestbp-angle)
print("    angular dispersion: ",angulardispersion)
print("    Resolution for 7.4 deg: ",np.interp(7.4*u.deg,gamma_arr,magnification)*Rbase)

n1 = n_6600

beta0= littrow
beta1= np.arcsin(np.sin(beta0)/n1)
gamma_arr = np.arange(95)/10.*u.deg

betaprime = np.arcsin(n1*np.sin(beta1+gamma_arr))

magnification = n1*np.tan(beta1)/np.tan(beta0)*np.cos(beta1+gamma_arr)/np.cos(betaprime)

angle=np.interp(resolution/Rbase, magnification, gamma_arr)
bestbp = np.interp(resolution/Rbase, magnification,betaprime).to(u.deg)

lineden_red = 2*np.sin(beta0)/redwave
angulardispersion = (np.cos(beta1+angle)/np.cos(bestbp)/np.cos(beta1)*u.rad).to(u.deg)*lineden_red

print(' Wavelength %4i A:' % (redwave.value*10))
print("    Prism angle:",angle)
print("    alphaprime:",(beta1+angle).to(u.deg))
print("    betaprime: ",bestbp)
print("    overall grating tilt: ",bestbp-angle)
print("    angular dispersion: ",angulardispersion)
print("    Resolution for 7.4 deg: ",np.interp(7.4*u.deg,gamma_arr,magnification)*Rbase)
fixangle = 7.0*u.deg
bp = np.interp(fixangle, gamma_arr, betaprime)
print("    Angular dispersion for for 7.0 deg: ",(np.cos(beta1+fixangle)/np.cos(bp)/np.cos(beta1)*u.rad).to(u.deg)*lineden_red)

