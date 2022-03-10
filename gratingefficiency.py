import os
import math
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

delta_n = 0.2
thickness = 0.562/delta_n # micron  ; 0.565 for red R=16948; 0.427 for blue R=16948; 0.602 for blue R=8500. 0.573 for blue centered at 4700A. ; 0.562 for blue centered at 4610A
wavelength = 0.461 #0.4943 # micron
theta = 35.26*math.pi/180 # 54.7*math.pi/180.
n=1.4563  # 1.4563 for Quartz
spacing = wavelength/(2*n*np.sin(theta)) # 1000/3732. or 1000/4966.
# 3732 l/mm for red, 4966 l/mm for blue.

dwave = np.arange(1001)/10000. # 0A to 1000A in microns.
nu = math.pi*delta_n*thickness/wavelength/np.cos(theta)
zeta = dwave*math.pi*thickness/(2*n*spacing**2*np.cos(theta))
eta_denom = 1+(dwave**2*wavelength**2/4/n**2/spacing**4/delta_n**2) 
eta_numer = np.sin(math.pi*thickness/np.cos(theta)
            *np.sqrt((delta_n/wavelength)**2 + (dwave/2/n/spacing**2)**2)
            /eta_denom)**2

plt.plot(dwave,eta_numer/eta_denom)
plt.show()
