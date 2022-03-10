import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import sys
from astropy.io import ascii
from astropy import units as u


alpha = 60.*u.deg
offset = 9. #mm
fcam = 58. #mm
beta = np.arctan(offset/fcam)/2.*u.rad
gamma = np.arcsin(np.sin(beta)/np.sin(alpha))
alphaprime = np.arcsin(np.sin(alpha)*np.cos(gamma)/np.sqrt(1-np.sin(alpha)**2*np.sin(gamma)**2))

delta = np.arcsin(np.cos(alpha)*np.sin(gamma)/np.sqrt(1-np.sin(alpha)**2*np.sin(gamma)**2))

print("Spectral Tilt:",delta.to(u.deg))
print("Fringe rotation:",gamma.to(u.deg))


