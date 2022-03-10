import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import scipy.integrate as int
import math

h = 6.63e-34
cc = 2.9979245e8
k = 1.38e-23
T=3200 #3200
#Bnu = 2*h*nu**3/cc**2 /(np.exp(h*nu/k/T)-1)

lam = np.arange(9900)+100 # in unit of nm, 100nm to 10 micron
Blam = 2*math.pi*h*cc**2/(lam/1.e9)**5/(np.exp(h*cc/k/T/(lam/1.e9))-1)
photopic=ascii.read('photopic.txt')
lumfun = np.interp(lam, photopic['wave'],photopic['photopic'])

lm=683.*int.trapz(Blam*lumfun,x=lam/1.e9)/int.trapz(Blam,x=lam/1.e9) # lm/W.
int.trapz(Blam,x=lam/1.e9)

#The Sunlite series 5W bulb at 3200k produces 25lm. Thus it is equivalent to
#a 0.92W ideal blackbody. 
#The labeling is messy.  According to https://www.any-lamp.com/lumen-to-watt, 20W halogen should give about 250lumens, which would be a 9.2 W ideal blackbody
#The luminosity density (W/A) would be 
power = 9.2
flam = Blam/int.trapz(Blam,x=lam/1.e9)*power/1.e10 #W/Angstrom. (1.e10 is to convert meter to Angstrom)

maskouter = 1.671*2.54 #1.54*2.54 
maskinner = 0.713*2.54
nlamp=4
maskarea = (maskouter**2-maskinner**2)*math.pi/nlamp #cm^2
lamp_solidangle = maskarea/5**2 *1.8 #1.8=1+0.8 where 0.8 is reflection from the mirror below.
projector_efficiency=0.96**6
screen_efficiency = 0.7
fiber_solidangle = math.pi*((1./3600.)/180.*math.pi)**2
fiber_anglefraction = fiber_solidangle/math.pi
instrumentefficiencypeak=0.4
flux = nlamp*flam*lamp_solidangle/(4*math.pi)*screen_efficiency*fiber_anglefraction*projector_efficiency*instrumentefficiencypeak  # Watt/Angstrom
#glassfilter = ascii.read()
photonflux = flux/(h*cc/(lam/1.e9))  #photons/Angstrom

maskouter_k = 11/2./15.*20 # outer radius in cm
maskinner_k = 2.5/2./15.*20 # inner radius in cm
nlamp_k = 6
maskarea_k = (maskouter_k**2-maskinner_k**2)*math.pi/nlamp_k #cm^2
lamp_solidangle_k = maskarea_k/5**2 *1.8 #1.8=1+0.8 where 0.8 is reflection from the mirror below.
projector_efficiency=0.96**6
screen_efficiency = 0.7
fiber_solidangle_k = math.pi*((0.4/3600.)/180.*math.pi)**2
fiber_anglefraction_k = fiber_solidangle_k/math.pi
instrumentefficiencypeak=0.4
flux_k = nlamp_k*flam*lamp_solidangle_k/(4*math.pi)*screen_efficiency*fiber_anglefraction_k*projector_efficiency*instrumentefficiencypeak  # Watt/Angstrom
#glassfilter = ascii.read()
photonflux_k = flux_k/(h*cc/(lam/1.e9))  #photons/Angstrom

#Onboard
focallength = 1500 # mm
fiberdiameter = 0.175  #mm
fiber_solidangle_on = math.pi*(fiberdiameter/2./focallength)**2
fiber_anglefraction_on=fiber_solidangle_on/math.pi
# Assuming 30cm projector to project to 1.5m. The diaphragm will be 2cm in diameter.  Assuming condensing lens has a focal length of 5cm.
lamp_solidangle_on = math.pi/5**2

