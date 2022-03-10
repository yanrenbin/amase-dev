import numpy as np
import matplotlib.pyplot as plt
import math

#wavelength = 15450*1.e-4
#wavelength = 6604.*1.e-4
wavelength = 4903.*1.e-4
#wavelength = 4610*1.e-4 #4943.*1.e-4
delta_n = 0.1
alpha = np.array(range(800),dtype=float)/10.*math.pi/180.
#const = 1.85 # math.pi*delta_n*thickness/wavelength
darr = np.array(range(2000),dtype=float)/100.+0.5
alpha2 = np.tile(alpha,[len(darr),1])
darr2 = np.transpose(np.tile(darr,[len(alpha),1]))
const = math.pi*delta_n*darr2/wavelength
eta_s = np.sin(const/np.cos(alpha2))**2
eta_p = np.sin(const/np.cos(alpha2)*np.cos(2*alpha2))**2

levels = np.array(range(11))/10
levels[10] = 0.999

plt.contour(alpha*180/math.pi,darr,eta_s,levels=levels)
plt.contour(alpha*180/math.pi,darr,eta_p,levels=levels)
plt.xlabel('Incidence Angle within Grating Medium (deg)')
plt.ylabel('Grating Thickness (assuming delta_n='+str(delta_n)+')')
plt.title('Efficiencies for the two polarizations ('+str(int(wavelength*1.e4))+'A)')
plt.show()


plt.contour(alpha*180/math.pi,darr,(eta_p+eta_s)/2.,levels=levels)
plt.xlabel('Incidence Angle within Grating Medium (deg)')
plt.ylabel('Grating Thickness (assuming delta_n='+str(delta_n)+')')
plt.title('Average efficiency for the two polarizations ('+str(int(wavelength*1.e4))+'A)')
plt.show()
