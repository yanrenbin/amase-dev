import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt


dir = '/Users/renbin/astro/amase/labdata/qhy163m/fiber50/'
cendata= ascii.read(dir+'cendiag_tilt_b.txt')
offdata= ascii.read(dir+'offdiag_tilt_b.txt')
cenres = np.polyfit(cendata['X']+80,cendata['Y']+60,1)
offres = np.polyfit(offdata['X']+80,offdata['Y']+60,1)

pixsize=0.0038 #mm
offset = 7.0 #mm
microunit = 2.5 #mm

x0 = 4656/2.
y0_cen = x0*cenres[0]+cenres[1]
y0_off = x0*offres[0]+offres[1]

signcen = (cendata['X']+80-x0)/abs(cendata['X']+80-x0)
signoff = (offdata['X']+80-x0)/abs(offdata['X']+80-x0) 
distcen = np.sqrt((cendata['X']+80-x0)**2+(cendata['Y']+60-y0_cen)**2)*signcen
distoff = np.sqrt((offdata['X']+80-x0)**2+(offdata['Y']+60-y0_off)**2)*signoff

coeff = np.polyfit(distcen,cendata['Micro'],2)
coeffoff = np.polyfit(distoff+offset/pixsize,offdata['Micro'],2)

xx = np.arange(4300)-1300.
ycen = coeff[0]*xx**2+coeff[1]*xx+coeff[2]
yoff = coeffoff[0]*xx**2+coeffoff[1]*xx+coeffoff[2]
tiltcoe = np.polyfit(xx,ycen-yoff,1)
newoffmicro = offdata['Micro']+tiltcoe[0]*(distoff+offset/pixsize)+tiltcoe[1]

combxx = np.concatenate([distcen,distoff+offset/pixsize])
combmicro = np.concatenate([cendata['Micro'],newoffmicro])

combcoe = np.polyfit(combxx,combmicro,2)

corrmicro = combmicro-combxx*combcoe[1]

xx = (np.arange(37)-18.)/pixsize
plt.plot(distcen*pixsize,(cendata['Micro']-distcen*combcoe[1])*microunit,'x')
plt.plot((distoff+offset/pixsize)*pixsize,(newoffmicro-(distoff+offset/pixsize)*combcoe[1])*microunit,'x')
plt.plot(xx*pixsize,(combcoe[0]*xx**2+combcoe[2])*microunit)
plt.xlabel('Distance from center (mm)')
plt.ylabel('Slit piston (mm)')
plt.show()





