from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt


dir = '/Users/renbin/astro/amase/optical_model/'
flatslit = ascii.read(dir+'vignettingmodel_flatslit_default.txt')
#curvslit = ascii.read(dir+'vignettingmodel_curvedslit.txt')
curvslit = ascii.read(dir+'vignettingmodel_flatslit_6040.txt')

flatthrupt=np.vstack([flatslit['col2'],flatslit['col3'],flatslit['col4'],flatslit['col5']])#,flatslit['col6'],flatslit['col7'],flatslit['col8']])
curvthrupt=np.vstack([curvslit['col2'],curvslit['col3'],curvslit['col4'],curvslit['col5']])#,curvslit['col6'],curvslit['col7'],curvslit['col8']])

flatrelative = flatthrupt/np.max(flatthrupt)
curvrelative = curvthrupt/np.max(curvthrupt)

waveflat = flatslit['col1']/10
wavecurv = curvslit['col1']/10

#xarr=[0,10,20,25,30,35,40]
xarr = [0,10,20,30]
labelflat = [str("de %3d nm" % p) for p in waveflat]
labelcurv = [str("64 %3d nm" % p) for p in wavecurv]
colorarr = ['C0','C1','C2','C3','C4','C5','C6','C7']

lineall = plt.plot(xarr,flatrelative[:,0],label=labelflat[0],color=colorarr[0])
for i in np.arange(np.shape(flatrelative)[1]-1)+1:
   line=plt.plot(xarr,flatrelative[:,i],label=labelflat[i],color=colorarr[i])
   lineall = lineall+line

#plt.xlabel('Distance from slit center')
#plt.ylabel('Throughtput relative to the center of the detector')
#plt.title('Flat slit')
#plt.legend()
#plt.show()

lineall = plt.plot(xarr,curvrelative[:,0],'--',label=labelcurv[0],color=colorarr[0])
for i in np.arange(np.shape(curvrelative)[1]-1)+1:
   line=plt.plot(xarr,curvrelative[:,i],'--',label=labelcurv[i],color=colorarr[i])
   lineall = lineall+line

plt.xlabel('Distance from slit center')
plt.ylabel('Throughtput relative to the center of the detector')
#plt.title('Curved slit')
plt.legend()
plt.show()


