import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import pdb
import scipy.signal as signal

#dir='/Users/renbin/astro/specdragon/sigma105/vignetting/fits/'
#dir='/Users/renbin/astro/mase/GFX50R/vignetting/'
dir ='/Users/renbin/astro/amase/labdata/Z7/'

#hdu=fits.open(dir+'LB9X0745.fit')
#hdu =fits.open(dir+'DSCF0054.fit')  # 46mm filter , effectively f/2.4
#hdu =fits.open(dir+'DSCF0059.fit')  # 37mm filter , effectively f/2.97
#hdu =fits.open(dir+'sevencm/DSCF0226.fit')  # 37mm filter at 73mm away , effectively f/2.97
#hdu =fits.open(dir+'sevencm/DSCF0232.fit')  # 46mm filter at 73mm away , effectively f/2.97
#hdu = fits.open(dir+'DSC_0352.fits') # Nikon 58mm lens wide open at f/0.95, 77mm pupil at 61.5mm away
#hdu = fits.open(dir+'DSC_0098.fits') # Canon 200mm as collimator, Nikon 58mm lens as camera, front of lens separated by 251.8mm.
#hdu = fits.open(dir+'DSC_0099.fits')# Canon 200mm as collimator, Nikon 58mm lens as camera, front of lens separated by 277.2mm. 
hdu = fits.open(dir+'DSC_0359.fits')# Nikon 58mm lens wide open at f/0.95, 77mm pupil at 0.5mm away from the front edge.
data = hdu[0].data
hdu.close()

#camera = 'fuji110'
camera = 'Nikon_Z7'

ys,xs = data.shape
if camera == 'sigma105':
    image = data[2:ys-100,126:xs].copy()
    overscan = data[2:ys,0:126].copy()
    bias = np.median(overscan)
    image = image-bias
    pixel = 6.94 # this was tested using the EOS-1DX
if camera == 'fuji110':
    image = data[2:ys,0:8280]
    biashdu = fits.open(dir+'/../bias/avebias.fits')
    bias = biashdu[0].data
    image = image-bias
    pixel = 5.0 # this was tested using the Fuji GFX50R
if camera == 'Nikon_Z7':
    image = data
    biashdu = fits.open(dir+'bias/avebias.fits')
    bias = biashdu[0].data
    image = image-bias
    pixel = 4.35

ysize,xsize = image.shape
xsize = int(xsize)
ysize = int(ysize)
xodd = np.array(range(int(xsize/2)))*2+1
xeven = np.array(range(int(xsize/2)))*2
yodd = np.array(range(int(ysize/2)))*2+1
yeven = np.array(range(int(ysize/2)))*2
flag = np.zeros([ysize,xsize],dtype=int)
# red is 0, green is 1, blue is 2.
# Nikon Z7 follows the above rule.
for y in yeven:
    flag[y,xeven] = 1
    flag[y,xodd] = 2
for y in yodd:
    flag[y,xodd] = 1
    flag[y,xeven] = 0

xcoords = np.array(range(xsize))
ycoords = np.array(range(ysize))
x2 = np.tile(xcoords,[ysize,1])
y2 = np.transpose(np.tile(ycoords,[xsize,1]))
#dist = np.sqrt((x2-(xsize/2.-0.5))**2 + (y2-(ysize/2.-0.5))**2)

midx = int(xsize/2)
midy = int(ysize/2)
imgcutouts = image[midy-10:midy+10,midx-10:midx+10]
flagcutouts = flag[midy-10:midy+10,midx-10:midx+10] 
maxval = np.median(imgcutouts[flagcutouts == 1])

gindex = np.where(flag[:,midx] ==1)[0]
uppart=np.where((image[gindex,midx] < maxval*0.8)
    & (image[gindex,midx] > maxval*0.7) & (ycoords[gindex] < midy))[0]
coeff=np.polyfit(ycoords[gindex[uppart]],image[gindex[uppart],midx],1)
uphalf = (0.75*maxval-coeff[1])/coeff[0]

downpart=np.where((image[gindex,midx] < maxval*0.8) 
    & (image[gindex,midx] > maxval*0.7) & (ycoords[gindex] > midy))[0]
coeff = np.polyfit(ycoords[gindex[downpart]],image[gindex[downpart],midx],1)    
downhalf = (0.75*maxval-coeff[1])/coeff[0]
ypos = (uphalf+downhalf)/2.

gindex = np.where(flag[midy,:] ==1)[0]
leftpart = np.where((image[midy,gindex] < maxval*0.8) 
    & (image[midy,gindex] > maxval*0.7) & (xcoords[gindex] < midx))[0]
coeff = np.polyfit(xcoords[gindex[leftpart]],image[midy,gindex[leftpart]],1)
lefthalf = (0.75*maxval-coeff[1])/coeff[0]
rightpart = np.where((image[midy,gindex] < maxval*0.8) 
    & (image[midy,gindex] > maxval*0.7) & (xcoords[gindex] > midx))[0]
coeff = np.polyfit(xcoords[gindex[rightpart]],image[midy,gindex[rightpart]],1)
righthalf = (0.75*maxval-coeff[1])/coeff[0]
xpos = (lefthalf+righthalf)/2.


dist = np.sqrt((x2-xpos)**2 + (y2-ypos)**2)
#plt.plot(dist[flag ==1],image[flag ==1])
x2b=np.reshape(x2[flag==2],[int(ysize/2),int(xsize/2)])
y2b=np.reshape(y2[flag==2],[int(ysize/2),int(xsize/2)])
x2r=np.reshape(x2[flag==0],[int(ysize/2),int(xsize/2)])
y2r=np.reshape(y2[flag==0],[int(ysize/2),int(xsize/2)])
x2g=np.reshape(x2[flag==1],[int(ysize),int(xsize/2)])
y2g=np.reshape(y2[flag==1],[int(ysize),int(xsize/2)])

imgr = np.reshape(image[flag==0],[int(ysize/2),int(xsize/2)])
imgb = np.reshape(image[flag==2],[int(ysize/2),int(xsize/2)])
imgg = np.reshape(image[flag==1],[int(ysize),int(xsize/2)]) 
smimgr = signal.medfilt2d(imgr,[7,7])
smimgb = signal.medfilt2d(imgb,[7,7])
smimgg = signal.medfilt2d(imgg,[7,7])

xquarter = int(xsize/4)
yquarter = int(ysize/4)
xhalf = int(xsize/2)
yhalf = int(ysize/2)
maxr = np.median(smimgr[yquarter-10:yquarter+10,xquarter-10:xquarter+10])
maxb = np.median(smimgb[yquarter-10:yquarter+10,xquarter-10:xquarter+10])
maxg = np.median(smimgg[yhalf-10:yhalf+10,xquarter-10:xquarter+10])

level = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.98])
#plt.contour((x2b-xpos)*pixel/1000.,(y2b-ypos)*pixel/1000.,smimgb,levels=maxb*level)
#plt.contour((x2r-xpos)*pixel/1000.,(y2r-ypos)*pixel/1000.,smimgr,levels=maxr*level)
plt.contour((x2g-xpos)*pixel/1000.,(y2g-ypos)*pixel/1000.,smimgg,levels=maxg*level)
#plt.plot([12,12,-12,-12,12],[-12,12,12,-12,-12],':')
plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')
#plt.title('Sigma 105mm f/1.4 ART Vignetting Profile')
#plt.title(camera+' at f/2.4 (pupil 73mm away)')
#plt.title(camera+' at f/0.95 (pupil 61.5mm away)')
#plt.title('Canon 200mm F/2 + Nikon 58mm F/0.95, separated by 251.8mm')
plt.title(camera+'at f/0.95 (77mm pupil 0.5mm away)')
plt.grid(linestyle='-')
plt.show()



