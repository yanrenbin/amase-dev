import numpy as np
import math
from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib

# 1 Rayleigh is defined as a column emission rate of 10^6 photons per cm^2 (column) per second emission in all directions (4pi).
#  corresponds to a surface brightness of 
#               2.41e-7 erg/cm^2/s/steridian at Halpha
#            = 5.665e-18 erg/cm^2/s/arcsec^2

cc = 299792458. 
hh = 6.626e-34
halpha_photon_energy_in_J = hh*cc/(6563.*1.e-10)
sdss_mirror = 3.5814e4 # in cm^2
manga_fiberarea = math.pi # in arcsec^2
rayleigh = 5.665e-18 # erg/cm^2/s/arcssec^2

exposuretime = 900 # sec

transmission_air= 0.88 # Atmosphere ransmission at 1.4 Airmass at Halpha wavelength
n_sb=41
sb_array = 10**(np.arange(n_sb)/10-1)
#fraction_fwhm = 0.76  # 76% of the flux is within the FWHM for a Gaussian.
fraction_fwhm = 0.60  # For a line with a 20km/s intrinsic FWHM. After convolving with our resolution. then within the central 3.86 pixel (instrument FWHM), there are 60% of the flux.

efficiency = 0.4
photons_halpha = sb_array*rayleigh*sdss_mirror*manga_fiberarea/(halpha_photon_energy_in_J*1.e7)*exposuretime*transmission_air*efficiency

mag = np.arange(33)/2. # AB magnitude of star in r-band
star_fnu = 10**(-0.4*mag)*3631*1.e-23 # erg/s/Hz/cm^2
star_flambda = star_fnu*3.e8*1.e10/(6563.**2)

fcoll = 200.0 # mm
fcam = 58.0 # mm
fibersize = 50 # mm
pixsize=3.76 # micron
incidenceangle = 54.7*math.pi/180.
refindex = 1.46 # index of the prism glass
readnoise = 1.68 #/2.
darkcurrent = 0.0022/2. #0.02/16 #0.002
fratio = 3.7 #3.7 is the effective f-ratio after considering the average vignetting. #2.8 #3.5 is the pre-FRD focal ratio for a spectrograph that only accept f/3.33
pixsampling = fibersize/fcoll*fcam/pixsize
psfsize = pixsampling
resolution = 15000. #2*fcoll/(fibersize/1000.)*np.tan(incidenceangle)*refindex
dispersion = 6563./resolution/pixsampling # Angstrom/pixel
skyfactor = 1 #2.5 

photons_star_amase = star_flambda*(40/fratio)**2/4.*math.pi*dispersion*psfsize/(halpha_photon_energy_in_J*1.e7)*exposuretime*transmission_air*efficiency # number of photons per resolution element.



readnoise_manga = 2.75
fratio_manga = 5
fibersize_manga = 120.
resolution_manga = 1875.
pixsampling_manga= 2.5
psfsize_manga = pixsampling_manga
dispersion_manga = 6563./resolution_manga/pixsampling_manga

photons_star_manga = star_flambda*sdss_mirror*dispersion_manga*psfsize_manga/(halpha_photon_energy_in_J*1.e7)*exposuretime*transmission_air*efficiency # number of photons per resolution element

background_manga = 90./psfsize_manga # 80 cts/exposure/pixel ; 900 sec exposure in reduced 1-d spectrum

n_sky = 21
skyarray = 10**(np.arange(21)/10.) 
sky2d = np.transpose(np.tile(skyarray,[n_sb,1]))
photons_halpha_2d = np.tile(photons_halpha,[n_sky,1])

# For an unresolved emission line at the limit of detection
signalratio =  (fibersize/fibersize_manga)**2*(fratio_manga/fratio)**2 # integrated flux in an unresolved line.
noisesq_perpixel = sky2d*signalratio*background_manga/dispersion_manga*dispersion+readnoise**2+darkcurrent*900.
noisesq_perpixel_manga = sky2d*background_manga + readnoise_manga**2
noisesq_ratio = (noisesq_perpixel*psfsize**2)/(noisesq_perpixel_manga*psfsize_manga**2)
snratio = signalratio/np.sqrt(noisesq_ratio)

signal_manga = photons_halpha_2d*fraction_fwhm
signal_amase = signal_manga*signalratio
noisesq_manga = signal_manga+noisesq_perpixel_manga*psfsize_manga**2
noisesq_amase = signal_amase+noisesq_perpixel*psfsize**2

sn_manga = signal_manga/np.sqrt(noisesq_manga)
sn_amase = signal_amase/np.sqrt(noisesq_amase)

noisesq_star_amase = photons_star_amase + noisesq_perpixel[0,0]*psfsize**2
sn_star_amase = photons_star_amase/np.sqrt(noisesq_star_amase)

print("MaNGA snratio    :",1.0)
print("      signalratio:",1.0)
print("           sky   :",skyfactor*background_manga*psfsize_manga**2)
print("           RN**2 :",readnoise_manga**2*psfsize_manga**2)
print("           DC    :",0)
print("AMASE:")
#print("AMASE: snratio    :",snratio)
print("      signalratio:",signalratio)
print("           sky   :",skyfactor*signalratio*background_manga/dispersion_manga*dispersion*psfsize**2)
print("           RN**2 :",readnoise**2*psfsize**2)
print("           DC    :",darkcurrent*900.*psfsize**2)
print("       Fiber size:",fibersize)
print("       Resolution:",resolution)
print("       Pixel FWHM:",pixsampling)

readnoise_lvm = 2.5 #1.8
fratio_lvm = 3.7
fibersize_lvm = 107.
resolution_lvm = 4160.
pixsampling_lvm = 3.37
psfsize_lvm = pixsampling_lvm
dispersion_lvm = 6563./resolution_lvm/pixsampling_lvm

signalratio_lvm = (fibersize_lvm/fibersize_manga)**2*(fratio_manga/fratio_lvm)**2 
noisesq_perpixel_lvm = sky2d*signalratio_lvm*background_manga/dispersion_manga*dispersion_lvm+readnoise_lvm**2
noisesq_ratio_lvm = (noisesq_perpixel_lvm*psfsize_lvm**2)/(noisesq_perpixel_manga*psfsize_manga**2)
snratio_lvm = signalratio_lvm/np.sqrt(noisesq_ratio_lvm)

signal_lvm = signal_manga*signalratio_lvm
noisesq_lvm = signal_lvm+noisesq_perpixel_lvm*psfsize_lvm**2
sn_lvm = signal_lvm/np.sqrt(noisesq_lvm)

print("LVM : ")
#print("LVM : snratio    :",snratio_lvm)
print("      signalratio:",signalratio_lvm)
print("           sky   :",skyfactor*signalratio_lvm*background_manga/dispersion_manga*dispersion_lvm*psfsize_lvm**2)
print("           RN**2 :",readnoise_lvm**2*psfsize_lvm**2)
print("           DC    :",0)
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 13}
matplotlib.rc('font', **font)

darktime=plt.plot(sb_array,sn_amase[0,:],label='Under Dark Time Sky at APO')
darkx10 =plt.plot(sb_array,sn_amase[10,:],label='Under 10x brighter sky')
darkx100=plt.plot(sb_array,sn_amase[20,:],label='Under 100x brighter sky')
#plt.plot(sb_array,sn_manga[0,:],'--')
#plt.plot(sb_array,sn_manga[10,:],'--')
#plt.plot(sb_array,sn_manga[20,:],'--')


plt.xscale("log")
plt.yscale("log")
plt.xlabel('Halpha Surface Brightness (Rayleigh)')
plt.ylabel('S/N within FWHM in 15 min')
plt.legend()
plt.show()

darktime=plt.plot(sb_array,sn_amase[0,:]*2,label='AMASE-P DarkTime (S/N x 2)')
darkx10 =plt.plot(sb_array,sn_amase[10,:]*2,label='AMASE-P 10x brighter sky (S/N x 2)')
darkx100=plt.plot(sb_array,sn_amase[20,:]*2,label='AMASE-P 100x brighter sky (S/N x 2)')
darktime_lvm=plt.plot(sb_array,sn_lvm[0,:],'--',label='DESI DarkTime')
darkx10_lvm=plt.plot(sb_array,sn_lvm[10,:],'--',label='DESI 10x brighter sky')
darkx100_vlm=plt.plot(sb_array,sn_lvm[20,:],'--',label='DESI 100x brighter sky')

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'H$\alpha$ Surface Brightness (Rayleigh)')
plt.ylabel(r'S/N per fiber of H$\alpha$ in 15 min')
plt.xlim([0.5,3000])
plt.ylim([0.2,300])
plt.legend()
plt.show()
