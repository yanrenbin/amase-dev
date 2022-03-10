import numpy as np
import matplotlib.pyplot as plt
import math

arcinrad = 206265.
amasefocal = 400 #mm
fiberdia = 0.05 #mm
amase_angres = 206265/amasefocal*fiberdia
nfibers = 547.
amase_angcov = np.sqrt(nfibers)*amase_angres #in arcsec. This is the effective diameter of the actual coverage of all the fibers. It is the same as first calculating the total covered area of all the fibers, then converting that to the diameter of the circle with the same area. 

manga_angres = 2.5
manga_angcov = np.sqrt(127.)*2.0
lvm_angres = arcinrad/(161*11.42)*0.330
lvm_angcov = np.sqrt(1801)*lvm_angres

kcwi_angres_s = np.sqrt(0.35**2+0.7**2)
kcwi_angcov_s = 2*np.sqrt(8*20./math.pi)

muse_angres = np.sqrt(0.6**2+0.35**2)
muse_angcov = np.sqrt(60.*60./math.pi)*2.

wham_angres = 3600.
wham_angcov = 2*3600.

mw_low = 0.4  #kpc
mw_hi = 3.0  #kpc

lmc = 50. #kpc
m31 = 780.

tenmpc = 10000

virgo = 16000

hundredmpc = 100000

amase_res = 15000
kcwi_res = 18000
muse_res = 3000
lvm_res = 3500
manga_res = 1800
wham_res = 30000

amase_dist = np.array([mw_low,mw_hi,lmc,m31,tenmpc,virgo])
amase_phyres = amase_angres*amase_dist*1000./arcinrad
amase_phycov = amase_angcov*amase_dist*1000./arcinrad

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([1000,35000])
for i in range(len(amase_dist)-1,-1,-1):
    xlow = amase_phyres[i]
    xhi = amase_phycov[i]
    ylow = amase_res*0.9
    ycen = amase_res
    yhi = amase_res*1.11
    polyx = [xlow,xhi,xlow,xlow]
    polyy = [ylow,ycen,yhi,ylow]
    ax.plot(polyx,polyy)
#    ax.plot([xlow,xlow],[ylow,yhi])

manga_dist = 100.*1000.
manga_phyres = manga_angres*manga_dist*1000./arcinrad
manga_phycov = manga_angcov*manga_dist*1000./arcinrad
xlow = manga_phyres
xhi = manga_phycov
ylow = manga_res*0.9
yhi = manga_res*1.11
ycen = manga_res
polyx = [xlow,xhi,xlow,xlow]
polyy = [ylow,ycen,yhi,ylow]
ax.plot(polyx,polyy)

lvm_phyres = amase_dist*lvm_angres*1000./arcinrad
lvm_phycov = amase_dist*lvm_angcov*1000./arcinrad
for i in range(len(amase_dist)-1,-1,-1):
    xlow = lvm_phyres[i]
    xhi = lvm_phycov[i]
    ylow = lvm_res*0.9
    ycen = lvm_res
    yhi = lvm_res*1.11
    polyx = [xlow,xhi,xlow,xlow]
    polyy = [ylow,ycen,yhi,ylow]
    ax.plot(polyx,polyy)

kcwi_dist = np.array([mw_low,mw_hi,tenmpc,hundredmpc])
kcwi_phyres = kcwi_dist*kcwi_angres_s*1000./arcinrad
kcwi_phycov = kcwi_dist*kcwi_angcov_s*1000./arcinrad
for i in range(len(kcwi_dist)-1,-1,-1):
    xlow = kcwi_phyres[i]
    xhi = kcwi_phycov[i]
    ylow = kcwi_res*0.9
    ycen = kcwi_res
    yhi = kcwi_res*1.11
    polyx = [xlow,xhi,xlow,xlow]
    polyy = [ylow,ycen,yhi,ylow]
    ax.plot(polyx,polyy)
    
wham_dist = np.array([mw_low, mw_hi])
wham_phyres = wham_dist*wham_angres*1000./arcinrad
wham_phycov = wham_dist*wham_angcov*1000./arcinrad
for i in range(len(wham_dist)-1,-1,-1):
    xlow = wham_phyres[i]
    xhi = wham_phycov[i]
    ylow = wham_res*0.9
    ycen = wham_res
    yhi = wham_res*1.11
    polyx = [xlow,xhi,xlow,xlow]
    polyy = [ylow,ycen,yhi,ylow]
    ax.plot(polyx,polyy)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim([1000,35000])
ax2.set_xlabel('Spatial Resolution and Coverage (arcsec)')
ax2.set_ylabel('Spectral Resolution')
ax2.plot([amase_angres,amase_angcov],[amase_res,amase_res],linewidth=6)
ax2.plot([lvm_angres,lvm_angcov],[lvm_res,lvm_res],linewidth=6)
ax2.plot([manga_angres,manga_angcov],[manga_res,manga_res],linewidth=6)
ax2.plot([muse_angres,muse_angcov],[muse_res,muse_res],linewidth=6)
#ax2.plot([ppak_angres,ppak_angcov],[ppak_res,ppak_res])
ax2.plot([kcwi_angres_s,kcwi_angcov_s],[kcwi_res,kcwi_res],linewidth=6)
ax2.plot([hetdex_angred,hetdex_angcov],[hetdex_res,hetdex_res],linewidth=6)

