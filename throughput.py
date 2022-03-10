import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import matplotlib.pyplot as plt
import pdb
import numpy.ma as ma


#dir = '/Users/renbin/ASICAP/May25/'
#dir = '/Users/renbin/ASICAP/Jun2/'
#dir = '/Users/renbin/ASICAP/Jun8/'
dir = '/Users/renbin/astro/amase/labdata/qhy163m/2019-12-20/'

hdu = fits.open(dir+'arcsolution.fits')
arcsol= hdu[1].data
hdu.close()

grpct = 0
for kk in range(0,len(arcsol)): 
    grpname = arcsol['grpname'][kk]
    grpname = grpname.replace('wave_','')
    filename='median_'+grpname+'.fits'
    hdu = fits.open(dir+filename)
    data = hdu[0].data
#    exptime = float(hdu[0].header['EXPOINUS'])/1000000.
    exptime = float(hdu[0].header['EXPTIME'])
    ysize,xsize = data.shape
    spatial=np.sum(data,axis=1)
    indmax = np.argmax(spatial)
    upperpart = spatial[indmax:len(spatial)]
    lowerpart = spatial[0:indmax]
    uppdiff = upperpart-np.roll(upperpart,-1)
    bound_up = indmax+np.where(uppdiff < 0)[0][0]
    lowdiff = lowerpart-np.roll(lowerpart,1)
    bound_low = np.where(lowdiff<0)[0][-1]
    cutouts = data[bound_low:bound_up+1,:]
    indymax = np.argmax(cutouts,axis=0)
    indymin = np.argmin(cutouts,axis=0)
    maxval = cutouts[indymax,np.array(range(0,xsize))]
    minval = cutouts[indymin,np.array(range(0,xsize))]
    val70 = minval+(maxval-minval)*0.7
    above = cutouts > np.tile(val70,[bound_up-bound_low+1,1])
    y2d = np.transpose(np.tile(np.arange(bound_up-bound_low+1),[xsize,1]))
    my2d = ma.masked_array(y2d,mask=(1-above))
    center = np.mean(my2d,axis=0)
    
    p = np.where(maxval> 200)[0]

    coeff = np.polyfit(p,center[p],1)
    xx = np.array(range(0,xsize))
    yy = xx*coeff[0]+coeff[1]

    plt.plot(p,center[p],'+')
    plt.plot(xx,yy)
    plt.show()
    redoflag = str(input('Enter Y or y if you want to redo the fit.'))
    while (redoflag == 'y') | (redoflag == 'Y'):
        xlow =  int(input('Enter the lower limit for X'))
        xhi = int(input('Enter the upper limit for X'))
        order = int(input('Enter the order of the polynomial'))
        good = np.where( (p > xlow) & (p < xhi))[0]
        coeff = np.polyfit(p[good],center[p[good]],order)
        if order == 1: 
            yy = xx*coeff[0]+coeff[1]
        if order == 2:
            yy = xx*xx*coeff[0]+xx*coeff[1]+coeff[2]
        plt.plot(p,center[p],'+')
        plt.plot(xx,yy)
        plt.show()
        redoflag = 'N'
        redoflag = str(input('Enter Y or y if you want to redo the fit.'))
    cenrow = np.around(yy+bound_low+0.5).astype(int)
    shift = cenrow-cenrow[0]
    accu = np.zeros([21,xsize]) # was 30 for sigma105
    flux1 = np.zeros(xsize)
    maxv = np.zeros(xsize)
    minv = np.zeros(xsize)
    for j in range(0,xsize):
# Sigma105
#        med = np.median(data[cenrow[j]+50:cenrow[j]+80,j])
##        xp = np.concatenate([np.array(range(cenrow[j]-30,cenrow[j]-20)),np.array(range(cenrow[j]+21,cenrow[j]+31))])
##        cont = np.polyfit(xp,data[xp,j],1)
##        xin = np.array(range(cenrow[j]-20,cenrow[j]+21))
##        flux1[j] = np.sum(data[cenrow[j]-20:cenrow[j]+21,j]-(xin*cont[0]+cont[1]))
        xp = np.concatenate([np.array(range(cenrow[j]-400,cenrow[j]-200)),np.array(range(cenrow[j]+200,cenrow[j]+300))])
        cont = np.polyfit(xp,data[xp,j],1)
        xin = np.array(range(cenrow[j]-200,cenrow[j]+201))
#        flux1[j] = np.sum(data[cenrow[j]-200:cenrow[j]+201,j]-(xin*cont[0]+cont[1]))
        for i in range(0,21):
# Sigma105
#            accu[i,j] = np.sum(data[cenrow[j]-i*20:cenrow[j]+i*20+1,j]-med)
            accu[i,j] = np.sum(data[cenrow[j]-i*10:cenrow[j]+i*10+1,j]-(xin[200-i*10:200+i*10+1]*cont[0]+cont[1]))
        maxv[j] = np.max(accu[:,j])
        minv[j] = np.min(accu[10:21,j])
# Sigma105
#    flux = accu[25,:]
    flux = accu[10,:]
    err = np.amax(np.vstack((maxv-flux,flux-minv)),axis=0)
# convert to per Angstrom per second
#    pdb.set_trace()
    flux = flux/abs(arcsol['slope'][kk])/exptime 
    wave = xx*arcsol['slope'][kk]+arcsol['intercept'][kk]
    plt.plot(wave,flux)
    plt.title(grpname)
    plt.show()

#    if grpname=='grp3_quartz' or grpname =='grp3_sigma105':
#       pdb.set_trace()
    
    res = np.array([(grpname, arcsol['arcy1'][kk],arcsol['arcy2'][kk],arcsol['slope'][kk],arcsol['intercept'][kk],wave,flux,err)],dtype=[('grpname','U20'),('arcy1','i8'),('arcy2','i8'),('slope','f8'),('intercept','f8'),('wave','4656f4'),('flux','4656f4'),('err','4656f4')])
    if kk == 0: 
       totres = res 
    else:
       totres = np.append(totres,res,axis=0)
    grpct +=1
    hdu.close()
    

restable = Table(totres)
restable.write(dir+'throughput.fits',format='fits',overwrite='True')

