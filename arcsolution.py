import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import matplotlib.pyplot as plt
import pdb


#dir = '/Users/renbin/ASICAP/May25/'
#dir = '/Users/renbin/ASICAP/Jun2/' # Jun 1, 2 for Canon85mm
#dir = '/Users/renbin/ASICAP/Jun8/' # Jun 8 for Fujinon 110mm
dir = '/Users/renbin/astro/amase/labdata/qhy163m/2019-12-20/' # 2019-12-20 for Nikkor 58mm f/0.95. Log file in /Users/renbin/amase/labdata/qhy163m/log2019_12_18.txt

grpct = 0
f = open(dir+'parameters.txt','r')
for line in f:
    line = line.strip()
    columns = line.split()
    grpname = columns[0]
    if grpname == '#': continue
    if len(columns) == 1:
        addf = open(dir+'addparameters.txt','a')
        file = 'median_'+grpname+'.fits'
        hdu = fits.open(dir+file)
        data = hdu[0].data
        ysize,xsize = data.shape
        spatial = np.sum(data,axis=1)
        plt.plot(spatial)
        plt.yscale('log')
        plt.title(grpname)
        plt.show()
        arcy1 = int(input('Please input the Y1 for the arc region')) #1846
        arcy2 = int(input('Please input the Y2 for the arc region')) #2203
        arc = np.sum(data[arcy1:arcy2,:],axis=0)
        plt.plot(arc)
        plt.title(grpname)
        plt.show()
        nline = int(input('Please input the number of arc lines:'))
        line_x = np.zeros(nline)
        wave = np.zeros(nline)
        for i in range(0,nline):
            line_x[i] = float(input('Please input the x-position'))
            wave[i] = float(input('Please input the wavelength'))
        print(grpname,arcy1,arcy2,end=' ',file=addf)
        for x in line_x: print(x,end=' ',file=addf)
        for w in wave:  print(w,end=' ',file=addf)
        print('',file=addf)
        addf.close()
    else:
        xsize=4656
        ysize=3520
        arcy1 = columns[1]
        arcy2 = columns[2]
        nline = (len(columns)-3)/2
        if nline-int(nline) != 0: print('Error!')
        else: nline=int(nline)
        line_x = np.array(columns[3:3+nline])
        wave = np.array(columns[3+nline:len(columns)])
        line_x = line_x.astype(float)
        wave = wave.astype(float)


    #print('Please input the Y1 and Y2 for the arc region.')
    #spectral = np.sum(data,axis=0)

#    print('Please input the x position of arc lines.')
#    line_x = np.array(input()) #np.array([1319.,1973.,2045.,2959.])
#    print('Please input the wavelength of the arc lines in the same order.')
#    wave  = np.array(input()) #np.array([4358.337, 4077.838, 4046.572,3663.279])

    coeff = np.polyfit(line_x,wave,1)
    x0 = np.array([0,xsize-1])
    y0 = x0*coeff[0]+coeff[1]

    plt.plot(line_x,wave,'+')
    plt.plot(x0,y0)
    plt.title(grpname)
    plt.show()
#    pdb.set_trace()

    res = np.array([(grpname, arcy1,arcy2,coeff[0],coeff[1])],dtype=[('grpname','U20'),('arcy1','i8'),('arcy2','i8'),('slope','f8'),('intercept','f8')])
    if grpct == 0: 
       totres = res 
    else:
       totres = np.append(totres,res,axis=0)
    grpct +=1
f.close()

restable = Table(totres)
restable.write(dir+'arcsolution.fits',format='fits',overwrite=True)
