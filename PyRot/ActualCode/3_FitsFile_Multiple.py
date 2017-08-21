# -*- coding: utf-8 -*-

from specutils.wcs import specwcs
from astropy.io import fits
from astropy import units as u
from specutils.io import read_fits
import os
import glob
import numpy as np
import fnmatch

Halpha=6562.760


for file in os.listdir('.'):
    if fnmatch.fnmatch(file, 'tar*.dat'):
        os.remove(file)   
        
        
            
            


myPath="/Users/Serena/Desktop/DEF_TRES"
fitsCounter = len(glob.glob1(myPath,"*.fits"))


cc=[[]]
ww=[[]]
for j in range(0, fitsCounter, 1):
    name='sp'+str(j)+'.fits'
    spectra_list = read_fits.read_fits_spectrum1d(name)

    

    for i in range(0, len(spectra_list)-1, 1):
        arr=spectra_list[i].dispersion.value #value is to convert object quantity
        w_min=min(arr)
        w_max=max(arr)
        if (Halpha > w_min) & (Halpha < w_max):
            level=i
        
    cc.append(spectra_list[level].flux.value)
    ww.append(spectra_list[level].dispersion.value)
        
ww_new=np.delete(ww,0)
cc_new=np.delete(cc,0)     

for k in range(len(ww_new)):
    print len(ww_new)
    f=open('target'+str(k)+'.dat', 'w');
    for n in range(0, len(ww_new[k]), 1):
        f.write('%f\t%f\n' % (ww_new[k][n], cc_new[k][n]))  
    
    
        