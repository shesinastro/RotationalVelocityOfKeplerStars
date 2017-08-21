from astropy.io import fits


hdulist = fits.open('HD39587_482732_55586_UVB+VIS.fits')
prihdr = hdulist[0].header
prihdr['CRVAL1']=3.47612
prihdr['CDELT1']=0.000011155 
prihdr['CD1_1']=0.000011155
hdulist.writeto('G0_new.fits')


