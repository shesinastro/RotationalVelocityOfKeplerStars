import numpy as np
from scipy import interpolate
import os
import fnmatch
import glob
#import matplotlib.pyplot as plt


    
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*norm.dat'):
        os.remove(file)     
            
def region_around_line(w, flux, cont):
    '''cut out and normalize flux around a line

    Parameters
    ----------
    w : 1 dim np.ndarray
    array of wavelengths
    flux : np.ndarray of shape (N, len(w))
    array of flux values for different spectra in the series
    cont : list of lists
    wavelengths for continuum normalization [[low1,up1],[low2, up2]]
    that described two areas on both sides of the line
    '''
    #index is true in the region where we fit the polynomial
    indcont = ((w > cont[0][0]) & (w < cont[0][1])) |((w > cont[1][0]) & (w < cont[1][1]))
    indrange = len(w)
    f = np.zeros(indrange)
    
    # fit polynomial of second order to the continuum region
    linecoeff = np.polyfit(w[indcont], flux[indcont], 2)
    # divide the flux by the polynomial and put the result in our
    # new flux array
    f_norm = flux/np.polyval(linecoeff, w)
    return w, f_norm
  
myPath="/Users/Serena/Desktop/DEF_TRES"
Counter = len(glob.glob1(myPath,"target*.dat"))

target=np.genfromtxt('target0.dat', dtype=None, names=('l_target', 'N_target'))
ll = target['l_target']    
mat=np.zeros((len(ll), Counter))   

    
for k in range(Counter):  
    target=np.genfromtxt('target'+str(k)+'.dat', dtype=None, names=('l_target', 'N_target'))
    ww_temp = target['l_target']
    cc_temp = target['N_target'] 
    #ww=ww_temp[(ww_temp>6540) & (ww_temp< 6580)]
    #cc = cc_temp[(ww_temp>6540) & (ww_temp< 6580)]    
    ww=ww_temp
    cc=cc_temp
    ww_norm, cc_norm = region_around_line(ww, cc, [[6540, 6545],[6576, 6580]])    




    
    mat[:,k]=cc_norm

        
  
cc_norm_ave=[]    
for i in range(len(mat)):
    cc_norm_ave.append(sum(mat[i,:])/Counter)
    
    
f=open('target_norm.dat', 'a');
for i in range(0, len(cc_norm), 1):
    f.write('%f\t%f\n' % (ww_norm[i], cc_norm_ave[i]))
    
    
#plt.plot(ww_norm, cc_norm_ave)
#plt.show()



    
    