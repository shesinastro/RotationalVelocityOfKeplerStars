# -*- coding: utf-8 -*-

#from __future__ import print_function
#from astropy import stats
#import emcee
#import corner
import numpy as np
#import scipy.optimize as op
import matplotlib.pyplot as pl
#from matplotlib.ticker import MaxNLocator
#from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl
#import scipy.integrate as sci
#import scipy.stats
#from scipy.optimize import curve_fit
#from PyAstronomy.pyasl import binningx0dt
import os
import sys

VSINI=float(sys.argv[1])
sigma_noise=float(sys.argv[2])

filename = "BroadenedSpectrum.dat" 
if os.path.exists(filename):
    os.remove(filename)
    
filename='SynthSpectrum_selected_norm_withErrors.dat'   
if os.path.exists(filename):
    os.remove(filename)   
    
    
            
f=open('BroadenedSpectrum.dat', 'a');
g=open('SynthSpectrum_selected_norm_withErrors.dat', 'a');
template=np.genfromtxt('SynthSpectrum_selected_norm.dat', dtype=None, names=('l_template', 'N_template', 'error_template'))

ll_template_temp = template['l_template']
counts_template_temp = template['N_template'] 
ll_template=ll_template_temp[(ll_template_temp>6540) & (ll_template_temp< 6580)] #mask all the spectrum but small region around line of interest
counts_template = counts_template_temp[(ll_template_temp>6540) & (ll_template_temp< 6580)]

ll_template_array = np.asarray(ll_template)
counts_template_array = np.asarray(counts_template)
counts_template_array_selected = counts_template_array[(ll_template_array > 6525.) & (ll_template_array < 6550)]
mean_count = np.mean(counts_template_array_selected)
rms = np.sqrt(np.sum((counts_template_array_selected-mean_count)**2)/len(counts_template_array_selected))
errors_template = np.repeat(rms, len(ll_template)) #error is equal to rms of the flat part of the spectrum

noise = np.random.normal(size=len(counts_template),scale=sigma_noise)
rflux_with_sigma = pyasl.rotBroad(ll_template, counts_template, 0.5, VSINI) + noise
errors_target= errors_template


pl.close()
pl.plot(ll_template, rflux_with_sigma, 'b')
pl.plot(ll_template, counts_template, 'm')
pl.show()


for i in range(0, len(rflux_with_sigma)-1, 1):
    f.write('%f\t%f\t%f\n' % (ll_template[i], rflux_with_sigma[i], errors_target[i]))
    
    
for i in range(0, len(ll_template)-1, 1):
    g.write('%f\t%f\t%f\n' % (ll_template[i], counts_template[i], errors_template[i]))    