# -*- coding: utf-8 -*-



from __future__ import print_function

from astropy import stats

import emcee
import corner
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl
import scipy.integrate as sci
import scipy.stats
from scipy.optimize import curve_fit
from PyAstronomy.pyasl import binningx0dt
import sys

Halpha = 6562.760

temp = np.loadtxt("vsini_chi_sqr.txt")
vsini_guess = float(temp)
if vsini_guess == 0.0:
    vsini_guess = 0.1

# Reproducible results!
np.random.seed(123)

mine = '#008B8B'





template=np.genfromtxt('SynthSpectrum_selected_norm_withErrors.dat', dtype=None, names=('l_template', 'N_template', 'error_template'))
target=np.genfromtxt('BroadenedSpectrum.dat', dtype=None, names=('l_target', 'N_target', 'error_target'))




template=np.genfromtxt('SynthSpectrum_selected_norm_withErrors.dat', dtype=None, names=('l_template', 'N_template', 'error_template'))
target=np.genfromtxt('BroadenedSpectrum.dat', dtype=None, names=('l_target', 'N_target', 'error_target'))

ll_template_temp = template['l_template']
counts_template_temp = template['N_template']-1 #normalize to 0 to avoid MCMC forcing model to reach continuum at 0
ll_template=ll_template_temp[(ll_template_temp>6540) & (ll_template_temp< 6580)] #mask all the spectrum but small region around line of interest
counts_template = counts_template_temp[(ll_template_temp>6540) & (ll_template_temp< 6580)]


ll_target_temp = target['l_target']
counts_target_temp = target['N_target']-1 #normalize to 0 to avoid MCMC forcing model to reach continuum at 0
errors_target_temp = target['error_target']
ll_target=ll_target_temp[(ll_target_temp>6540) & (ll_target_temp< 6580)]
counts_target = counts_target_temp[(ll_target_temp>6540) & (ll_target_temp< 6580)]
errors_target = errors_target_temp[(ll_target_temp>6540) & (ll_target_temp< 6580)] #error is equal to rms of the flat part of the spectrum





# interpolate the target spectrum, as the spectrum has to be evaluated at evenly spaced wavelenghts
def interpol_template(x):
    return np.interp(x, ll_template, counts_template)
    
    
def interpol_target(x):
    return np.interp(x, ll_target, counts_target)   
    
def interpol_error(x):
    return np.interp(x, ll_target, errors_target) 
    

l_min_template=min(ll_template)
l_max_template=max(ll_template)
l_min_target=min(ll_target)
l_max_target=max(ll_target)


sum_template = 0

for i in range(0, len(ll_template)-1, 1):
    sum_template += ll_template[i+1]-ll_template[i]
    
delta_l_template =  sum_template/len(ll_template)   

sum_target = 0

for i in range(0, len(ll_target)-1, 1):
    sum_target += ll_target[i+1]-ll_target[i]

delta_l_target=  sum_target/len(ll_target)
       

#rebin to same wavelenth scale, using as standard the lower-res spectrum
delta_l_def = max(delta_l_template, delta_l_target)  


new_wave = np.arange(l_min_target, l_max_target, delta_l_def)
new_counts_template = interpol_template(new_wave)
new_counts_target = interpol_target(new_wave)
new_errors_target = interpol_error(new_wave)    
 
 
x = np.asarray(new_wave)
y = np.asarray(new_counts_target)
yerr = np.asarray(new_errors_target)
 





#cross-correlation
#rv, cc = pyasl.crosscorrRV(x, y, x, new_counts_template, -100., 100., 1.,skipedge=20)
spec = np.column_stack((x, y))
bootres=stats.bootstrap(spec, bootnum=100)

RV_list = []

#use full template for crosscorrelation
for i in range(0, 99, 1):
    rv, cc = pyasl.crosscorrRV(bootres[i][:,0], bootres[i][:,1], ll_template_temp, counts_template_temp, -5., 5., 1.,skipedge=20)
    maxind = np.argmax(cc) 
    RV_list.append(rv[maxind])
    
RV_mean = np.mean(RV_list)
RV_std = np.std(RV_list)

print("RV mean %f" % RV_mean)
print("RV std %f" % RV_std)


c_template_shifted, w_template_shifted = pyasl.dopplerShift(x, new_counts_template,RV_mean/2., edgeHandling=None, fillValue=None)
c_template_shifted_new = c_template_shifted[np.logical_not(np.isnan(c_template_shifted))]
w_template_shifted_new = w_template_shifted[np.logical_not(np.isnan(c_template_shifted))]

c_template_shifted_rebinned = [] 
for ll in x: 
    c_template_shifted_rebinned.append(np.interp(ll, w_template_shifted_new, c_template_shifted_new))
c_template_shifted_rebinned_arr = np.asarray(c_template_shifted_rebinned)



#plt.plot(x, c_template_shifted_rebinned_arr)
#plt.plot(x, y)
#plt.show()

RV_mean_new=RV_mean/2.

l0_guess = (RV_mean_new/300000.)*Halpha + Halpha
delta_l0_guess = (RV_std/300000.)*l0_guess
delta_l0_guess_prior = 50.0*delta_l0_guess #multiply by 10 to be generous


#prior function
def lnprior(theta):
    vsini, l0, lnf = theta #ldark=limb darkening #lnf, log_e(f), f how much errors are underestimated
    #if 0 < vsini < 437 and -10.0 < lnf < 1.0 and Halpha-delta_l0_guess_prior < l0 < Halpha+delta_l0_guess_prior: #437 km/s is break up speed Sun like star
    if 0 < vsini < 200 and -10.0 < lnf < 1.0 and l0_guess-delta_l0_guess_prior < l0 < l0_guess+delta_l0_guess_prior: #437 km/s is break up speed Sun like star    
        return 0.0
    return -np.inf
    
#likelihood function  
def lnlike(theta, x, y, yerr):
    vsini, l0, lnf = theta
    model = pyasl.fastRotBroad(x, c_template_shifted_rebinned_arr, 0.5, vsini, l0) #assume limb darkening of 0.5
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
     
     
        
    
#posterior probability function    
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
    
print(lnprob([vsini_guess, l0_guess, -7], x, y, yerr)) #f=-0.5

ndim, nwalkers = 3, 100

#ldark_guess = 0.5
#l0_guess = Halpha
lnf_guess = -0.7
pos = [[vsini_guess, Halpha, lnf_guess]+1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))   

# Clear and run the production chain.
print("Running MCMC..."); sampler.run_mcmc(pos, 500, rstate0=np.random.get_state()); print("Done.") #pos, 500


pl.clf()
fig, axes = pl.subplots(3, 1, sharex=True, figsize=(8, 9))
axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
#axes[0].axhline(m_true, color="#888888", lw=2)
axes[0].set_ylabel("$V_\mathrm{sini}$ [km/s]")
fig.tight_layout(h_pad=0.0)


axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].set_ylabel("l0")
axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
axes[2].yaxis.set_major_locator(MaxNLocator(5))
axes[2].set_ylabel("ln(f)")
axes[2].set_xlabel("step number")

#fig.show()
fig.savefig("line-time.pdf")



#burnin = 50  #to disgregard initial 50 steps
#samples=sampler.chain.reshape(-1) #to flatten array, if one parameter
#samples = sampler.chain[:, burnin:, :].reshape(-1) #to flatten array if more than one param
burnin = 50
samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
#samples = sampler.chain.reshape((-1, ndim))


fig = corner.corner(samples, labels=["$V_\mathrm{sini}$ [km/s]", "$\lambda_0$", "ln(f)"])
fig.savefig("line-triangle.pdf")

vsini_mcmc, l0_mcmc, lnf_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
#print(np.percentile(samples, [16, 50, 84],axis=0))
#vsini_mcmc, ldark_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],  axis=0)))
print("""MCMC result: vsini = {0[0]} +{0[1]} -{0[2]} l0 = {1[0]} +{1[1]} -{1[2]} lnf = {2[0]} +{2[1]} -{2[2]}""".format(vsini_mcmc, l0_mcmc, lnf_mcmc))



pl.figure()
for vsini, l0, lnf in samples[np.random.randint(len(samples), size=100)]: #size = 100
    mod = pyasl.fastRotBroad(x, c_template_shifted_rebinned_arr, 0.5, vsini, l0)
    pl.plot(x, mod, color=mine, alpha=0.1)
pl.plot(x, y, color="k", lw=1.2, alpha=0.8)
#pl.plot(x, c_template_shifted_rebinned_arr, 'b')
#pl.errorbar(x, y, yerr=yerr, fmt=".k")
pl.axis([6539, 6580, -1, 0.2])
pl.xlabel("$\lambda$ [A]")
pl.ylabel("Normalized counts")
pl.tight_layout()
pl.savefig("target_with_mcmcfit.pdf")



#find value of vsini where the histogram has a maximum
n, b, patches = pl.hist(samples[:,0], 100, range= [0,200], normed=1, histtype='stepfilled')
bin_max = np.where(n == n.max())
vsini_DEF = b[bin_max][0]
print("optimal vsini is: %f" % vsini_DEF)

f=open('vsini.dat', 'a')
f.write('%f\n' % vsini_DEF)




'''
#find value of vsini where the histogram has a maximum
n, b, patches = pl.hist(samples[:,1], 100, range= [0,1], normed=1, histtype='stepfilled')
bin_max = np.where(n == n.max())
ldark_DEF = b[bin_max][0]
print("limb darkening is: %f" % ldark_DEF)



#find value of vsini where the histogram has a maximum
n, b, patches = pl.hist(samples[:,2], 100, range= [CentralLambda_P1,CentralLambda_P2], normed=1, histtype='stepfilled')
bin_max = np.where(n == n.max())
l0_DEF = b[bin_max][0]
RV = ((l0_DEF-Halpha)/Halpha)*300000.
print("Central lambda is: %f" % l0_DEF)  #convert shift into velocity
print("RV shift is : %f" % RV)  #convert shift into velocity




    
'''








    
    
    
    



