from __future__ import print_function, division
import numpy as np
import matplotlib.pylab as plt
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl
import scipy.integrate as sci
import scipy.stats
from scipy.optimize import curve_fit


f = open("vsini_chi_sqr.txt", "w")


template=np.genfromtxt('SynthSpectrum_selected_norm_withErrors.dat', dtype=None, names=('l_template', 'N_template', 'error_template'))
target=np.genfromtxt('BroadenedSpectrum.dat', dtype=None, names=('l_target', 'N_target', 'error_target'))


midwave = 6562.76 #Halpha wavelength


ll_template_temp = template['l_template']
counts_template_temp = template['N_template']
ll_template=ll_template_temp[(ll_template_temp>6540) & (ll_template_temp< 6580)] #mask all the spectrum but small region around line of interest
counts_template = counts_template_temp[(ll_template_temp>6540) & (ll_template_temp< 6580)]


ll_target_temp = target['l_target']
counts_target_temp = target['N_target']
ll_target=ll_target_temp[(ll_target_temp>6540) & (ll_target_temp< 6580)]
counts_target = counts_target_temp[(ll_target_temp>6540) & (ll_target_temp< 6580)]
errors_target_temp = np.sqrt(abs(counts_target_temp)) #create new errors, dato che gli errori erano riferiti ancora ai count, non allo spettro normalizzato
errors_target = errors_target_temp[(ll_target_temp>6540) & (ll_target_temp< 6580)]




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


N_elements_template = (l_max_template-l_min_template)/delta_l_def
N_elements_target = (l_max_target-l_min_target)/delta_l_def



new_wave = np.arange(l_min_target, l_max_target, delta_l_def)


new_counts_template = interpol_template(new_wave)
new_counts_target = interpol_target(new_wave)
new_errors_target = interpol_error(new_wave)

# Obtain the broadened spectrum using
# 4th param = vsini; 3rd param = limb darkening

chi_squared = []
vsini_list = []


vsini_list = []
chi_squared = []
    
# Scipy chi^2
for vsini in np.logspace(-2, 2, 100):
    vsini_list.append(vsini)
    BrdSpectrum = pyasl.rotBroad(new_wave, new_counts_template, 0.5, vsini) #assume limb darkening of 0.5 as typical for Sun-like star
    chi_squared.append(scipy.stats.chisquare(new_counts_target,BrdSpectrum, axis=0)[-2])
        

min_chi_squared = min(chi_squared)   
index = chi_squared.index(min_chi_squared)  
vsini_index = vsini_list[index]


  
#vsini_array_sel = vsini_array[(vsini_array> vsini_index -10) & (vsini_array< vsini_index+10)] #define range where to do polynomial fit    
#res=curve_fit(thirdpol, vsini_list, chi_squared, bounds=(min(vsini_array_sel), max(vsini_array_sel)))
#res=curve_fit(thirdpol, vsini_list, chi_squared)


vsini_array = np.asarray(vsini_list)
chi_squared_array=np.asarray(chi_squared)
vsini_array_sel=vsini_array[(vsini_array>vsini_index-30) & (vsini_array<vsini_index+30)]
chi_squared_array_sel=chi_squared_array[(vsini_array>vsini_index-30) & (vsini_array<vsini_index+30)]

    
#fit_fun = thirdpol(vsini_array, 6.07754638e-05)
z=np.polyfit(vsini_array_sel, chi_squared_array_sel, 3)
p = np.poly1d(z)
fitfun = np.poly1d(z)
xp = np.linspace(0, 400, 100)
min_chi_test1 = scipy.optimize.minimize_scalar(p)

print("manual vsini: %d, vsini from cubic fitting: %f" % (vsini_index, min_chi_test1["x"]))
f.write(str(min_chi_test1["x"])) 

#iterate procedure, in interval around vsini
#guess1 = min_chi_test1["x"]
# polynimial fit only in small region around 



plt.plot(vsini_list, chi_squared, '.', xp, p(xp), '-')
plt.xlabel("vsini [km/s]")
plt.ylabel("Chi^2")
plt.savefig('test2.pdf')
plt.show() 
