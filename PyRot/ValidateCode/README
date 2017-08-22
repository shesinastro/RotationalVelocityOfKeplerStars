Execute as:<br />

$BroadenSpectrum.py $vsini $sigma_noise

where $vsini is the input vsini<br />
where $sigma_noise is the sigma of the Gaussian to simulate noise<br />

Input: SynthSpectrum_selected_norm.dat

Output: Broadened+noised synthetic spectrum (synthetic spectrum with errors as rms)



$python RotBroad_ChiSquared_Iterative.py <br />
Get vsini_guess through chi^2 minimization<br />
Input: SynthSpectrum_selected_norm_withErrors.dat<br />
       BroadenedSpectrum.dat<br />
Output: a file in which the vsini that mimimizes the chi^2 is stored <br />


$python TestRealTemplate_LimbDark_CentralLambda_ToZero_BootstrapPrior_NewErrors_FixedLimb_UnderError.py<br />
emcee routine<br />
Output: EW of some of the models<br />
        cloud plot<br />
        time evolution of the walkers plot<br />
        triangle plot<br />


