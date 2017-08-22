Execute with:
$ /.exec.sh

——
1_Delete.py 
Clean directory of previous output files

2_DownloadSpectra_Multiple.py
Download spectra from a certain spectrograph for all the input KOIs

3_FitsFile_Multiple.py
Convert fits files into text files

4_Normalize_Multiple.py
Normalize each spectrum and take average

5_ChiSquared.py
Get vsini_guess through classical chi^2 minimization
Output: a file in which the vsini that mimimizes the chi^2 is stored


6_MeasureRotBroad.py
emcee routine to measure vsini
Output: EW of some of the models
        cloud plot
        time evolution of the walkers plot
        triangle plot

----
