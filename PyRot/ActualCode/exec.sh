#!/bin/sh

export PYTHONPATH=/usr/local/bin/

alias python="/usr/local/bin/python"

cd /Users/Serena/


source profile
source .bash_profile

cd /Users/Serena/Desktop/DEF_TRES


for KOI in $(cat KOI_TRES_third_part.dat); 
do 
echo $KOI
python 1_Delete.py
python -W ignore 2_DownloadSpectra_Multiple.py $KOI
if [ -f sp0.fits ]; then
echo $KOI >> existing.dat	
python -W ignore 3_FitsFile_Multiple.py
python 4_Normalize_Multiple.py
python 5_ChiSquared.py
python 6_MeasureRotBroad.py $KOI
fi

done