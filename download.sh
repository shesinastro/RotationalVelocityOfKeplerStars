#!/bin/bash

cd /Users/Serena/Desktop/ScriptsForEugene

for KOI in $(cat Kepler_KOI_new.txt); 
do 
    

curl -u repetto:R7XiXkn https://exofop.ipac.caltech.edu/kepler/edit_target.php?id=$KOI > kep.txt
tid=`grep -E -o ".{0}tid.{3,5}" kep.txt | head -1 | tr -dc '0-9'`
 
wget --user repetto --password R7XiXkn "https://exofop.ipac.caltech.edu/kepler/download_table.php?table=stellar&koi=$KOI&tid=$tid"
cat "download_table.p
id" >> "paramsKOI=$KOI.dat"
rm -r "download_table.php?table=stellar&koi=$KOI&tid=$tid"
   
done
