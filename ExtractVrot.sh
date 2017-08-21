#!/bin/bash

cd /Users/Serena/Desktop/ScriptsForEugene

rm -r vrot_DEF.dat


for KOI in $(cat Kepler_KOI.txt);  
do
	
	cp paramsKOI=$KOI.dat params.dat
	python charachter.py >> tempfile.txt
	var=`cat tempfile.txt`
     cat paramsKOI=$KOI.dat | awk -v a=$var '{print substr($0,'a',15)}' >> vrot.dat
	 python charachter1.py $KOI
     rm -r tempfile.txt
	 rm -r vrot.dat


done
