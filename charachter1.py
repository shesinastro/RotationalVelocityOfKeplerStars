from string import split
import re
import numpy as np
import sys

f=open('vrot.dat', 'r')
g=open('vrot_DEF.dat', 'a')

KOI=sys.argv[1]
listone=[]
listerror=[]

b=0

for lines in f.readlines():
    
    if 'Vsini' not in lines:
        
        if '+/-' in lines:
            a=split(lines, '+/-')

            
            
            val = a[0]
            val2 = a[1]
            
            if val != '  ':
                listone.append(float(a[0]))
                try:
                    listerror.append(float(a[1]))  
                except ValueError:
                    listerror.append(1) 
                    
                    
                    
                    #appendo un valore tipico dell'incertezza su vsini quando il valore e vuoto
                
               
               



if len(listone) != 0:
    average=sum(listone)/len(listone)
    error_on_average=np.sqrt(sum(np.square(listerror)))/len(listerror)
    g.write('%s\t%s\t%s\n' % (average, error_on_average, KOI))
                








