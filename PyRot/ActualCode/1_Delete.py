import glob
import os


dir = "/Users/Serena/Desktop/DEF_TRES"

for fitspath in glob.iglob(os.path.join(dir, '*.fits')):
    os.remove(fitspath)    
    
    
filename='target.dat'
if os.path.exists(filename):
    os.remove(filename)    
    
filename='target_norm.dat'
if os.path.exists(filename):
    os.remove(filename)        