from string import split
f=open('params.dat', 'r')

for lines in f.readlines():
    
    if 'Vsini' in lines:
        a=len(split(lines, 'Vsini')[0])
        print a
