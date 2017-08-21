# -*- coding: utf-8 -*-

import urllib2
from bs4 import BeautifulSoup
import os
import requests
import sys
from urlparse import urlparse


KOI=int(sys.argv[1]) 

# create a password manager
password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()

# Add the username and password.
# If we knew the realm, we could use it instead of None.
top_level_url = "https://exofop.ipac.caltech.edu/kepler/"
password_mgr.add_password(None, top_level_url, "repetto", "R7XiXkn")

handler = urllib2.HTTPBasicAuthHandler(password_mgr)

# create "opener" (OpenerDirector instance)
opener = urllib2.build_opener(handler)

# use the opener to fetch a URL
opener.open("https://exofop.ipac.caltech.edu/kepler/edit_target.php?id="+str(KOI))


#soup=BeautifulSoup(opener.open("https://exofop.ipac.caltech.edu/kepler/edit_target.php?id="+str(KOI)).read(), "lxml")
soup=BeautifulSoup(opener.open("https://exofop.ipac.caltech.edu/kepler/edit_target.php?id="+str(KOI)).read(), "lxml", from_encoding="utf-8")

     
        

rows = soup.findAll('tr')
#cols=rows.findAll('td')
j=0
for i in range(0, len(rows),1):
    cols=rows[i].findAll('td') 
    if len(cols)==6: 
        x1,x2,x3,x4,x5,x6= [c.text for c in cols]
        #if "FIES spectrum" in x4.replace(u'\xa0', u' ') or "TRES spectrum" in x4.replace(u'\xa0', u' '):
        if "TRES spectrum" in x4.replace(u'\xa0', u' ') or "TRES extracted spectra" in x4.replace(u'\xa0', u' '):
            name=x3
            username = 'repetto'
            password = 'R7XiXkn'
            print name
             
            if name.endswith("fits"): 

                url = 'https://exofop.ipac.caltech.edu/kepler/files/target'+str(KOI)+'/Spectra/'+str(name)
#filename = os.path.basename(urlparse(url).path)
                filename="sp"+str(j)+".fits"
                j=j+1

                r = requests.get(url, auth=(username,password))

                if r.status_code == 200:
                    with open(filename, 'wb') as out:
                        for bits in r.iter_content():
                            out.write(bits)    
    



            
            
            
    #print a
    #if "FIES spectrum" in a:
     #  print "ok"
    