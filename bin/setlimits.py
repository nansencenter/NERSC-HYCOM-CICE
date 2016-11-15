#!/usr/bin/python -E
##!/usr/bin/env python
import forday 
import re
import sys
import os


# Get reference year and start/stop times
f = open("infile.in", "r")
text = f.readlines()
refyear=int(text[2].split()[0])
tmp=text[3].split() ; fday=int(tmp[0]); fhour=int(tmp[1]);
tmp=text[4].split() ; lday=int(tmp[0]); lhour=int(tmp[1]);
f.close()

# Get hycom yearflag
f = open("blkdat.input", "r")
p = re.compile(".*'yrflag'.*")
for line in f.readlines(): 
   if p.match(line) :
      yrflag=int(line.split()[0])
f.close()

#Test if file INITIALIZE exists
restart=1
if os.path.exists(r'INITIALIZE')  :
   restart=-1

# NB Add 1 to fday, ldat - they start from 0 but hycom
# starts from 1
fdtime=forday.dayfor(yrflag,refyear,fday+1,fhour)
ldtime=forday.dayfor(yrflag,refyear,lday+1,lhour)
#print "forday(fdtime) : ", forday.forday(yrflag,fdtime)
#print "forday(ldtime) : ", forday.forday(yrflag,ldtime)

f = open("limits", "w")
f.write( "%14.5f %14.5f" % ( (restart*fdtime), ldtime ))
f.close

## Dump number of segments
#f = open("NUMSEGMENTS", "w") 
#f.write(str(numseg))
#f.close


