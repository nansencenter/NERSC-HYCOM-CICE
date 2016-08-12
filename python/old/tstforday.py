#!/usr/bin/python
import forday 

print forday.forday(3,36264)
print forday.forday(3,3)
print forday.forday(3,365)
print forday.forday(3,366)
print forday.forday(3,367)
print forday.forday(3,368)


if (0==1) :
   yrflag=3
   for yrflag in range(0,3) :
      for i in range(1,1000) :
         dtime=float(i)
         ( iyear, iday, ihour) =  forday.forday(yrflag,dtime)
         dtime2 =  forday.dayfor(yrflag, iyear, iday, ihour)
         if (dtime != dtime2 ) :
            print dtime, dtime2
         
         


