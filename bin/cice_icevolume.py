#!/usr/bin/env python
import netCDF4
import numpy
import matplotlib.pyplot
import abfile
from mpl_toolkits.basemap import Basemap
import sys
import argparse
import glob
import cfunits
import datetime



def plot_fig(m,x,y,fld,filename,title,cmap=None,clim=None) :
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   if cmap is not None :
      P=m.pcolormesh(x,y,fld,cmap=cmap)
   else :
      P=m.pcolormesh(x,y,fld)
   m.drawcoastlines()
   m.fillcontinents(color='coral',lake_color='aqua')
   # draw parallels and meridians.
   m.drawparallels(numpy.arange(-80.,81.,20.))
   m.drawmeridians(numpy.arange(-180.,181.,20.))
   m.drawmapboundary(fill_color='aqua')
   ax.figure.colorbar(P)
   ax.set_title(title)
   if clim is not None : P.set_clim(clim)
   figure.canvas.print_figure(filename)



def main(myglobs) :

   # Read cice files, plot volume
   reg = abfile.ABFileGrid("regional.grid","r")
   plon=reg.read_field("plon")
   plat=reg.read_field("plat")
   scpx=reg.read_field("scpx")
   scpy=reg.read_field("scpy")

   figure1 = matplotlib.pyplot.figure(figsize=(8,8))
   ax1=figure1.add_subplot(111)
   ax1.set_title("Total Ice Area [ 1.000.000 km^2 ]")
   ax1.grid(True)

   figure2 = matplotlib.pyplot.figure(figsize=(8,8))
   ax2=figure2.add_subplot(111)
   ax2.set_title("Total Ice Volume [ km^3 ]")
   ax2.grid(True)


   for myglob in myglobs :
      print myglob
      files =glob.glob(myglob)


      l_area=[]
      l_vol =[]
      l_time=[]

      for file in files :

         print file

         nc = netCDF4.Dataset(file,"r") 
         aice=nc.variables["aice"][0,:,:]
         hice=nc.variables["hi"][0,:,:]
         newt=nc.variables["time"][0]

         vol=numpy.sum(aice*hice*scpx*scpy)
         area=numpy.sum(aice*scpx*scpy)

         t_unit = cfunits.Units(nc.variables["time"].units)
         my_t_unit = cfunits.Units('days since 1900-1-1')
         newt = cfunits.Units.conform(newt,t_unit,my_t_unit)
         newt = int(newt*86400.)
         newdt= datetime.datetime(1900,1,1,0,0,0) + datetime.timedelta(seconds=newt)
         #print file, "%14.6gm**3    %14.6gm**2" % (vol,area)
         nc.close()

         l_area.append(area)
         l_vol.append(vol)
         l_time.append(newdt)


      #  Sort
      I = sorted(range(len(l_time)),key=lambda x:l_time[x])
      l_area = [l_area[i] for i in I]
      l_vol  = [l_vol [i] for i in I]
      l_time = [l_time[i] for i in I]

      ax1.plot(l_time,numpy.array(l_area)*1e-12,label=myglob,lw=3)
      ax2.plot(l_time,numpy.array(l_vol )*1e-12,label=myglob,lw=3)


   ax1.legend()
   figure1.canvas.print_figure("icearea.png")

   ax2.legend()
   figure2.canvas.print_figure("icevolume.png")

if __name__ == "__main__" :
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('globs',     nargs='+')
   args = parser.parse_args()
   main(args.globs)
