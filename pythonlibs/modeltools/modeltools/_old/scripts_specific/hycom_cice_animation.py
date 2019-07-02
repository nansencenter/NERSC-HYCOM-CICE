import sys

def main(files) :
   import netCDF4
   import datetime
   import numpy
   import matplotlib
   import cfunits
   matplotlib.use('Agg')

   from mpl_toolkits.basemap import Basemap
   import matplotlib.pyplot
   from matplotlib import cm


   # lo/la from first datafile
   ncid = netCDF4.Dataset(files[0])
   lo = ncid.variables["TLON"][:]
   la = ncid.variables["TLAT"][:]


   #m = Basemap(projection='stere',width=6000000,height=6000000,lon_0=0,lat_0=80,resolution="i")
   m = Basemap(projection='stere',width=8000000,height=7500000,lon_0=-45,lat_0=85,resolution="i")
   x,y=m(lo,la)

   ncid.close()


   icnt=0
   for file in files :
      icnt+=1
      ncid = netCDF4.Dataset(file)
      t  = ncid.variables["time"]

      t_unit = cfunits.Units(t.units)
      my_t_unit = cfunits.Units('days since 1900-1-1')
      newt = cfunits.Units.conform(t,t_unit,my_t_unit)
      newt = [int(elem*86400.) for elem in newt]
      #print newt[0]
      newdt= [ datetime.datetime(1900,1,1,0,0,0) + datetime.timedelta(seconds=elem) for elem in newt]
      print newdt[0]


      sst  = ncid.variables["sst"][0,:,:]
      cice = ncid.variables["hi"][0,:,:]
      cice = numpy.ma.masked_where(cice<.05,cice)

      fig = matplotlib.pyplot.figure(figsize=(16,12))
      ax=fig.add_axes([0,0,1,1])

      m.drawcoastlines()
      m.fillcontinents(color='.3',lake_color='aqua')
      m.drawparallels(numpy.arange(-80.,81.,20.))
      m.drawmeridians(numpy.arange(-180.,181.,20.))


      v=numpy.linspace(-2.,18,60,endpoint=True)
      v2=numpy.linspace(-2.,18,10,endpoint=True)
      CF=m.contourf(x,y,sst,v,extend="both")
      CF.set_clim(0,18)
      CB=matplotlib.pyplot.colorbar(ticks=v2)

      ax.text(0.02,0.02,str(newdt[0]),fontsize=36,
         transform=ax.transAxes,
         bbox={'facecolor':'.7', 'alpha':0.9, 'pad':10})

      cmap = cm.get_cmap("YlGn")
      matplotlib.pyplot.hold(True)
      v3=numpy.linspace(.1,3,6,endpoint=True)
      CF=m.contourf(x,y,cice,v3,cmap=cmap,extend="both")
      CF.set_clim(0,3)
      CB2=matplotlib.pyplot.colorbar(ticks=v3)

      matplotlib.pyplot.title(str(newdt[0]))
      matplotlib.pyplot.savefig("tst%04d.png"%icnt)

      matplotlib.pyplot.close()







if __name__ == "__main__" : 
   main(sys.argv[1:])
