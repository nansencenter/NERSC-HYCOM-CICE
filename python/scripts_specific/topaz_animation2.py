def main() :
   import netCDF4
   import datetime
   import numpy
   import matplotlib
   matplotlib.use('Agg')

   from mpl_toolkits.basemap import Basemap
   import matplotlib.pyplot
   from matplotlib import cm


   url="http://thredds.met.no/thredds/dodsC/topaz/dataset-topaz4-arc-myoceanv2-be"
   print url
   ncid = netCDF4.Dataset(url)
   lo = ncid.variables["longitude"][:]
   la = ncid.variables["latitude"][:]
   t  = ncid.variables["time"]
   print  ncid.variables["time"].__dict__

   dt0 = datetime.datetime(1950,1,1,0,0,0)
   t = [ dt0 + datetime.timedelta(hours=elem) for elem in t]
   print t


   rangelo = datetime.datetime.utcnow() - datetime.timedelta(days=10)
   rangehi = datetime.datetime.utcnow() + datetime.timedelta(days=20)

   I=numpy.where(numpy.logical_and(numpy.array(t)>rangelo,numpy.array(t)<rangehi))

   I=I[0]

   #m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='i')
   #m = Basemap(projection='mill',llcrnrlat=50,llcrnrlon=-60,urcrnrlat=80,urcrnrlon=40)
   m = Basemap(projection='stere',width=6000000,height=6000000,lon_0=0,lat_0=80,resolution="i")

   x,y=m(lo,la)

   for icnt,i in enumerate(I) :
   #for icnt,i in enumerate(I[:4]) :
      #print ncid.variables["temperature"]
      print t[i]

      sst  = ncid.variables["temperature"][i,0,:,:]
      cice = ncid.variables["fice"][i,:,:]
      cice = numpy.ma.masked_where(cice<.15,cice)

      fig = matplotlib.pyplot.figure(figsize=(16,12))
      ax=fig.add_axes([0,0,1,1])

      m.drawcoastlines()
      m.fillcontinents(color='.3',lake_color='aqua')
      m.drawparallels(numpy.arange(-80.,81.,20.))
      m.drawmeridians(numpy.arange(-180.,181.,20.))


      v=numpy.linspace(-2.,18,60,endpoint=True)
      v2=numpy.linspace(-2.,18,10,endpoint=True)
      CF=m.contourf(x,y,sst,v,extend="both")
      CF.set_clim(-2,18)
      CB=matplotlib.pyplot.colorbar(ticks=v2)

      ax.text(0.02,0.02,str(t[i]),fontsize=36,
         transform=ax.transAxes,
         bbox={'facecolor':'.7', 'alpha':0.9, 'pad':10})

      cmap = cm.get_cmap("gray")
      matplotlib.pyplot.hold(True)
      CF=m.contourf(x,y,cice,cmap=cmap)
      CF.set_clim(0,1)

      matplotlib.pyplot.title(str(t[i]))
      matplotlib.pyplot.savefig("tst%04d.png"%icnt)







if __name__ == "__main__" : 
   main()
