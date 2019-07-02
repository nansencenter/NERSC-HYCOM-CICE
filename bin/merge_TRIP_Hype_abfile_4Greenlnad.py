from netCDF4 import Dataset
import argparse
import numpy as np
import abfile
import logging
import matplotlib
#

def greenlandmask(fld,lats,lons):
   #--------
   maskd=np.ones(lons.shape)
   maskd1= np.ma.masked_where(lats<57,maskd)
   maskd1= np.ma.masked_where(lats>68,maskd1)
   maskd1= np.ma.masked_where(lons<-57,maskd1)
   maskd1= np.ma.masked_where(lons>-26,maskd1)
   #
   maskd2= np.ma.masked_where(lats<68,maskd)
   maskd2= np.ma.masked_where(lats>75,maskd2)
   maskd2= np.ma.masked_where(lons<-62,maskd2)
   maskd2= np.ma.masked_where(lons>-5,maskd2)
   #
   maskd3= np.ma.masked_where(lats<75,maskd)
   maskd3= np.ma.masked_where(lats>83,maskd3)
   maskd3= np.ma.masked_where(lons<-75,maskd3)
   maskd3= np.ma.masked_where(lons>-5,maskd3)
   # mask for Greenl
   tmp_full=fld
   GR_mask1=tmp_full*maskd1
   GR_mask2=tmp_full*maskd2
   GR_mask3=tmp_full*maskd3
   GR_mask1[GR_mask1.mask]=0.0
   GR_mask2[GR_mask2.mask]=0.0
   GR_mask3[GR_mask3.mask]=0.0
   tot=GR_mask1 + GR_mask2 +GR_mask3
   return tot


# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False

if __name__ == "__main__" :
   class ClimParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values.split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       setattr(args, self.dest, tmp)
   class WindowParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values.split(",")
       tmp = [int(elem) for elem in tmp[0:4]]
       setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='')
  # parser.add_argument('--clim',       action=ClimParseAction,default=None)
  # parser.add_argument('--masklim',    action=ClimParseAction,default=None,help="mask limits")
  # parser.add_argument('--window',     action=WindowParseAction, help='firsti,firstj,lasti,lastj', default=None)
   parser.add_argument('Trip_river_afile',  help='Trip_river_afile')
   parser.add_argument('Hype_river_afile',  help='Hype_river_afile')
   parser.add_argument('griver_afile',  help='greenland_glacier_river_afile')

   #parser.add_argument('filename', help="")
   #parser.add_argument('records',  nargs="+",    type=int)
   #examplie
   # The first input is from ETRIP and the second one is from HYPE
   # we have added the third input corresponding to griver (Green land Glacier from /cluster/projects/nn2993k/TRIP/triver_Roshin/)
   #python ./merge_TRIP_Hype_abfile_4Greenlnad.py ./SCRATCH_ERAI-TRIP/rivers.a ./SCRATCH_hype_rev2/rev2_rivers.a

   args = parser.parse_args()


   gfile = abfile.ABFileGrid("regional.grid","r")
   scpx=gfile.read_field("scpx")
   scpy=gfile.read_field("scpy")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   #
   #
   #read lons from rivers.nc form A-Ehype
   #file_AEhype_river_nc="./SCRATCH_Ahype/rivers.nc"
   #fh_ahype = Dataset(file_AEhype_river_nc, mode='r')
   lons=plon
   lats=plat
   #lons = fh_ahype.variables['modlon'][:]
   #lats = fh_ahype.variables['modlat'][:]
   #fh_ahype.close()
   #----------------- read from ab file 
   #---- ERAI TRIP
   #file_TRIP_river="./SCRATCH_ERAI-TRIP/rivers.a"
   file_TRIP_river=args.Trip_river_afile
   afile_trip = abfile.AFile(800,760,file_TRIP_river,"r")
    # compute monthly flux:
   #tripafil_discharge_month = np.NaN * np.ones([12,])
   Greenland_trip_discharge_month = np.NaN * np.ones([12,])
   Greenland_trip_riverflux=np.zeros((12,760,800))
   tripafil_total_riverflux=np.zeros(lons.shape)
   countr=0
   for record in range(12) :
      trp_fld = afile_trip.read_record(record)
      Gr_trp_fld=greenlandmask(trp_fld,lats,lons)
      Greenland_trip_riverflux[record,:,:]=Gr_trp_fld
      ## convert from m^3/s km^3/year
      trp_fld=trp_fld*(365*86400.)*1e-9/12
      Gr_trp_fld=Gr_trp_fld*(365*86400.)*1e-9/12
      Greenland_trip_discharge_month[record] = np.nansum(Gr_trp_fld*scpx*scpy)
      #uns_discharge_month[record] = np.nansum(uns_riverflux[record,:,:]*scpx*scpy)
      tripafil_total_riverflux=tripafil_total_riverflux+trp_fld*scpx*scpy
      countr=countr+1
   
   
   # apply mask for Greenland
   tripafil_total_riverflux=greenlandmask(tripafil_total_riverflux,lats,lons)
   Maskout_only_Gr=np.ones(lons.shape)
   Maskout_only_Gr= np.ma.masked_where(tripafil_total_riverflux>0,Maskout_only_Gr)
   Maskout_only_Gr[Maskout_only_Gr.mask]=0.0
   
   width=190/25.4
   height = width/1.618
   print 'shapeafil of trip_total_riverflux = {}'.format(tripafil_total_riverflux.shape)
   print 'Tot annu. discharge from ETRIP-Greenland climatology is \
         {:.3f}'.format(np.sum(Greenland_trip_discharge_month))
   GrETRIP_sum="{:.3f}".format(np.sum(Greenland_trip_discharge_month))     
   #--------
   
   #-------------- Read from Ahype-Ehype
   #-- Ahype
   #file_AEhype_river="./SCRATCH_hype_rev2/rev2_rivers.a"
   file_AEhype_river=args.Hype_river_afile
   #afile_ahype = abfile.AFile(800,760,"./SCRATCH_Ahype/rivers.a","r")
   #afile_ahy = abfile.AFile(800,760,"./SCRATCH_Ahype/rivers.a","r")
   afile_ahy = abfile.AFile(800,760,file_AEhype_river,"r")
    # compute monthly flux:
   ahy_afil_total_riverflux=np.zeros(lons.shape)
   tot_hyp_trip_riverflux=np.zeros((12,760,800))
   countr=0
   for record in range(12) :
      ahyfld = afile_ahy.read_record(record)
      tot_hyp_trip_riverflux[record,:,:]=ahyfld*Maskout_only_Gr + Greenland_trip_riverflux[record,:,:] 
      ## convert from m^3/s km^3/year
      ahyfld=ahyfld*(365*86400.)*1e-9/12
      #
      ahy_afil_total_riverflux=ahy_afil_total_riverflux+ahyfld*scpx*scpy
      countr=countr+1
   
   #-------------- Read from griver
   
   file_griver=args.griver_afile
   afile_gr = abfile.AFile(800,760,file_griver,"r")
    # compute monthly flux:
   tot_griverflux=np.zeros((12,760,800))
   for record in range(12) :
      grfld = afile_gr.read_record(record)
      tot_griverflux[record,:,:]=grfld



  
 
   
   print "Number of records=", countr
   
   import matplotlib.pyplot as plt
   from mpl_toolkits.basemap import Basemap
   
   # Not necessary if coordinates are already in 2D arrays.
   cmap="jet"
   #-- create figure and axes instances
   dpi = 100
   figure = plt.figure(figsize=(1100/dpi, 1100/dpi), dpi=dpi)
   ax  = figure.add_axes([0.1,0.1,0.8,0.9])
   #-- create map
   # plot Ahype
   from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
   bm = Basemap(width=7400000,height=7400000,
         resolution='i',projection='stere',\
         lat_ts=70,lat_0=85,lon_0=-40.)
   
   #lon, lat = np.meshgrid(lons, lats)
   lon=lons
   lat=lats
   xi, yi = bm(lon, lat)
   print 'lon.shap=',lon.shape
   ahype=1
   if ahype :
      fig = plt.figure(figsize=(8,8))
      ax=fig.add_subplot(111)
      #cmap=plt.get_cmap("jet")
      #mld_avg=ahy_total_riverflu
      #for Greenland remove Unesco and add Trip data
      mld_avg=ahy_afil_total_riverflux*Maskout_only_Gr + tripafil_total_riverflux
      #mld_avg=ahy_afil_total_riverflux*Maskout_only_Gr 
      mld_avg = np.ma.masked_where(mld_avg<=0.0,mld_avg)
      mx=mld_avg.max()
      mn=mld_avg.min()
      clim=[mn, mx]
      clim=[0.0, 2.0]
      #clim=[0.2, 2.0]
      P=bm.pcolormesh(xi,yi,mld_avg,cmap=cmap)
      bm.drawcoastlines(linewidth=0.5)
      #bm.fillcontinents(color='.8',lake_color='white')
      #bm.drawparallels(np.arange(-80.,81.,20.))
      #bm.drawmeridians(np.arange(-180.,181.,20.))
      #bm.drawmapboundary() #fill_color='aqua')
      #
      aspect = 20
      pad_fraction = 0.5
      divider = make_axes_locatable(ax)
      width = axes_size.AxesY(ax, aspect=1./aspect)
      pad = axes_size.Fraction(pad_fraction, width)
      cax = divider.append_axes("right", size=width, pad=pad)
      cb=ax.figure.colorbar(P,cax=cax,extend='both')
      P.set_clim(clim)
      ax.set_title('Rev Ahype+EHYPE clim total river flux')
      fig.canvas.print_figure("Rev2_Green_Fig_Ahype_totalflux%03d.png"%(record+9))
      ax.clear()
      cb.remove()
   
   xticks = np.arange(0,12)
   xticklabels = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec') 
   #
   tripafil=1
   if tripafil :
      fig = plt.figure(figsize=(8,8))
      ax=fig.add_subplot(111)
      #cmap=plt.get_cmap("jet")
      mld_avg=tripafil_total_riverflux
      mld_avg = np.ma.masked_where(mld_avg<=0.0,mld_avg)
      mx=mld_avg.max()
      mn=mld_avg.min()
      clim=[mn, mx]
      clim=[0.0, 2.0]
      P=bm.pcolormesh(xi,yi,mld_avg,cmap=cmap)
      bm.drawcoastlines(linewidth=0.5)
      #bm.fillcontinents(color='.8',lake_color='white')
      #bm.drawparallels(np.arange(-80.,81.,20.))
      #bm.drawmeridians(np.arange(-180.,181.,20.))
      #bm.drawmapboundary() #fill_color='aqua')
      #
      #x, y = bm(lon, lat)
      #bm.plot(x, y, 'o-', markersize=5, linewidth=1)
      aspect = 20
      pad_fraction = 0.5
      divider = make_axes_locatable(ax)
      width = axes_size.AxesY(ax, aspect=1./aspect)
      pad = axes_size.Fraction(pad_fraction, width)
      cax = divider.append_axes("right", size=width, pad=pad)
      cb=ax.figure.colorbar(P,cax=cax,extend='both')
      P.set_clim(clim)
      ax.set_title('Rev2_Greenland-ETRIP total annual discharge='+GrETRIP_sum)
      fig.canvas.print_figure("GreenFig_trippeafil_totalflux%03d.png"%(record+9))
      ax.clear()
      cb.remove()
      plt.close()
   
   
   #Now write to ab files the tot river flux
   #print '-----------------'
   #write to ab file
   af = abfile.AFile(lons.shape[1],lons.shape[0],"merged_tot_rev2_rivers.a","w")
   bf = open("merged_tot_rev2_rivers.b","w")
   bf.write("Tot. river mass flux from Hype + ETRIP(Greenland)  \n")
   bf.write("\n")
   bf.write("\n")
   bf.write("\n")
   bf.write("i/jdm =  %5d %5d\n"%(lons.shape[1],lons.shape[0]))
   for recrd in range(12) : 
      dummy_sum=tot_griverflux[recrd,:,:]+tot_hyp_trip_riverflux[recrd,:,:]
      hmin,hmax = af.writerecord(dummy_sum,None,record=recrd)
      bf.write(" rivers:month,range = %2.2i%16.8e%16.8e\n"%((recrd+1),hmin,hmax))
   af.close()
   bf.close()
   
    ## convert from m^3/s km^3/year
   
  
  
  
