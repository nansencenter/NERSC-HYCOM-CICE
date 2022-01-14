#!/usr/bin/env python
import modeltools.hycom
import modeltools.tools
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import abfile
import numpy
from mpl_toolkits.basemap import Basemap
import netCDF4
import logging
import re
import os.path
import gridxsec
import myMOD 

# Set up logger
_loglevel=logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False

'''
usage: example
   python ./hycom_plot_section_TOPAZ_Avg.py --section_map --clim=30,35.6    --sectionid='Fram_Strait' -19 79 12 79 'salin' ./TP6archv.2021_{001..030}_12.a
   python ./hycom_plot_section_TOPAZ_Avg.py --section_map --clim=-1,6    --sectionid='Fram_Strait' -19 79 12 79 'temp' ./TP6archv.2021_120_{00..03}.a
'''


def main(lon1,lat1,lon2,lat2,variable,files,filetype="archive",clim=None,sectionid="",
      ijspace=False,xaxis="distance",section_map=False,ncfiles="",dpi=180) :
   #TP4Grd='/cluster/work/users/aal069/TP4a0.12/mfile/'
   logger.info("Filetype is %s"% filetype)
   jip='/cluster/work/users/xiejp/work_2017/TP2_EXP/'
   if 'TP2' in files[0]: 
      gfile = abfile.ABFileGrid(jip+"regional.grid","r")
   else:
      gfile = abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   qlon=gfile.read_field("qlon")
   qlat=gfile.read_field("qlat")

   # Set up section info
   if ijspace :
      sec = gridxsec.SectionIJSpace([lon1,lon2],[lat1,lat2],plon,plat)
   else  :
      sec = gridxsec.Section([lon1,lon2],[lat1,lat2],plon,plat)
   I,J=sec.grid_indexes
   dist=sec.distance
   print 'dit.shae=',dist.shape
   slon=sec.longitude
   slat=sec.latitude

   # In testing
   #J,I,slon,slat,case,dist=sec.find_intersection(qlon,qlat)
   #print I,J
   #raise NameError,"test"



   logger.info("Min max I-index (starts from 0):%d %d"%(I.min(),I.max()))
   logger.info("Min max J-index (starts from 0):%d %d"%(J.min(),J.max()))


   if section_map :
      ll_lon=slon.min()-10.
      ur_lon=slon.max()+10.
      ll_lat=numpy.maximum(-90.,slat.min()-10.)
      ur_lat=numpy.minimum(90. ,slat.max()+10.)
      m = Basemap(projection='mill', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, resolution='l')
      (x,y) = m(slon,slat)
      figure = matplotlib.pyplot.figure()
      ax=figure.add_subplot(111)
      m.drawcoastlines()
      #m.fillcontinents(color='coral',lake_color='aqua')
      m.drawparallels(numpy.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
      m.drawmeridians(numpy.arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
      m.drawmapboundary() # draw a line around the map region
      m.plot(x,y,"r",lw=3)
      m.etopo()
      #m.scatter(x,y,s=20,c=dist)
      pos = ax.get_position()
      #print pos
      asp=pos.height/pos.width
      #print asp
      w=figure.get_figwidth()
      #print w
      h=asp*w
      figure.set_figheight(h)
      if sectionid :
         figure.canvas.print_figure("map_%s.png"%sectionid,dpi=dpi)
      else :
         figure.canvas.print_figure("map.png",dpi=dpi)

   # Get layer thickness variable used in hycom
   dpname = modeltools.hycom.layer_thickness_variable[filetype]
   logger.info("Filetype %s: layer thickness variable is %s"%(filetype,dpname))


   if xaxis == "distance" :
      x=dist/1000.
      xlab="Distance along section[km]"
   elif xaxis == "i" :
      x=I
      xlab="i-index"
   elif xaxis == "j" :
      x=J
      xlab="j-index"
   elif xaxis == "lon" :
      x=slon
      xlab="longitude"
   elif xaxis == "lat" :
      x=slat
      xlab="latitude"
   else :
      logger.warning("xaxis must be i,j,lo,lat or distance")
      x=dist/1000.
      xlab="Distance along section[km]"

   # get kdm from the first file:
   # Remove [ab] ending if present
   print 'firstfilw', files[0]
   m=re.match("(.*)\.[ab]",files[0])
   print 'm=',m.group(1)
   myf=m.group(1)
   fi_abfile = abfile.ABFileArchv(myf,"r")
   kdm=max(fi_abfile.fieldlevels)

   # Loop over archive files
   figure = matplotlib.pyplot.figure()
   ax=figure.add_subplot(111)
   pos = ax.get_position()
   count_sum=0
   intfsec_sum=numpy.zeros((kdm+1,I.size))
   datasec_sum=numpy.zeros((kdm+1,I.size))
   for fcnt,myfile0 in enumerate(files) :
      count_sum=count_sum+1
      print 'count_sum==', count_sum
      print 'fcnt=', fcnt
      print 'mfile0=', myfile0
      # Remove [ab] ending if present
      m=re.match("(.*)\.[ab]",myfile0)
      if m :
         myfile=m.group(1)
      else :
         myfile=myfile0
      # Add more filetypes if needed. By def we assume archive
      if filetype == "archive" :
         i_abfile = abfile.ABFileArchv(myfile,"r")
      elif filetype == "restart" :
         i_abfile = abfile.ABFileRestart(myfile,"r",idm=gfile.idm,jdm=gfile.jdm)
      else :
         raise NotImplementedError,"Filetype %s not implemented"%filetype
      # kdm assumed to be max level in ab file
      kdm=max(i_abfile.fieldlevels)

      # Set up interface and daat arrays
      
      xx=numpy.zeros((kdm+1,I.size))
      intfsec=numpy.zeros((kdm+1,I.size))
      datasec=numpy.zeros((kdm+1,I.size))
      # Loop over layers in file. 
      logger.info("File %s"%(myfile))
      for k in range(kdm) :
         logger.debug("File %s, layer %03d/%03d"%(myfile,k,kdm))

         # Get 2D fields
         dp2d=i_abfile.read_field(dpname,k+1)
         data2d=i_abfile.read_field(variable,k+1)
         #print('---mn,mx  data=',  data2d.min(),data2d.max())
         if (k%kdm==49):
            print("---Reach bottom layer" )
         dp2d=numpy.ma.filled(dp2d,0.)/modeltools.hycom.onem
         data2d=numpy.ma.filled(data2d,1e30)
         # Place data into section arrays
         intfsec[k+1,:] = intfsec[k,:] + dp2d[J,I]
         if k==0 : datasec[k,:] = data2d[J,I]
         datasec[k+1,:] = data2d[J,I]
      

      intfsec_sum=intfsec_sum + intfsec
      datasec_sum=datasec_sum + datasec
      #print 'prs_intafce=', numpy.transpose(intfsec[:,15]) 
      i_abfile.close()

      # end loop over files
      
   intfsec_avg=intfsec_sum/count_sum
   datasec_avg=datasec_sum/count_sum

   if ncfiles :
      MLDGS_sum=numpy.zeros((1,I.size))
      count_sum=0
      for fcnt,ncfile in enumerate(ncfiles) :
         count_sum=count_sum+1
         print 'ncfile count_sum==', count_sum
         print 'ncfile fcnt=', fcnt
         print 'ncfilefile=', ncfile
         MLDGS=numpy.zeros((1,I.size))
         ncfile0 = netCDF4.Dataset(ncfile,'r')
         MLD_2D  = ncfile0.variables['GS_MLD'][:]
         #MLD_2D  = ncfile0.variables['mlp'][:]
         MLDGS[0,:]=MLD_2D[0,J,I]
         MLDGS_sum= MLDGS_sum + MLDGS
         ncfile0.close()
      # end loop over files
      MLDGS_avg=MLDGS_sum/count_sum
   #
   #-----------------------------------------------------------------
   # read from clim mld TP5netcdf
   if ncfiles :
      if 'TP2' in files[0]:
         fh=netCDF4.Dataset('mld_dr003_l3_modif_Interp_TP2grd.nc')
      else:
         fh=netCDF4.Dataset('mld_dr003_l3_modif_Interp_TP5grd.nc')
      fhmldintrp = fh.variables['TP5mld'][:]
      fh.close()
      #fhMLDintrp_sum=np.zeros((760,800))
      MLDclim_sum=numpy.zeros((1,I.size))
      cunt_sum=0
      for ii in range(12) :
          cunt_sum=cunt_sum +1
          MLDclim=numpy.zeros((1,I.size))
          MLDclim[0,:]=fhmldintrp[ii,J,I]

          MLDclim_sum= MLDclim_sum + MLDclim
          print 'clim count_sum==', cunt_sum
      MLDclim_avg=MLDclim_sum/cunt_sum
   #-----------------------------------------------------------------   
   i_maxd=numpy.argmax(numpy.abs(intfsec_avg[kdm,:]))
   #print i_maxd
   for k in range(kdm+1) :
      xx[k,:] = x[:]
   # Set up section plot
   #datasec = numpy.ma.masked_where(datasec==1e30,datasec)
   datasec_avg = numpy.ma.masked_where(datasec_avg>0.5*1e30,datasec_avg)
   #print datasec.min(),datasec.max()
   #P=ax.pcolormesh(dist/1000.,-intfsec,datasec)
   #print i_maxd
   for k in range(kdm+1) :
      xx[k,:] = x[:]
   
   if clim is not None : lvls = MaxNLocator(nbins=30).tick_values(clim[0], clim[1])
   #print 'levels=', lvls
   mf='sawtooth_0-1.txt'
   LinDic=myMOD.cmap_dict(mf)
   my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',LinDic)
   cmap=my_cmap
   #cmap = matplotlib.pyplot.get_cmap('gist_rainbow_r')
   norm = BoundaryNorm(lvls, ncolors=cmap.N, clip=True)
   print 'x.shape=' ,      x.shape
   print 'x.min,xmax=' ,  x.min(),x.max()
   print 'xx.shape=' ,      xx.shape
   print 'xx.min,xxmax=' ,  xx.min(),xx.max()
   print 'intfsec_avg.shape=', intfsec_avg.shape
   print 'datasec_avg.shape=', datasec_avg.shape
   #P=ax.pcolormesh(x,-intfsec,datasec,cmap=cmap)
   P=ax.contourf(xx,-intfsec_avg,datasec_avg,extend='both',cmap=cmap,levels=lvls)
   if 'sal' in variable:
      P1=ax.contour(xx,-intfsec_avg,datasec_avg,levels=[32.0,33.0,34.0,35.0,35.5],
          colors=('k',),linestyles=('-',),linewidths=(1.5,))
   else:
      P1=ax.contour(xx,-intfsec_avg,datasec_avg,levels=[-1,0.0,2.0],
          colors=('k',),linestyles=('-',),linewidths=(1.5,))
   matplotlib.pyplot.clabel(P1, fmt = '%2.1d', colors = 'k', fontsize=10) #contour line labels
   # Plot layer interfaces
   for k in range(1,kdm+1) :
      if k%100 == 0 : 
         PL=ax.plot(x,-intfsec_avg[k,:],"-",color="k")
      elif k%5 == 0 and k <= 10: 
         PL=ax.plot(x,-intfsec_avg[k,:],"--",color="k", linewidth=0.5)
         textx = x[i_maxd]
         texty = -0.5*(intfsec_avg[k-1,i_maxd] + intfsec_avg[k,i_maxd])
         ax.text(textx,texty,str(k),verticalalignment="center",horizontalalignment="center",fontsize=6)
      elif k%2 and k > 10 : 
         PL=ax.plot(x,-intfsec_avg[k,:],"--",color="k", linewidth=0.5)
         textx = x[i_maxd]
         texty = -0.5*(intfsec_avg[k-1,i_maxd] + intfsec_avg[k,i_maxd])
         ax.text(textx,texty,str(k),verticalalignment="center",horizontalalignment="center",fontsize=6)
   if ncfiles :
      PL=ax.plot(x,-MLDGS_avg[0,:],"-",color="w", linewidth=1.50)
      PL=ax.plot(x,-MLDclim_avg[0,:],"--",color="r", linewidth=1.50)
###    else :
###       PL=ax.plot(x,-intfsec_avg[k,:],"-",color=".5")
  # Print figure and remove wite space.
   aspect = 50
   pad_fraction = 0.25
   divider = make_axes_locatable(ax)
   width = axes_size.AxesY(ax, aspect=1./aspect)
   print 'width=',width
   pad = axes_size.Fraction(pad_fraction, width)
   cax = divider.append_axes("right", size=width, pad=pad)
   cb=ax.figure.colorbar(P,cax=cax,extend='both')
   #cb=ax.figure.colorbar(P,extend='both')
   if clim is not None : P.set_clim(clim)
   #cb=ax.figure.colorbar(P,extend='both')
   ax.set_title(variable+':'+myfile+'AVG-')
   ax.set_ylabel('Depth [m]')
   ax.set_xlabel(xlab)
   #ax.set_position(pos)
   #matplotlib.pyplot.tight_layout()

   # Print in different y-lims 
   suff=os.path.basename(myfile)
   if sectionid : suff=suff+"_"+sectionid
   figure.canvas.print_figure("sec_AVG_%s_full_%s.png"%(variable,suff),dpi=dpi)
   #ax.set_ylim(-1000,0)
   if 'Fram' in sectionid or 'Svin' in sectionid:
      print 'sectionid=', sectionid
      ax.set_ylim(-600,0)
      figure.canvas.print_figure("sec_AVG_%s_600m_%s.png"%(variable,suff),dpi=dpi)
   else:
      #ax.set_ylim(-2500,0)
      #figure.canvas.print_figure("sec_AVG_%s_2500m_%s.png"%(variable,suff),dpi=dpi)
      ax.set_ylim(-3000,0)
      figure.canvas.print_figure("sec_AVG_%s_3000m_%s.png"%(variable,suff),dpi=dpi)

   # Close input file
   #i_abfile.close()
   #
   ax.clear()
   cb.remove()


if __name__ == "__main__" :
   class ClimParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values.split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       setattr(args, self.dest, tmp)


   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--dpi',     type=int,default=180)
   parser.add_argument('--clim',     action=ClimParseAction,default=None)
   parser.add_argument('--filetype'    ,     type=str, help='',default="archive")
   parser.add_argument('--ij'      , action="store_true",default=False)
   parser.add_argument('--sectionid'      , type=str,default="") 
   parser.add_argument('--xaxis'          , type=str,default="distance") 
   parser.add_argument('--section_map'    , action="store_true",default=False,help="Produces a simple map of the section")
   parser.add_argument('--ncfiles'      , nargs="+")
   parser.add_argument('lon1',     type=int, help='')
   parser.add_argument('lat1',     type=int, help='')
   parser.add_argument('lon2',     type=int, help='')
   parser.add_argument('lat2',     type=int, help='')
   parser.add_argument('variable', type=str, help='')
   parser.add_argument('files',nargs="+")

   args = parser.parse_args()

   main(args.lon1,args.lat1,args.lon2,args.lat2,args.variable,args.files,filetype=args.filetype,clim=args.clim,sectionid=args.sectionid,ijspace=args.ij,
         xaxis=args.xaxis,section_map=args.section_map,ncfiles=args.ncfiles,
         dpi=args.dpi) 
