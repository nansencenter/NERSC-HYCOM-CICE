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
Usage
python ./hycom_plot_section_TOPAZ_Avg_relax.py --clim=-2.0,4    --sectionid='Fram_Strait_relax' -19 79 12 79 'tem' ../../relax/020/relax_tem.a
python ./hycom_plot_section_TOPAZ_Avg_relax.py --clim=30,35.6    --sectionid='Fram_Strait_TP5_relax' -19 79 12 79 'sal' ../../relax/020/relax_sal.a

'''


def main(lon1,lat1,lon2,lat2,variable,files,filetype="archive",clim=None,sectionid="",
      ijspace=False,xaxis="distance",section_map=False,dpi=180) :
   TP4Grd='/cluster/work/users/aal069/TP4a0.12/mfile/'
   logger.info("Filetype is %s"% filetype)
   if 'TP4' in files[0]:
      gfile = abfile.ABFileGrid(TP4Grd+"regional.grid","r")
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
   print 'dist.shae=',dist.shape
   print("I.shape=", I.size)
   print("J.shape=", J.shape)
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
   myfile0=files[0]
   print 'myfile0', myfile0

   m=re.match("(.*)\.[ab]",myfile0)
   print 'm=',m.group(1)
   if m :
      myfile=m.group(1)
   else :
      myfile=myfile0
   dta_afile = abfile.AFile(gfile.idm,gfile.jdm,myfile0,"r")
   #intfl='../../../relax/010/relax_int.a'
   intfl=myfile[:-3] + 'int.a'
   int_afile = abfile.AFile(gfile.idm,gfile.jdm,intfl,"r")
   #
   lyr=1
   #record_num,xmn,xmx=dta_afile.get_record_number(variable,lyr)
   #print 'record_num, variable, layer ===', record_num-1, variable, lyr
   #print 'mn,mx=',xmn,xmx
   # for record in record_num :
   record_num=1
   record_var=record_num-1 
   fld = dta_afile.read_record(record_var)
   print 'mn,mx  data=',fld.min(),fld.max()
   #numpy.testing.assert_approx_equal(xmn,fld.min(),significant=8)
   #numpy.testing.assert_approx_equal(xmx,fld.max(),significant=8)
   # presure interfc
   #record_num,xmn,xmx=int_afile.get_record_number('int',lyr)
   #print 'record_num, variable, layer ===', record_num-1, 'int', lyr
   #print 'mn,mx=',xmn,xmx
   # for record in record_num :
   record_prs=record_num-1 
   fld = int_afile.read_record(record_prs)
   print 'mn,mx  intface=',fld.min(),fld.max()
   #numpy.testing.assert_approx_equal(xmn,fld.min(),significant=8)
   #numpy.testing.assert_approx_equal(xmx,fld.max(),significant=8)
   #
   #kdm=max(fi_abfile.fieldlevels)
   kdm=50
   # Loop over archive files
   figure = matplotlib.pyplot.figure()
   ax=figure.add_subplot(111)
   pos = ax.get_position()
   count_sum=0
   intfsec_sum=numpy.zeros((kdm+1,I.size))
   datasec_sum=numpy.zeros((kdm+1,I.size))
   #
   for mnth in range(12) :
      count_sum=count_sum+1

      logger.info("Reading data and presure interface for month %s"%(mnth+1))
      record_var_pnt=record_var + mnth*kdm
      record_prs_pnt=record_prs + mnth*kdm
      print('pointing at record num:  record_var_pnt', record_var_pnt)
      print('pointing at record num:  record_prs_pnt', record_prs_pnt)
      # Set up interface and daat arrays
      xx=numpy.zeros((kdm+1,I.size))
      intfsec=numpy.zeros((kdm+1,I.size))
      datasec=numpy.zeros((kdm+1,I.size))
      # Loop over layers in file. 
      logger.info("File %s"%(myfile))
      logger.info("intfac_File %s"%(intfl))
      #loop over 50 records
      for k in range(kdm) :
         logger.debug("File %s, layer %03d/%03d"%(myfile,k,kdm))
         r_var=record_var_pnt+k
         r_prs=record_prs_pnt+k
         # Get 2D fields
         dp2d = int_afile.read_record(r_prs)
         data2d=dta_afile.read_record(r_var)
         #print('reading rec num: r_var,r_dp=', r_var, r_prs,' ','data.min,data.max=',data2d.min(),data2d.max())
         print('data: month,layer, range', (mnth+1),'',(r_var)%kdm,data2d.min(),data2d.max())
         #dp2d=i_abfile.read_field(dpname,k+1)
         #data2d=i_abfile.read_field(variable,k+1)
         if ((r_var)%kdm==49):
            print('mn,mx  intface=',dp2d.min(),dp2d.max())
            print('mn,mx  data=',  data2d.min(),data2d.max())
            print( "Reach bottom layer" )
         dp2d=numpy.ma.filled(dp2d,0.)/modeltools.hycom.onem
         data2d=numpy.ma.filled(data2d,1e30)

         # Place data into section arrays
         #intfsec[k+1,:] = intfsec[k,:] + dp2d[J,I]
         #print("data2d.shape=",data2d.shape)
         #print("data2d[J,I].size=",data2d[J,I].size)
         intfsec[k+1,:] = dp2d[J,I]
         if k==0 : datasec[k,:] = data2d[J,I]
         datasec[k+1,:] = data2d[J,I]
      

      intfsec_sum=intfsec_sum + intfsec
      datasec_sum=datasec_sum + datasec
      #print 'prs_intafce=', numpy.transpose(intfsec[:,15]) 
     
   dta_afile.close()
   int_afile.close()
      # end loop over files

   print 'count_sum=',count_sum 
   intfsec_avg=intfsec_sum/count_sum
   datasec_avg=datasec_sum/count_sum
   #
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
   if 'sal' in variable :
      P1=ax.contour(xx,-intfsec_avg,datasec_avg,levels=[32.0,33.0,34.0,35.0,35.5],
            colors=('k',),linestyles=('-',),linewidths=(1.5,))
   else :
      P1=ax.contour(xx,-intfsec_avg,datasec_avg,levels=[-1.0,0.0,2.0],
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
      elif k%2 == 0 and k > 10: 
         PL=ax.plot(x,-intfsec_avg[k,:],"--",color="k", linewidth=0.5)
         textx = x[i_maxd]
         texty = -0.5*(intfsec_avg[k-1,i_maxd] + intfsec_avg[k,i_maxd])
         ax.text(textx,texty,str(k),verticalalignment="center",horizontalalignment="center",fontsize=6)
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
   if 'Fram' in sectionid or 'Svin' in sectionid:
      ax.set_ylim(-600,0)
      figure.canvas.print_figure("sec_AVG_%s_600m_%s.png"%(variable,suff),dpi=dpi)
   else :
      ax.set_ylim(-3000,0)
      figure.canvas.print_figure("sec_AVG_%s_3000m_%s.png"%(variable,suff),dpi=dpi)
      #ax.set_ylim(-600,0)
      #figure.canvas.print_figure("sec_AVG_%s_600m_%s.png"%(variable,suff),dpi=dpi)

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
   parser.add_argument('lon1',     type=int, help='')
   parser.add_argument('lat1',     type=int, help='')
   parser.add_argument('lon2',     type=int, help='')
   parser.add_argument('lat2',     type=int, help='')
   parser.add_argument('variable', type=str, help='')
   parser.add_argument('files',nargs="+")

   args = parser.parse_args()

   main(args.lon1,args.lat1,args.lon2,args.lat2,args.variable,args.files,filetype=args.filetype,clim=args.clim,sectionid=args.sectionid,ijspace=args.ij,
         xaxis=args.xaxis,section_map=args.section_map,
         dpi=args.dpi) 
