#!/usr/bin/env python
import warnings
warnings.filterwarnings("ignore")
import modeltools.hycom
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
#import modeltools.hycom.io
import abfile
import modeltools.cice.io
import numpy
##AO from mpl_toolkits.basemap import Basemap
import netCDF4
import logging
import scipy.ndimage.measurements
import re
#  Modified June 2018, Mostafa Bakhoday-Paskyabi
#  Modified May 2019, Mostafa Bakhoday-Paskyabi


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


def port_setup(kdport,in_depth_m) :
   labels=numpy.zeros(in_depth_m.shape) 
   if kdport == 2 :
      jfport=1 # starts from 0. Test on first ocean cell
      jlport=1 # starts from 0. Test on first ocean cell
      my_label, my_nf = scipy.ndimage.measurements.label(~in_depth_m.mask[jfport,:])
      labels[jfport,my_label>0] = my_label[my_label>0]
      edge="south"
   elif kdport==1 :
      jfport=in_depth_m.shape[0]-2 # starts from 0. Test on first ocean cell
      jlport =jfport             # starts from 0. Test on first ocean cell
      my_label, my_nf = scipy.ndimage.measurements.label(~in_depth_m.mask[jfport,:])
      labels[jfport,my_label>0] = my_label[my_label>0]
      edge="north"
   elif kdport==3 :
      ifport=in_depth_m.shape[1]-2 # starts from 0. Test on first ocean cell
      ilport=ifport                # starts from 0. Test on first ocean cell
      my_label, my_nf = scipy.ndimage.measurements.label(~in_depth_m.mask[:,ifport])
      labels[my_label>0,ifport] = my_label[my_label>0] + labels.max()
      edge="east"
   elif kdport==4 :
      ifport=1  # Starts from 0. Test on first ocean cell
      ilport=1  # Starts from 0. Test on first ocean cell
      my_label, my_nf = scipy.ndimage.measurements.label(~in_depth_m.mask[:,ifport])
      labels[my_label>0,ifport] = my_label[my_label>0] + labels.max()
      edge="west"

   tmp=[str(elem) for elem in my_label]
   tmp="".join(tmp)
   logger.info("%s edge segments:%s"%(edge,tmp))
   ifports=[]
   ilports=[]
   jfports=[]
   jlports=[]
   kdports=[]
   for i in range(my_nf ) :
      numport=len(ifports)+1
      I,=numpy.where(my_label==i+1)
      kdports.append(kdport)
      if kdport == 2 :
         ifport=I[0]
         ilport=I[-1]
         logger.info("Port %5d on %s edge: First and last i: %5d %5d"%(numport,edge, ifport+1,ilport+1))
         ifports.append(ifport)
         ilports.append(ilport)
         jfports.append(jfport)
         jlports.append(jlport)
      elif kdport == 1 :
         ifport=I[0]
         ilport=I[-1]
         logger.info("Port %5d on %s edge: First and last i: %5d %5d"%(numport,edge, ifport+1,ilport+1))
         ifports.append(ifport)
         ilports.append(ilport)
         jfports.append(jfport+1) # Move towards edge, accounting for iv mask
         jlports.append(jlport+1) # Move towards edge, accounting for iv mask
      elif kdport == 3 :
         jfport=I[0]
         jlport=I[-1]
         logger.info("Port %5d on %s  edge: First and last j: %5d %5d"%( numport, edge, jfport+1,jlport+1))
         ifports.append(ifport+1) # Move towards edge, accounting for iu mask
         ilports.append(ilport+1) # Move towards edge, accounting for iu mask
         jfports.append(jfport)
         jlports.append(jlport)
      elif kdport == 4 :
         jfport=I[0]
         jlport=I[-1]
         logger.info("Port %5d on %s  edge: First and last j: %5d %5d"%( numport, edge, jfport+1,jlport+1))
         ifports.append(ifport)
         ilports.append(ilport)
         jfports.append(jfport)
         jlports.append(jlport)
      else :
         raise ValueError,"Unknown kdport value %d"%kdport

   return ifports,ilports,jfports,jlports,kdports,labels





def relaxation_mask(rmumask,ifport,ilport,jfport,jlport,kdport,rmu_width) :

   allj,alli=numpy.where(numpy.ones(rmumask.shape) == 1)

   # Create relaxation mask
   # Start point -Tag points in rmu_width distance
   dist = numpy.sqrt((ifport-alli)**2 + (jfport-allj)**2) 
   tmp  = numpy.where(dist < rmu_width,1.-dist/rmu_width,0.)
   rmumask = numpy.maximum(rmumask,tmp.reshape(rmumask.shape))

   dist = numpy.sqrt((ilport-alli)**2 + (jlport-allj)**2) 
   tmp2 = numpy.where(dist < rmu_width,1.-dist/rmu_width,0.)
   rmumask = numpy.maximum(rmumask,tmp2.reshape(rmumask.shape))

   distedge=numpy.zeros(rmumask.shape)
   tmp = numpy.arange(rmu_width,dtype=float) ; tmp=1.-tmp/rmu_width ; tmp.shape=(tmp.size,1)
   if kdport == 2 : # South 
      #print "south"
      #print  distedge[jfport:rmu_width+jfport,ifport:ilport+1].shape, tmp.shape
      distedge[jfport:rmu_width+jfport,ifport:ilport+1] = tmp
   elif kdport == 1 : # North
      #print "north"
      #distedge[jfport-rmu_width:jfport,ifport:ilport+1] = tmp[::-1]
      #print distedge[jfport-rmu_width-1:jfport-1,ifport:ilport+1].shape, tmp[::-1].shape
      distedge[jfport-rmu_width:jfport,ifport:ilport+1] = tmp[::-1]
   elif kdport == 3 : # East
      #print "east"
      #distedge[jfport:jlport+1,ifport-rmu_width:ifport] = tmp[::-1]
      #distedge[jfport:jlport+1,ifport-rmu_width-1:ifport-1] = tmp[::-1]
      distedge[jfport:jlport+1,ifport-rmu_width-1:ifport-1] = tmp[::-1].transpose()
   elif kdport == 4 : # West
      #print "west"
      #print distedge[jfport:jlport+1,ifport:ifport+rmu_width].shape, tmp.shape
      distedge[jfport:jlport+1,ifport:ifport+rmu_width] = tmp.transpose()
   rmumask = numpy.maximum(rmumask,distedge)

   return rmumask



def write_port_location(fid,kdport,ifport,ilport,jfport,jlport) :
   #c ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)
   #c ---     'ifport' = first i-index
   #c ---     'ilport' = last  i-index (=ifport for N or S orientation)
   #c ---     'jfport' = first j-index
   #c ---     'jlport' = last  j-index (=jfport for E or W orientation)
   #c ---     'lnport' = port length (calculated, not input)
   fid.write("%6d  'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)\n"%kdport)
   fid.write("%6d  'ifport' = first i-index\n"%ifport)
   fid.write("%6d  'ilport' = last  i-index (=ifport for N or S orientation)\n"%ilport)
   fid.write("%6d  'jfport' = first j-index\n"%jfport)
   fid.write("%6d  'jlport' = last  j-index (=jfport for E or W orientation)\n"%jlport)
   return


def check_consistency(ifport,ilport,jfport,jlport,kdport,iu,iv,port_number) :
   fatal = False
   if kdport==2 :
      for  i2 in range(ifport,ilport+1) :
         if not iv[jfport+1,i2] or not iv[jfport+2,i2]  :
            fatal=True
            #fatal = False
            # Write using fortran indexes
            logger.error("Port number %d(kdport=%d)  : Point  (i=%5d,j=%5d) has land too close to port"%(port_number,kdport,i2+1,jfport+1))
   elif kdport == 1 :
      for  i2 in range(ifport,ilport+1) :
         if not iv[jfport-1,i2] or not iv[jfport-2,i2]  :
            fatal=True
            # Write using fortran indexes
            logger.error("Port number %d(kdport=%d)  : Point  (i=%5d,j=%5d) has land too close to port"%(port_number,kdport,i2+1,jfport+1))
   elif kdport == 3 :
      for  j2 in range(jfport,jlport+1) :
         if not iu[j2,ifport-1] or not iu[j2,ifport-2]  :
            fatal=True
            # Write using fortran indexes
            logger.error("Port number %d(kdport=%d)  : Point  (i=%5d,j=%5d) has land too close to port"%(port_number,kdport,ifport+1,j2+1))
   elif kdport == 4 :
      for  j2 in range(jfport,jlport+1) :
         if not iu[j2,ifport+1] or not iu[j2,ifport+2]  :
            fatal=True
            # Write using fortran indexes
            logger.error("Port number %d(kdport=%d)  : Point  (i=%5d,j=%5d) has land too close to port"%(port_number,kdport,ifport+1,j2+1))
   return fatal


def main(infile,rmu_width,rmu_efold,dpi=180):

   bathy_threshold=0. # TODO

   # Read plon,plat
   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()

   # Read input bathymetri
   m=re.match( "^(.*)(\.[ab])", infile)
   if m : infile=m.group(1)
   bfile=abfile.ABFileBathy(infile,"r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   in_depth_m=bfile.read_field("depth")
   bfile.close()
   in_depth=numpy.ma.filled(in_depth_m,bathy_threshold)

   #print in_depth.min(),in_depth.max()
   ip=~in_depth_m.mask
   ip=numpy.copy(ip)
   iu=numpy.copy(ip)
   iv=numpy.copy(ip)
   iu[:,1:] = numpy.logical_and(ip[:,1:],ip[:,:-1])
   iv[1:,:] = numpy.logical_and(ip[1:,:],ip[:-1,:])


   ifports=[]
   ilports=[]
   jfports=[]
   jlports=[]
   kdports=[]

   process_south=True
   process_north=True
   process_west=True
   process_east=True

   fatal=False

   rmumask=numpy.zeros(in_depth.shape)
   labels=numpy.zeros(in_depth.shape)

   # Test ocean mask in 2nd grid cell from edge. If ocean, mark as nesting boundary.
   # When written to ports.input, we write 
   # NB: All diag output is "Fortran" indexes (Starting from 1) - thats why we add 1 here and there


   for kdport in [1,2,3,4] :
      t_ifports,t_ilports,t_jfports,t_jlports,t_kdports,t_labels, = port_setup(kdport,in_depth_m)
      labels[t_labels>0] = t_labels[t_labels>0] + labels.max()
      ifports.extend(t_ifports)
      ilports.extend(t_ilports)
      jfports.extend(t_jfports)
      jlports.extend(t_jlports)
      kdports.extend(t_kdports)



   # Build mask
   for i in range(len(ifports)) :
      rmumask = relaxation_mask(rmumask,ifports[i],ilports[i],jfports[i],jlports[i],kdports[i],rmu_width) 
   #print rmumask.min(),rmumask.max()
   rmumask = numpy.minimum(rmumask,1.) * 1./ (rmu_efold * 86400.)
   rmumask_m = numpy.ma.masked_where(in_depth_m.mask,rmumask)

   # Check consistency
   fatal = False
   for i in range(len(ifports)) :
      fatal = fatal or check_consistency(ifports[i],ilports[i],jfports[i],jlports[i],kdports[i],iu,iv,i+1) 

   # Open port output file
   logger.info("Writing to ports.input.tmp")
   fid=open("ports.input.tmp","w")
   fid.write("%6d  'nports' = Number of ports \n"%len(kdports))
   for i in range(len(kdports)) :
      write_port_location(fid,kdports[i],ifports[i]+1,ilports[i]+1,jfports[i]+1,jlports[i]+1)
   fid.close()

   # Write rmu file
   rmufile=abfile.ABFileRmu("rmu","w",
         cline1="Relaxation mask",
         cline2="Relaxation mask created by topo_ports.py. rel zone width=%d, efold time=%d days"%(rmu_width,rmu_efold),
         mask=True)
   rmufile.write_field(rmumask,in_depth_m.mask,"rmu")
   rmufile.close()




   # Plot rmu with pcolormesh
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   cmap=matplotlib.pyplot.get_cmap("Greys_r")
   cmap2=matplotlib.pyplot.get_cmap("jet")
   ax.add_patch(matplotlib.patches.Rectangle((1,1),in_depth.shape[1],in_depth.shape[0],color=".5",alpha=.5))
   P=ax.pcolormesh(in_depth_m,cmap=cmap)
   P=ax.pcolormesh(rmumask_m,cmap=cmap2)
   CB=ax.figure.colorbar(P)
   figure.canvas.print_figure("rmu.png",dpi=dpi)


   # Plot ports with pcolormesh
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   cmap=matplotlib.pyplot.get_cmap("Greys_r")
   ax.add_patch(matplotlib.patches.Rectangle((1,1),in_depth.shape[1],in_depth.shape[0],color=".5",alpha=.5))
   P=ax.pcolormesh(in_depth_m,cmap=cmap)
   I,J=numpy.where(labels>0)
   S=ax.scatter(J,I,50,labels[I,J],edgecolor='none')
   CB=ax.figure.colorbar(S)
   ax.set_xlim(0,in_depth.shape[1])
   ax.set_ylim(0,in_depth.shape[0])
   CB.ax.set_title("Port number")
   logger.info("Writing to ports_all.png")
   figure.canvas.print_figure("ports_all.png",dpi=dpi)



   # Port diagnostics plot
   figure2 = matplotlib.pyplot.figure(figsize=(8,8))
   ax2=figure2.add_subplot(111)
   cmap=matplotlib.pyplot.get_cmap("Greys_r")
   P=ax2.pcolormesh(in_depth_m,cmap=cmap,edgecolor=".4",alpha=.5,linewidth=.05)
   ax2.hold()
   Ps=[]; Ls=[]
   for i in range(len(kdports)) :


      iwidth = ilports[i] - ifports[i]+1
      jwidth = jlports[i] - jfports[i]+1
      #print ifports[i],jfports[i],iwidth,jwidth
      d=1
      if kdports[i] == 1 :
         xy = (ifports[i],jfports[i])
         jwidth=d
         c="r"
      elif kdports[i] == 2 :
         xy = (ifports[i],jfports[i])
         jwidth=d
         c="g"
      elif kdports[i] == 3 :
         xy = (ifports[i],jfports[i])
         iwidth=d
         c="b"
      elif kdports[i] == 4 :
         xy = (ifports[i],jfports[i])
         iwidth=d
         c="m"


      figure.clf()
      ax=figure.add_subplot(111)
      P=ax.pcolormesh(in_depth_m,cmap=cmap,edgecolor=".4",alpha=.5,linewidth=.05)
      ax.add_patch(matplotlib.patches.Rectangle(xy,iwidth,jwidth, color=c,alpha=.5))
      ax.grid()
      ax.set_xlim(xy[0]-20,xy[0]+iwidth+20)
      ax.set_ylim(xy[1]-20,xy[1]+jwidth+20)
      ax.set_title("Port number %d - kdport=%d" %(i+1,kdports[i]))

      R=ax2.add_patch(matplotlib.patches.Rectangle(xy,iwidth,jwidth, color=c,alpha=.5))
      Ps.append(R) ; Ls.append("Port %d"%(i+1))


      fname="port_%03d.png"%(i+1)
      logger.info("Writing Diagnostics to %s"%fname)
      figure.canvas.print_figure(fname,bbox_inches='tight',dpi=dpi)

   fname="ports_all_2.png"
   logger.info("Writing Diagnostics to %s"%fname)
   ax2.legend(Ps,Ls)
   figure2.canvas.print_figure(fname,bbox_inches='tight',dpi=dpi)


   if fatal :
      logger.error("Errors were encountered - see errors above, and consult diag files. You may need to modify your topo file")
      raise NameError,"fatal exit"
   return rmumask,rmumask_m


if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Find nesting port locations along edge, set up relaxation mask')
   parser.add_argument('--dpi',    type=int,help="dpi of output plots")
   parser.add_argument('depthfile',    type=str,help="hycom bathymetri file")
   parser.add_argument('rmu_width', type=int,help="width of relaxation zone")
   parser.add_argument('rmu_efold', type=int,help="relaxation zone weights ")

   args = parser.parse_args()

   main(args.depthfile,args.rmu_width,args.rmu_efold,
         dpi=args.dpi)
