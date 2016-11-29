#!/usr/bin/env python
import modeltools.hycom
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile
import numpy
import logging
import datetime
import re
import scipy.interpolate
import os.path

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


def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    # https://ocefpaf.github.io/python4oceanographers/blog/2014/03/10/gearth/
    aspect = numpy.cos(numpy.mean([llcrnrlat, urcrnrlat]) * numpy.pi/180.0)
    xsize = numpy.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = numpy.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    fig = matplotlib.pyplot.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    ax.set_axis_off()
    return fig, ax


def interpolate_to_latlon(lon,lat,data,res=0.1) :
   lon2 = numpy.mod(lon+360.,360.)
   # New grid
   minlon=numpy.floor((numpy.min(lon2)/res))*res
   minlat=max(-90.,numpy.floor((numpy.min(lat)/res))*res)
   maxlon=numpy.ceil((numpy.max(lon2)/res))*res
   maxlat=min(90.,numpy.ceil((numpy.max(lat)/res))*res)
   #maxlat=90.
   lo1d = numpy.arange(minlon,maxlon+res,res)
   la1d = numpy.arange(minlat,maxlat,res)
   lo2d,la2d=numpy.meshgrid(lo1d,la1d)
   print minlon,maxlon,minlat,maxlat

   if os.path.exists("grid.info") :
      import modelgrid

      grid=modelgrid.ConformalGrid.init_from_file("grid.info")
      map=grid.mapping

      # Index into model data, using grid info
      I,J=map.ll2gind(la2d,lo2d)

      # Location of model p-cell corner 
      I=I-0.5
      J=J-0.5

      # Mask out points not on grid
      K=J<data.shape[0]-1
      K=numpy.logical_and(K,I<data.shape[1]-1)
      K=numpy.logical_and(K,J>=0)
      K=numpy.logical_and(K,I>=0)

      # Pivot point 
      Ii=I.astype('i')
      Ji=J.astype('i')

      # Takes into account data mask
      tmp =numpy.logical_and(K[K],~data.mask[Ji[K],Ii[K]])
      K[K]=tmp
      
      tmp=data[Ji[K],Ii[K]]
      a,b=numpy.where(K) 
      z=numpy.zeros(K.shape)
      z[a,b] = tmp
      z=numpy.ma.masked_where(~K,z)

   # Brute force ...
   else  :
      K=numpy.where(~data.mask)
      z=scipy.interpolate.griddata( (lon2[K],lat[K]),data[K],(lo2d,la2d),'cubic')
      z=numpy.ma.masked_invalid(z)
      z2=scipy.interpolate.griddata( (lon2.flatten(),lat.flatten()),data.mask.flatten(),(lo2d,la2d),'nearest')
      z2=z2>.1
      z=numpy.ma.masked_where(z2,z)

   return lo2d,la2d,z


# Simple routine to create a kml file from field

def to_kml(lon,lat,data,fieldname,cmap,res=0.1,clim=None) :
   import simplekml
   
   lo2d,la2d,z= interpolate_to_latlon(lon,lat,data,res=res)
   urcrnrlon = lo2d[-1,-1]
   urcrnrlat = la2d[-1,-1]
   llcrnrlon = lo2d[0,0]
   llcrnrlat = la2d[0,0]

   fig,ax= gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024)
   figname="overlay.png"
   P=ax.pcolormesh(lo2d,la2d,z,cmap=cmap)
   if clim is not None : P.set_clim(clim)
   fig.canvas.print_figure(figname,Transparent=True)

   kw={}
   kml=simplekml.Kml()
   draworder = 0
   draworder += 1
   ground = kml.newgroundoverlay(name='GroundOverlay')
   ground.draworder = draworder
   ground.visibility = kw.pop('visibility', 1)
   ground.name = kw.pop('name', fieldname)
   ground.color = kw.pop('color', 'ddffffff') # First hex gives transparency
   ground.atomauthor = kw.pop('author', 'NERSC')
   ground.latlonbox.rotation = kw.pop('rotation', 0)
   ground.description = kw.pop('description', 'Matplotlib figure')
   ground.gxaltitudemode = kw.pop('gxaltitudemode', 'clampToGround')
   ground.icon.href = figname
   ground.latlonbox.east = llcrnrlon
   ground.latlonbox.south = llcrnrlat
   ground.latlonbox.north = urcrnrlat
   ground.latlonbox.west = urcrnrlon
   kmzfile="overlay.kmz"
   kml.savekmz(kmzfile)


def open_file(myfile0,filetype,fieldname,fieldlevel,datetime1=None,datetime2=None,vector="",
      idm=None,
      jdm=None) :

   logger.info("Now processing  %s"%myfile0)
   m=re.match("(.*)\.[ab]",myfile0)
   if m :
      myfile=m.group(1)
   else :
      myfile=myfile0

   ab2=None
   rdtimes=[]
   if filetype == "archive" :
      ab = abfile.ABFileArchv(myfile,"r")
      n_intloop=1
   #elif filetype == "restart" :
   #   tmp = abfile.ABFileRestart(myfile,"r",idm=gfile.idm,jdm=gfile.jdm)
   elif filetype == "regional.depth" :
      ab = abfile.ABFileBathy(myfile,"r",idm=idm,jdm=jdm)
      n_intloop=1
   elif filetype == "forcing" :
      ab = abfile.ABFileForcing(myfile,"r",idm=idm,jdm=jdm)
      if vector :
         file2=myfile.replace(fieldname,vector)
         logger.info("Opening file %s for vector component nr 2"%file2)
         ab2=abfile.ABFileForcing(file2,"r",idm=idm,jdm=jdm)
      if datetime1 is None or datetime2 is None :
         raise NameError,"datetime1 and datetime2 must be specified when plotting forcing files"
      else :
         iday1,ihour1,isec1 = modeltools.hycom.datetime_to_ordinal(datetime1,3)
         rdtime1 = modeltools.hycom.dayfor(datetime1.year,iday1,ihour1,3)
         #
         iday2,ihour2,isec2 = modeltools.hycom.datetime_to_ordinal(datetime2,3)
         rdtime2 = modeltools.hycom.dayfor(datetime2.year,iday2,ihour2,3)
         rdtimes=sorted([elem for elem in ab.field_times if elem >rdtime1 and elem < rdtime2])
         n_intloop=len(rdtimes)
   else :
      raise NotImplementedError,"Filetype %s not implemented"%filetype
   # Check that fieldname is actually in file
   if fieldname not in  ab.fieldnames :
      logger.error("Unknown field %s at level %d"%(fieldname,fieldlevel))
      logger.error("Available fields : %s"%(" ".join(ab.fieldnames)))
      raise ValueError,"Unknown field %s at level %d"%(fieldname,fieldlevel)

   return n_intloop,ab,ab2,rdtimes


def main(myfiles,fieldname,fieldlevel,
      idm=None,
      jdm=None,
      clim=None,
      filetype="archive",
      window=None,
      cmap="jet",
      datetime1=None,
      datetime2=None,
      vector="",
      tokml=False,
      masklim=None,
      dpi=180) :


   #cmap=matplotlib.pyplot.get_cmap("jet")
   cmap=matplotlib.pyplot.get_cmap(cmap)

   if tokml :
      ab = abfile.ABFileGrid("regional.grid","r")
      plon=ab.read_field("plon")
      plat=ab.read_field("plat")
      ab.close()

   # Prelim support for projections. import basemap only if needed since its usually not needed
   # aaaand installation can sometimes be a bit painful .... bmn is None now, define it if needed
   bm=None
   #from mpl_toolkits.basemap import Basemap
   #ab = abfile.ABFileGrid("regional.grid","r")
   #plon=ab.read_field("plon")
   #plat=ab.read_field("plat")
   #bm = Basemap(width=9000000,height=9000000,
   #      resolution='i',projection='stere',\
   #      lat_ts=70,lat_0=90,lon_0=-40.)
   #x,y=bm(plon,plat)

   if vector :
      logger.info("Vector component 1:%s"%fieldname)
      logger.info("Vector component 2:%s"%vector) 

   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)

   counter=0
   for myfile0 in myfiles :

      # Open files, and return some useful stuff.
      # ab2 i used in case of vector
      # rdtimes is used for plotting forcing fields
      n_intloop,ab,ab2,rdtimes = open_file(myfile0,filetype,fieldname,fieldlevel,datetime1=datetime1,datetime2=datetime2,vector=vector,idm=idm,jdm=jdm)


      # Intloop used to read more fields in one file. Only for forcing for now
      for i_intloop in range(n_intloop) :

         # Read ab file of different types
         if filetype == "archive" :
            fld1 = ab.read_field(fieldname,fieldlevel)
            if vector :fld2=ab.read_field(vector,fieldlevel)
         elif filetype == "regional.depth" :
            fld1 = ab.read_field(fieldname)
         elif filetype == "forcing" :
            fld1 = ab.read_field(fieldname,rdtimes[i_intloop])
            if vector :fld2=ab2.read_field(vector,rdtimes[i_intloop])
            logger.info("Processing time %.2f"%rdtimes[i_intloop])
         else :
            raise NotImplementedError,"Filetype %s not implemented"%filetype

         if not window :
            J,I=numpy.meshgrid(numpy.arange(fld1.shape[0]),numpy.arange(fld1.shape[1]))
         else :
            J,I=numpy.meshgrid( numpy.arange(window[1],window[3]),numpy.arange(window[0],window[2]))

         # Create scalar field for vectors
         if vector : 
            fld = numpy.sqrt(fld1**2+fld2**2)
         else :
            fld=fld1

         # Apply mask if requested
         if masklim :
            fld = numpy.ma.masked_where(fld<=masklim[0],fld)
            fld = numpy.ma.masked_where(fld>=masklim[1],fld)

         if tokml :
            #NB: Window not used
            if vector :
               logger.warning("Vector plots in kml not implemented, only plottign vector speed")
            to_kml(plon,plat,fld,fieldname,cmap)
         else :

            if bm is not None :
               P=bm.pcolormesh(x[J,I],y[J,I],fld[J,I],cmap=cmap)
               bm.drawcoastlines()
               bm.fillcontinents(color='.5',lake_color='aqua')
               # draw parallels and meridians.
               bm.drawparallels(numpy.arange(-80.,81.,20.))
               bm.drawmeridians(numpy.arange(-180.,181.,20.))
               bm.drawmapboundary() #fill_color='aqua')
            else :
               P=ax.pcolormesh(I,J,fld[J,I],cmap=cmap)
            cb=ax.figure.colorbar(P)
            if clim is not None : P.set_clim(clim)
            ax.set_title("%s:%s(%d)"%(myfile0,fieldname,fieldlevel))
            if vector: 
               skip=10
               I2=I[::skip,::skip]
               J2=J[::skip,::skip]
               ax.quiver(I2,J2,fld1[J2,I2],fld2[J2,I2])


            # Print figure.
            fnamepng_template="%s_%d_%03d.png"
            fnamepng=fnamepng_template%(fieldname,fieldlevel,counter)
            logger.info("output in  %s"%fnamepng)
            figure.canvas.print_figure(fnamepng,dpi=dpi)
            ax.clear()
            cb.remove()

         counter=counter+1

      # End i_intloop



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
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--dpi',        type=int, default=180)
   parser.add_argument('--clim',       action=ClimParseAction,default=None,help="range of colormap")
   parser.add_argument('--masklim',    action=ClimParseAction,default=None,help="mask limits")
   parser.add_argument('--cmap',       type=str,default="jet",help="matplotlib colormap to use")
   parser.add_argument('--filetype',   type=str, default="archive",help="filetype to plot (archive by default)")
   parser.add_argument('--window',     action=WindowParseAction, help='firsti,firstj,lasti,lastj', default=None)
   parser.add_argument('--idm',        type=int, help='Grid dimension 1st index []', default=None)
   parser.add_argument('--jdm',        type=int, help='Grid dimension 2nd index []', default=None)
   parser.add_argument('--datetime1',  action=DateTimeParseAction)
   parser.add_argument('--datetime2',  action=DateTimeParseAction)
   parser.add_argument('--vector',     type=str,default="",help="Denotes second vector component")
   parser.add_argument('--tokml',      action="store_true", default=False)
   parser.add_argument('fieldname',  type=str)
   parser.add_argument('fieldlevel', type=int)
   parser.add_argument('filename',   help="",nargs="+")
   args = parser.parse_args()

   main(args.filename,args.fieldname,args.fieldlevel,
         idm=args.idm,jdm=args.jdm,clim=args.clim,filetype=args.filetype,
         window=args.window,cmap=args.cmap,
         datetime1=args.datetime1,datetime2=args.datetime2,
         vector=args.vector,
         tokml=args.tokml,
         masklim=args.masklim,
         dpi=args.dpi)


