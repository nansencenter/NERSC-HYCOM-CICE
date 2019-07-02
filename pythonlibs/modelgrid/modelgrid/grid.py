""" Module for generating grids, mainly for use by models.  """

import pyproj
import numpy
import logging
import re 
import confmap


# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False # Dont propagate to parent in hierarchy (determined by "." in __name__)

class Grid(object) :
   """ Class inherited by others . not meant to be instantiated"""
   def _grid(self,deltax,deltay,extended=False) : 
      raise NotImplementedError,""


   def proj_is_latlong(self) :
      raise NotImplementedError,""


   def p_azimuth(self) :
      plon,plat = self.pgrid(extended=True)
      #return fwd_azimuth( plon[1:-1,1:-1],plat[1:-1,1:-1], 

      tmp = fwd_azimuth( plon[1:-1,1:-1],plat[1:-1,1:-1],plon[1:-1,2:],plat[1:-1,2:])

      # Compass angle -> angle rel lat line
      tmp = numpy.radians(90-tmp)

      # Put in range [-pi,pi)
      tmp = numpy.where(tmp >= numpy.pi,tmp-2*numpy.pi,tmp) # Safest option
      tmp = numpy.where(tmp < -numpy.pi,tmp+2*numpy.pi,tmp) # Safest option
      #mult = max(-tmp.min()/numpy.pi,1)
      #mult = 2*(mult/2) + 1                                # use any odd number as long as tmp+mult*numpy.pi > 0...
      #tmp = numpy.fmod(tmp+9*numpy.pi,2*numpy.pi)-numpy.pi # use any odd number as long as tmp+mult*numpy.pi > 0...
      return tmp
      
   def corio(self) :
      qlon,qlat = self.qgrid()
      #return numpy.sin(numpy.radians(qlat)) * 4. * numpy.pi / 86400.0
      return numpy.sin(numpy.radians(qlat)) * 4. * numpy.pi / 86164.0 # Sidereal day


   def aspect_ratio(self) :
      scpx=self.scpx()
      scpy=self.scpy()
      asp = numpy.where(scpy==0.,99.0,scpx/scpy)
      return asp
   

   def plotgrid(self,fac=1.) :
      return plotgrid(*self.pgrid(),width=self.width*fac,height=self.height*fac)

      

# --- ------------------------------------------------------------------
# --- Calc scuy and scvx from qlat, qlon
# ---               *--*--*
# ---               |  |  |
# --- stencil:      u--p--*
# ---               |  |  |
# --- We have q   ->q--v--*
# --- NB: Python uses reversed indexing rel fortran
# --- ------------------------------------------------------------------

   def scuy(self) :
      qlon,qlat = self.qgrid(extended=True)
      I   = numpy.s_[1:-1]
      J   = numpy.s_[1:-1]
      Jp1 = numpy.s_[2:]
      return actual_grid_spacing(qlon[J,I],qlat[J,I], qlon[Jp1,I], qlat[Jp1,I])

   def cice_hte(self) :
      # length of easternmost grid cell wall
      qlon,qlat = self.qgrid(extended=True)
      Ip1  = numpy.s_[2:]
      Jp1  = numpy.s_[2:]
      J    = numpy.s_[1:-1]
      return actual_grid_spacing(qlon[Jp1,Ip1],qlat[Jp1,Ip1], qlon[J,Ip1], qlat[J,Ip1])

   def scvx(self) :
      qlon,qlat = self.qgrid(extended=True)
      I   = numpy.s_[1:-1]
      J   = numpy.s_[1:-1]
      Ip1 = numpy.s_[2:]
      return actual_grid_spacing(qlon[J,I],qlat[J,I], qlon[J,Ip1],qlat[J,Ip1])

   def cice_htn(self) :
      # length of northernmost grid cell wall
      qlon,qlat = self.qgrid(extended=True)
      I   = numpy.s_[1:-1]
      Ip1 = numpy.s_[2:]
      Jp1 = numpy.s_[2:]
      return actual_grid_spacing(qlon[Jp1,I],qlat[Jp1,I], qlon[Jp1,Ip1],qlat[Jp1,Ip1])

# --- ------------------------------------------------------------------
# --- Calc scvy and scux from plat, plon
# ---               *--*--*
# ---               |  |  |
# --- stencil:      u--p--*
# ---               |  |  |
# --- We have p   ->q--v--*
# --- NB: Python uses reversed indexing rel fortran
# --- ------------------------------------------------------------------

   def scvy(self) :
      plon,plat = self.pgrid(extended=True)
      I  =numpy.s_[1:-1]
      J  =numpy.s_[1:-1]
      Jm1=numpy.s_[:-2]
      return actual_grid_spacing(plon[J,I],plat[J,I], plon[Jm1,I], plat[Jm1,I])

   def cice_huw(self) :
      plon,plat = self.pgrid(extended=True)
      I  =numpy.s_[1:-1]
      Jp1=numpy.s_[2:]
      J  =numpy.s_[1:-1]
      return actual_grid_spacing(plon[J,I],plat[J,I], plon[Jp1,I], plat[Jp1,I])

   def scux(self) :
      plon,plat = self.pgrid(extended=True)
      I  =numpy.s_[1:-1]
      J  =numpy.s_[1:-1]
      Im1=numpy.s_[:-2]
      return actual_grid_spacing(plon[J,I],plat[J,I], plon[J,Im1],plat[J,Im1])

   def cice_hus(self) :
      plon,plat = self.pgrid(extended=True)
      I  =numpy.s_[1:-1]
      J  =numpy.s_[1:-1]
      Ip1=numpy.s_[:-2]
      return actual_grid_spacing(plon[J,I],plat[J,I], plon[J,Ip1],plat[J,Ip1])

# --- ------------------------------------------------------------------
# --- Calc scpx and scqy from ulat, ulon
# ---               *--*--*
# ---               |  |  |
# --- stencil:      u--p--*
# ---               |  |  |
# --- We have u   ->q--v--*
# --- ------------------------------------------------------------------

   def scpx(self) :
      ulon,ulat = self.ugrid(extended=True)
      I  =numpy.s_[1:-1]
      J  =numpy.s_[1:-1]
      Ip1=numpy.s_[2:]
      return actual_grid_spacing(ulon[J,I],ulat[J,I], ulon[J,Ip1], ulat[J,Ip1])

   def scqy(self) :
      ulon,ulat = self.ugrid(extended=True)
      I  =numpy.s_[1:-1]
      J  =numpy.s_[1:-1]
      Jm1=numpy.s_[:-2]
      return actual_grid_spacing(ulon[J,I],ulat[J,I], ulon[Jm1,I],ulat[Jm1,I])


# --- ------------------------------------------------------------------
# --- Calc scpy and scqx from vlat,vlon
# ---               *--*--*
# ---               |  |  |
# --- stencil:      u--p--*
# ---               |  |  |
# --- We have v   ->q--v--*
# --- ------------------------------------------------------------------

   def scpy(self) :
      vlon,vlat = self.vgrid(extended=True)
      I  =numpy.s_[1:-1]
      J  =numpy.s_[1:-1]
      Jp1=numpy.s_[2:]
      return actual_grid_spacing(vlon[J,I],vlat[J,I], vlon[Jp1,I], vlat[Jp1,I])

   def scqx(self) :
      vlon,vlat = self.vgrid(extended=True)
      I  =numpy.s_[1:-1]
      J  =numpy.s_[1:-1]
      Im1=numpy.s_[:-2]
      return actual_grid_spacing(vlon[J,I],vlat[J,I], vlon[J,Im1],vlat[J,Im1])

# --- End of grid-related stuff

   def save_to_scrip(self,filename,mask=None) :
      import matplotlib
      #import scipy.io.netcdf
      import netCDF4
      #nc = scipy.io.netcdf.netcdf_file("tst.nc","w")
      nc = netCDF4.Dataset(filename,"w")
      plon,plat=self.pgrid(extended=False)
      qlon,qlat=self.qgrid(extended=True)

      nc.createDimension("grid_size",plon.size)
      nc.createDimension("grid_corners",4)
      nc.createDimension("grid_rank",2)

      nc.createVariable("grid_dims","i8",("grid_rank",))
      nc.createVariable("grid_center_lon","d",("grid_size",))
      nc.createVariable("grid_center_lat","d",("grid_size",))
      nc.createVariable("grid_imask","i8",("grid_size",))
      nc.createVariable("grid_corner_lon","d",("grid_size","grid_corners",))
      nc.createVariable("grid_corner_lat","d",("grid_size","grid_corners",))

      nc["grid_center_lon"].setncattr("units","degrees")
      nc["grid_corner_lon"].setncattr("units","degrees")
      nc["grid_center_lat"].setncattr("units","degrees")
      nc["grid_corner_lat"].setncattr("units","degrees")

      nc["grid_center_lon"][:]=plon.flatten()
      nc["grid_center_lat"][:]=plat.flatten()

      nc["grid_dims"][0]=plon.shape[1]
      nc["grid_dims"][1]=plon.shape[0]

      if mask is None :
         nc["grid_imask"][:]=1
      else :
         nc["grid_imask"][:]=mask[:]


      nc["grid_corner_lon"][:,0]=qlon[1:-1,1:-1].flatten()
      nc["grid_corner_lon"][:,1]=qlon[1:-1,2:  ].flatten()
      nc["grid_corner_lon"][:,2]=qlon[2:  ,2:  ].flatten()
      nc["grid_corner_lon"][:,3]=qlon[2:  ,1:-1].flatten()

      nc["grid_corner_lat"][:,0]=qlat[1:-1,1:-1].flatten()
      nc["grid_corner_lat"][:,1]=qlat[1:-1,2:  ].flatten()
      nc["grid_corner_lat"][:,2]=qlat[2:  ,2:  ].flatten()
      nc["grid_corner_lat"][:,3]=qlat[2:  ,1:-1].flatten()


      itest=0
      #print nc["grid_center_lon"][itest]
      #print nc["grid_corner_lon"][itest]
      #print nc["grid_center_lat"][itest]
      #print nc["grid_corner_lat"][itest]
      figure = matplotlib.pyplot.figure(figsize=(8,8))
      ax=figure.add_subplot(111)
      ax.hold(True)
      ax.set_xlim(numpy.amin(nc["grid_corner_lon"][itest,:])-.1,numpy.amax(nc["grid_corner_lon"][itest,:])+.1)
      ax.set_ylim(numpy.amin(nc["grid_corner_lat"][itest,:])-.1,numpy.amax(nc["grid_corner_lat"][itest,:])+.1)
      for i in range(4) :
         x=nc["grid_corner_lon"][itest,i]
         y=nc["grid_corner_lat"][itest,i]
         dx=nc["grid_corner_lon"][itest,(i+1)%4]-nc["grid_corner_lon"][itest,i]
         dy=nc["grid_corner_lat"][itest,(i+1)%4]-nc["grid_corner_lat"][itest,i]
         ax.plot([x,x+dx],[y,y+dy],color=".5",lw=3)
         ax.plot(x,y,"*",label="corner " + str(i),markersize=20)
      ax.plot(nc["grid_center_lon"][itest],nc["grid_center_lat"][itest],"s",label="center")
      ax.legend()


      nc.close()

   def create_datadict_hycom(self) :
      """ Used when writing regional.grid files for hycom"""
      plon,plat=self.pgrid()
      ulon,ulat=self.ugrid()
      vlon,vlat=self.vgrid()
      qlon,qlat=self.qgrid()
      datadict = {}
      datadict["plon"]=plon
      datadict["plat"]=plat
      datadict["qlon"]=qlon
      datadict["qlat"]=qlat
      datadict["ulon"]=ulon
      datadict["ulat"]=ulat
      datadict["vlon"]=vlon
      datadict["vlat"]=vlat
      datadict["pang"]=self.p_azimuth()
      datadict["scpx"]=self.scpx()
      datadict["scpy"]=self.scpy()
      datadict["scqx"]=self.scqx()
      datadict["scqy"]=self.scqy()
      datadict["scux"]=self.scux()
      datadict["scuy"]=self.scuy()
      datadict["scvx"]=self.scvx()
      datadict["scvy"]=self.scvy()
      datadict["cori"]=self.corio()
      datadict["pasp"]=self.aspect_ratio()
      return datadict

   
   @property
   def width(self) : 
      if self.proj_is_latlong() :
         return self._Nx*self.dx*111000.
      else :
         return self._Nx*self.dx

   @property
   def height(self) :
      if self.proj_is_latlong() :
         return self._Ny*self.dy*111000.
      else :
         return self._Ny*self.dy

      
   @property
   def Nx(self) : return self._Nx

   @property
   def Ny(self) : return self._Ny

class ConformalGrid(Grid) :
   """ Grid generator based on Bentsen et al conformal mapping """

   def __init__(self,
         lat_a,lon_a,lat_b,lon_b,
         wlim,elim,ires,
         slim,nlim,jres,
         mercator,
         mercfac,lold) :

      self._Nx = ires
      self._Ny = jres

      self._conformal_mapping = confmap.ConformalMapping(
         lat_a,lon_a,lat_b,lon_b,
         wlim,elim,ires,
         slim,nlim,jres,
         mercator,
         mercfac,lold) 


   @classmethod
   def init_from_file(cls,filename) :
      fid=open(filename,"r")

      logger.debug("Contents of %s:\n%s"%(filename,open(filename,"r").read()))

      #-40.0 140.0     ! Position of pole N (lattitude, longitude):
      tmp=fid.readline().split("!")[0].strip().split()
      lat_a,lon_a = [float(elem) for elem in tmp ]

      #-50.0 140.0      ! Position of pole S (lattitude, longitude):
      tmp=fid.readline().split("!")[0].strip().split()
      lat_b,lon_b = [float(elem) for elem in tmp ]

      #178.2 181.42 800 ! Longitude interval West lon,   East lon, idim
      tmp=fid.readline().split("!")[0].strip().split()
      wlim,elim = [float(elem) for elem in tmp[0:2]]
      ires      = int(tmp[2])

      # 0.0  80.0 760 ! Lattitude interval south limit, north limit, jdim
      tmp=fid.readline().split("!")[0].strip().split()
      slim,nlim = [float(elem) for elem in tmp[0:2]]
      jres      = int(tmp[2])

      #.true.            ! Generate topography
      fid.readline() # Not needed for grid generation

      #.true.            ! dump teclatlon.dat
      fid.readline() # Not needed for grid generation

      #.true.            ! dump micom latlon.dat
      fid.readline() # Not needed for grid generation

      #.true.            ! mercator grid (true, false)
      tmp=fid.readline().split("!")[0].strip().split()
      if tmp[0] == ".true." : 
         mercator=True
      elif tmp[0] == ".false." : 
         mercator=False
      else :
         raise ValueError,"Unable to safely resolve value of mercator flag for confmap"

      #   0.365 .false.   ! merc fac
      tmp=fid.readline().split("!")[0].strip().split()
      mercfac = float(tmp[0])
      if tmp[1] == ".true." : 
         lold=True
      elif tmp[1] == ".false." : 
         lold=False
      else :
         raise ValueError,"Unable to safely resolve value of lold flag for confmap"

      #.false.           ! Smoothing, Shapiro filter
      fid.readline() # Not needed for grid generation

      # 8   2            ! Order of Shapiro filter,  number of passes
      fid.readline() # Not needed for grid generation

      fid.close()

      return cls(lat_a,lon_a,lat_b,lon_b,
         wlim,elim,ires,
         slim,nlim,jres,
         mercator,
         mercfac,lold) 





      
   def _grid(self,deltax,deltay,extended=False) : 
      """ TODO: delta x in grid distance  - make sure same across subclasses"""
      tmp= self._conformal_mapping.get_grid(shifti=deltax,shiftj=deltay,extended=extended)

      # Order reversed
      return tmp[1],tmp[0]


   # Shift is shift in new coordinate spacing
   def pgrid(self,extended=False) : return self._grid(0.,0.,extended)
   def ugrid(self,extended=False) : return self._grid(-0.5,0.,extended)
   def vgrid(self,extended=False) : return self._grid(0.,-0.5,extended)
   def qgrid(self,extended=False) : return self._grid(-0.5,-0.5,extended)
   

   def cice_ugrid(self,extended=False) :
      # TODO: Check out why its like this ....
      #return self._grid(+0.5*self._dx,+0.5*self._dy,extended)
      #return self._grid(-0.5,-0.5,extended)
      return self._grid(0.5,0.5,extended)

   @property
   def dx(self) : 
      return numpy.median(self.scpx())

   @property
   def dy(self) : 
      return numpy.median(self.scpy())
   
   
   def proj_is_latlong(self) :
      return False

   @property
   def mapping(self) : return self._conformal_mapping


class Proj4Grid(Grid) :
   """ Grid generator based on proj4 projections """

   def __init__(self,proj4string,ll_lon,ll_lat,dx,dy,Nx,Ny) :
      self._proj4string=proj4string
      self._proj=pyproj.Proj(proj4string)
      self._initgrid(ll_lon,ll_lat,dx,dy,Nx,Ny)


   def _initgrid(self,ll_lon,ll_lat,dx,dy,Nx,Ny) :
      self._dx=dx
      self._dy=dy
      self._Nx=Nx
      self._Ny=Ny
      self._ll_lon=ll_lon
      self._ll_lat=ll_lat

      # Calculate 0,0 (LL=Lower Left) in projection coordinates. 
      if self.proj_is_latlong() :
         self._ll_x,self._ll_y = ll_lon,ll_lat
      else :
         self._ll_x,self._ll_y = self._proj(ll_lon,ll_lat)

      # Create grid. This is the P-grid. Note increase of stencil - used to calculate grid sizes
      self._x=self._ll_x + numpy.linspace(-dx,dx*(Nx),Nx+2)
      self._y=self._ll_y + numpy.linspace(-dy,dy*(Ny),Ny+2)
      self._X,self._Y = numpy.meshgrid(self._x,self._y)
      #print self._X.shape,self._x.shape,self._y.shape

      if self.proj_is_latlong() :
         tmp= self._X,self._Y
         #print self._X
      else :
         tmp= self._proj(self._X,self._Y,inverse=True)
      logger.debug("Initialized P-grid using projection %s"%self._proj4string)
      logger.debug("Lower left corner lon/lat of grid: (%.3g,%.3g)" % (ll_lon,ll_lat))
      logger.debug("Grid spacing in projection coords: (%.3g,%.3g)" % (self._dx,self._dy))
      logger.debug("Number of grid Nodes in x/y      : (%5d,%5d)"   % (self._Nx,self._Ny))

      logger.debug("Min   x projection coordinate = %.3g"%self._x.min())
      logger.debug("Max   x projection coordinate = %.3g"%self._x.max())
      logger.debug("Min   y projection coordinate = %.3g"%self._y.min())
      logger.debug("Max   y projection coordinate = %.3g"%self._y.max())

      logger.debug("Min lon = %.3g"%tmp[0].min())
      logger.debug("Max lon = %.3g"%tmp[0].max())
      logger.debug("Min lat = %.3g"%tmp[1].min())
      logger.debug("Max lat = %.3g"%tmp[1].max())

      # Should insist that the ellipse is spherical
      #if bla bla bla :
      #   msg = "The geoid used in the projection is not spherical. Aborting"
      #   logger.error(msg)
      #   raise ValueError,msg

   def _grid(self,deltax,deltay,extended=False) : 
      """ TODO: delta x in ACTUAL distance  - make sure same across subclasses"""
      if extended :
         tmp = (self._X+deltax,self._Y+deltay)
         #return self._proj(self._X+deltax,self._Y+deltay,inverse=True)
      else :
         #return self._proj(self._X[1:-1,1:-1]-deltax,self._Y[1:-1,1:-1]-deltay,inverse=True)
         tmp = (self._X[1:-1,1:-1]-deltax,self._Y[1:-1,1:-1]-deltay)

      if self.proj_is_latlong() :
         return tmp
      else :
         return self._proj(*tmp,inverse=True)

   # Shift is shift in model projection
   def pgrid(self,extended=False) : return self._grid(0.,0.,extended)
   def ugrid(self,extended=False) : return self._grid(-0.5*self._dx,0.,extended)
   def vgrid(self,extended=False) : return self._grid(0.,-0.5*self._dy,extended)
   def qgrid(self,extended=False) : return self._grid(-0.5*self._dx,-0.5*self._dy,extended)

   def cice_ugrid(self,extended=False) :
      # TODO: Check out why its like this ....
      #return self._grid(+0.5*self._dx,+0.5*self._dy,extended)
      raise NotImplementedEerror,"check cice_ugrid return values"
      return self._grid(-0.5*self._dx,-0.5*self._dy,extended)

   def proj_is_latlong(self) :
      return self._proj.is_latlong() 


   def write_my_projection_info(self) :
      self.write_projection_info(self._proj4string,self._ll_lon,self._ll_lat,self._dx,self._dy,
         self._Nx,self._Ny)

   @classmethod
   def write_projection_info(cls,p4s,ll_lon,ll_lat,dx,dy,nx,ny) :
      fid=open("proj.info","w")
      fid.write("File Version         : %.1f\n"%1.0)
      fid.write("Proj4 String         : %s\n"%p4s)
      fid.write("Lower Left Longitude : %16.10f\n"%ll_lon)
      fid.write("Lower Left Latitude  : %16.10f\n"%ll_lat)
      fid.write("Projection Delta X   : %16.10g\n"%dx)
      fid.write("Projection Delta Y   : %16.10g\n"%dy)
      fid.write("Projection Nodes X   : %6i\n"%nx)
      fid.write("Projection Nodes Y   : %6i\n"%ny)
      fid.close()
      
   @classmethod
   def grid_from_file(cls,filename="proj.info") :
      p4s,ll_lon,ll_lat,dx,dy,nx,ny = cls.read_projection_info(filename)
      return cls(p4s,ll_lon,ll_lat,dx,dy,nx,ny)

   
   @classmethod
   def read_projection_info(cls,filename="proj.info") :
      fid=open(filename,"r")

      res = True

      line=fid.readline().strip()
      m = re.search("File Version[ ]*:(.*)",line)
      if m :
         #print m.group(1)
         v = float(m.group(1))
      else :
         msg = "Failed to read version number from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg


      line=fid.readline()
      m = re.search("Proj4 String[ ]*:(.*)",line)
      if m :
         p4s = m.group(1)
      else :
         msg = "Failed to read proj 4 string from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Lower Left Longitude[ ]*:(.*)",line)
      if m :
         ll_lon = float(m.group(1))
      else :
         msg = "Failed to read lower left longitude from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Lower Left Latitude[ ]*:(.*)",line)
      if m :
         ll_lat = float(m.group(1))
      else :
         msg = "Failed to read lower left latitude number from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Projection Delta X[ ]*:(.*)",line)
      if m :
         dx = float(m.group(1))
      else :
         msg = "Failed to read delta x from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Projection Delta Y[ ]*:(.*)",line)
      if m :
         dy = float(m.group(1))
      else :
         msg = "Failed to read delta y from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Projection Nodes X[ ]*:(.*)",line)
      if m :
         nx = int(m.group(1))
      else :
         msg = "Failed to read number of x nodes from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg

      line=fid.readline()
      m = re.search("Projection Nodes Y[ ]*:(.*)",line)
      if m :
         ny = int(m.group(1))
      else :
         msg = "Failed to read number of y nodes from projection file %s"%filename
         logger.error(msg)
         raise ValueError,msg
      fid.close()

      return  p4s,ll_lon,ll_lat,dx,dy,nx,ny



   def raster(self) :
      pass

   @property
   def dx(self) : return self._dx

   @property
   def dy(self) : return self._dy
   







def actual_grid_spacing(lon1,lat1,lon2,lat2) : 
   geod=pyproj.Geod(ellps="sphere")
   az1,az2,dist = geod.inv(lon1,lat1,lon2,lat2)
   return dist

def fwd_azimuth(lon1,lat1,lon2,lat2) : 
   geod=pyproj.Geod(ellps="sphere")
   az1,az2,dist = geod.inv(lon1,lat1,lon2,lat2)
   return az1

   


def plotgrid(lon,lat,width=3000000,height=3000000) :
   import matplotlib
   import matplotlib.pyplot
   from matplotlib.figure import Figure
   from matplotlib.backends.backend_agg import FigureCanvasAgg
   from mpl_toolkits.basemap import Basemap

   #print "testlon, testlat:",lon[1,1],lat[1,1]

   #figure = Figure()
   #ax     = figure.add_subplot(111)
   #canvas = FigureCanvasAgg(figure)

   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)

   #Pick center longitude
   ix,iy=[elem/2 for elem in lon.shape]
   clon=lon[ix,iy]
   clat=lat[ix,iy]

   # Probably a way of estimating the width here...
   #print width,height
   #print clon,clat
   m = Basemap(projection='stere',lon_0=clon,lat_0=clat,resolution='l',width=width,height=height,ax=ax)
   x,y = m(lon,lat)

   # Pick a suitable set of grid lines
   nlines=10
   stepx,stepy=[elem/nlines for elem in lon.shape]
   x2=numpy.zeros((nlines+1,nlines+1))
   y2=numpy.zeros((nlines+1,nlines+1))

   #print x2.shape,x[::stepx,::stepy].shape
   x2[:-1,:-1]=x[::stepx,::stepy]
   x2[-1,:-1]=x[-1,::stepy]
   x2[:-1,-1]=x[::stepx,-1]
   x2[-1,-1]=x[-1,-1]
   y2[:-1,:-1]=y[::stepx,::stepy]
   y2[:-1,-1]=y[::stepx,-1]
   y2[-1,:-1]=y[-1,::stepy]
   y2[-1,-1]=y[-1,-1]

   m.drawcoastlines()
   m.drawmapboundary() # draw a line around the map region
   m.drawparallels(numpy.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
   m.drawmeridians(numpy.arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
   v=numpy.zeros(x.shape)
   col=".8"
   cmap=matplotlib.colors.ListedColormap([col,col])
   #m.pcolormesh(x,y,v,ax=ax,edgecolor="k",cmap=cmap)
   m.pcolormesh(x,y,v,ax=ax,edgecolor="none",cmap=cmap)
   for j in range(y2.shape[1]) :
      m.plot(x2[:,j],y2[:,j],color="b",lw=2)
   for i in range(y2.shape[0]) :
      m.plot(x2[i,:],y2[i,:],color="r",lw=2)
   ax.set_title("Every %d in x(blue) and every %d in y(red) shown"%(stepy,stepx))

   return figure


