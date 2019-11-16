#import scipy.io.netcdf
import netCDF4 
import numpy
import logging
import re  
import xml.etree.ElementTree
import datetime
import cfunits
import netcdftime
import scipy

class FieldInterpolatorError(object) :
    """Base class for exceptions in this module."""
    pass


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



class FieldInterpolator(object) :
   def __init__(self,x,y,targetx,targety,x_is_longitude=True) :
      self._x=x
      self._y=y


      dx=self._x[1:]-self._x[:-1]
      absdiffx= numpy.abs((numpy.max(dx)-numpy.min(dx))/numpy.abs(numpy.max(dx)))
      dy=self._y[1:]-self._y[:-1]
      absdiffy= numpy.abs((numpy.max(dy)-numpy.min(dy))/numpy.abs(numpy.max(dy)))

      if len(self._x.shape) <> 1 or len(self._x.shape) <> 1 :
         raise FieldInterpolatorError,"Field coordinates (x) must be 1D and regular"
      elif absdiffx>1e-5 :
         raise FieldInterpolatorError,"Field coordinates (x) must be 1D and regular"

      if len(self._y.shape) <> 1 or len(self._y.shape) <> 1 :
         raise FieldInterpolatorError,"Field coordinates (y) must be 1D and regular"
      elif absdiffy>1e-5 :
         raise FieldInterpolatorError,"Field coordinates (y) must be 1D and regular"


      # Simplify internal calculations by flipping coordinate directions to be increasing
      self._dx=numpy.abs(dx[0])
      self._nx=self._x.size
      if dx[0] < 0 :
         self._x0=x[-1]
         self._flipx=True
         self._x=-self._x
      else :
         self._x0=x[0]
         self._flipx=False
      self._dy=numpy.abs(dy[0])
      self._ny=self._y.size
      if dy[0] < 0 :
         self._y0=y[-1]
         self._flipy=True
         self._y=-self._y
      else :
         self._y0=y[0]
         self._flipy=False


      if x_is_longitude and abs(self._dx*self._nx - 360.0) < 1e-4:
         self._wrapmode=("raise","wrap")
      self._targetx=targetx
      self._targety=targety





class FieldInterpolatorBilinear(FieldInterpolator)  :
   def __init__(self,x,y,targetx,targety,x_is_longitude=True) :
      super(FieldInterpolatorBilinear,self).__init__(x,y,targetx,targety,x_is_longitude=x_is_longitude)
      self._calc_bilinear_weights()


   def _calc_pivot_points(self) :

      float_ip = (self._targetx-self._x0)/self._dx
      float_jp = (self._targety-self._y0)/self._dy

      # Pivot points
      self._ip=numpy.floor(float_ip).astype("int")
      self._jp=numpy.floor(float_jp).astype("int")

      # Remainder used in weights
      self._irem = numpy.remainder(float_ip,1.)
      self._jrem = numpy.remainder(float_jp,1.)

      # Pivot points displaced
      self._ip1=self._ip+1
      self._jp1=self._jp+1

      # Indexes into flattened array
      self._Ill = numpy.ravel_multi_index([self._jp ,self._ip  ],(self._ny,self._nx),mode=self._wrapmode)
      self._Ilr = numpy.ravel_multi_index([self._jp ,self._ip1 ],(self._ny,self._nx),mode=self._wrapmode)
      self._Iur = numpy.ravel_multi_index([self._jp1,self._ip1 ],(self._ny,self._nx),mode=self._wrapmode)
      self._Iul = numpy.ravel_multi_index([self._jp1,self._ip  ],(self._ny,self._nx),mode=self._wrapmode)


   def _calc_bilinear_weights(self) :
      self._calc_pivot_points()
      # Weights of points used in bilinear calculation
      self._wll = (1.-self._irem)*(1.-self._jrem) # Lower left
      self._wlr = (   self._irem)*(1.-self._jrem) # Lower right
      self._wur = (   self._irem)*(   self._jrem) # Upper right
      self._wul = (1.-self._irem)*(   self._jrem) # Upper left


   def interpolate(self,infld) :
      tmp = numpy.array(infld)

      # Flip array if needed for the interpolation (see init)
      if self._flipx :
         tmp = numpy.fliplr(tmp)
      if self._flipy :
         tmp = numpy.flipud(tmp)

      tmp=tmp.flatten()

      newfld = tmp[self._Ill]*self._wll 
      newfld+= tmp[self._Ilr]*self._wlr 
      newfld+= tmp[self._Iur]*self._wur
      newfld+= tmp[self._Iul]*self._wul

      newfld.shape=self._targetx.shape
      return newfld




class FieldInterpolatorRectBivariateSpline(FieldInterpolator)  :
   def __init__(self,x,y,targetx,targety,x_is_longitude=True) :
      super(FieldInterpolatorRectBivariateSpline,self).__init__(x,y,targetx,targety,x_is_longitude=x_is_longitude)
      tmp=numpy.array(field)
      if self._flipx :
         tmp=numpy.fliplr(tmp)
      if self._flipy :
         tmp=numpy.flipud(tmp)
      self._intobj = scipy.interpolate.RectBivariateSpline(self._x,self._y,tmp.transpose())

   def interpolate(self,fld) :
      return self._intobj(self._targetx,self._targety,grid=False)


def extrapolate_data(infld,method) :
   outfld=numpy.ma.MaskedArray.copy(infld)
   mask=numpy.copy(outfld.mask)

   # Find points that are defined, and has undefined neighbours
   #logger.info("1")
   tmp=numpy.zeros(outfld.shape)
   tmp[1:-1,1:-1] = tmp[1:-1,1:-1] + mask[1:-1,0:-2]
   tmp[1:-1,1:-1] = tmp[1:-1,1:-1] + mask[1:-1,2:]
   tmp[1:-1,1:-1] = tmp[1:-1,1:-1] + mask[0:-2,1:-1]
   tmp[1:-1,1:-1] = tmp[1:-1,1:-1] + mask[2:,1:-1]
   #logger.debug("unmask_data 2")
   tmp1=numpy.where(tmp>=1,True,False)              # More than one undefined neighbour
   tmp2=numpy.where(mask,False,True)                # Ocean point (mask=True is land)
   I,J=numpy.where(numpy.logical_and(tmp1,tmp2))    # Combination
   #logger.debug("unmask_data 3")

   # Found no points , return unmodified field
   #print len(I)
   if len(I) == 0  :
      pass
   else :
      gd=outfld[I,J]                                   # Points input to griddata
      # diag
      #tmp=numpy.zeros(outfld.shape)
      #tmp[I,J]=gd
      #plot_test(tmp,"s_coast.png")
      # Interpolate using griddata
      #logger.debug("unmask_data 3.1")
      grid_x,grid_y=numpy.meshgrid(range(tmp1.shape[0]),range(tmp1.shape[1]))
      #logger.debug("unmask_data 3.2")
      new = scipy.interpolate.griddata((I,J,),gd,(grid_x.transpose(),grid_y.transpose()),method)
      #logger.debug("unmask_data 3.3")
      outfld[outfld.mask] = new[outfld.mask]
      #plot_test(infld,"s_out.png")
      #plot_test(new,"s_griddata.png")
      #plot_test(outfld,"s_final.png")
   #logger.debug("unmask_data 4")
   outfld=numpy.ma.masked_invalid(outfld)
   return outfld








