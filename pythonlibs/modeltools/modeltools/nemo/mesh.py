# Various utility routines for the NEMO grid (from Mercator)
import modeltools.tools
import logging
import numpy
import netCDF4

#  NEMO stagger with same i,j index:
#  -------------
#  x----V----F
#  |    |    |
#  |    |    |
#  x----P----U
#  |    |    |
#  |    |    |
#  x----|----x

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



class NemoMesh(object) :

   def __init__(self,meshfile,first_i=None,last_i=None,first_j=None,last_j=None)  :

      # Limited subset capability for now (only first_j)
      if last_j is not None or first_i is not None or last_i is not None :
         msg="only first_j may be specified for subregion for now"
         logger.error(msg)
         raise ValueError(msg)


      self._meshfile=meshfile
      self._ncid=netCDF4.Dataset(self._meshfile,"r")

      self._is_periodic_i   = self.periodic_i  ()
      self._is_arctic_patch = self.arctic_patch()
      logger.info("Grid periodic in i  :%s"%str(self._is_periodic_i))
      logger.info("Arctic Patch        :%s"%str(self._is_arctic_patch))

      self._slices_set(first_j=first_j)


   def periodic_i(self) :
      #1) Check if periodic in i. Done by checking if there is a overlap in mesh. Currently there is a two point  overlap 
      plon=self._ncid.variables["glamt"][0,:,:]
      plat=self._ncid.variables["gphit"][0,:,:]
      per_i=True
      for j in range(plon.shape[0]) :
         per_i=per_i and (plon[j,0] == plon[j,-2]  and plat[j,0] == plat[j,-2])
      return per_i


   def arctic_patch(self):
      plon=self._ncid.variables["glamt"][0,:,:]
      plat=self._ncid.variables["gphit"][0,:,:]
      jtst=-2
      #Operate on 2nd last row. point 1 must match point idm-1, point 2 must match i2m -2. etc etc
      # For testing on plon, plat.
      #print plon[jtst,:5] ,plat[jtst,:5]
      #print plon[jtst,-5:],plat[jtst,-5:]
      # For testing on vlon, vlat. v-points move 1 step down as we move from 1 edge to another (Consider flow into p-cell)
      #print plon[-2,:5] ,plat[-2,:5]      
      #print plon[-3,-5:],plat[-3,-5:]
      idm=plon.shape[1]
      maxdist=0.
      for i in range(1,idm): # NB: 1st point is skipped. Compare 2nd with last, 3rd with 2nd last, etc etc
         i2 = idm- i 
         d=modeltools.tools.spherdist_haversine( plon[jtst,i], plat[jtst,i],plon[jtst,i2],plat[jtst,i2])
         maxdist=max(maxdist,d)
      return maxdist == 0.


   # The Nemo mesh contains repeated points - this returns the unique ones for each variable, taking into account periodic and arctic patch properties
   # NB: Arctic patch is "repeated" in a sence
   def _slices_set(self,first_i=1,last_i=-1,first_j=0,last_j=-2) :
      if not (self._is_periodic_i and self._is_arctic_patch) :
         raise ValueError("nemo mesh is not periodic with arctic patch - fixme") 
      if last_j is not -2 or first_i is not 1 or last_i is not -1 :
         msg="only first_j may be specified for subregion for now"
         logger.error(msg)
         raise ValueError(msg)

      # Slices for selected domain 
      self._sj_sel=slice(first_j,last_j)
      self._si_sel=slice(first_i,last_i)


   @property
   def variables(self) :
      return self._ncid.variables


   def __getitem__(self,item) :
      return self._ncid.variables[item]

   # nemo U at i,j -> HYCOM U at i+1,j. 
   # NB: field.shape=(jdm,idm)
   def u_to_hycom_u(self,field2d)  :
      if not self._is_periodic_i :
         raise ValueError("nemo mesh is not periodic in i - fixme") 
      return numpy.roll(field2d,1,axis=1)


   # nemo V at i,j -> HYCOM U at i,j+1. 
   # NB: field.shape=(jdm,idm)
   def v_to_hycom_v(self,field2d,extrapolate="none")  :
      myfield=numpy.copy(field2d)
      if self._is_arctic_patch :

         # For now the bottom row is just replicated
         myfield[1:,:] = myfield[:-1,:]

         ## Bottom row must be extrapolated... Some approaches for extrapolation:
         ## 1) All values on bottom row the same?
         #if extrapolate=="longitude" :
         #   dlon = numpy.mod(myfield[0,:] - numpy.roll(myfield[0,:],1) + 360.,360.)
         #   if dlon.min() == dlon.max() :
         #      logger.debug("Extrapolating longitude by replicating bottom row")
         ## 2) 
         #elif extrapolate=="latitude" :
         #   dlat = myfield[2,:]-myfield[1,:]
         #   if dlat.min() == dlat.max() :
         #      logger.debug("Extrapolating latitude ")
         #      myfield[0,:] = myfield[0,:] - dlat
         #elif extrapolate=="" :

      else :
         raise ValueError("nemo mesh is not periodic in i - fixme") 
      return myfield



   # nemo f at i,j -> HYCOM Q at i+1,j+1. 
   # NB: field.shape=(jdm,idm)
   def f_to_hycom_q(self,field,extrapolate="none")  :
      myfield=self.u_to_hycom_u(field)
      myfield=self.v_to_hycom_v(myfield)
      return myfield



   def arctic_patch_shift_up(self,field,jstep) :
      # shift field down
      if jstep != 1 :
         raise NameError("Arctic_patch_shift only with jstep=1 for now")
      field2 = numpy.copy(field)
      field2[1:,:] = field2[0:-1,:] # Shift up
      # NB:  row 0 is same as row 1 (extrapolated
      return field2

   def arctic_patch_shift_down(self,field,jstep) :
      # shift field down
      if jstep != 1 :
         raise NameError("Arctic_patch_shift only with jstep=1 for now")
      field2 = numpy.copy(field)
      field2[0:-1,:] = field2[1:,:] # Shift down
      tmp=field2[-1,:]              # Top row as top ...
      field2[-1,:] = tmp[::-1]      # .. but reversed direction
      return field2

   def periodic_i_shift_right(self,field,istep) :
      # shift field left by istep steps
      field2  = numpy.roll(field,istep,axis=1)
      return field2

      
   # Estimate partial steps and layer in u-points
   def depth_u_points(self) :
      depth     = self["hdepw" ][0,:,:] 
      mbathy    = self["mbathy"][0,:,:] 
      gdepw     = self["gdepw_0"][0,:]       # Depth of w points
      depthip1  = self.periodic_i_shift_right(depth ,-1)    # nemo values at cell i+1
      mbathyip1 = self.periodic_i_shift_right(mbathy,-1)    # nemo values at cell i+1
      depthu  =numpy.minimum(depth , depthip1)
      mbathy_u=numpy.minimum(mbathy,mbathyip1)
      e3u_ps  =numpy.zeros(depthu.shape)
      J,I= numpy.where(mbathy_u>-1)
      e3u_ps[J,I] = depthu[J,I] - gdepw[mbathy_u[J,I]]
      return mbathy_u,e3u_ps,depthu


      
   # Estimate partial steps and layer in v-points
   def depth_v_points(self) :
      depth     = self["hdepw" ][0,:,:] 
      mbathy    = self["mbathy"][0,:,:] 
      gdepw     = self["gdepw_0"][0,:]       # Depth of w points
      depthjp1  = self.arctic_patch_shift_down(depth ,1)    # nemo values at cell j+1
      mbathyjp1 = self.arctic_patch_shift_down(mbathy,1)    # nemo values at cell j+1 
      depthv  =numpy.minimum(depth ,depthjp1)
      mbathy_v=numpy.minimum(mbathy,mbathyjp1)
      e3v_ps  =numpy.zeros(depthv.shape)
      I= numpy.where(mbathy_v>-1)
      e3v_ps[I] = depthv[I] - gdepw[mbathy_v[I]]
      return mbathy_v,e3v_ps,depthv


#   @property
#   def slicei(self) : return self._slicei
#
#   @property
#   def slicej(self) : return self._slicej


   def sliced(self,field2d) :
      #print field2d.shape
      #print self._sj_sel,self._si_sel
      return field2d[self._sj_sel,self._si_sel]


