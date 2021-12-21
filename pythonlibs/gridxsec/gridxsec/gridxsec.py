#!/usr/bin/env python
import numpy

#TODO: Add Sections going along lines in projections? Or along lines in xy-space?
#TODO: Gather section methods and make subclasses

class SectionError(Exception):
   pass


class SectionBase(object) :

   @property 
   def mask(self) :
      """ Mask used to determine section points and transport masks """
      return self._mask

   @property 
   def mask2(self) :
      """ Mask used to determine points within section wedge """
      return self._mask2

   @property 
   def flagu(self) :
      """ Transport mask in u-direction """
      return self._flagu

   @property 
   def flagv(self) :
      """ Transport mask in v-direction """
      return self._flagv

   @property 
   def grid_indexes(self) :
      """ indices along section """
      return self._section_i,self._section_j

   @property 
   def distance(self) :
      """ Distance along section """
      return self._distance_along_section

   @property 
   def longitude(self) :
      """ Longitudes along section """
      return self._section_longitudes

   @property 
   def latitude(self) :
      """ Latitudes along section """
      return self._section_latitudes

   def plot_section_arrows(self,plon=None,plat=None) :
      """ Plots mask and transport directions """
      import matplotlib
      fig = matplotlib.pyplot.figure(figsize=(12,12))
      ax=fig.add_subplot(111)
      ax.hold(True)
      #
      i=numpy.arange(self._flagu.shape[1])-0.5
      j=numpy.arange(self._flagu.shape[0])-0.5
      x,y=numpy.meshgrid(i,j)
      ax.pcolormesh(i,j,self._mask,edgecolors=".5")
      #
      i=numpy.arange(self._flagu.shape[1]) - 0.5
      j=numpy.arange(self._flagu.shape[0])
      x,y=numpy.meshgrid(i,j)
      I=numpy.where(self._flagu!=0)
      ax.quiver(x[I],y[I],self._flagu[I],numpy.zeros(self._flagu.shape)[I],width=.005,color="y")
      #
      i=numpy.arange(self._flagu.shape[1])
      j=numpy.arange(self._flagu.shape[0]) - 0.5
      x,y=numpy.meshgrid(i,j)
      I=numpy.where(self._flagv!=0)
      ax.quiver(x[I],y[I],numpy.zeros(self._flagv.shape)[I], self._flagv[I],width=.005,color="y")
      #
      ax.set_xlim(I[1].min(),I[1].max())
      ax.set_ylim(I[0].min(),I[0].max())
      ax.set_title("positive direction across section")
      self._add_gridlines(ax,plon,plat) 
      return fig


   def plot_section_mask(self,plon=None,plat=None) :
      import matplotlib
      fig = matplotlib.pyplot.figure(figsize=(12,24))
      ax=fig.add_subplot(211)
      J,I=numpy.where(self._flagu!=0)
      ax.scatter(I,J,30,self._flagu[J,I],edgecolors="face")
      self._add_gridlines(ax,plon,plat) 
      ax.set_title("u-flag; negative values mean negative\n grid direction is treated as positive direction")
      ax=fig.add_subplot(212)
      J,I=numpy.where(self._flagv!=0)
      ax.scatter(I,J,30,self._flagv[J,I],edgecolors="face")
      self._add_gridlines(ax,plon,plat) 
      ax.set_title("v-flag; negative values mean negative\n grid direction is treated as positive direction")
      return fig


   def plot_section_1d(self,plon=None,plat=None) :
      import matplotlib
      fig = matplotlib.pyplot.figure(figsize=(12,18))
#
      ax=fig.add_subplot(321)
      ax.plot(self._distance_along_section,self._section_longitudes)
      ax.set_title("longitude along section")
#
      ax=fig.add_subplot(322)
      ax.plot(self._distance_along_section,self._section_latitudes)
      ax.set_title("latitude along section")
#
      ax=fig.add_subplot(323)
      ax.plot(self._distance_along_section,self._section_i)
      ax.set_title("i pivot along section")
#
      ax=fig.add_subplot(324)
      ax.plot(self._distance_along_section,self._section_j)
      ax.set_title("j pivot along section")
#
      ax=fig.add_subplot(325)
      ax.plot(self._distance_along_section,self._distance_along_section_1)
      ax.set_title("distance measure 1")
#
      ax=fig.add_subplot(326)
      ax.plot(self._distance_along_section,self._distance_along_section_2)
      ax.set_title("distance measure 2")
#
      return fig


   def _add_gridlines(self,ax,plon,plat) :
      if plon is not None :
         CS=ax.contour(plon,numpy.arange(-180,180,10),colors="k")
         ax.clabel(CS, inline=1, fontsize=10,fmt="%1.1f")
      if plat is not None :
         CS = ax.contour(plat,numpy.arange(-90,90,10),colors="k")
         ax.clabel(CS, inline=1, fontsize=10,fmt="%1.1f")


class Section(SectionBase) :
   """Class describes the grid points, longitudes, latitues and distance when going along a section on a 2D grid. It
   also estimates flags that can be used to calculate transport estimates across a section.

   This version uses great circles to estimate the section on the model grid.
   
   """


   def __init__(self,waypoints_lon,waypoints_lat,grid_c_lon,grid_c_lat) :
      """Initializes instance. Input :
         waypoints_lon : waypoints along section (Only 2 points for now)
         waypoints_lat : waypoints along section (Only 2 points for now)
         grid_c_lon    : longitude of cell centers
         grid_c_lat    : latitude of  cell centers

      __init__ calculates :
         self._mask  : The mask of points in the grid that are involved in transport calculations
         self._flagu : Denotes if positive u-velocity contribute positively or negatively to transport across section
         self._flagv : Denotes if positive v-velocity contribute positively or negatively to transport across section
      By calling find_section_pivots, __init__ also calculates the section indexes, longitudes, latitude and distance along section
       """

      self._grid_c_lon=grid_c_lon
      self._grid_c_lat=grid_c_lat
      self._jdm,self._idm = self._grid_c_lon.shape


      if len(waypoints_lon) !=2 or len(waypoints_lat) != 2 :
         raise SectionError("Only two waypoints all2oed")

      ########### The following assumes crossections go along great circles ###############

      # Radius vectors for start and end points
      self._waypoints_lon=waypoints_lon
      self._waypoints_lat=waypoints_lat
      self._waypoint_vectors=geo2cart(self._waypoints_lon,self._waypoints_lat)

      # Normal vector of plane defined by cross product of
      # 1) Vector from earth center to start of section
      # 2) Vector from earth center to end   of section
      self._normal_vector=numpy.cross(self._waypoint_vectors[0,:],self._waypoint_vectors[1,:])

      # The section plane is indeterminate if waypoints and origo is on the same line. Abort
      if numpy.sqrt(numpy.sum(self._normal_vector**2)) < 1e-2:
         raise SectionError("Section is on opposite sides of the earth")

      # Make normal vector a unit vector
      self._normal_vector= self._normal_vector / numpy.sqrt(numpy.sum(self._normal_vector**2))

      # Radius vector to all grid points
      self._grid_c_rvec=geo2cart(self._grid_c_lon,self._grid_c_lat)

      # Now go through grid and mark all points on one side of the great circle
      tmp=dot_product_last( self._normal_vector,self._grid_c_rvec)
      self._mask = numpy.where(tmp < 0,1,0)
      self._flagu=numpy.zeros(self._grid_c_lon.shape)
      self._flagv=numpy.zeros(self._grid_c_lon.shape)
      J,I=numpy.where(self._mask==1)

      # TODO: handleperiodic grids
      Ip1 = numpy.minimum(I+1,self._idm-1)
      Jp1 = numpy.minimum(J+1,self._jdm-1)

      # Now calculate the node points along the hemisphere line by  using a ``telescopic'' sum
      self._flagu[J,I  ] = self._flagu[J,I  ] + 1
      self._flagu[J,Ip1] = self._flagu[J,Ip1] - 1
      self._flagv[J,I  ] = self._flagv[J,I  ] + 1
      self._flagv[Jp1,I] = self._flagv[Jp1,I] - 1

      # Reduce the number of points to those between section endpoints
      # Cross product first section point X grid vectors
      cp1=cross_product_last(self._waypoint_vectors[0],self._grid_c_rvec)
      cp2=cross_product_last(self._grid_c_rvec,self._waypoint_vectors[1])
      #self._mask2 = dot_product_last(cp1,cp2)
      tmp  = dot_product_last(cp1,self._normal_vector)
      tmp2 = dot_product_last(cp2,self._normal_vector)
      self._mask2 = numpy.minimum(tmp,tmp2)
    
      # Modify mask and u/v flags
      self._flagu[self._mask2<0.] = 0.
      self._flagv[self._mask2<0.] = 0.
      self._mask [self._mask2<0.] = 0.

      # Remove boundary points
      # TODO: handleperiodic grids
      self._flagu[0,:]=0
      self._flagu[-1,:]=0
      self._flagu[:,0]=0
      self._flagu[:,-1]=0
      self._flagv[0,:]=0
      self._flagv[-1,:]=0
      self._flagv[:,0]=0
      self._flagv[:,-1]=0

      # Locate points along section. This is the simplest approach
      self.find_section_pivots()


   def find_section_pivots(self) :
      """Sets the points defining the section.  The following attributes of this instance is set :
       _section_i          : i-index defining the section
       _section_j          : j-index defining the section
       _section_longitudes : longitudes along the section (extracted from :grid_c_lon using i and j indexes)
       _section_latitudes  : latitudes  along the section (extracted from _grid_c_lat using i and j indexes)
       _section_distance   : Cumulative distance along section
       These attributes are sorted using the section distance.
       """


      # Find pivot points along section. This is the simplest approach where we use cell centers  to plot section
      J,I = numpy.where(numpy.logical_or(self._flagu != 0,self._flagv != 0))
      self._section_i=I
      self._section_j=J
      self._section_longitudes = self._grid_c_lon  [J,I]
      self._section_latitudes  = self._grid_c_lat  [J,I]
      self._section_vectors    = self._grid_c_rvec [J,I,:]

      # Approach 1) Haversine for distance. TODO: Will not work when section > 180 degrees...
      self._distance_along_section_1 = [ haversine(e[0],e[1],self._waypoints_lon[0],self._waypoints_lat[0]) for e in 
            zip(self._section_longitudes,self._section_latitudes) ]
      self._distance_along_section_1 = numpy.array(self._distance_along_section_1)

      # Approach 2) Angle in projected plane
      #x-coord in projected plane:
      xvec = self._waypoint_vectors[0]

      #y-coord in projected plane:
      yvec = numpy.cross(self._normal_vector,xvec) # TODO: check if on same line with origo!

      # Projection of section vectors onto plane
      xcomp = dot_product_last(self._section_vectors,xvec)
      ycomp = dot_product_last(self._section_vectors,yvec)

      # Angle in plane
      angle = numpy.arctan2(ycomp,xcomp)

      # Let angle be in 0 to 2pi
      angle[angle<0] = angle[angle<0] + 2.*numpy.pi
      self._distance_along_section_2 = angle * 6371000
      
      self._distance_along_section = numpy.copy(self._distance_along_section_2)
      I=numpy.argsort(self._distance_along_section)

      self._section_i         =self._section_i[I]
      self._section_j         =self._section_j[I]
      self._section_longitudes=self._section_longitudes[I]
      self._section_latitudes =self._section_latitudes[I]
      self._distance_along_section_1 =self._distance_along_section_1[I]
      self._distance_along_section_2 =self._distance_along_section_2[I]
      self._distance_along_section   =self._distance_along_section[I]

   
   def _find_intersection_helper(self,grid_ll_lon,grid_ll_lat,case) :

      rvec=geo2cart(grid_ll_lon,grid_ll_lat)
      fac=dot_product_last( self._normal_vector,rvec)

      # Case 1 : great circle intersects between lower left and upper left corner of grid cell
      if case == 1 :
         J,I=numpy.where(fac * numpy.roll(fac,-1,axis=0) < 0)
         K=numpy.where(numpy.logical_and(J != fac.shape[0]-1,I != fac.shape[1]-1))
         I=I[K]
         J=J[K]

         #Plane defined by ll_corner, upper right corner and origo
         rvec1=geo2cart(grid_ll_lon[J,I]  ,grid_ll_lat[J,I])
         rvec2=geo2cart(grid_ll_lon[J+1,I],grid_ll_lat[J+1,I])

      # Case 2 : great circle intersects between lower left and lower right corner of grid cell
      elif case == 2 :
         J,I=numpy.where(fac * numpy.roll(fac,-1,axis=1) < 0)
         K=numpy.where(numpy.logical_and(J != fac.shape[0]-1,I != fac.shape[1]-1))
         I=I[K]
         J=J[K]

         #Plane defined by ll_corner, upper right corner and origo
         rvec1=geo2cart(grid_ll_lon[J,I]  ,grid_ll_lat[J,I])
         rvec2=geo2cart(grid_ll_lon[J,I+1],grid_ll_lat[J,I+1])

      nvec=cross_product_last(rvec1,rvec2)

      # Cross product between section normal and above normal yields cartesian vector
      tmp = cross_product_last(nvec,self._normal_vector)

      # Not 100% sure at this stage what side of the earth is on. Make sure we are correctly located...
      tmp = numpy.transpose( 
            numpy.sign(dot_product_last(tmp,rvec[J,I]))*tmp.transpose()
            )
      
      # Convert cartesian vector to longitude and latitude
      lo,la = cart2geo(tmp)

      # Reduce the number of points to those between section endpoints.
      #TODO: Use waypoints?
      K=numpy.where(self._mask2[J,I]>0.)
      lo=lo[K]
      la=la[K]
      I=I[K]
      J=J[K]

      return I,J,lo,la,numpy.ones(lo.shape)*case



   def find_intersection(self,grid_ll_lon,grid_ll_lat) :
      """Sets the points defining the section. 
         grid_ll_lon    : longitude of cell ll corner
         grid_ll_lat    : latitude of  cell ll corner
      The following is returned
         section_i          : i-index defining the section
         section_j          : j-index defining the section
         section_longitudes : longitudes along the section (extracted from :grid_c_lon using i and j indexes)
         section_latitudes  : latitudes  along the section (extracted from _grid_c_lat using i and j indexes)
         section_case       : 
         section_distance   : Cumulative distance along section
       These attributes are sorted using the section distance.

       This method is more precise than find_pivots, and calculates crossing between grid edges and section line.
       It will give a section that appears straighter on a map if you plot section_longitudes and section_latitudes.
       """

      # TODO: Find weights  for interpolation along section

      # Case 1 : great circle intersects between lower left and upper left corner of grid cell
      I1,J1,lo1,la1,case1 = self._find_intersection_helper(grid_ll_lon,grid_ll_lat,1)

      # Case 2 : great circle intersects between lower left and lower right corner of grid cell
      I2,J2,lo2,la2,case2 = self._find_intersection_helper(grid_ll_lon,grid_ll_lat,2)


      # Create return vectors
      I=numpy.append(I1,I2)
      J=numpy.append(J1,J2)
      lo=numpy.append(lo1,lo2)
      la=numpy.append(la1,la2)
      case=numpy.append(numpy.ones(lo1.shape),numpy.ones(lo2.shape)*2)

      # Sort according to distance
      distance = [ haversine(e[0],e[1],self._waypoints_lon[0],self._waypoints_lat[0]) for e in 
            zip(lo,la) ]
      distance = numpy.array(distance)
      K=numpy.argsort(distance)
      lo=lo[K]
      la=la[K]
      I=I[K]
      J=J[K]
      case=case[K]
      distance=distance[K]

      #return I1,J1,I2,J2
      return I,J,lo,la,case,distance







class SectionIJSpace(SectionBase) :
   """Class describes the grid points, longitudes, latitues and distance when going along a section on a 2D grid. It
   also estimates flags that can be used to calculate transport estimates across a section.

   This version uses index space to estimate the section on the model grid.
   
   """
   def __init__(self,waypoints_i,waypoints_j,grid_c_lon,grid_c_lat) :
      """Initializes instance. Input :
         waypoints_i : waypoints along section (Only 2 points for now)
         waypoints_j : waypoints along section (Only 2 points for now)
         grid_c_lon  : longitude of cell centers
         grid_c_lat  : latitude of  cell centers

      __init__ calculates :
         self._mask  : The mask of points in the grid that are involved in transport calculations
         self._flagu : Denotes if positive u-velocity contribute positively or negatively to transport across section
         self._flagv : Denotes if positive v-velocity contribute positively or negatively to transport across section
      __init__ also calculates the section indexes, longitudes, latitude and distance along section
       """

      self._grid_c_lon = grid_c_lon
      self._grid_c_lat = grid_c_lat
      jishape = self._grid_c_lon.shape
      self._jdm,self._idm = jishape

      if len(waypoints_i) != 2 or len(waypoints_j) != 2 :
         raise SectionError("Only two waypoints all2oed")

      self._waypoints_x=waypoints_i
      self._waypoints_y=waypoints_j

      # Waypoint vector in IJ-plane
      v1=numpy.array([self._waypoints_x[1]-self._waypoints_x[0],self._waypoints_y[1]-self._waypoints_y[0],0])

      # Vector from first point to all I,J points
      Jg,Ig=numpy.meshgrid(range(jishape[1]),range(jishape[0]))
      Jg=Jg.flatten()
      Ig=Ig.flatten()
      v2=numpy.zeros((Ig.size,3))
      v3=numpy.zeros((Ig.size,3))
      v2[:,0]=Ig-self._waypoints_x[0]
      v2[:,1]=Jg-self._waypoints_y[0]
      v2[:,2]=0.
      v3[:,0]=Ig-self._waypoints_x[1]
      v3[:,1]=Jg-self._waypoints_y[1]
      v3[:,2]=0.

      # Angle all points in IJ space and vector in 
      self._normal_vector=numpy.cross(v1,v2)

      # Now go through grid and mark all points on one side of waypoint vector
      #self._mask = numpy.where(self._normal_vector[:,2] < 0,1,0)
      self._mask = numpy.where(self._normal_vector[:,2] <= 0,True,False)
      self._mask.shape=tuple(jishape)


      # TODO: handleperiodic grids
      self._flagu=numpy.zeros(self._mask.shape)
      self._flagv=numpy.zeros(self._mask.shape)
      J,I=numpy.where(self._mask)
      Ip1 = numpy.minimum(I+1,self._idm-1)
      Jp1 = numpy.minimum(J+1,self._jdm-1)

      ## Both of I,J and Ip1,Jp1 locations must be inside mask?
      #tmp=numpy.where(numpy.logical_and(self._mask[J,I],self._mask[Jp1,Ip1]))
      #J=J[tmp]
      #I=I[tmp]
      #Jp1=Jp1[tmp]
      #Ip1=Ip1[tmp]

      self._flagu[J,I  ] = self._flagu[J,I  ] + 1
      self._flagu[J,Ip1] = self._flagu[J,Ip1] - 1
      self._flagv[J,I  ] = self._flagv[J,I  ] + 1
      self._flagv[Jp1,I] = self._flagv[Jp1,I] - 1

      # Remove boundary points
      # TODO: handleperiodic grids
      self._flagu[0,:]=0
      self._flagu[-1,:]=0
      self._flagu[:,0]=0
      self._flagu[:,-1]=0
      self._flagv[0,:]=0
      self._flagv[-1,:]=0
      self._flagv[:,0]=0
      self._flagv[:,-1]=0

      # Reduce the number of points to those between section endpoints
      tmp  = dot_product_last(v1,v2)
      tmp2 = dot_product_last(v3,-v1)
      self._mask2 = tmp*tmp2 > 0.
      self._mask2.shape=tuple(jishape)

      # Modify mask and u/v flags
      self._flagu[~self._mask2] = 0.
      self._flagv[~self._mask2] = 0.
      #self._mask [~self._mask2] = 0.

      # Find pivot points along section. This is the simplest approach where we use cell centers 
      # TODO: find actual crossing points of grid ?
      I,J = numpy.where(numpy.logical_or(self._flagu!=0,self._flagv!=0))
      self._section_i=I
      self._section_j=J

      # Approach 1) Haversine for distance. TODO: Will not work when section > 180 degrees...
      self._section_longitudes     = self._grid_c_lon[J,I]
      self._section_latitudes      = self._grid_c_lat[J,I]
      #
      self._distance_along_section = numpy.sqrt((self._section_i-self._waypoints_x[0])**2 + (self._section_j-self._waypoints_y[0])**2) # Index space
      I=numpy.argsort(self._distance_along_section)
      self._section_i              = self._section_i[I]
      self._section_j              = self._section_j[I]
      self._distance_along_section = self._distance_along_section[I]
      self._section_longitudes     = self._section_longitudes[I]
      self._section_latitudes      = self._section_latitudes[I]


      # Cumulative distance (nb : includes zig-zag)
      self._distance_along_section = numpy.zeros(I.size)
      for k in range(1,I.size) :
         self._distance_along_section[k] = self._distance_along_section[k-1] + haversine(
               self._section_longitudes[k]  ,self._section_latitudes[k]  ,
               self._section_longitudes[k-1],self._section_latitudes[k-1]
               )



   def _find_intersection_helper(self,grid_ll_lon,grid_ll_lat,case) :

      x3= self._waypoints_x[0]
      x4= self._waypoints_x[1]
      y3= self._waypoints_y[0]
      y4= self._waypoints_y[1]

      # Case 1 : Section intersects lines where i is constant
      if case == 1 :
         x1=numpy.arange(self._idm)-.5
         x2=numpy.arange(self._idm)-.5
         y1=numpy.zeros(x1.shape)
         y2=numpy.ones(x1.shape)*self._jdm-1

      # Case 2 : great circle intersects between lower left and lower right corner of grid cell
      elif case == 2 :
         y1=numpy.arange(self._jdm) -.5
         y2=numpy.arange(self._jdm) -.5
         x1=numpy.zeros(y1.shape)
         x2=numpy.ones (y1.shape)*self._idm-1

      den = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)


      # TODO: Weight a tolerance by grid dimensions here...
      K   = numpy.abs(den)  == 0
      # The way the above is set up - either all are zero or none are zero :
      if numpy.count_nonzero(K) > 0 :
         return [],[],[],[],[]
      else :
         pass
      py = (x1*y2-y1*x2)*(x3-x4) - (x1 - x2) * (x3*y4-y3*x4) 
      py = py / den
      px = (x1*y2-y1*x2)*(y3-y4) - (y1 - y2) * (x3*y4-y3*x4) 
      px = px / den

      # Only points on grid allowed
      I=(px+.5).astype(int) # corner cell fraction -> center cell index
      J=(py+.5).astype(int) # corner cell fraction -> center cell index
      K=(I>=0) * (I<self._idm) * (J>=0) * (J<self._jdm)
      I=I[K]
      J=J[K]
      px=px[K]
      py=py[K]

      # Only points within section allowed
      #TODO: Use waypoints?
      K=numpy.where(self._mask2[J,I]>0.)
      I=I[K]
      J=J[K]
      px=px[K]
      py=py[K]

      return I,J,px,py,numpy.ones(px.shape)*case




   def find_intersection(self,grid_ll_lon,grid_ll_lat) :
      """Sets the points defining the section. 
         grid_ll_lon    : longitude of cell ll corner
         grid_ll_lat    : latitude of  cell ll corner
      The following is calculated and returned :
         section_i          : i-index defining the section
         section_j          : j-index defining the section
         section_px         : coordinates along the section 
         section_py         : coordinates along the section 
         section_case       : ...
         section_distance   : Cumulative distance along section
       These attributes are sorted using the section distance.

       This method is more precise than find_pivots, and calculates crossing between grid edges and section line.
       It will give a section that appears straighter on a map if you plot section_longitudes and section_latitudes.
       """

      # Case 1 : great circle intersects between lower left and upper left corner of grid cell
      I1,J1,px1,py1,case1 = self._find_intersection_helper(grid_ll_lon,grid_ll_lat,1)

      # Case 2 : great circle intersects between lower left and lower right corner of grid cell
      I2,J2,px2,py2,case2 = self._find_intersection_helper(grid_ll_lon,grid_ll_lat,2)

      # Create return vectors
      I=numpy.append(I1,I2)
      J=numpy.append(J1,J2)
      px=numpy.append(px1,px2)
      py=numpy.append(py1,py2)
      case=numpy.append(case1,case2)

      # Sort according to distance
      distance = [ numpy.sqrt( (e[0]-self._waypoints_x[0])**2 + (e[1]-self._waypoints_y[0])**2) 
         for e in zip(px,py) ]
      distance = numpy.array(distance)
      K=numpy.argsort(distance)
      px=px[K]
      py=py[K]
      I=I[K]
      J=J[K]
      case=case[K]
      distance=distance[K]

      #TODO: Find bilin weights 

      return I,J,px,py,case,distance


def geo2cart(lon,lat) :
   deg2rad = numpy.pi / 180.
   lmbda=numpy.array(lat)*deg2rad
   theta=numpy.array(lon)*deg2rad
   tmp1=numpy.cos(lmbda)*numpy.cos(theta)
   tmp2=numpy.cos(lmbda)*numpy.sin(theta)
   tmp3=numpy.sin(lmbda)
   tmp1.shape = tuple(list(tmp1.shape)+[1])
   tmp2.shape = tuple(list(tmp2.shape)+[1])
   tmp3.shape = tuple(list(tmp3.shape)+[1])
   t=numpy.concatenate((tmp1,tmp2,tmp3),axis=-1)
   return t

def cart2geo(cvec) :
   R = numpy.sqrt(numpy.sum(cvec**2,axis=1))
   #print cvec[:,2]
   lat=numpy.rad2deg(
         numpy.arcsin (cvec[:,2]/R)  
         )
   lon=numpy.rad2deg(
         numpy.arctan2(cvec[:,1],cvec[:,0])
         )
   return lon,lat



def cross_product_last(v1,v2) :
   tmp1 =  v1[...,1] * v2[...,2] - v1[...,2]*v2[...,1]
   tmp2 =  v1[...,2] * v2[...,0] - v1[...,0]*v2[...,2]
   tmp3 =  v1[...,0] * v2[...,1] - v1[...,1]*v2[...,0]
   tmp1.shape = tuple(list(tmp1.shape)+[1])
   tmp2.shape = tuple(list(tmp2.shape)+[1])
   tmp3.shape = tuple(list(tmp3.shape)+[1])
   t=numpy.concatenate((tmp1,tmp2,tmp3),axis=-1)
   #print "cross_product_last :",v1.shape,v2.shape,t.shape
   return t



def dot_product_last(v1,v2) :
   #print "dot_product_last :",v1.shape,v2.shape
   tmp  =  v1[...,0] * v2[...,0]
   tmp +=  v1[...,1] * v2[...,1]
   tmp +=  v1[...,2] * v2[...,2]
   return tmp

def haversine(lon1,lat1,lon2,lat2) :
   deg2rad = numpy.pi / 180.
   dlon=lon2-lon1
   dlat=lat2-lat1
   dlon=dlon*deg2rad
   dlat=dlat*deg2rad
   #a=numpy.sin(dlat/2.)**2 + numpy.cos(lon1*deg2rad) * numpy.cos(lon2*deg2rad) * numpy.sin(dlon/2.)**2
   a=numpy.sin(dlat/2.)**2 + numpy.cos(lat1*deg2rad) * numpy.cos(lat2*deg2rad) * numpy.sin(dlon/2.)**2
   #c=2.*numpy.arctan2(numpy.sqrt(a),numpy.sqrt(1.-a))
   c=2.*numpy.arcsin(numpy.sqrt(a))
   d=6371000.*c
   return d




