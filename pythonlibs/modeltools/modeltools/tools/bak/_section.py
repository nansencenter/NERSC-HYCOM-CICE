#!/usr/bin/env python
import numpy

#   ! This troutine sets up the node points (p-cell) of the grid along the section
#   subroutine section_nodepoints(lon,lat,idm,jdm,periodic)
#   use mod_sphere_tools
#   !use m_sort
#   use m_handle_err
#   use netcdf
#   implicit none
#   logical, intent(in) :: periodic
#   integer, intent(in) :: idm,jdm
#   real, dimension(idm,jdm), intent(in) :: lon,lat
#
#   real   , dimension(idm,jdm) :: rmask
#   logical, dimension(idm,jdm) :: mask
#   integer, dimension(idm,jdm) :: flagu, flagv
#   integer :: isec,i,j,ib,jb, ipnt,n
#   integer, allocatable :: indx(:)
#   real,dimension(3) :: nvec, rvec, rvec1, rvec2, cp1, cp2
#
#   character(len=3) :: css
#   integer :: idmid, jdmid, ncid, varid
#
#
#   ! Sanity check
#   if (.not.(allocated(seclon) .and. allocated(seclat) .and. &
#             allocated(secname) )) then
#      print *,'Error -- section arrays are unallocated'
#      print *, '(transport_init)'
#      call exit(1)
#   end if
#
#   allocate(ndeipiv (max_sdm,nsec))
#   allocate(ndejpiv (max_sdm,nsec))
#   allocate(ndelon  (max_sdm,nsec))
#   allocate(ndelat  (max_sdm,nsec))
#   allocate(ndeflagu(max_sdm,nsec))
#   allocate(ndeflagv(max_sdm,nsec))
#   allocate(ndedist (max_sdm,nsec))
#   allocate(sdm     (max_sec))
#   allocate(indx    (max_sdm))
#
#   ! Loop over all sections
#   do isec=1,nsec
#
#      do i=1,max_sdm
#         indx(i)=i
#      end do
#
#      ! Radius vectors for start and end points
#      rvec1=geo2cart(seclon(isec,1),seclat(isec,1))
#      rvec2=geo2cart(seclon(isec,2),seclat(isec,2))
#
#      ! Normal vector of plane defined by cross product of
#      ! 1) Vector from earth center to start of section
#      ! 2) Vector from earth center to end   of section
#      nvec=cross_product(rvec1,rvec2)
#
#      ! Now go through grid and mark all points on one side of the sphere
#      ! (i.e. all points whose (radius vecor X normal vector) is negative)
#      do j=1,jdm
#      do i=1,idm
#         
#         ! radius vector : earth center -> point
#         rvec=geo2cart(lon(i,j),lat(i,j))
#
#         ! dot product of radius vector with normal vector sets the mask
#         mask(i,j)=dot_product(nvec,rvec)<0.
#
#         ! init flags
#         flagu(i,j)=0
#         flagv(i,j)=0
#
#      end do
#      end do
#
#
#
#      ! Now calculate the node points along the hemisphere line by 
#      ! using a ``telescopic'' sum
#      do j=1,jdm
#      do i=1,idm
#
#         if (periodic) then
#            ib=mod(i,idm)+1
#         else
#            ib=min(i+1,idm)
#         end if
#         jb=min(j+1,jdm)
#
#         if (mask(i,j)) then
#            flagu(i,j) = flagu(i,j)+1
#            flagv(i,j) = flagv(i,j)+1
#            flagu(ib,j) = flagu(ib,j)-1
#            flagv(i,jb) = flagv(i,jb)-1
#         end if
#      end do
#      end do
#
#
#      ! Remove points along boundary
#      do i=1,idm
#         flagu(i,  1)=0
#         flagv(i,  1)=0
#         flagu(i,jdm)=0
#         flagv(i,jdm)=0
#      end do
#
#      if (.not.periodic)  then
#      do j=1,jdm
#         flagu(1,  j)=0
#         flagv(1,  j)=0
#         flagu(idm,j)=0
#         flagv(idm,j)=0
#      end do
#      end if
#
#
#      ! Finally reduce the number of points to those that are between
#      ! start and end points of section
#      do j=1,jdm
#      do i=1,idm
#
#         ! radius vector : earth center -> point
#         rvec=geo2cart(lon(i,j),lat(i,j))
#
#         ! Cross product start/end and this point
#         cp1=cross_product(rvec1,rvec)
#         cp2=cross_product(rvec,rvec2)
#
#         ! These must have the same sign for rvec to be between 
#         ! rvec1 and rvec2
#         if (dot_product(cp1,cp2)<0.) then
#            flagu(i,j)=0
#            flagv(i,j)=0
#         end if
#
#      end do
#      end do
#
#
#      ! Number of node points
#      sdm(isec)= count(flagu /=0 .or. flagv /=0 )
#
#
#      ! security check on length (max_sdm)
#      if (sdm(isec)>max_sdm) then
#         print *,'Section length  of '//trim(secname(isec))//' exceeds max_sdm ',sdm(isec)
#         call exit(1)
#      else if (sdm(isec)==0) then
#         print *,'Section length  of '//trim(secname(isec))//' is zero (on this grid). Remove section'
#         call exit(1)
#      end if
#
#      ! Unsorted node points (for now)
#      ipnt=1
#      do j=1,jdm
#      do i=1,idm
#         if (flagu(i,j)/=0 .or. flagv(i,j)/=0) then
#            ndelon (ipnt,isec)=lon(i,j)
#            ndelat (ipnt,isec)=lat(i,j)
#            ndeipiv(ipnt,isec)=i
#            ndejpiv(ipnt,isec)=j
#            ndeflagu(ipnt,isec)=flagu(i,j)
#            ndeflagv(ipnt,isec)=flagv(i,j)
#            ndedist(ipnt,isec)=spherdist(lon(i,j),lat(i,j), &
#               seclon(isec,1),seclat(isec,1))
#            ipnt =ipnt+1
#         end if
#      end do
#      end do
#
#      ! Sort by distance ...
#      n=sdm(isec)
#      call sort(n,ndedist(1:n,isec),indx(1:n))
#      ndelon(1:n  ,isec)=ndelon(indx(1:n)  ,isec)
#      ndelat(1:n  ,isec)=ndelat(indx(1:n)  ,isec)
#      ndeipiv(1:n ,isec)=ndeipiv(indx(1:n) ,isec)
#      ndejpiv(1:n ,isec)=ndejpiv(indx(1:n) ,isec)
#      ndeflagu(1:n,isec)=ndeflagu(indx(1:n),isec)
#      ndeflagv(1:n,isec)=ndeflagv(indx(1:n),isec)
#
#
#      write(css,'(i3.3)') isec
#      call handle_err(NF90_create('tst'//css//'.nc',NF90_CLOBBER,ncid))
#      call handle_err(NF90_DEF_DIM(ncid,'idm',idm,idmid))
#      call handle_err(NF90_DEF_DIM(ncid,'jdm',jdm,jdmid))
#      call handle_err(NF90_DEF_VAR(ncid,'mask',NF90_Float,(/idmid,jdmid/),varid))
#      call handle_err(NF90_ENDDEF(ncid))
#      rmask=0. ; where(mask) rmask=1.
#      call handle_err(NF90_PUT_VAR(ncid,varid,rmask))
#      call handle_err(NF90_CLOSE(ncid))
#
#
#
#   end do
#   end subroutine section_nodepoints


class SectionError(Exception):
   pass


class Section(object) :
   def __init__(self,waypoints_lon,waypoints_lat,grid_lon,grid_lat) :

      self._jdm,self._idm = grid_lon.shape

      if len(waypoints_lon) <> 2 or len(waypoints_lat) <> 2 :
         raise SectionError,"Only two waypoints all2oed"

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
         raise SectionError,"Section is on opposite sides of the earth"

      # Make normal vector a unit vector
      self._normal_vector= self._normal_vector / numpy.sqrt(numpy.sum(self._normal_vector**2))


      # Radius vector to all grid points
      rvec=geo2cart(grid_lon,grid_lat)

      # Now go through grid and mark all points on one side of the great circle
      tmp=dot_product_last( self._normal_vector,rvec)
      self._mask = numpy.where(tmp < 0,1,0)

      self._flagu=numpy.zeros(grid_lon.shape)
      self._flagv=numpy.zeros(grid_lon.shape)
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
      cp1=cross_product_last(self._waypoint_vectors[0],rvec)
      cp2=cross_product_last(rvec,self._waypoint_vectors[1])
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

      # Find pivot points along section. This is the simplest approach where we use cell centers 
      # TODO: find actual crossing points of grid ?
      J,I = numpy.where(numpy.logical_or(self._flagu<>0,self._flagv<>0))
      self._section_i=I
      self._section_j=J
      self._section_longitudes = grid_lon[J,I]
      self._section_latitudes  = grid_lat[J,I]
      self._section_vectors    = rvec    [J,I,:]

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
      #print self._distance_along_section 
      #print self._distance_along_section_1
      #print I
      self._section_i         =self._section_i[I]
      self._section_j         =self._section_j[I]
      self._section_longitudes=self._section_longitudes[I]
      self._section_latitudes =self._section_latitudes[I]
      self._distance_along_section_1 =self._distance_along_section_1[I]
      self._distance_along_section_2 =self._distance_along_section_2[I]
      self._distance_along_section   =self._distance_along_section[I]


   @property 
   def grid_indexes(self) :
      return self._section_i,self._section_j

   @property 
   def longitude(self) :
      return self._section_longitudes

   @property 
   def latitude(self) :
      return self._section_latitudes

   @property 
   def distance(self) :
      return self._distance_along_section

   @property 
   def mask(self) :
      return self._mask

   @property 
   def mask2(self) :
      return self._mask2

   @property 
   def flagu(self) :
      return self._flagu

   @property 
   def flagv(self) :
      return self._flagv


   def plot_section_arrows(self,plon=None,plat=None) :
      import matplotlib
      fig = matplotlib.pyplot.figure(figsize=(12,12))
      ax=fig.add_subplot(111)
      ax.hold(True)

      ax.pcolor(self._mask)

      i=numpy.arange(self._flagu.shape[1]) - 0.5
      j=numpy.arange(self._flagu.shape[0])
      x,y=numpy.meshgrid(i,j)
      I=numpy.where(self._flagu<>0)
      ax.quiver(x[I],y[I],self._flagu[I],numpy.zeros(self._flagu.shape)[I],width=.002)

      i=numpy.arange(self._flagu.shape[1])
      j=numpy.arange(self._flagu.shape[0]) - 0.5
      x,y=numpy.meshgrid(i,j)
      I=numpy.where(self._flagv<>0)
      ax.quiver(x[I],y[I],numpy.zeros(self._flagv.shape)[I], self._flagv[I],width=.002)
      ax.set_title("positive direction across section")


      self._add_gridlines(ax,plon,plat) 
      return fig


   def plot_section_mask(self,plon=None,plat=None) :
      import matplotlib
      fig = matplotlib.pyplot.figure(figsize=(12,24))

      ax=fig.add_subplot(211)
      J,I=numpy.where(self._flagu<>0)
      ax.scatter(I,J,30,self._flagu[J,I],edgecolors="face")
      self._add_gridlines(ax,plon,plat) 
      ax.set_title("u-flag; negative values mean negative\n grid direction is treated as positive direction")

      ax=fig.add_subplot(212)
      J,I=numpy.where(self._flagv<>0)
      ax.scatter(I,J,30,self._flagv[J,I],edgecolors="face")
      self._add_gridlines(ax,plon,plat) 
      ax.set_title("v-flag; negative values mean negative\n grid direction is treated as positive direction")

      return fig


   def plot_section_1d(self,plon=None,plat=None) :
      import matplotlib
      fig = matplotlib.pyplot.figure(figsize=(12,18))

      ax=fig.add_subplot(321)
      ax.plot(self._distance_along_section,self._section_longitudes)
      ax.set_title("longitude along section")

      ax=fig.add_subplot(322)
      ax.plot(self._distance_along_section,self._section_latitudes)
      ax.set_title("latitude along section")

      ax=fig.add_subplot(323)
      ax.plot(self._distance_along_section,self._section_i)
      ax.set_title("i pivot along section")

      ax=fig.add_subplot(324)
      ax.plot(self._distance_along_section,self._section_j)
      ax.set_title("j pivot along section")

      ax=fig.add_subplot(325)
      ax.plot(self._distance_along_section,self._distance_along_section_1)
      ax.set_title("distance measure 1")

      ax=fig.add_subplot(326)
      ax.plot(self._distance_along_section,self._distance_along_section_2)
      ax.set_title("distance measure 2")

      return fig


   def _add_gridlines(self,ax,plon,plat) :
      if plon is not None :
         CS=ax.contour(plon,numpy.arange(-180,180,10),colors="k")
         ax.clabel(CS, inline=1, fontsize=10,fmt="%1.1f")
      if plat is not None :
         CS = ax.contour(plat,numpy.arange(-90,90,10),colors="k")
         ax.clabel(CS, inline=1, fontsize=10,fmt="%1.1f")


#
#
#
#P=m.pcolor(x,y,sec.mask)
#m.drawcoastlines()
#pl=m.drawparallels(numpy.arange(-90.,120.,10.),labels=[1,0,0,0]) # draw parallels
#mer=m.drawmeridians(numpy.arange(0.,420.,10.),labels=[0,0,0,1]) # draw meridians



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


