module mod_sphere_tools

contains
  
   ! Routine to get cartesian coordinates from geographical coordinates
   function geo2cart(lon,lat)
      real, parameter :: rad=1.7453292519943295E-02,deg=57.29577951308232
      real, intent(in) :: lon,lat
      real, dimension(3) :: geo2cart

      real :: lambda

      lambda=lat*rad
      theta=lon*rad
      geo2cart(1)=cos(lambda)*cos(theta)
      geo2cart(2)=cos(lambda)*sin(theta)
      geo2cart(3)=sin(lambda)
   end function geo2cart


   ! Routine to calculate cross product of two 3D vectors.
   function cross_product(v1,v2)
      implicit none
      real, intent(in), dimension(3) :: v1,v2
      real, dimension(3) :: cross_product

      cross_product(1) = v1(2)*v2(3) - v1(3)*v2(2)
      cross_product(2) = v1(3)*v2(1) - v1(1)*v2(3)
      cross_product(3) = v1(1)*v2(2) - v1(2)*v2(1)
   end function cross_product


   ! Routine to calculate vector norm (more precise the 2-norm)
   function norm2(vector,p)
      implicit none
      integer, intent(in) :: p
      real,    intent(in) :: vector(p)
      integer :: i
      real    :: norm2

      norm2=0.
      do i=1,p
         norm2 = norm2 + vector(i)**2
      end do
      norm2=sqrt(norm2)
   end function norm2






   ! Routine to calculate wether the point (plon,plat)  is in the box
   ! defined by crnlon,crnlat. Cnrlon/crnlat must be traversed so that 
   ! they form a convex polygon in 3D coordinates when following indices.
   ! It should work for for all regions defined by crnlon/crnlat...
   function inbox(crnlon,crnlat,npt,plon,plat)
      implicit none
      real, parameter :: rad=1.7453292519943295E-02,deg=57.29577951308232

      integer, intent(in) :: npt
      real, dimension(npt), intent(in) :: crnlon,crnlat
      real, intent(in) :: plon,plat
      logical :: inbox

      real, dimension(npt,3) :: cvec
      real, dimension(3)     :: pvec
      real, dimension(3,3)   :: rvec
      real, dimension(3)     :: nvec, nvec_prev, cprod
      integer :: i,im1,ip1
      logical :: lsign
      real    :: rotsign,old_rotsign


      ! point vector from origo
      pvec = geo2cart(plon,plat)

      !print *,crnlon
      !print *,crnlat

      ! vector to rectangle corner
      do i=1,npt
         cvec(i,:) = geo2cart(crnlon(i),crnlat(i))
      end do

      ! Traverse box boundaries -- Check that traversion is
      ! consistent and that point is in box
      lsign=.true.
      i=1
      old_rotsign=0.
      do while (i<npt+1 .and. lsign)
         
         ip1= mod(i      ,npt)+1
         im1= mod(i-2+npt,npt)+1


         ! Vectors used to span planes
         rvec(3,:) = cvec(ip1,:)
         rvec(2,:) = cvec(i  ,:)
         rvec(1,:) = cvec(im1,:)

         ! Normal vector to two spanning planes
         nvec      = cross_product(rvec(2,:),rvec(3,:))
         nvec_prev = cross_product(rvec(1,:),rvec(2,:))

         ! As we move to new planes, the cross product rotates in 
         ! a certain direction
         cprod = cross_product(nvec_prev,nvec)

         ! For anticlockwise rotation, this should be positive
         rotsign=sign(1.,dot_product(cprod,rvec(2,:)))

         ! Check that box is consistently traversed
         if (i>1 .and. rotsign * old_rotsign < 0) then
            print *,'Grid cell not consistently traversed'
            print *,'or polygon is not convex'
            stop '(inbox2)'
         end if
         old_rotsign=rotsign

         ! If this is true for all four planes, we are in grid box
         lsign = lsign .and. (dot_product(nvec,pvec)*rotsign)>0
         i=i+1

      end do
      inbox = lsign
   end function inbox








   
   ! Routine to get angle between two vectors defined in
   ! geographical coordinates. 
   function secangle(lon1,lat1,lon2,lat2)
      implicit none

      real, intent(in), dimension(2) :: lon1,lat1,lon2,lat2
      real, dimension(3) :: nx, n2, ny
      real :: cos1, cos2
      real :: secangle

      ! Normal of the planes defined by positions and origo
      nx = cross_product(geo2cart(lon1(1),lat1(1)),geo2cart(lon1(2),lat1(2)))
      n2 = cross_product(geo2cart(lon2(1),lat2(1)),geo2cart(lon2(2),lat2(2)))

      ! Normal to position 1 and vector 1 (x) -- forms rh system
      ny = cross_product(geo2cart(lon1(1),lat1(1)), nx)
      ny = ny / norm2(ny,3)

      ! Angle info 1 -- cosine of angle between planes nx and n2
      cos1 = dot_product(nx,n2)

      ! Angle info 2 -- Cosine of angle between planes ny and n2
      cos2 = dot_product(ny,n2)

      ! Angle between vectors 1 and 2
      secangle = atan2(cos1,cos2)

   end function secangle



         

   ! Intersection routine by Mats Bentsen.
   ! --- this routine computes the lat/lon coordinates for the intersection
   ! --- of the two geodesic lines which connects the lat/lon pairs a1,a2
   ! --- and b1,b2
   logical function intersect(lat_a1,lon_a1,lat_a2,lon_a2, &
                              lat_b1,lon_b1,lat_b2,lon_b2, &
                              lat_i,lon_i)
      implicit none


      real lat_a1,lon_a1,lat_a2,lon_a2, &
           lat_b1,lon_b1,lat_b2,lon_b2, &
           lat_i,lon_i
 
      real lambda,theta,                  &
           x_a1,y_a1,z_a1,x_a2,y_a2,z_a2, &
           x_b1,y_b1,z_b1,x_b2,y_b2,z_b2, &
           x_na,y_na,z_na,x_nb,y_nb,z_nb, &
           x_i,y_i,z_i,l_i,               &
           x_a,y_a,z_a,x_b,y_b,z_b,l_a,l_b,l_l,a_a,a_b
 
      real rad,deg
      parameter(rad=1.7453292519943295E-02,deg=57.29577951308232)
 
! --- transforming from spherical to cartesian coordinates
      lambda=lat_a1*rad
      theta=lon_a1*rad
      x_a1=cos(lambda)*cos(theta)
      y_a1=cos(lambda)*sin(theta)
      z_a1=sin(lambda)
 
      lambda=lat_a2*rad
      theta=lon_a2*rad
      x_a2=cos(lambda)*cos(theta)
      y_a2=cos(lambda)*sin(theta)
      z_a2=sin(lambda)
 
      lambda=lat_b1*rad
      theta=lon_b1*rad
      x_b1=cos(lambda)*cos(theta)
      y_b1=cos(lambda)*sin(theta)
      z_b1=sin(lambda)
 
      lambda=lat_b2*rad
      theta=lon_b2*rad
      x_b2=cos(lambda)*cos(theta)
      y_b2=cos(lambda)*sin(theta)
      z_b2=sin(lambda)
 
      x_na=y_a1*z_a2-y_a2*z_a1
      y_na=z_a1*x_a2-z_a2*x_a1
      z_na=x_a1*y_a2-x_a2*y_a1
 
      x_nb=y_b1*z_b2-y_b2*z_b1
      y_nb=z_b1*x_b2-z_b2*x_b1
      z_nb=x_b1*y_b2-x_b2*y_b1
 
! --- Let a1 be the vector from the center of the sphere to the point
! --- (lat_a1,lon_a1) on the sphere. Similar with vectors a2, b1 and b2.
! --- Then we compute the components and length of a vector i pointing
! --- along the intersection of the two planes spanned out by the
! --- vectors a1, a2 and b1, b2 respectively.
      x_i=y_na*z_nb-y_nb*z_na
      y_i=z_na*x_nb-z_nb*x_na
      z_i=x_na*y_nb-x_nb*y_na
 
      l_i=sqrt(x_i*x_i+y_i*y_i+z_i*z_i)
!
! --- check if i lies between a1 and a2
!
      intersect=.true.
!
! --- first find the vector a between a1 and a2 and its angle a_a to a1
! --- and a2
!
      x_a=x_a1+x_a2
      y_a=y_a1+y_a2
      z_a=z_a1+z_a2
 
      l_a=sqrt(x_a*x_a+y_a*y_a+z_a*z_a)
 
      l_l=sign(1.,l_i*l_a)*max(1.e-9,abs(l_i*l_a))
      a_a=acos(max(-1.,min(1.,x_a1*x_a2+y_a1*y_a2+z_a1*z_a2)))*0.5



! --- if the angle between i and a is greater than
! --- a_a, then intersect=.false.
      if (acos(max(-1.,min(1.,(x_i*x_a+y_i*y_a+z_i*z_a)/l_l))) &
          .gt.a_a) then

         ! --- - test the opposite directed intersection vector
         x_i=-x_i
         y_i=-y_i
         z_i=-z_i
  
         if (acos(max(-1.,min(1.,(x_i*x_a+y_i*y_a+z_i*z_a)/l_l))) &
            .gt.a_a) intersect=.false.
      endif
 


      ! do similar test for b1 and b2
      if (intersect) then
 
        x_b=x_b1+x_b2
        y_b=y_b1+y_b2
        z_b=z_b1+z_b2
 
        l_b=sqrt(x_b*x_b+y_b*y_b+z_b*z_b)
 
        l_l=sign(1.,l_i*l_b)*max(1.e-9,abs(l_i*l_b))
        a_b=acos(max(-1.,min(1.,x_b1*x_b2+y_b1*y_b2+z_b1*z_b2)))*0.5
 
        if (acos(max(-1.,min(1.,(x_i*x_b+y_i*y_b+z_i*z_b)/l_l))) &
            .gt.a_b) intersect=.false.
      endif
 
      ! represent the intersection in lat,lon coordinates
      lat_i=atan2(z_i,sqrt(x_i*x_i+y_i*y_i))*deg
      lon_i=atan2(y_i,x_i)*deg
 
   end function






elemental real function spherdist(lon1,lat1,lon2,lat2)
! --- -----------------------------------------
! --- Computes the distance between geo. pos.
! --- lon1,lat1 and lon2,lat2. 
! --- INPUT is in degrees.
! --- -----------------------------------------

   implicit none
   REAL, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees

   real, parameter :: invradian=0.017453292
   real, parameter :: rearth=6371001.0     ! Radius of earth

   real  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
   real  x1,y1,z1,x2,y2,z2                 ! Cartesian position
   real  dx,dy,dz,dr                       ! Cartesian distances


   rlon1=lon1*invradian             !lon1 in rad
   rlat1=(90.-lat1)*invradian       !90-lat1 in rad 

   rlon2=lon2*invradian             !lon2 in rad
   rlat2=(90.-lat2)*invradian       !90-lat2 in rad 

   x1= SIN(rlat1)*COS(rlon1)        !x,y,z of pos 1.
   y1= SIN(rlat1)*SIN(rlon1)
   z1= COS(rlat1) 

   x2= SIN(rlat2)*COS(rlon2)        !x,y,z of pos 2.
   y2= SIN(rlat2)*SIN(rlon2)
   z2= COS(rlat2) 

   dx=x2-x1                         !distances in x, y, z 
   dy=y2-y1
   dz=z2-z1

   dr=SQRT(dx*dx+dy*dy+dz*dz)       !distance pytagaros
   dr=acos(x1*x2+y1*y2+z1*z2)       ! Acr length

   spherdist=dr*rearth

end function spherdist

end module mod_sphere_tools
