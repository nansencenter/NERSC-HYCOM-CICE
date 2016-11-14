module m_spherdist

contains
real function spherdist(lon1,lat1,lon2,lat2)
! --- -----------------------------------------
! --- Computes the distance between geo. pos.
! --- lon1,lat1 and lon2,lat2. 
! --- INPUT is in degrees.
! --- -----------------------------------------

   implicit none
   REAL, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees

   real, parameter :: invradian=0.017453292
   real, parameter :: rearth=6371001.0     ! Radius of earth

   double precision  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
   double precision  x1,y1,z1,x2,y2,z2                 ! Cartesian position
   double precision  dx,dy,dz,dr                       ! Cartesian distances
   double precision  old


   rlon1=lon1*invradian             !lon1 in rad
   rlat1=(90.-lat1)*invradian       !90-lat1 in rad 

   rlon2=lon2*invradian             !lon2 in rad
   rlat2=(90.-lat2)*invradian       !90-lat2 in rad 

   x1= DSIN(rlat1)*DCOS(rlon1)        !x,y,z of pos 1.
   y1= DSIN(rlat1)*DSIN(rlon1)
   z1= DCOS(rlat1) 

   x2= DSIN(rlat2)*DCOS(rlon2)        !x,y,z of pos 2.
   y2= DSIN(rlat2)*DSIN(rlon2)
   z2= DCOS(rlat2) 

   dx=x2-x1                         !distances in x, y, z 
   dy=y2-y1
   dz=z2-z1

   dr=DSQRT(dx*dx+dy*dy+dz*dz)       !distance
   old=dr*rearth

   dr=Dacos(x1*x2+y1*y2+z1*z2)
   spherdist=dr*rearth
   !spherdist=old
end function spherdist

real function spherdist8(lon1,lat1,lon2,lat2)
! --- -----------------------------------------
! --- Computes the distance between geo. pos.
! --- lon1,lat1 and lon2,lat2. 
! --- INPUT is in degrees. Input is real 8 in this version
! --- -----------------------------------------

   implicit none
   real*8, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees

   real*8, parameter :: invradian=0.017453292
   real*8, parameter :: rearth=6371001.0     ! Radius of earth

   double precision  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
   double precision  x1,y1,z1,x2,y2,z2                 ! Cartesian position
   double precision  dx,dy,dz,dr                       ! Cartesian distances
   double precision  old


   rlon1=lon1*invradian             !lon1 in rad
   rlat1=(90.-lat1)*invradian       !90-lat1 in rad 

   rlon2=lon2*invradian             !lon2 in rad
   rlat2=(90.-lat2)*invradian       !90-lat2 in rad 

   x1= DSIN(rlat1)*DCOS(rlon1)        !x,y,z of pos 1.
   y1= DSIN(rlat1)*DSIN(rlon1)
   z1= DCOS(rlat1) 

   x2= DSIN(rlat2)*DCOS(rlon2)        !x,y,z of pos 2.
   y2= DSIN(rlat2)*DSIN(rlon2)
   z2= DCOS(rlat2) 

   dx=x2-x1                         !distances in x, y, z 
   dy=y2-y1
   dz=z2-z1

   dr=DSQRT(dx*dx+dy*dy+dz*dz)       !distance
   old=dr*rearth

   dr=Dacos(x1*x2+y1*y2+z1*z2)
   spherdist8=dr*rearth
   !spherdist=old
end function spherdist8

end module m_spherdist
