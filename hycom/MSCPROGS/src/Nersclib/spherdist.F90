elemental function spherdist(lon1,lat1,lon2,lat2)
! --- -----------------------------------------
! --- Computes the distance between geo. pos.
! --- lon1,lat1 and lon2,lat2. 
! --- INPUT is in degrees.
! --- -----------------------------------------

   implicit none
   REAL, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees

   real*8, parameter :: invradian=0.017453292
   real*8, parameter :: rearth=6371001.0     ! Radius of earth

   real*8  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
   real*8  x1,y1,z1,x2,y2,z2                 ! Cartesian position
   real*8  dx,dy,dz,dr,dot                   ! Cartesian distances
   real :: spherdist


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
   
   !dr=acos(x1*x2+y1*y2+z1*z2)       ! Acr length
   dot=min(max(-1.,x1*x2+y1*y2+z1*z2),1.)
   dr=acos(dot)       ! Acr length
   spherdist=dr*rearth
end function spherdist
