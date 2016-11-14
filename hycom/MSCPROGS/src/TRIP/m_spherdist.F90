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

   !dr=SQRT(dx*dx+dy*dy+dz*dz)       !distance pytagaros
   !dr=acos(x1*x2+y1*y2+z1*z2)       ! Acr length


   ! Hrmph ... laws of small numbers. In theory this shouldnt be necessary
   dr=acos(min(max(-1.,x1*x2+y1*y2+z1*z2),1.))     ! Acr length

   spherdist=dr*rearth

end function spherdist
end module m_spherdist
