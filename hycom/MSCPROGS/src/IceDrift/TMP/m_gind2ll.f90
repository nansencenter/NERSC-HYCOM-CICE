module m_gind2ll
contains
subroutine gind2ll(ipiv,jpiv,lon,lat)
! This subroutine computes the pivot point of each of the observations
! in the temporary array tmpobs of type observation. The pivot point
! is the biggest i and the biggest j, (i,j) is the computation points/
! the grid, that is less than the position of the observation.

   use mod_confmap
   use m_newtoold
   implicit none

   real, intent(out) ::  lon,lat
   real, intent(in) :: ipiv,jpiv

   real tmptan,tmp
   real lontmp
   real lon_n,lat_n

   !print '(a,2i5,4f10.2)','pivotp:',ipiv,jpiv,lontmp,lon,wlim,di
   !ipiv=int((lontmp-wlim)/di)+1
   lon_n = (ipiv-1)*di + wlim

!   if (mercator) then
!      if (abs(lat) < 89.999) then
!         tmptan=tan(0.5*rad*lat+0.25*pi_1)
!         jpiv=int( (log(tmptan)-slim*rad)/(rad*dj) ) +1
!      else
!         jpiv=-999
!      endif 
!   else
!      jpiv=int((lat-slim)/dj)+1
!   endif

   if (mercator) then
      tmptan = (jpiv-1) * (rad*dj) + slim*rad
      tmptan = exp(tmptan)
      lat_n=(atan(tmptan)-0.25*pi_1)*2/rad
   else
      lat_n = (jpiv-1) *dj + slim
   endif

   call newtoold(lat_n,lon_n,lat,lon)


! Inverse transformation to check pivot point jpiv
!   tmp=slim+(jpiv-1)*dj
!   tmp=(2.*atan(exp(tmp*rad))-pi_1*.5)*deg

end subroutine gind2ll
end module m_gind2ll
