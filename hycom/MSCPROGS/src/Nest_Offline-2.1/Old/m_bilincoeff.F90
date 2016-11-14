module m_bilincoeff
contains
subroutine bilincoeff(glon,glat,nx,ny,lon,lat,ipiv,jpiv,a1,a2,a3,a4,lperiodic)
! This subroutine uses bilinear interpolation to interpolate the field
! computed by the model (MICOM) to the position defined by lon, lat
! The output is the interpolation coeffisients a[1-4]
! NB  NO locations on land.

   use mod_confmap
   use m_oldtonew
   implicit none

   integer, intent(in)  :: nx,ny
   real,    intent(in)  :: glon(nx,ny),glat(nx,ny)
   real,    intent(in)  :: lon,lat
   integer, intent(in)  :: ipiv,jpiv
   real,    intent(out) :: a1,a2,a3,a4
   logical, intent(in ) :: lperiodic

   real t,u
   real lat_1,lon_1,lat_2,lon_2,lat_n,lon_n,lat_t,lon_t
   integer :: ipib

   ! Security check
   if (lperiodic .and. ipiv==nx) then
      print *,'Bilincoeff handles periodic grids - but needs testing!'
      stop
   end if

   if (lperiodic) then
      ipib=mod(ipiv,nx)+1
   else
      ipib=ipiv+1
   end if

   call oldtonew(glat(ipiv,jpiv),glon(ipiv,jpiv),lat_1,lon_1)
   call oldtonew(glat(ipib,jpiv+1),glon(ipib,jpiv+1),lat_2,lon_2)
   call oldtonew(lat,lon,lat_n,lon_n)

!   print *,lat_1,lon_1,lat_2,lon_2,lat_n,lon_n
!   call oldtonew(glat(ipib,jpiv),glon(ipib,jpiv),lat_t,lon_t)
!   print *,lat_t,lon_t,lat_t,lon_t
!
!   call oldtonew(glat(ipiv,jpiv+1),glon(ipiv,jpiv+1),lat_t,lon_t)
!   print *,lat_t,lon_t,lat_t,lon_t
!   print *,'--------------'


   t=(lon_n-lon_1)/(lon_2-lon_1)
   u=(lat_n-lat_1)/(lat_2-lat_1)

   a1=(1-t)*(1-u)
   a2=t*(1-u)
   a3=t*u
   a4=(1-t)*u

end subroutine bilincoeff
end module m_bilincoeff
