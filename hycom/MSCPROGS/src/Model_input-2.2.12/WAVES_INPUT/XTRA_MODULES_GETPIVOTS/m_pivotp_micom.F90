module m_pivotp_micom
implicit none

contains
  ! F. Counillon (adapted from an algorithm of Mats Bentsen) 
  ! This subroutine search the pivot point for a given observations
  ! The search is linear moving toward the neigboring grid cell that minimize
  ! the distance to the obs. This search is in the worst case in O(n). The input
  ! ipiv, jpiv corresponds to the pivot point from the previous search. If there
  ! is a kind of order in the way the observation are given the search will be
  ! very fast.
  !
  subroutine pivotp_micom(lon, lat,modlon,modlat, ipiv, jpiv, nx, ny,min_d,r_bdy)
! use m_spherdist
!  use mod_grid
   real, intent(in) ::  lon, lat
   integer, intent(in) :: nx, ny
   integer, intent(out) :: r_bdy
   real, intent(in), dimension(nx,ny) :: modlon,modlat
   integer, intent(inout) :: ipiv, jpiv
   real*8, intent(out) :: min_d
   real*8 :: d
   integer :: i, j, ito, jto
   logical :: piv_on_bdy
   real, external :: spherdist


      min_d = spherdist(modlon(ipiv,jpiv), modlat(ipiv,jpiv), lon, lat)

      ito = ipiv+1 
      jto = jpiv+1 

do while ( .not. (abs(ipiv-ito)<.001 .and. abs(jpiv - jto)<.001))
!do while (ipiv==ito .and. jpiv==jto)
!      print *,'new loop',ipiv,ito,jpiv,jto,abs(ipiv-ito),abs(jpiv - jto)

      ito = ipiv
      jto = jpiv
   !   i = max(ipiv-1,1)
      i = ipiv-1
      if (i <1) i=nx
      j = jpiv
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
!      print * ,'i-1', d,min_d,ipiv,jpiv
        ipiv = i
        jpiv = j
        min_d = d
      endif
!      i = min(ipiv+1,nx)
      i=ipiv+1
      if(i>nx) i=1
      j = jpiv
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
!     print * ,'i+1', d,min_d,ipiv,jpiv
        ipiv = i
        jpiv = j
        min_d = d
      endif
      j = max(jpiv-1,1)
      i = ipiv
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
!      print * ,'j-1', d,min_d,ipiv,jpiv
        ipiv = i
        jpiv = j
        min_d = d
      endif
      i = ipiv
      j = min(jpiv+1,ny)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
!      print * ,'j+1', d,min_d,ipiv,jpiv
        ipiv = i
        jpiv = j
        min_d = d
      endif

enddo

piv_on_bdy  = (ipiv==1 .or. jpiv==1 .or. ipiv==nx .or. jpiv==ny )
if (piv_on_bdy) then
   r_bdy = bdy_index(ipiv,jpiv,nx,ny)
else
   r_bdy = 0
end if


end subroutine pivotp_micom

integer function bdy_index(i,j,nlon,nlat)
!! this function gives points on the boundary of a 2d array
!! an index of 1 for (1,1), 2 for (2,1)
!! and continuing anticlockwise
implicit none
integer, intent(in) :: i,j,nlon,nlat
integer  :: r,r1,r2,r3,r4  !r = 1 : (i,j) = (1,1)

r1 = nlon                  !r = r1: (i,j) = (nlon,1)
r2 = nlon+(nlat-2)         !r = r2: (i,j) = (nlon,nlat-1) ok
r3 = 2*nlon+(nlat-2)       !r = r3: (i,j) = (1,nlat)
r4 = 2*nlon+2*(nlat-2)     !r = r4: (i,j) = (1,2)


if (j==1) then
   r  = i
elseif (i==nlon) then
   r  = nlon+j-1
elseif (j==nlat) then
   r  = r3+1-i
else!i==1
   r  = r4+2-j
end if
bdy_index   = r

end function bdy_index

end module m_pivotp_micom
