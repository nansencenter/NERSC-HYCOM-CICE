module m_avecell
contains
subroutine avecell(deep,modlon,modlat,nocean,nland,&
                         obs,lat,lon,nxd,nyd,ilatf,ilatl,nx,ny)
! Maps a data grid to the model grid based on an averaging procedure.
use mod_confmap
!use m_oldtonew
!use m_pivotp
implicit none
integer, intent(in) :: nx,ny
real,    intent(out)   :: deep(nx,ny)       ! Computed depth
real,    intent(in)    :: modlon(nx,ny)     ! model longitudes
real,    intent(in)    :: modlat(nx,ny)     ! model latitudes
integer, intent(out)   :: nocean(nx,ny)     ! Computed depth
integer, intent(out)   :: nland(nx,ny)      ! Computed depth
integer, intent(in)    :: nxd               ! data grid i-dimension
integer, intent(in)    :: nyd               ! data grid j-dimension
real,    intent(in)    :: obs(nxd,nyd)      ! Observed depth
real,    intent(in)    :: lat(nxd,nyd)      ! data longitudes
real,    intent(in)    :: lon(nxd,nyd)      ! data latitudes
integer, intent(in)    :: ilatf             ! first lat point in data to use
integer, intent(in)    :: ilatl             ! last lat point in data to use

integer :: tmpdeep(nx,ny)

integer i,j,im,jm
real lat_n,lon_n  !,pi_4,di_inv,dj_inv,dm_inv
integer icount

! Parameters used in interpolation
   !pi_4=atan(1.0)
   !rad=pi_4/45.0
   !deg=45.0/pi_4
   !di_inv=1.0/di
   !dj_inv=1.0/dj
   !dm_inv=1.0/dm

   nland=0
   nocean=0
   tmpdeep=0.0
   do j=ilatf,ilatl
   do i=1,nxd
      call oldtonew(lat(i,j),lon(i,j),lat_n,lon_n)
      call pivotp(lon_n,lat_n,im,jm,0.5)

      if (1 <= im .and. im <= nx .and. 1 <= jm .and. jm <= ny) then
         if (obs(i,j) > 0) then
            nocean(im,jm)=nocean(im,jm)+1
            tmpdeep(im,jm)=tmpdeep(im,jm)+obs(i,j)
         else
            nland(im,jm)=nland(im,jm)+1
         endif
      endif

   enddo
   enddo


! If there are more ocean points than land points in a gridcell,
! and if the total number of data points in that cell is >= 16
! average over the ocean points rather than staying with the 

   icount=0
   do j=1,ny
      do i=1,nx
       if (nland(i,j)+nocean(i,j) >= 16) then
          if (nocean(i,j) > nland(i,j)) then
             icount=icount+1
             deep(i,j)=tmpdeep(i,j)/float(nocean(i,j))
          endif
       endif
     enddo
   enddo
   print *,'    avecell: updated ',icount,' number of grid points'

end subroutine avecell
end module m_avecell
