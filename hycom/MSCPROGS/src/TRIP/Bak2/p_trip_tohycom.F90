! --- -------------------------------------------------------------------
! --- River routine trip_tohycom
! --- -------------------------------------------------------------------
! --- Program to convert netcdf files produced by "trip_riverflow" to
! --- hycom forcing files
! ---
! --- Output from this routine is:
! ---   forcing.rivers.[ab]  
! --- -------------------------------------------------------------------
! --- Prerequisites:
! --- You must have run trip_riverweight + trip_riverflow before this routine
! --- For now only the climatology is produced
! --- -------------------------------------------------------------------

program trip_tohycom
   use mod_xc
   use mod_za
   use mod_grid
   use mod_confmap
   use netcdf
   use m_handle_err
   implicit none
#if defined (TRIP05)
   integer, parameter :: nxd=720, nyd=360 ! grid dimensions 0.5 X 0.5 grid
   real,    parameter :: dx=0.5, dy=0.5 ! grid spacing
   character(len=*), parameter :: tfile='trip05.txt'
   character(len=720) txtline
#elif defined (TRIP)
   integer, parameter :: nxd=360, nyd=180 ! grid dimensions 1 X 1 grid
   real,    parameter :: dx=1., dy=1. ! grid spacing
   character(len=*), parameter :: tfile='trip.txt'
   character(len=360) txtline
#endif
   real, dimension(nxd) :: lon
   real, dimension(nyd) :: lat
   real :: riv_flux(nxd,nyd,12)
   real, allocatable :: mod_riv_flux(:,:,:)

   integer :: i,j, ncid, varid
   real    :: ri, rj


   ! trip05 center data points
   do j=1,nyd
      ! Flip it when we are working with it here
      !lat(j) = 90-dx/2 - (j-1)*dy 
      lat(nyd-j+1) = 90-dx/2 - (j-1)*dy 
   end do
   do i=1,nxd
      lon(i) = (i-1)*dx + dx/2
   end do

   ! Read climatology (pointwise)
   call handle_err(nf90_open('trip_era40_clim.nc',NF90_CLOBBER,ncid))
   call handle_err(NF90_INQ_VARID(ncid,'river',varid))
   call handle_err(NF90_GET_VAR(ncid,varid,riv_flux))


   ! Init HYCOM fields
   call xcspmd()
   call zaiost()
   call get_grid()

   allocate(mod_riv_flux(idm,jdm,12))

   print *,minval(depths),maxval(depths)
   print *,nxd,nyd


   ! Go through TRIP points, and place rivers on map
   call initconfmap(idm,jdm)
   do j=1,nyd
   do i=1,nxd

      call ll2gind(lon(i),lat(j),ri,rj)
      if (floor(ri)>=1.and.floor(ri)<idm.and.floor(rj)>=1.and.floor(rj)<jdm) then
         print *,lon(i),lat(j),ri,rj
         mod_riv_flux(floor(ri),floor(rj),:)=riv_flux(i,j,:)
      end if


   end do
   end do

   end program


         



