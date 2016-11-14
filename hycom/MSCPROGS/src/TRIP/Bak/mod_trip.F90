! --- -------------------------------------------------------------------
! --- River routine trip_riverweights
! --- -------------------------------------------------------------------
! --- Program to map ERA40 grid cells onto TRIP grid cells. Required
! --- by trip_riverflow.
! ---
! --- For now this routine uses ERA40 data, but it can easily be changed 
! --- to other runoff products.
! ---
! --- Output from this routine is:
! --- unformatted file containing mapping from ERA40 runoff grid -> TRIP grid
! --- -------------------------------------------------------------------
! --- Prerequisites:
! --- 1) ERA40 landmask must be available in either the directory "./ERA40/",
! ---    or in the directory set in env variable ERA40_PATH
! --- 2) TRIP data base must be available in either the directory "./Data/",
! ---    or in the directory set in env variable TRIP_PATH
! --- -------------------------------------------------------------------

module mod_trip
   implicit none
#if defined (TRIP05)
   integer, parameter :: nxd=720, nyd=360 ! grid dimensions 0.5 X 0.5 grid
   real,    parameter :: dx=0.5, dy=0.5 ! grid spacing
   character(len=*), parameter :: tfile='trip05.txt'
   character(len=720), save :: txtline
#elif defined (TRIP)
   integer, parameter :: nxd=360, nyd=180 ! grid dimensions 1 X 1 grid
   real,    parameter :: dx=1., dy=1. ! grid spacing
   character(len=*), parameter :: tfile='trip.txt'
   character(len=360), save :: txtline
#endif
   real, dimension(nxd), save :: lon
   real, dimension(nyd), save :: lat
   character(len=200), save, private :: trip_path0
   integer, dimension(nxd,nyd), save :: direction


contains

   subroutine init_trip()
   implicit none
   integer :: i,j
   character(len=200):: cenv
   logical :: ex
   integer :: ios
   ! trip05 center data points
   do j=1,nyd
      ! Flip it when we are working with it here
      !lat(j) = 90-dx/2 - (j-1)*dy 
      lat(nyd-j+1) = 90-dx/2 - (j-1)*dy 
   end do
   do i=1,nxd
      lon(i) = (i-1)*dx + dx/2
   end do


   ! Look for environment variable pointing to TRIP data dir
   call getenv('TRIP_PATH',cenv)
   if (trim(cenv)=='') then
      trip_path0='./Data/'
   else
      trip_path0=trim(cenv)//'/'
   end if
      
   inquire(exist=ex,file=trim(trip_path0)//tfile)
   if (.not. ex) then 
      print *,'TRIP database not found in '//trim(trip_path0)
      print *,'Set the enviroment variable TRIP_PATH to the location of the database'
      stop
   end if

   ! Read trip05 database, put in grid (txt file)
   print '(a)','Reading TRIP database'
   open(10,file=trim(trip_path0)//tfile,status='old',action='read',iostat=ios)
   do j=1,nyd
      read(10,*,iostat=ios) txtline

      if (ios/=0) then
         print *,'error readint trip file '//trim(trip_path0)//tfile
         stop
      end if

      do i=1,nxd
         if (txtline(i:i)=='.') then
            direction(i,j)=0
         elseif (txtline(i:i)=='+') then
            direction(i,j)=9
         else
            read(txtline(i:i),'(i1)') direction(i,j)
         end if
      end do
   end do
   close(10)
   end subroutine

end module


         



