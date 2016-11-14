module mod_etopo5
! ETOP5 variables
   integer, public, parameter :: etopo5_nx=4320
   integer, public, parameter :: etopo5_ny=2160
   type ETOPO5_data
      real deep(etopo5_nx,etopo5_ny)
      real lat(etopo5_nx,etopo5_ny)
      real lon(etopo5_nx,etopo5_ny)
      real minlon,maxlon,minlat,maxlat
   end type ETOPO5_data
   type(ETOPO5_data), public :: etopo5
   real, public, parameter :: five_min=1.0/12.0
   integer shift
   integer, allocatable :: itmp(:,:)


contains

subroutine read_etopo5()
! The 5 minute elevation data are stored in full 360 degree longitude
! circles for each 5 minutes of latitude, starting at the North Pole
! and stepping southward.  Each logical record contains $360*12=4320$
! 5-minute gridpoints.  The first record contains 4320 values for the
! elevation at 90 degrees north (all the same).  Each record starts at
! longitude 0 minutes east and ends at 359 degrees 55 minutes east. 
! there are  $180*12=2160$ records in the file.   Data values are positive
! above sealevel and negative below sealevel.

   implicit none
   character(len=200) :: path0
   integer i,j
   logical :: ex


   ! Retrieve etopo5 path
   call getenv('ETOPO5_PATH',path0)
   if (trim(path0)=='') then
      print *,'No path set for ETOPO5. Set the environment variable ETOPO5_PATH'
      print *,'To point to the location of the ETOPO5 bathymetry'
      call exit(1)
   end if
   path0=trim(path0)//'/'


   allocate(itmp(etopo5_nx,etopo5_ny))
   inquire(exist=ex,file=trim(path0)//'/DS759.2.uf')
   if (.not.ex) then
      print *,'Can not find ETOPO 5 file '//trim(path0)//'/DS759.2.uf'
      print *,'path0=',trim(path0)
      stop '(read_etopo5)'
   end if
   open(10,file=trim(path0)//'/DS759.2.uf',form='unformatted')
      read(10)itmp
   close(10)

! Copy to regular lon-lat grid (-180:180, -90:90))
   do j=1,etopo5_ny
   do i=1,etopo5_nx
      etopo5%deep(i,j)=float(itmp(i,etopo5_ny-j+1))
      etopo5%lon(i,j)=float(i-1)*five_min
      etopo5%lat(i,j)=-90.0+float(j)*five_min
      if (etopo5%deep(i,j) < 0.0) then
         etopo5%deep(i,j)=-etopo5%deep(i,j)
      else
         etopo5%deep(i,j)=0.0
      endif
   enddo
   enddo
   where (etopo5%lon >= 180.0) etopo5%lon=etopo5%lon-360.0
   shift=etopo5_nx/2
   etopo5%deep=cshift(etopo5%deep,-shift,1)
   etopo5%lon=cshift(etopo5%lon,-shift,1)
   etopo5%lat=cshift(etopo5%lat,-shift,1)
   etopo5%minlat=etopo5%lat(1,1)
   etopo5%maxlat=etopo5%lat(1,etopo5_ny)
   etopo5%minlon=etopo5%lon(1,1)
   etopo5%maxlon=etopo5%lon(etopo5_nx,1)
   print *,'ETOPO5:  bounds= ',etopo5%minlat,etopo5%maxlat,etopo5%minlon,etopo5%maxlon
   deallocate(itmp)
end subroutine read_etopo5
end module mod_etopo5
