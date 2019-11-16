module mod_bodc
! BODC variables
   integer, parameter :: bodc_nx=1321
   integer, parameter :: bodc_ny=1501
   type BODC_data
      real deep(bodc_nx,bodc_ny)
      real lat(bodc_nx,bodc_ny)
      real lon(bodc_nx,bodc_ny)
      real weight(bodc_nx,bodc_ny)
      real minlon,maxlon,minlat,maxlat
   end type BODC_data
   type(BODC_data)   :: bodc
   real, allocatable :: tmp(:,:)
   real :: dlonb,dlatb
contains

subroutine read_bodc()
! Read BODC dataset
   implicit none
   integer i,j
   logical ex
   character(len=200) :: path0


   ! Retrieve etopo5 path
   call getenv('BODC_PATH',path0)
   if (trim(path0)=='') then
      print *,'No path set for BODC. Set the environment variable BODC_PATH'
      print *,'To point to the location of the BODC bathymetry'
      call exit(1)
   end if
   path0=trim(path0)//'/'

   inquire(file=trim(path0)//'bodc.uf',exist=ex)
   if (ex) then
      print *,'BODC: reading ',trim(path0),'bodc.uf'
      open(10,file=trim(path0)//'bodc.uf',form='unformatted')
         read(10)bodc%deep,bodc%minlon,bodc%maxlon,bodc%minlat,bodc%maxlat
      close(10)
   else
      inquire(file=trim(path0)//'ammp_bathy.dat',exist=ex)
      if (.not.ex) then
         print *,'Can not find BODC file '//trim(path0)//'/ammp_bathy.dat'
         print *,'path0=',trim(path0)
         stop '(mod_bodc)'
      end if

      print *,'BODC: reading ',trim(path0),'ammp_bathy.dat'
      open(10,file=trim(path0)//'ammp_bathy.dat')
         read(10,'(4f6.1,2i5)')bodc%minlon,bodc%maxlon,bodc%minlat,bodc%maxlat,i,j
         print *,bodc%minlon,bodc%maxlon,bodc%minlat,bodc%maxlat,i,j
         allocate(tmp(bodc_nx,bodc_ny))
         read(10,'(f6.1,10f7.1)')tmp(1:bodc_nx,1:bodc_ny)
         do j=1,bodc_ny
         do i=1,bodc_nx
            bodc%deep(i,j)=tmp(i,bodc_ny-j+1)
         enddo
         enddo
      close(10)
      deallocate (tmp)
      open(10,file=trim(path0)//'bodc.uf',form='unformatted')
         write(10)bodc%deep,bodc%minlon,bodc%maxlon,bodc%minlat,bodc%maxlat
      close(10)
   endif

   dlonb=(bodc%maxlon-bodc%minlon)/float(bodc_nx-1); print '(a,2f8.5)','dlonb=',dlonb,dlonb*60.0
   dlatb=(bodc%maxlat-bodc%minlat)/float(bodc_ny-1); print '(a,2f8.5)','dlatb=',dlatb,dlatb*60.0
   do j=1,bodc_ny
   do i=1,bodc_nx
      bodc%lon(i,j)=bodc%minlon+float(i-1)*dlonb
      bodc%lat(i,j)=bodc%minlat+float(j-1)*dlatb
      bodc%weight(i,j)=min( min(float(i)/20.0, 1.0),        &
                            min(float(j)/20.0, 1.0),        &
                            min(float(bodc_nx-i)/20.0, 1.0),&
                            min(float(bodc_ny-j)/20.0, 1.0) )
      if (bodc%deep(i,j) > 0.0) then
      else
         bodc%deep(i,j)=0.0
      endif
   enddo
   enddo
   print *,'BODC min/max weight:',minval(bodc%weight),maxval(bodc%weight)
end subroutine read_bodc
end module mod_bodc
