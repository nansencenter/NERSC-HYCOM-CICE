module mod_ibcao
! IBCAO variables
   integer, parameter :: ibcao_nx=2323
   integer, parameter :: ibcao_ny=2323
   type IBCAO_data
      real deep(ibcao_nx,ibcao_ny)
      real lat(ibcao_nx,ibcao_ny)
      real lon(ibcao_nx,ibcao_ny)
      real weight(ibcao_nx,ibcao_ny)
      real minlon,maxlon,minlat,maxlat
   end type IBCAO_data
   type (IBCAO_data) ibcao
   real, allocatable :: tmp3(:,:,:)

contains
subroutine read_ibcao()
! Read IBCAO dataset http://www.ngdc.noaa.gov/mgg/bathymetry/arctic/ibcaobetagrid.html
! Grid type: Cartesian (XYZ)
! Projection: Polarstereographic, true scale at 75 °N.
! Horizontal datum: WGS 84
! Grid spacing: 2.5 x 2.5 km
! Origin: The Cartesian origin (0,0) is located at the center of the grid.
!                       This corresponds to the North Pole
! Data organization: Grid registration, row major with first coordinate at Upper left
! Grid corner coordinates: corresponding geographic coordinates in parentheses(lon lat)
! UL -2902500,2902500 (-135, 53:49:1.4687)
! UR 2902500, 2902500 (135, 53:49:1.4687)
! LL -2902500,-2902500 (-45, 53:49:1.4687)
! LR 2902500, -2902500 (45, 53:49:1.4687)
   implicit none
   integer i,j
   logical ex
   logical :: exx=.false.
   character(len=200) :: path0



   ! Retrieve etopo5 path
   call getenv('IBCAO_PATH',path0)
   if (trim(path0)=='') then
      print *,'No path set for IBCAO. Set the environment variable IBCAO_PATH'
      print *,'To point to the location of the IBCAO bathymetry'
      call exit(1)
   end if
   path0=trim(path0)//'/'

   inquire(file=trim(path0)//'ibcao.uf',exist=ex)
   if (ex) then
      open(10,file=trim(path0)//'ibcao.uf',form='unformatted')
      read(10)ibcao%deep,ibcao%lon,ibcao%lat
      close(10)
   else
      inquire(file=trim(path0)//'IBCAO_beta_lonlat',exist=ex)
      if (.not.ex) then
         print *,'Can not find IBCAO file '//trim(path0)//'/IBCAO_beta_lonlat'
         print *,'path0=',trim(path0)
         stop '(mod_ibcao)'
      end if
      open(10,file=trim(path0)//'IBCAO_beta_lonlat')
         allocate(tmp3(3,ibcao_nx,ibcao_ny))
         read(10,*)tmp3(1:3,1:ibcao_nx,1:ibcao_ny)
      close(10)
      do j=1,ibcao_ny
      do i=1,ibcao_nx
         ibcao%lon(i,j) =tmp3(1,j,i)
         ibcao%lat(i,j) =tmp3(2,j,i)
         ibcao%deep(i,j)=tmp3(3,j,i)
      enddo
      enddo
      deallocate(tmp3)
      open(10,file=trim(path0)//'ibcao.uf',form='unformatted')
         write(10)ibcao%deep,ibcao%lon,ibcao%lat
      close(10)
   endif
   ibcao%minlon=minval(ibcao%lon)
   ibcao%maxlon=maxval(ibcao%lon)
   ibcao%minlat=minval(ibcao%lat)
   ibcao%maxlat=maxval(ibcao%lat)
   do j=1,ibcao_ny
   do i=1,ibcao_nx
      ibcao%weight(i,j)=min(min(float(i)/20.0, 1.0),         &
                            min(float(j)/20.0, 1.0),         &
                            min(float(ibcao_nx-i)/20.0, 1.0),&
                            min(float(ibcao_ny-j)/20.0, 1.0) )
      if (ibcao%deep(i,j) < 0.0) then
         ibcao%deep(i,j)=-ibcao%deep(i,j)
      else
         ibcao%deep(i,j)=0.0
      endif
   enddo
   enddo
   print *,'IBCAO min/max weight:',minval(ibcao%weight),maxval(ibcao%weight)
end subroutine read_ibcao
end module mod_ibcao
