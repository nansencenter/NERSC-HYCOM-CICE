module mod_ibcao
! IBCAO variables
   integer :: ibcao_nx=11617
   integer :: ibcao_ny=11617
   real, allocatable ::  ibcao_depth(:,:) ! reduce mem usage
   real*4, allocatable ::  tmp(:,:) ! File from IBCAO is in real*4
   real, allocatable ::  ibcao_weight(:,:)
   real, allocatable ::  ibcao_lat  (:,:)
   real, allocatable ::  ibcao_lon  (:,:)
   real ibcao_minlon,ibcao_maxlon,ibcao_minlat,ibcao_maxlat


!   type (IBCAO_data) ibcao
!   real, allocatable :: tmp3(:,:,:)

contains

subroutine read_ibcao()
   use netcdf
   implicit none
! Read IBCAO dataset http://www.ngdc.noaa.gov/mgg/bathymetry/arctic/grids/version3_0/
! Grid type: Cartesian (XYZ)
! Projection: Polarstereographic, true scale at 75 °N.
! Horizontal datum: WGS 84
! Grid spacing: .5 x .5 km
! Origin: The Cartesian origin (0,0) is located at the center of the grid.
!                       This corresponds to the North Pole
! Data organization: Grid registration, row major with first coordinate at Upper left
! Grid corner coordinates: corresponding geographic coordinates in parentheses(lon lat)
! FC -290400,2904000 (-135, 53:49:1.4687)
   integer i,j
   integer xmin,xmax,ymin,ymax
   real*4 dx,dy
   logical ex
   logical :: exx=.false.
   character(len=200) :: path0
   character(len=220) :: ibcao_grid_file
   integer :: ncid
   integer :: nvars,ndims, nc_stat
   integer :: depthvarid
   integer :: err


   call getenv('IBCAO_PATH',path0)
   if (trim(path0)=='') then
      print *,'No path set for IBCAO. Set the environment variable IBCAO_PATH'
      print *,'To point to the location of the IBCAO bathymetry'
      call exit(1)
   end if
   path0=trim(path0)//'/'
   ibcao_grid_file=trim(path0)//trim('IBCAO_V3_500m_RR.grd')

   allocate(ibcao_depth(ibcao_nx,ibcao_ny))
   allocate(tmp(ibcao_nx,ibcao_ny))
   inquire(file=trim(ibcao_grid_file),exist=ex)
   if (.not.ex) then
      print *,'Can not find IBCAO file '//ibcao_grid_file
      stop '(mod_IBCAO)'
   end if
   ! Get grid dims  here
   nc_stat=0
   err=NF90_OPEN(ibcao_grid_file,NF90_NOCLOBBER,ncid)
   if (err/=0) then
      print *,'Could not open '//trim(ibcao_grid_file)
      stop '(read_ibcao)'
   end if

   call nf90_handle_err('var-depth',nf90_inq_varid(ncid,'z',depthvarid))
   call nf90_handle_err('var-depth',NF90_GET_VAR  (ncid,depthvarid  ,tmp))

   ibcao_depth(1:ibcao_nx,1:ibcao_ny)=tmp(1:ibcao_nx,1:ibcao_ny)
   !deallocate(tmp)
   allocate(ibcao_lon(ibcao_nx,ibcao_ny))
   allocate(ibcao_lat(ibcao_nx,ibcao_ny))
   inquire(file=trim(path0)//'Lat_ibcao.uf',exist=ex)
   if (ex) then
      print *,'Openning file Lat_ibcao.uf'
      open(10,file=trim(path0)//'Lat_ibcao.uf',form='unformatted')
      read(10) ibcao_lat
      close(10)
   endif
   inquire(file=trim(path0)//'Lon_ibcao.uf',exist=ex)
   if (ex) then
      print *,'Openning file Lon_ibcao.uf'
      open(10,file=trim(path0)//'Lon_ibcao.uf',form='unformatted')
      read(10) ibcao_lon
      close(10)
   endif
   allocate(ibcao_weight(ibcao_nx,ibcao_ny))
   do j=1,ibcao_ny
   do i=1,ibcao_nx
!set the weight on the first 20 grid point gradually smoothed
      ibcao_weight(i,j)=min(min(float(i)/20.0, 1.0),         &
                            min(float(j)/20.0, 1.0),         &
                            min(float(ibcao_nx-i)/20.0, 1.0),&
                            min(float(ibcao_ny-j)/20.0, 1.0) )
      if (ibcao_depth(i,j) < 0.0) then
         ibcao_depth(i,j)=-ibcao_depth(i,j)
      else
         ibcao_depth(i,j)=0.0
      endif
   enddo
   enddo
    ibcao_minlon=minval(ibcao_lon)
    ibcao_maxlon=maxval(ibcao_lon)
    ibcao_minlat=minval(ibcao_lat)
    ibcao_maxlat=maxval(ibcao_lat)

  !print *,'Final summary depth min/max', minval(ibcao_depth),maxval(ibcao_depth)
  !print *,'Final summary lon   min/max', minval(ibcao_lon),maxval(ibcao_lon)
  !print *,'Final summary lat   min/max', minval(ibcao_lat),maxval(ibcao_lat)
  !print *,'Final summary weight min/max', minval(ibcao_weight),maxval(ibcao_weight)
end subroutine read_ibcao
subroutine nf90_handle_err(cinfo,error_code)
     use netcdf
     implicit none
      integer :: error_code
      character(len=*) :: cinfo

      if (error_code/=nf90_noerr) then
         write(6,'(a)') 'Netcdf error for '//cinfo//' error string is '// &
            nf90_strerror(error_code)
      end if
end subroutine nf90_handle_err


end module mod_ibcao
