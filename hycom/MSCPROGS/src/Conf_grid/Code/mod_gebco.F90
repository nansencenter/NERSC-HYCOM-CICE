module mod_gebco

   ! GEBCO variables - use 2 minute to avoid filling memory
   ! TODO: use real*4 variables
#define GEBCO_2_MIN
#if defined (GEBCO_2_MIN)
   integer :: gebco_nx=10800
   integer :: gebco_ny=5400
#else
   integer :: gebco_nx=21600
   integer :: gebco_ny=10800
#endif
   real, allocatable ::  gebco_depth(:,:) ! reduce mem usage
   real, allocatable ::  gebco_lat  (:,:)
   real, allocatable ::  gebco_lon  (:,:)
   real gebco_minlon,gebco_maxlon,gebco_minlat,gebco_maxlat
   real, public, parameter :: one_min=1.0/60.0
   integer shift

   real :: gebco_flat, gebco_flon
   real :: gebco_dlat, gebco_dlon
   integer :: readchunk

   character(len=*), parameter :: gebco_grid_file1= 'gridone.grd'



contains

   subroutine read_gebco()
      use netcdf
      implicit none
      integer i,j,j2
      character(len=200):: path0

      integer :: ncid
      integer :: nvars,ndims
      integer :: nc_stat
      integer :: side_dimid, xysize_dimid
      integer :: side_dimsize, xysize_dimsize
      integer :: dimvarid, rangevarid, depthvarid, spacingvarid

      integer,allocatable, dimension(:) :: actual_dimensions, itmp
      real, allocatable, dimension(:) :: actual_xrange, actual_yrange, actual_spacing
      character(len=220) :: gebco_grid_file
      integer :: indx, nx2,err
      logical :: ex



      ! Retrieve etopo5 path
      call getenv('GEBCO_PATH',path0)
      if (trim(path0)=='') then
         print *,'No path set for GEBCO. Set the environment variable GEBCO_PATH'
         print *,'To point to the location of the GEBCO bathymetry'
         call exit(1)
      end if
      path0=trim(path0)//'/'

      gebco_grid_file=trim(path0)//trim(gebco_grid_file1)

      inquire(file=trim(gebco_grid_file),exist=ex)
      if (.not.ex) then
         print *,'Can not find GEBCO file '//trim(gebco_grid_file)
         print *,'path0=',trim(path0)
         stop '(mod_GEBCO)'
      end if


! --- Get Gebco grid dimensions and variables, allocate gebco variables
      ! Data producers decided to put bathymetry into one long matrix
      ! for some obscure reason....


      ! Get grid dims  here
      nc_stat=0
      err=NF90_OPEN(gebco_grid_file,NF90_NOCLOBBER,ncid)
      if (err/=0) then
         print *,'Could not open '//trim(gebco_grid_file)
         stop '(read_gebco)'
      end if
      call nf90_handle_err('file',NF90_INQUIRE(ncid,ndims,nvars))
      print '(a,i2)','Reading GEBCO File '//gebco_grid_file

      ! Inquire on dimensions (id, size and name)
      call nf90_handle_err('dims',NF90_INQ_DIMID(ncid,'side',side_dimid))
      call nf90_handle_err('dims',NF90_INQ_DIMID(ncid,'xysize',xysize_dimid) )
      call nf90_handle_err('dims',NF90_INQUIRE_DIMENSION(ncid,  side_dimid,len=  side_dimsize))
      call nf90_handle_err('dims',NF90_INQUIRE_DIMENSION(ncid,xysize_dimid,len=xysize_dimsize))
      print '(a,2i10)','Read GEBCO netcdf Dimensions : ',side_dimsize,xysize_dimsize

      ! Real (2D) dimensions are actually in a variable (?!)
      write(6,'(a22)',advance='no') 'Reading dimension var'
      allocate(actual_dimensions(side_dimsize))
      call nf90_handle_err('var-dim',nf90_inq_varid(ncid,'dimension',dimvarid       ))
      call nf90_handle_err('var-dim',NF90_GET_VAR  (ncid,dimvarid  ,actual_dimensions))
      write(6,'(a,2i8)',advance='yes') '..done. Actual dims are ',actual_dimensions

      ! Real (2D) ranges are actually in a variable (?!)
      write(6,'(a22)',advance='no') 'Reading xrange var'
      allocate(actual_xrange(side_dimsize))
      call nf90_handle_err('var-xrange',nf90_inq_varid(ncid,'x_range',rangevarid))
      call nf90_handle_err('var-xrange',NF90_GET_VAR  (ncid,rangevarid  ,actual_xrange))
      write(6,'(a,2f10.2)',advance='yes') '..done. Xrange is ',actual_xrange

      ! Real (2D) ranges are actually in a variable (?!)
      write(6,'(a22)',advance='no') 'Reading yrange var'
      allocate(actual_yrange(side_dimsize))
      call nf90_handle_err('var-yrange',nf90_inq_varid(ncid,'y_range',rangevarid))
      call nf90_handle_err('var-yrange',NF90_GET_VAR  (ncid,rangevarid  ,actual_yrange))
      write(6,'(a,2f10.2)',advance='yes') '..done. Yrange is ',actual_yrange

      ! Real (2D) spacings are actually in a variable (?!)
      write(6,'(a22)',advance='no') 'Reading spacing var'
      allocate(actual_spacing(side_dimsize))
      call nf90_handle_err('var-spacing',nf90_inq_varid(ncid,'spacing',spacingvarid))
      call nf90_handle_err('var-spacing',NF90_GET_VAR  (ncid,spacingvarid  ,actual_spacing))
      write(6,'(a,2f10.4)',advance='yes') '..done. Spacing is ',actual_spacing
      print *,gebco_nx,gebco_ny,actual_dimensions
      !stop

#if defined (GEBCO_2_MIN)


      print  * ,'Using 2 minute gebco subsampling'


      ! Low-memory version.... The data is reduced to 2-minute resolution

      ! NB: the gebco data wraps around the globe (points on 180 longitude
      ! are covered twice)
      nx2=actual_dimensions(1)
      gebco_nx=actual_dimensions(1)/2
      gebco_ny=actual_dimensions(2)/2
      gebco_dlon=actual_spacing(1)*2.
      gebco_dlat=actual_spacing(2)*2.
      gebco_flon=actual_xrange(1)+gebco_dlon/2.
      gebco_flat=actual_yrange(1)+gebco_dlat/2.
      readchunk=2*nx2
      print *,gebco_flon,gebco_flat,gebco_dlon,gebco_dlat


      allocate(gebco_depth  (gebco_nx,gebco_ny))
      allocate(gebco_lon    (gebco_nx,gebco_ny))
      allocate(gebco_lat    (gebco_nx,gebco_ny))
      allocate(itmp(readchunk))
      write(6,'(a22)',advance='no') 'Reading depth var'
      print *
      call nf90_handle_err('var-depth',nf90_inq_varid(ncid,'z',depthvarid))
      !open(99,file='gebco.out',form='formatted',status='replace')

      ! Real (2D) bathymetry is in a vector (?!)
      do j=1,gebco_ny
         !call nf90_handle_err('var-depth',
         err=NF90_GET_VAR  (ncid,depthvarid  ,itmp, start=(/2*(j-1)*nx2+1/),count=(/readchunk/))
         if (mod(j,100)==0) print *,j,gebco_ny,gebco_nx
         !print *,itmp(1)
         !call nf90_handle_err('var-depth',NF90_GET_VAR  (ncid,depthvarid ,itmp))
         !write(6,'(a,2i7)',advance='yes') '..done  Max/min depth is ', maxval(itmp),minval(itmp)

         ! Copy to regular lon-lat grid (-180:180, -90:90))
         do i=1,gebco_nx
               gebco_lon  (i,j) = (i-1)*gebco_dlon + gebco_flon
               gebco_lat  (i,j) = (j-1)*gebco_dlat + gebco_flat

               !print *,i,j,  gebco_lon  (i,j) ,   gebco_lat  (i,j) 

               ! Hmmm  
               !indx = (j-1)*(gebco_nx+1) + i 
               !indx = (gebco_ny-j)*(nx2) + i 
               j2=gebco_ny-j+1
               gebco_depth(i,j2) =                     float(itmp(    2*(i-1)+1))
               gebco_depth(i,j2) = gebco_depth(i,j2) + float(itmp(    2*(i-1)+2))
               gebco_depth(i,j2) = gebco_depth(i,j2) + float(itmp(nx2+2*(i-1)+1))
               gebco_depth(i,j2) = gebco_depth(i,j2) + float(itmp(nx2+2*(i-1)+2))
               gebco_depth(i,j2) = gebco_depth(i,j2) / 4.

               ! gebco uses negative values for depth -- 
               ! set land regions to depth zero
               gebco_depth(i,j2) = max(0.,-gebco_depth(i,j2))

               !write (99,*) gebco_lon(i,j),gebco_lat(i,j),gebco_depth(i,j)
         end do
      end do
      !close(99)
      deallocate(itmp)
#else
      print  * ,'Using Full gebco dataset'

      ! Real (2D) bathymetry is in a vector (?!)
      write(6,'(a22)',advance='no') 'Reading depth var'
      allocate(itmp(xysize_dimsize))
      call nf90_handle_err('var-depth',nf90_inq_varid(ncid,'z',depthvarid))
      call nf90_handle_err('var-depth',NF90_GET_VAR  (ncid,depthvarid  ,itmp))
      write(6,'(a,2i7)',advance='yes') '..done  Max/min depth is ', maxval(itmp),minval(itmp)
      ! NB: the gebco data wraps around the globe (points on 180 longitude
      ! are covered twice)
      nx2=actual_dimensions(1)
      !gebco_nx=actual_dimensions(1)-1
      gebco_nx=actual_dimensions(1)
      gebco_ny=actual_dimensions(2) 
      gebco_dlon=actual_spacing(1)
      gebco_dlat=actual_spacing(2)
      gebco_flon=actual_xrange(1)
      gebco_flat=actual_yrange(1)
      !print *,gebco_flon,gebco_flat

      allocate(gebco_depth  (gebco_nx,gebco_ny))
      allocate(gebco_lon    (gebco_nx,gebco_ny))
      allocate(gebco_lat    (gebco_nx,gebco_ny))

      ! Copy to regular lon-lat grid (-180:180, -90:90))
      do j=1,gebco_ny
      do i=1,gebco_nx
            gebco_lon  (i,j) = (i-1)*gebco_dlon + gebco_flon
            !gebco_lat  (i,j) = (j-1)*gebco_dlat + gebco_flat
            gebco_lat  (i,j) = (j-1)*gebco_dlat + gebco_flat

            ! Hmmm  
            !indx = (j-1)*(gebco_nx+1) + i 
            indx = (gebco_ny-j)*(nx2) + i 
            gebco_depth(i,j) = float(itmp( indx ))

            ! gebco uses negative values for depth -- 
            ! set land regions to depth zero
            gebco_depth(i,j) = max(0.,-gebco_depth(i,j))
      end do
      end do
      deallocate(itmp)
#endif

      where (gebco_lon >= 180.0) gebco_lon=gebco_lon-360.0
      !shift=etopo5_nx/2
      !etopo5%deep=cshift(etopo5%deep,-shift,1)
      !etopo5%lon=cshift(etopo5%lon,-shift,1)
      !etopo5%lat=cshift(etopo5%lat,-shift,1)
      gebco_minlat=minval(gebco_lat)
      gebco_maxlat=maxval(gebco_lat)
      gebco_minlon=minval(gebco_lon)
      gebco_maxlon=maxval(gebco_lon)
      print *,'GEBCO:  bounds= ',gebco_minlat,gebco_maxlat,gebco_minlon,gebco_maxlon
   end subroutine read_gebco


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





!KAL -- Changes:
!KAL -- 20070425 - fixed using integer as logical in ifs (bndflag)

subroutine gebco_fix65N(glon,glat,gdepths,gnx,gny,edepths,eflon,eflat,five_min,enx,eny)
use mod_sphere_tools
implicit none
integer, intent(in) :: gnx,gny,enx,eny
real, intent(in) :: eflon,eflat,five_min
real, dimension(gnx,gny), intent(in)    :: glon,glat
real, dimension(gnx,gny), intent(inout) :: gdepths
real, dimension(enx,eny), intent(in   ) :: edepths

real, dimension(4) :: crnlon, crnlat
integer :: i,j,k,l,i1,j1,i2,j2,im1,ip1
integer :: ninbox
integer, allocatable, dimension(:) :: boxi,boxj,bndflag
real   , allocatable, dimension(:) :: bnddist,boxdepths,boxweights
real :: dist, mindist

real :: aa,bb,a1,a2,a3,a4,x
integer :: ei,ej,eip1,ejp1

! Define problematic region
crnlat=(/66.0, 66.0,63.5 , 63.5/)
crnlon=(/-5.5,  1.0,-1.5 , -8.0/)

! Find points inside that box
ninbox=0
do j=2,gny-1
do i=1,gnx
   if(inbox(crnlon,crnlat,4,glon(i,j),glat(i,j))) then
      !print *,i,j
      !gdepths(i,j)=100.
      ninbox=ninbox+1
   end if
end do
end do
print *,'GEBCO Points in 65N box:',ninbox

allocate(boxi(ninbox))
allocate(boxj(ninbox))
allocate(bndflag(ninbox))
allocate(bnddist(ninbox))
allocate(boxdepths(ninbox))
allocate(boxweights(ninbox))

! Find box points
ninbox=0
do j=2,gny-1
do i=1,gnx
   if(inbox(crnlon,crnlat,4,glon(i,j),glat(i,j))) then
      !print *,i,j
      !gdepths(i,j)=100.
      ninbox=ninbox+1
      boxi(ninbox)=i
      boxj(ninbox)=j

      im1=mod(i-2+gnx,gnx)+1
      ip1=mod(i      ,gnx)+1

      ! Boundary point?
      if     (.not.inbox(crnlon,crnlat,4,glon(im1,j),glat(im1,j))) then
         bndflag(ninbox)=1
      else if(.not.inbox(crnlon,crnlat,4,glon(ip1,j),glat(ip1,j))) then
         bndflag(ninbox)=1
      else if(.not.inbox(crnlon,crnlat,4,glon(i,j-1),glat(i,j-1))) then
         bndflag(ninbox)=1
      else if(.not.inbox(crnlon,crnlat,4,glon(i,j+1),glat(i,j+1))) then
         bndflag(ninbox)=1
      else
         bndflag(ninbox)=0
      end if
   end if
end do
end do
print *,'Boundary points of 65N box:',sum(bndflag)

! Find minimum distance from every box point to boundary
do k=1,ninbox
   if (bndflag(k)==1) then
      bnddist(k)=0.
   else
      i1=boxi(k)
      j1=boxj(k)
      mindist=1e10
      do l=1,ninbox
         i2=boxi(l)
         j2=boxj(l)
         if (bndflag(l)==1) then
            dist=spherdist(glon(i1,j1),glat(i1,j1), &
                           glon(i2,j2),glat(i2,j2))
            mindist=min(dist,mindist)
         end if
      end do
      bnddist(k)=mindist
   end if
end do

! weights
do k=1,ninbox
   i=boxi(k)
   j=boxj(k)

   ! weights go as x squared

   boxweights(k)=bnddist(k)/maxval(bnddist)
   boxweights(k)=max(0.,boxweights(k)*5.)
   boxweights(k)=min(1.,boxweights(k))
   boxweights(k)= boxweights(k)**2
end do

! Find etopo5 pivot points for each point in box.
! Interpolate to temp variable "boxdepths"
do k=1,ninbox
   i=boxi(k)
   j=boxj(k)


   ! NB - mod fix here
   x=(glon(i,j)-eflon)/five_min + 1
   ei=floor(x)
   aa=x-ei

   x=(glat(i,j)-eflat)/five_min + 1
   ej=floor(x)
   bb=x-ej

   eip1=mod(ei      ,enx)+1
   ejp1=min(enx,ej+1)

   a1=(1.0-aa)*(1.0-bb)
   a2=aa*(1.0-bb)
   a3=aa*bb
   a4=(1.0-aa)*bb

   
   ! Interpolated etopo 5 weights in box:
   boxdepths(k) = a1*edepths(ei  ,ej  )+a2*edepths(eip1,ej  )&
                 +a3*edepths(eip1,ejp1)+a4*edepths(ei  ,ejp1)

   ! Weighted average of gebco and etopo5:
   x=boxweights(k)
   gdepths(i,j)=x*boxdepths(k) + (1.-x)*gdepths(i,j)
   !gdepths(i,j)=boxweights(k)

   ! For help
   !if (bndflag(k)) gdepths(i,j)=0.

end do




end subroutine gebco_fix65N


   end module mod_gebco
