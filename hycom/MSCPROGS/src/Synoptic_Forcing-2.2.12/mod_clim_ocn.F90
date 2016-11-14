module mod_clim_ocn
! -- sss_sst_nodc   - Read nodc sea surface climatologies. Does not really
!                     belong here either.
!
! NB - All these routines require that mod_common and mod_forcing is initialized 

private :: nf90_handle_error
character(len=200),save  :: path_woa2005    ='./WOA2005/'
contains


! #######################################################################
! Reading - manipulating sss_nodc and sst_nodc files 
! #######################################################################
! KAL -- This one is deprecated - Strange features in the Arctic..
      subroutine sss_sst_nodc(sstfld,sssfld,month)
      use mod_xc
      use mod_za
      use mod_grid
      implicit none
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: sstfld,sssfld
      integer, intent(in) :: month

      logical :: sssex,sssexa,sstex,sstexa,sstexb,sssexb
      integer :: i,j,k
      real    :: hmin, hmax, hmin2, hmax2
      real, dimension(idm,jdm) :: gfld
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: lfld



      ! Look for existing style files
      inquire(file='./Data/sst_nodc.dat',exist=sstex)
      inquire(file='./Data/sss_nodc.dat',exist=sssex)

      if (.not.sstex) then
         if (mnproc==1) then
            write(lp,*) 'Can not find file sst_nodc.dat'
            call flush(lp)
         end if
         print *, 'sss_sst_nodc)'
         call exit(1)
      end if
      open(14,file='./Data/sst_nodc.dat',form='formatted')
      do k=1,month
         read(14,'(10f9.4)')sstfld
      end do

      if (.not.sssex) then
         if (mnproc==1) then
            write(lp,*) 'Can not find file sss_nodc.{a,b,dat}'
            call flush(lp)
         end if
         print *, '(sss_sst_nodc)'
         call exit(1)
      end if
      open(14,file='./Data/sss_nodc.dat',form='formatted')
      do k=1,month
         read(14,'(10f9.4)')sssfld
      end do

      ! Strange arctic feature....
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
      do j=1,jdm
      do i=1,idm
         if (plat(i,j) > 71.0) sssfld(i,j)=min(35.0,sssfld(i,j))
      enddo
      enddo
!$OMP END PARALLEL DO

      end subroutine sss_sst_nodc
! ------- ------------ ----------- -------------- ------------ ---
! ------  Change notes:
! ------- ------------ ----------- -------------- ------------ ---
! --- Cleaned up main routine - created many subroutines            KAL -- 06022005
! --- Got rid of global arrays for climatology (reduces mem usage)  KAL -- 06022005
! --- Introduced preprocessing of clim. - placed in ".a" files      KAL -- 06022005
! --- Clim files now placed in .ab-files -- more robust             KAL -- 22022005
! ------- ------------ ----------- -------------- ------------ ---


      subroutine sss_sst_WOA2005(sstfld,sssfld,month)
      use netcdf
      use mod_xc
      use mod_za
      use mod_grid
      use mod_forcing_nersc
      use m_bilin_ecmwf2
      use m_nearestpoint
      implicit none
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: sstfld,sssfld
      integer, intent(in) :: month
      character(len=*), parameter :: fsal='s0112an1.nc'
      character(len=*), parameter :: ftem='t0112an1.nc'
      character(len=*), parameter :: vsal='s0112an1'
      character(len=*), parameter :: vtem='t0112an1'

      character(len=200) :: cenv
      integer :: ncid, varid
      integer i,j,ipiv, jpiv
      integer, parameter :: nlon=360
      integer, parameter :: nlat=180
      real, dimension(nlon,nlat) :: woasal, &
         woatem,woasal2,woatem2,woalon2,woalat2
      real, dimension(nlon) :: woalon
      real, dimension(nlat) :: woalat
      real :: flon, flat, dlon, dlat, a1, a2, a3, a4, hmin, hmax
      logical :: ass
      real,parameter    :: fill_sal=35.
      real,parameter    :: fill_tem=10.

      logical, save :: first=.true.
      integer, save :: ipivs (nlon,nlat)
      integer, save :: jpivs (nlon,nlat)
      integer, save :: npmask(nlon,nlat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! Default path to WOA2005 ocean climatology
      call getenv('WOA2005_PATH',cenv)
      if (trim(cenv)/='') then
         path_woa2005=trim(cenv)//'/'
      else
         print *
         print '(a)','Need path to WOA2005 data set - set path in '
         print '(a)','environment variable WOA2005_PATH'
         call exit(1)
      end if

      ! We know how the data is organized and we only need the surface 
      if (first) write(lp,'(a)') '---Look for WOA2005 data in '//trim(path_woa2005)

      ! Open NetCDF file
      call nf90_handle_err( NF90_OPEN(trim(path_woa2005)//fsal,NF90_NOCLOBBER,ncid),&
           'Can not find '//trim(path_woa2005)//fsal)

      call nf90_handle_err( nf90_inq_varid(ncid, 'lon', varid),   &
           fsal//' does not contain variable '//'lon')
      ! Get surface values
      call nf90_handle_err( NF90_GET_VAR(ncid,varid,woalon),&
              'Error in retrieving '//'lon'//' values') 

      call nf90_handle_err( nf90_inq_varid(ncid, 'lat', varid),   &
           fsal//' does not contain variable lat')
      ! Get surface values
      call nf90_handle_err( NF90_GET_VAR(ncid,varid,woalat),&
              'Error in retrieving '//'lat'//' values')

      call nf90_handle_err( nf90_inq_varid(ncid, vsal, varid),   &
           fsal//' does not contain variable '//vsal)
      ! Get surface values
      call nf90_handle_err( NF90_GET_VAR(ncid,varid,woasal,start=(/1,1,1,month/)),&
              'Error in retrieving '//vsal//' values')

      ! Open NetCDF file
      call nf90_handle_err( NF90_CLOSE(ncid),&
           'Can not close '//fsal)

      ! Open NetCDF file
      call nf90_handle_err( NF90_OPEN(trim(path_woa2005)//ftem,NF90_NOCLOBBER,ncid),&
           'Can not find '//trim(path_woa2005)//ftem)

      call nf90_handle_err( nf90_inq_varid(ncid, vtem, varid),   &
           ftem//' does not contain variable '//vtem)
      ! Get surface values
      call nf90_handle_err( NF90_GET_VAR(ncid,varid,woatem,start=(/1,1,1,month/)),&
              'Error in retrieving '//vtem//' values')

! --  Extrapolate dataset to land points
      !print *,count(woasal<1e10)
      woasal2=woasal
      woatem2=woatem
      do i=1,nlon
         woalat2(i,:)=woalat
      end do
      do i=1,nlat
         woalat2(:,i)=woalon
      end do
      if (first) then
         npmask=0
         print '(a)','---WOA2005 Extrapolation on first pass '
      end if
      do j=1,nlat
      do i=1,nlon
         if (woasal(i,j)>1e10) then
            if (first) then
               ipiv=i
               jpiv=j
               call  nearestpoint(woalon2,woalat2,nlon,nlat, &
                                  woalon2(i,j),woalat2(i,j),ipiv,jpiv,  &
                                  a1,a2,a3,a4,woasal<1e10,ass,.true.)

               ! extrapolation succeeded
               if(ass) then
                  woasal2(i,j)=woasal(ipiv,jpiv)
                  woatem2(i,j)=woatem(ipiv,jpiv)
                  ipivs(i,j)  = ipiv
                  jpivs(i,j)  = jpiv
                  npmask (i,j)=1
               ! out of  extrapolation range
               else
                  woasal2(i,j)=fill_sal
                  woatem2(i,j)=fill_tem
                  npmask (i,j)=2
               end if
            else
               ! earlier extrapolation succeeded
               if (npmask(i,j)==1) then
                  ipiv=ipivs(i,j)
                  jpiv=jpivs(i,j) 
                  woasal2(i,j)=woasal(ipiv,jpiv)
                  woatem2(i,j)=woatem(ipiv,jpiv)
               ! earlier extrapolation failed (out of range)
               elseif (npmask(i,j)==2) then
                  woasal2(i,j)=fill_sal
                  woatem2(i,j)=fill_tem
               ! This indicates a mask mismatch between months or tem/sal
               else
                  print *,'mask mismatch'
                  stop '(sss_sst_woa2005)'
               end if
            end if
         end if
      end do
      end do
      if (first) then
         print '(a)','--- ... WOA2005 Extrapolation done '
      end if
      woasal=woasal2
      woatem=woatem2
      !print *,count(woasal/=fill_sal)

! -- Interpolate dataset -- Dataset checked to be uniform
      flon=woalon(1)
      flat=woalat(1)
      dlon=woalon(2)-woalon(1)
      dlat=woalat(2)-woalat(1)

      call bilin_ecmwf2(woasal,nlon,nlat,flon,flat,dlon,dlat, &
                       sssfld,plon,plat,depths)   
      call bilin_ecmwf2(woatem,nlon,nlat,flon,flat,dlon,dlat, &
                       sstfld,plon,plat,depths)   
      first=.false.
      !call zaiopf('./tst.a','replace',24)
      !call zaiowr(sstfld,ip,.false.,hmin,hmax,24,.true.)
      !call zaiocl(24)
      !stop
      end subroutine


   subroutine nf90_handle_err(errcode,info)
      use mod_xc , only: xcstop, mnproc
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      character(len=*), intent(in) :: info

      if (errcode/=NF90_NOERR) then
         if (mnproc==1) then
            write(6,'(a)') NF90_STRERROR(errcode)
            write(6,'(a)') info
         end if
         stop '(ncvar_read)'
         call xcstop('(ncvar_read)')
      end if
   end subroutine

end module mod_clim_ocn
