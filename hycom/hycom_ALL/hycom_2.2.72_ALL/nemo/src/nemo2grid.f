      program nemo2grid
      use mod_za  ! HYCOM array I/O interface
      use netcdf  ! NetCDF fortran 90 interface
c
      implicit none
c
c --- extract nemo grid and bathymetry in HYCOM format.
c
c --- grid is just plon,plat (use grid_lonlat_2d to form regional.grid)
c
c --- Alan J. Wallcraft,  Naval Research Laboratory, August 2005.
c --- Mostafa Bakhoday-Paskyabi, NERSC, Sep. 2017.
c
      character*79  preambl(5)
      character*240 cfiler
      integer       lon,lat,i,j
      integer       ncFIDr,ncDIDx,ncVIDx
      real*4        plat_min,plat_max,plon_min,plon_max,
     &              xmax,xmin
c
      double precision, allocatable :: 
     &        bathy(:,:), !Final bathymetry at RHO-points
     &  lon_rho(:,:), !longitude of RHO-points
     &  lat_rho(:,:), !latitude of RHO-points
     & mask_rho(:,:)  !mask on RHO-points
c
      real*4, allocatable ::
     &   depth(:,:), !hycom bathymetry
     &    plon(:,:), !hycom longitude
     &    plat(:,:)  !hycom latitude
      integer, allocatable ::
     &     msk(:,:)  !hycom mask (not used)
c
      real*4, parameter :: spval = 2.0**100
c
c     array sizes
c
      CALL GETENV('CDF_NEMO',cfiler)
      write(6,*)
      write(6,*) 'CDF_NEMO = ',trim(cfiler)
      call zhflsh(6)
      ! open NetCDF file
      call ncheck(nf90_open(trim(cfiler), nf90_nowrite, ncFIDr))
      ! get longitude
      write(6,*) 'nf90_inq_dimid - ', 'lon'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFIDr, 'lon',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFIDr,ncDIDx,len=lon))
      ! get latitude
      write(6,*) 'nf90_inq_dimid - ','lat'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFIDr,'lat',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFIDr,ncDIDx,len=lat))
      write(6,*) 'longitude,latitude = ',lon,lat
      write(6,*) 
      call zhflsh(6)
      idm =  lon
      jdm = lat
c
      allocate(        bathy(idm,jdm),
     &           lon_rho(idm,jdm),
     &           lat_rho(idm,jdm),
     &          mask_rho(idm,jdm) )
c
      allocate(    depth(idm,jdm),
     &              plon(idm,jdm),
     &              plat(idm,jdm),
     &               msk(idm,jdm) )
c
c --- read nemo variables
c
      write(6,*) 'nf90_inq_varid - ','Bathymetry'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr,       'Bathymetry',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                          bathy(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','lon'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr, 'lon',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                    lon_rho(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','lat'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr, 'lat',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                    lat_rho(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','mask'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr,'mask',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                   mask_rho(:,:)    ))
      ! close NetCDF file
      call ncheck(nf90_close(ncFIDr))
c
c --- form and write HYCOM fields
c
      do j= 1,jdm
        do i= 1,idm
          plon(i,j) = lon_rho(i,j)
          plat(i,j) = lat_rho(i,j)
          if     (mask_rho(i,j).eq.1) then
            depth(i,j) = bathy(i,j)
          else
            depth(i,j) = spval
          endif
        enddo !i
      enddo !j
c
      call zaiost
c
      call zhopen(      11,'formatted','new',0)
      call zaiopn('new',11)
      call zaiowr( plon,msk,.false., xmin,xmax, 11, .false.)
      write(11,6100) 'plon',xmin,xmax
      write( 6,6100) 'plon',xmin,xmax
      call zhflsh(6)
      plon_min = xmin
      plon_max = xmax
      call zaiowr( plat,msk,.false., xmin,xmax, 11, .false.)
      write(11,6100) 'plat',xmin,xmax
      write( 6,6100) 'plat',xmin,xmax
      call zhflsh(6)
      call zaiocl(11)
      close( unit=11)
      plat_min = xmin
      plat_max = xmax
c
      call zaiopn('new',12)
      call zaiowr(depth,msk,.false., xmin,xmax, 12, .false.)
c
      call zhopen(      12,'formatted','new',0)
      preambl(1) = ' '
      preambl(2) = ' '
c     +  'bathymetery from nemo'

      write(12,'(i5,a)')
     +       idm,   "    'idm   ' = longitudinal array size"
      write(12,'(i5,a)')
     +       jdm,   "    'jdm   ' = latitudinal array size"
      write(12,'(a,2i5)')
     +        'i/jdm =',
     +       idm,jdm
      write(12,'(a,4f12.5)')
     +        'plon,plat range =',
     +       plon_min,plon_max,plat_min,plat_max
c      preambl(4) = ' '
c      preambl(5) = ' '
      write(12,'(A79)') ' '
      write(12,6200)    xmin,xmax
      write(6, '(A79)') ' '
      write( 6,6200)    xmin,xmax
      call zhflsh(6)
c
      call zaiocl(12)
      close( unit=12)
c
 6100 format(a,':  min,max = ', 2f15.5)
 6200 format('min,max depth = ',2f10.3)
      end

      subroutine ncheck(status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer, intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error from NetCDF library'
        write(6,'(a/)')   trim(nf90_strerror(status))
        call zhflsh(6)
        stop
      end if
      end subroutine ncheck
