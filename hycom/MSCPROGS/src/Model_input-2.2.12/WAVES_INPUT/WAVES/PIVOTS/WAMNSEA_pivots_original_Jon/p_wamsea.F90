program wamsea
   use mod_xc
   use mod_za
   use mod_grid, only : get_grid, plon, plat, qlat, qlon, scpx
   use mod_confmap
   use m_pivotp_micom

!   use mod_hycom_nersc

   use m_ncvar_read
   use netcdf

! --- TODO: add input argument for the wam data and grid data, read glon,glat,nlon,nlat

   implicit none

      integer*4, external :: iargc
      character(len=50) filenamea,filenameb
! , ncname, modgrid
      
! --- Hardcoded grid sizes in WAM North Sea nc files 
      integer, parameter :: nlon=248, nlat=400
! --- WAMNSEA dimension and lon,lat
      real :: glat(nlon,nlat), glon(nlon,nlat), gmask(nlon,nlat)
      integer :: ip(nlon,nlat)
! --- MODEL dimesion, lon,lat and related pivot points and weights
      integer,  allocatable,dimension(:,:) :: ipiv_o, jpiv_o,dist_mod

! --- local
      real :: lon, lat, lon_n, lat_n, a1, a2, a3, a4
      real*8 :: dist
      integer :: i, j, ipiv, jpiv

      logical :: ass
      real, external :: spherdist

      real cx,cy,cz,theta_c,phi_c,hmin,hmax
      complex c,w

!      if (iargc()/=1) then
!         print *,'Routine calculates ipiv,jpiv points in the WAMNSEA'
!         print *,'grid related to longitude/latitude in Model grid. '
!         print *,'Run program in /work/shared/nersc/msc/WAMNSEA/ .  '
!         print *,'!!!! Remember to copy grid files from region/topo/ .'
!         print *,'  Usage : wamnsea-2.2.12 '
!         call exit(1)
!      end if
      
!      filenamea='FR1a0.03_pivot_wamnsea.a'
!      filenameb='FR1a0.03_pivot_wamnsea.b'
      filenamea='BS1a0.045_pivot_wamnsea.a'
      filenameb='BS1a0.045_pivot_wamnsea.b'

      print *,'print to abfile: ',trim(filenamea)
      print *,'print to abfile: ',trim(filenameb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Initialize IO for .ab files
      CALL XCSPMD()  
      allocate(ipiv_o(idm,jdm),jpiv_o(idm,jdm),dist_mod(idm,jdm))
      CALL ZAIOST()

! --- Get info related the Model grid
      call get_grid()

! --- Get WAMNSEA lon/lat data
      print*, 'TW2'
      call ncvar_read('wam_nsea.an.20120228.nc','longitude',glon,nlon,nlat, 1,1,1)
      print*, 'TW3'
      call ncvar_read('wam_nsea.an.20120228.nc','latitude',glat,nlon,nlat, 1,1,1)
      print*, 'TW4'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




   ipiv=1
   jpiv=1
!   do i=1,idm
!      do j=1,jdm
   do i=1,idm
      !print *,'i    :',i
      do j=1,jdm
    !  do j=100,243
         call pivotp_micom(plon(i,j),plat(i,j),glon,glat, ipiv,jpiv, nlon,nlat,dist )
      !print *,'ipiv,jpiv    :',ipiv,jpiv
          
         ipiv_o(i,j)=ipiv
         jpiv_o(i,j)=jpiv
         dist_mod(i,j)=spherdist(plon(i,j),plat(i,j), glon(ipiv,jpiv), glat(ipiv,jpiv))
         if ((ipiv==1 .or. jpiv==1 .or. ipiv==nlon .or. jpiv==nlat ) ) then
         !if ((ipiv==1 .or. jpiv==1 .or. ipiv==idm .or. jpiv==jdm ) .and. dist > scpx(ipiv,jpiv)) then
            ipiv_o(i,j)=0.
            jpiv_o(i,j)=0.
         endif
      end do
   end do
!   print *,'test', maxval(ipiv_o),maxval(jpiv_o)
 call zaiopf(trim(filenamea), 'replace', 909)
 call zaiowr(real(ipiv_o(:,:)),ip,.false.,hmin,hmax,909,.true.) 
   open (unit=909, file=filenameb,status='replace', action='write')
    write(909,'(a)') 'Pivot positions in WAMNSEA related to model grid   '
    write(909,'(a)') ''
    write(909,'(a)') ''
    write(909,'(a)') ''
    write(909,'(a,2i5)') 'i/jdm = ',idm,jdm
    write(909,'(" ipiv:range = ",2e16.8)') hmin,hmax
 call zaiowr(real(jpiv_o(:,:)),ip,.false.,hmin,hmax,909,.true.)
    write(909,'(" jpiv:range = ",2e16.8)') hmin,hmax
 call zaiowr(real(dist_mod(:,:)),ip,.false.,hmin,hmax,909,.true.)
    write(909,'(" dist:range = ",2e16.8)') hmin,hmax
    close(909) 
 call zaiocl(909)

end program wamsea
