program getpivots
   !! pre-compiled modules;
   use mod_xc
   use mod_za
   use mod_grid, only : get_grid, plon, plat, qlat, qlon, scpx
   use mod_confmap


   !! in libhycnersc.a library;
   use netcdf

   !! extra
   use m_pivotp_micom
   use m_ncvar_read

! --- TODO: add input argument for the wam data and grid data, read glon,glat,nlon,nlat

   implicit none

      integer*4, external :: iargc
      character(len=80)   :: ncfil!! full name of input data file
                                  !! (determined using 'd_path', 'd_name' variables below)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Only need to change the parameters between the 2 lines of '!!!'
      !! cf README.txt for instructions on how to run the executable;

      !! Description of output in 'filenameb'
      character(len=*), parameter :: out_str    = 'Pivot positions in AMSR-E related to model grid   '

      !! Names of output files;
      character(len=*), parameter :: filenamea  = 'pivots_amsre.a'
      character(len=*), parameter :: filenameb  = 'pivots_amsre.b'

      !! Location and name of a file with lon/lat of input data points
      character(len=*), parameter :: d_path  = '/work/shared/nersc/msc/AMSRE6.25km/netcdf/'!!path of input data file
      character(len=*), parameter :: d_name  = 'LongitudeLatitudeGrid-n6250-Arctic.nc'!!local name of any input data file containing lon/lat

      !! Hardcoded grid sizes in data nc files 
      integer, parameter :: nlon=1216, nlat=1792

      !! Names of lon/lat variables in data nc files 
      character(len=*), parameter :: lon_name = 'Longitudes'
      character(len=*), parameter :: lat_name = 'Latitudes'

      !!some parameters for labels of output file names
!      character(len=*), parameter :: data_name  = 'ww3'
!      character(len=*), parameter :: model_name = 'BS1a0.045'
!      character(len=*), parameter :: model_name  = 'FR1a0.03'
!      character(len=*), parameter :: model_name  = 'TP4a0.12'

! , ncname, modgrid
      
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- input data dimension and lon,lat
      real :: glat(nlon,nlat), glon(nlon,nlat)!matrices of lon/lat;
      !real :: gmask(nlon,nlat)
      integer :: ip(nlon,nlat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!TW: add a brute force check on the boundary of the data set
      integer, dimension(2*nlon+2*(nlat-2))     :: i_bdy,j_bdy!!i,j values of data set's boundary
      real,    dimension(2*nlon+2*(nlat-2))     :: res_bdy    !!resolution of data set at its boundary
      integer, dimension(2*nlon+2*(nlat-2),2)   :: r_bdy2     !!index of neighbour on the boundary (to get 'res_bdy')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- MODEL dimension, lon,lat and related pivot points and weights
      integer,  allocatable,dimension(:,:) :: ipiv_o, jpiv_o,dist_mod

! --- local
      real :: lon, lat, lon_n, lat_n, a1, a2, a3, a4     &
     &         ,d1,d2,d_bdy,dmin_bdy,max_bdy,dtol
      real*8 :: dist
      integer :: i,j,ipiv,jpiv                           &
     &            ,r,r1,r2,r3,r4,r_bdy                   &
     &            ,r_1,r_2,rmin_bdy
     integer   :: dumloc(1)

      logical :: ass                                     &
     &           ,piv_on_bdy,chk_r
      real, external :: spherdist

      real cx,cy,cz,theta_c,phi_c,hmin,hmax
      complex c,w

      print *,'print to abfile: ',trim(filenamea)
      print *,'print to abfile: ',trim(filenameb)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Get lon/lat positions of input data from netcdf file:
      ncfil = d_path//d_name
      print *,'using nc file: ',trim(ncfil)
      call ncvar_read(trim(ncfil),lon_name,glon,nlon,nlat, 1,1,1)
      call ncvar_read(trim(ncfil),lat_name,glat,nlon,nlat, 1,1,1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Initialize IO for .ab files
      CALL XCSPMD()  
      allocate(ipiv_o(idm,jdm),jpiv_o(idm,jdm),dist_mod(idm,jdm))
      CALL ZAIOST()

! --- Get info related the Model grid
      call get_grid()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*, ' '
      print*, 'Getting resolution of data around source boundary'
      print*, '...'
      !!get indices of boundary points
      !!(anti-clockwise round boundary)
                               !r = 1 : (i,j) = (1,1)
      r1 = nlon                !r = r1: (i,j) = (nlon,1)
      r2 = nlon+(nlat-2)       !r = r2: (i,j) = (nlon,nlat-1) ok
      r3 = 2*nlon+(nlat-2)     !r = r3: (i,j) = (1,nlat)
      r4 = 2*nlon+2*(nlat-2)   !r = r4: (i,j) = (1,2)

      do i=1,nlon
         i_bdy(i)    = i
         j_bdy(i)    = 1
         i_bdy(r2+i) = nlon+1-i
         j_bdy(r2+i) = nlat
      end do
      do j=1,nlat-2
         i_bdy(r1+j) = nlon
         j_bdy(r1+j) = j+1
         i_bdy(r3+j) = 1
         j_bdy(r3+j) = nlat-j
      end do
      
      !!get indices of neighbours
      do r=2,r4-1
         r_bdy2(r,1)   = r-1
         r_bdy2(r,2)   = r+1
      end do
      r_bdy2(1,1)    = r4
      r_bdy2(1,2)    = 2
      r_bdy2(r4,1)   = r4-1
      r_bdy2(r4,2)   = 1

      !!get resolution on boundary
      do r=1,r4
         r_1         = r_bdy2(r,1)
         r_2         = r_bdy2(r,2)
         d1          = spherdist( glon(i_bdy(r),j_bdy(r)),glat(i_bdy(r),j_bdy(r)),  &
     &                   glon(i_bdy(r_1),j_bdy(r_1)),glat(i_bdy(r_1),j_bdy(r_1)) )
         d2          = spherdist( glon(i_bdy(r),j_bdy(r)),glat(i_bdy(r),j_bdy(r)),  &
     &                   glon(i_bdy(r_2),j_bdy(r_2)),glat(i_bdy(r_2),j_bdy(r_2)) )
         res_bdy(r)  = max(d1,d2)
!        chk_r = (r==r4).or.(r==1).or.(r==r1).or.(r==r2).or.(r==r3)                 &
         chk_r = .false.
         if (chk_r) then
            print*,' '
            print*, 'TWtest:', r,i_bdy(r),j_bdy(r)
            print*, 'TWtest:', r,res_bdy(r),glon(i_bdy(r),j_bdy(r)),glat(i_bdy(r),j_bdy(r))
            print*, 'TWtest:', r_1,i_bdy(r_1),j_bdy(r_1)
            print*, 'TWtest:', r_1,d1,glon(i_bdy(r_1),j_bdy(r_1)),glat(i_bdy(r_1),j_bdy(r_1))
            print*, 'TWtest:', r_2,i_bdy(r_2),j_bdy(r_2)
            print*, 'TWtest:', r_2,d2,glon(i_bdy(r_2),j_bdy(r_2)),glat(i_bdy(r_2),j_bdy(r_2))
            print*,' '
         end if
      end do

      max_bdy  = maxval(res_bdy)
      dumloc   = maxloc(res_bdy)
      !print*,'TW: ',dumloc,dumloc(1)
      r        = dumloc(1)
      r_1      = r_bdy2(r,1)
      r_2      = r_bdy2(r,2)
      print*,' '
      print*, 'TWtest max:', r,res_bdy(r),glon(i_bdy(r),j_bdy(r)),glat(i_bdy(r),j_bdy(r))
      print*, 'TWtest max:', r_1,d1,glon(i_bdy(r_1),j_bdy(r_1)),glat(i_bdy(r_1),j_bdy(r_1))
      print*, 'TWtest max:', r_2,d2,glon(i_bdy(r_2),j_bdy(r_2)),glat(i_bdy(r_2),j_bdy(r_2))
      print*,' '

      print*, 'Resolution of data around source boundary obtained'
      print*, 'Max res (m):',max_bdy
      print*, ' '
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- get the pivot points:
      ipiv  = 1
      jpiv  = 1
      do i=1,idm
         !print *,'i    :',i
         do j=1,jdm
       !  do j=100,243
            call pivotp_micom(plon(i,j),plat(i,j),glon,glat, ipiv,jpiv, nlon,nlat,dist,r_bdy )
         !print *,'ipiv,jpiv    :',ipiv,jpiv
             
            ipiv_o(i,j)    = ipiv
            jpiv_o(i,j)    = jpiv
            dist_mod(i,j)  = dist
            !dist_mod(i,j)  = spherdist(plon(i,j),plat(i,j), glon(ipiv,jpiv), glat(ipiv,jpiv))

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!double-check boundary pivot points
            !!(sometimes these are valid - may be important for low-res data)
            piv_on_bdy  = (ipiv==1 .or. jpiv==1 .or. ipiv==nlon .or. jpiv==nlat )

            if ( piv_on_bdy ) then
               !if ((jpiv==nlat).and.(dist.lt.1.5e4)) then
               !   print*,'pivot on bdy: ',ipiv,jpiv,r_bdy
               !   print*,'pivot on bdy: ',dist
               !end if
               dtol  = res_bdy(r_bdy)/2.0
               if (dist.gt.dtol) then
                  ipiv_o(i,j)    = 0
                  jpiv_o(i,j)    = 0
                  dist_mod(i,j)  = -1.0
               end if
            endif! piv_on_bdy

!           !! If pivot is on boundary, do a brute force check to see if it is valid
!           !! (sometimes this can be a valid point)
!           if ( piv_on_bdy ) then
!              if (dist.gt.max_bdy) then
!                 !!if it is very far away don't need to check
!                 ipiv_o(i,j)    = 0
!                 jpiv_o(i,j)    = 0
!                 dist_mod(i,j)  = -1.0
!              else
!                 dmin_bdy    = dist_mod(i,j)
!                 do r=1,r4
!                    d1 = spherdist( plon(i,j),plat(i,j),                              &
!       &                    glon(i_bdy(r),j_bdy(r)), glat(i_bdy(r),j_bdy(r)) )
!                    if (d1.le.dmin_bdy) then
!                       dmin_bdy = d1
!                       rmin_bdy = r
!                    end if
!                 end do

!                 r              = rmin_bdy
!                 ipiv_o(i,j)    = i_bdy(r)
!                 jpiv_o(i,j)    = j_bdy(r)
!                 dist_mod(i,j)  = dmin_bdy

!                 if (dmin_bdy.gt.0.5*res_bdy(r)) then
!                    ipiv_o(i,j)    = 0.0
!                    jpiv_o(i,j)    = 0.0
!                    dist_mod(i,j)  = -1.0
!                 end if

!                 print *,'dist, max_bdy, tol (m): ',dmin_bdy,max_bdy,res_bdy(r)/2.0
!                 print *,' '
!              end if

!           endif!piv_on_bdy
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         end do!j
      end do!i

      print*,'Max dist to pivot pt (m): ',maxval(dist_mod)

! --- write the pivot points into binary files:
      call zaiopf(trim(filenamea), 'replace', 909)

      call zaiowr(real(ipiv_o(:,:)),ip,.false.,hmin,hmax,909,.true.) 
      open (unit=909, file=filenameb,status='replace', action='write')
      write(909,'(a)') out_str
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

end program getpivots
