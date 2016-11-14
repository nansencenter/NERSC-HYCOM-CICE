!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! mod_toregugrid:
!
! A collection of routines and variables for converting
! the conformal mapping to a regular grid levels
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_toproj
use mod_parameters

! These names are "straightforward" for a regular grid. For other
! projections they have other meanings
   real, save :: firstx = -999. ! First x-point in projection coords
   real, save :: lastx  = -999. ! Last  x-point in projection coords
   real, save :: dx     = -999. ! x-step
   integer :: nxp               ! Number of points in x-dir

   real, save :: firsty = -999. ! First y-point in projection coords
   real, save :: lasty  = -999. ! Last  x-point in projection coords
   real, save :: dy     = -999. ! y-step
   integer :: nyp               ! Number of points in y-dir

   ! polar stereographic specific parameters
   real, save :: polar_stereographic_central_lon=-45;
   real, save :: polar_stereographic_central_lat= 90;
   real, save :: proj_conv_fac     = 1.

   logical, save :: ini_called = .false.

   real, save, allocatable :: lons(:,:) !Longitudes of xy-points
   real, save, allocatable :: lats(:,:) !Latitudes  of xy-points

   logical, save, dimension(:,:), allocatable :: ongrid             ! Mask
   real   , save, dimension(:,:), allocatable :: a1s, a2s, a3s, a4s ! Bilinear coeffs
   real   , save, dimension(:,:), allocatable :: xproj,yproj        ! Projection coords - used in rotation
   integer, save, dimension(:,:), allocatable :: ipivs, jpivs       ! Poivot points

   logical,save :: testflag=.false.
   logical,save :: gridrotate=.false. 
   character(len=40), save :: cprojection=''


   interface to_proj
      module procedure to_proj2d
      module procedure to_proj3d
   end interface

contains


subroutine proj_ini()
   use mod_xc, only : idm,jdm
   use mod_grid, only : plon, plat, ip
   use mod_confmap , only: bilincoeff, initconfmap, oldtonew, pivotp
   implicit none

   integer :: ipiv, jpiv, ipib, jpib,ipsum,ilon,jlat,i,j
   real    :: t1,t2,a1,a2,a3,a4,lon_n, lat_n, asum
   logical :: gmsk(idm,jdm),ass
   real :: tmp1,tmp2,e,x,y
   integer :: indx, indy


   ! Reads proj.in and does some error checks
   call get_proj()

   ! Check that the "native" projection option stays within model bounds
   ! (We couldn't check this until now)
   if (trim(cprojection)=='native') then

      ! -1 ivalues means idm, jdm for these two
      if ((nint(lastx)==-1)) lastx=idm
      if ((nint(lasty)==-1)) lasty=jdm

      if (firstx>idm .or. lastx>idm .or. &
          firsty>jdm .or. lasty>jdm ) then
          print *,'For native, bounds must be within grid'
          stop '(proj_ini)'
      end if
   end if

   ! Error check
   if (dx<0. .or. dy<0.) then
      print *,'negative dx or dy!'
      stop '(proj_ini)'
   end if

   if (firstx>lastx .or. firsty>lasty) then
      print *,'last lon/lat less than first lon/lat!'
      stop '(proj_ini)'
   end if

   nxp=floor((lastx-firstx)/dx) + 1
   nyp=floor((lasty-firsty)/dy) + 1
   ini_called=.true.

   allocate(ongrid(nxp,nyp))
   allocate(ipivs (nxp,nyp))
   allocate(jpivs (nxp,nyp))
   allocate(a1s   (nxp,nyp))
   allocate(a2s   (nxp,nyp))
   allocate(a3s   (nxp,nyp))
   allocate(a4s   (nxp,nyp))
   allocate(lons  (nxp,nyp))
   allocate(lats  (nxp,nyp))
   allocate(xproj (idm,jdm))
   allocate(yproj (idm,jdm))

   if (trim(cprojection)=='regular') then
      do ilon=1,nxp
         lons(ilon,:)=firstx+(ilon-1)*dx
      end do
      do jlat=1,nyp
         lats(:,jlat)=firsty+(jlat-1)*dy
      end do
   else if (trim(cprojection)=='polar_stereographic') then
      do jlat=1,nyp
      do ilon=1,nxp
         tmp1=firstx+(ilon-1)*dx
         tmp2=firsty+(jlat-1)*dy
         call polar_stereographic(lons(ilon,jlat),lats(ilon,jlat), &
                                  tmp1,tmp2,'xy2ll')
      end do
      end do

      ! Also fetch grid
      do j=1,jdm
      do i=1,idm
         call polar_stereographic(plon(i,j),plat(i,j), &
                                  xproj(i,j),yproj(i,j),'ll2xy')
      end do
      end do
   end if


   ! Initialize conformal mapping
   call initconfmap(idm,jdm)
   if (trim(cprojection)/='native') then

      do jlat=1,nyp
      do ilon=1,nxp

         ! Positions to interpolate to
         t1=lons(ilon,jlat)
         t2=lats(ilon,jlat)

         ! corresponding pivot points
         call oldtonew(t2,t1,lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)

         ! Interpolation to regugrid
         if ((ipiv < 1).or.(ipiv > idm-1).or.(jpiv < 1).or.(jpiv > jdm-1)) then
            ongrid(ilon,jlat)=.false.
         else
            ipib=min(ipiv+1,idm)
            JPib=min(jpiv+1,jdm)
            ipsum=sum(ip(ipiv:ipiv+1,jpiv:jpiv+1))

            ! Should have at least three neighbours
            if (ipsum >= 3) then
               call bilincoeff(plon,plat,idm,jdm,t1,t2,ipiv,jpiv,&
                                a1,a2,a3,a4)
               
               a1=a1*ip(ipiv,jpiv)
               a2=a2*ip(ipib,jpiv)
               a3=a3*ip(ipib,jpib)
               a4=a4*ip(ipiv,jpib)

               asum = a1+a2+a3+a4
               a1=a1/asum
               a2=a2/asum
               a3=a3/asum
               a4=a4/asum

               a1s(ilon,jlat)=a1
               a2s(ilon,jlat)=a2
               a3s(ilon,jlat)=a3
               a4s(ilon,jlat)=a4
               ipivs(ilon,jlat) = ipiv
               jpivs(ilon,jlat) = jpiv
               ongrid(ilon,jlat)=.true.

            else
               ongrid(ilon,jlat)=.false.
            endif
         endif
      end do
      end do
   else
      do jlat=1,nyp
      do ilon=1,nxp
         indx=nint(firstx)+(ilon-1)*nint(dx)
         indy=nint(firsty)+(jlat-1)*nint(dy)
         lons(ilon,jlat)=plon(indx,indy)
         lats(ilon,jlat)=plat(indx,indy)
         if (indx <1 .or. indx > idm .or. indy<1 .or. indy>jdm) then
            ongrid(ilon,jlat)=.false.
         else if (ip(indx,indy)==0) then
            ongrid(ilon,jlat)=.false.
         else
            a1s  (ilon,jlat)=1.
            a2s  (ilon,jlat)=0.
            a3s  (ilon,jlat)=0.
            a4s  (ilon,jlat)=0.
            ipivs(ilon,jlat)=indx
            jpivs(ilon,jlat)=indy
            ongrid(ilon,jlat)=.true.
         end if
      end do
      end do
   end if
end subroutine proj_ini



subroutine to_proj2d(varin,varout)
   use mod_xc, only: idm,jdm
   implicit none
   real, intent(in) :: varin (idm,jdm)
   real, intent(out) :: varout(nxp,nyp)

   logical :: tstundef
   integer :: ilon,jlat,ipiv,jpiv,jpib,ipib,indx,indy
   real a1,a2,a3,a4,asum

   real, parameter :: athresh=0.3


   if (trim(cprojection)/='native') then
      do jlat=1,nyp
      do ilon=1,nxp
      if (ongrid(ilon,jlat)) then
         a1=a1s(ilon,jlat)
         a2=a2s(ilon,jlat)
         a3=a3s(ilon,jlat)
         a4=a4s(ilon,jlat)
         ipiv=ipivs(ilon,jlat)
         jpiv=jpivs(ilon,jlat)
         ipib=min(ipiv+1,idm)
         JPib=min(jpiv+1,jdm)

         ! Test where the corner values are undefined
         if (varin(ipiv,jpiv)==undef) a1=0.
         if (varin(ipib,jpiv)==undef) a2=0.
         if (varin(ipib,jpib)==undef) a3=0.
         if (varin(ipiv,jpib)==undef) a4=0.

         ! Is the sum of "valid" a's enough to make a meaningful extrapolation?
         asum=a1+a2+a3+a4
         a1=a1/asum
         a2=a2/asum
         a3=a3/asum
         a4=a4/asum
         if (asum>athresh) then
            varout(ilon,jlat) = a1*varin(ipiv,jpiv) &
                              + a2*varin(ipib,jpiv) &
                              + a3*varin(ipib,jpib) &
                              + a4*varin(ipiv,jpib)
         else
            varout(ilon,jlat)=undef
         end if
      else
         varout(ilon,jlat)=undef
      end if
      end do
      end do

   ! Native (i.e. original confmap grid) option
   else
      do jlat=1,nyp
      do ilon=1,nxp
      if (ongrid(ilon,jlat)) then
         indx = nint(firstx)+(ilon-1)*nint(dx)
         indy = nint(firsty)+(jlat-1)*nint(dy)
         varout(ilon,jlat) = varin(indx,indy)
      else
         varout(ilon,jlat) = undef
      end if
      end do
      end do
   end if
end subroutine to_proj2d


subroutine to_proj3d(varin,varout,kdm)
   use mod_xc, only :idm,jdm
   implicit none
   integer, intent(in) :: kdm
   real, intent(in) :: varin (idm,jdm,kdm)
   real, intent(out) :: varout(nxp,nyp,kdm)
   integer :: k

   do k=1,kdm
      call to_proj2d(varin(:,:,k),varout(:,:,k))
   end do
end subroutine to_proj3d



subroutine get_proj()
   implicit none
   logical :: ex
   integer :: ios
   character(len=*), parameter :: infile_proj='proj.in'
   character(len=40)  :: tmpchar


   ! Get Regular grid specification
   inquire(exist=ex,file=infile_proj)
   if (.not. ex) then
      print '(a)','You must specify grid in'
      print '(a)','the file '//infile_proj
      stop '(get_proj)'
   endif

   ! Start reading
   open(10,file=infile_proj,form='formatted',status='old',action='read')
   read(10,*,iostat=ios) cprojection
   if (trim(cprojection)/='regular' .and. &
       trim(cprojection)/='polar_stereographic' .and. trim(cprojection)/='native') then
      print *,'unknown projection "'//trim(cprojection)//'"'
      stop 
   end if
   print '(a)','Projection:'//trim(cprojection)
   read(10,*,iostat=ios) firstx, tmpchar
   print '(a,f10.3)','SW Lon   :',firstx
   read(10,*,iostat=ios) lastx, tmpchar
   print '(a,f10.3)','NE lon   :',lastx
   read(10,*,iostat=ios) dx   , tmpchar
   print '(a,f10.3)','Delta lon:', dx
   read(10,*,iostat=ios) firsty, tmpchar
   print '(a,f10.3)','SW lat   :',firsty
   read(10,*,iostat=ios) lasty, tmpchar
   print '(a,f10.3)','NE lat   :',lasty
   read(10,*,iostat=ios) dy   , tmpchar
   print '(a,f10.3)','Delta lat:', dy
   if (trim(cprojection)=='polar_stereographic') then
      read(10,*,iostat=ios) polar_stereographic_central_lon,tmpchar
      print '(a,f10.3)','polar_stereographic central longitude:', &
         polar_stereographic_central_lon
      read(10,*,iostat=ios) polar_stereographic_central_lat,tmpchar
      print '(a,f10.3)','polar_stereographic central latitude:', &
         polar_stereographic_central_lat
      read(10,*,iostat=ios) proj_conv_fac,tmpchar
      print '(a,f10.3)','polar_stereographic conv factor:', proj_conv_fac

      firstx=firstx*proj_conv_fac
      lastx =lastx *proj_conv_fac
      firsty=firsty*proj_conv_fac
      lasty =lasty *proj_conv_fac
      dx    =dx    *proj_conv_fac
      dy    =dy    *proj_conv_fac

      read(10,*,iostat=ios) gridrotate,tmpchar
      print '(a,l5)','polar_stereographic rotate to output grid:', gridrotate


   end if
   close(10)

   ! Safet check for the "native" option
   if (trim(cprojection)=='native') then
      if (abs(dy-nint(dy))>1e-6 .or.  abs(dx-nint(dx))>1e-6) then
         print *,'For native - dellat and dellon should be integers'
         stop '(get_proj)'
      else if (dy<1. .or. dx<1. .or. firstx<1. .or. firsty<1.) then
         print *,'invalid spec of native output.'
         stop '(get_proj)'
      else if ( (nint(lastx)<1 .and. nint(lastx)/=-1) .or.  &
                (nint(lasty)<1 .and. nint(lasty)/=-1) ) then
         print *,'invalid spec of native output.'
         stop '(get_proj)'
      end if
   end if

   if (ios/=0) then
      print *,'ERROR: occured when reading '//infile_proj
      stop '(get_proj)'
   endif
end subroutine get_proj

   ! 
   subroutine polar_stereographic(lon,lat,x,y,direction)
      implicit none

      real, intent(inout) :: x,y
      real, intent(inout) :: lon,lat
      character(len=*), intent(in) :: direction
      real, parameter :: re  =6378.273      ! Reference radius of earth [km]
      real, parameter :: CDR=57.29577951    ! radian to degree conv factor

      real :: rho,c

      real :: x2,y2,k

      if (trim(direction)=='xy2ll') then

         x2=x/proj_conv_fac ! proj_conv_fac to get x/y in km
         y2=y/proj_conv_fac ! proj_conv_fac to get x/y in km

         rho=sqrt(x2**2 +y2**2)
         c  =2*atan(rho/(2*RE))

         !lat=asin( cos(c)  )
         !lon=atan2(x2*sin(c),-y2*sin(c))

         !lat=lat*CDR
         !lon=lon*CDR-45

         lon=polar_stereographic_central_lon + &
            atan2(x2*sin(c) ,&
                  rho*cos(polar_stereographic_central_lat/CDR)*cos(c) - &
                  y2*  sin(polar_stereographic_central_lat/CDR)*sin(c))*CDR

         lat=asin(cos(c)*sin(polar_stereographic_central_lat/CDR) +  &
                  y2*sin(c)*cos(polar_stereographic_central_lat/CDR)/rho)*CDR

      else if (trim(direction)=='ll2xy') then

         !k=2*Re/(1.+sin(lat/CDR))
         !x= k*cos(lat/CDR)*sin((lon+45)/CDR)*proj_conv_fac
         !y=-k*cos(lat/CDR)*cos((lon+45)/CDR)*proj_conv_fac
         k=2*Re/(1.+sin(polar_stereographic_central_lat/CDR)*sin(lat/CDR)  &
                   +cos(polar_stereographic_central_lat/CDR)*cos(lat/CDR) *&
                    cos((lon-polar_stereographic_central_lon)/CDR))

         y= k*( &
            cos(polar_stereographic_central_lat/CDR)*sin(lat/CDR) &
          - sin(polar_stereographic_central_lat/CDR)*cos(lat/CDR)* &
               cos((lon-polar_stereographic_central_lon)/CDR))*proj_conv_fac
         x= k*(cos(lat/CDR)*sin((lon-polar_stereographic_central_lon)/CDR))*proj_conv_fac

      else

         print *,'Invalid direction specified'
         stop '(polar_stereographic)'
      end if
   end subroutine polar_stereographic

end module mod_toproj
