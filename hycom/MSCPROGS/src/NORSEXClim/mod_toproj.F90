!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! mod_toregugrid:
!
! A collection of routines and variables for converting
! the conformal mapping to a regular grid levels
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_toproj


! These names are "straightforward" for a regular grid. For other
! projections they have other meanings
   real, save :: firstlon = -999.
   real, save :: lastlon  = -999.
   real, save :: dlon     = -999.

   real, save :: firstlat = -999.
   real, save :: lastlat  = -999.
   real, save :: dlat     = -999.

   ! polar stereographic specific parameters
   real, save :: polar_stereographic_central_lon=-45;
   real, save :: polar_stereographic_central_lat= 90;

   real, save :: proj_conv_fac     = 1.

   integer :: nlons
   integer :: nlats

   logical, save :: ini_called = .false.

   real, save, allocatable :: lons(:,:)
   real, save, allocatable :: lats(:,:)

   logical, save, dimension(:,:), allocatable :: ongrid
   real   , save, dimension(:,:), allocatable :: a1s, a2s, a3s, a4s
   real   , save, dimension(:,:), allocatable :: xproj,yproj
   integer, save, dimension(:,:), allocatable :: ipivs, jpivs

   logical,save :: testflag=.false.
   logical,save :: gridrotate=.false.

   real, save :: undef_regu

   character(len=40), save :: cprojection=''


   interface to_proj
      module procedure to_proj2d
      module procedure to_proj3d
   end interface



contains


subroutine proj_ini(undef,modlon,modlat,ip,idm,jdm)
#if defined (INITCONFMAP)
   use m_initconfmap
   use m_bilincoeff
   use m_oldtonew
   use m_pivotp
#endif
   implicit none
   integer, intent(in) :: idm,jdm
   integer, intent(in) :: ip(idm,jdm)
   real   , intent(in) :: modlon(idm,jdm)
   real   , intent(in) :: modlat(idm,jdm)
   real,    intent(in) :: undef


   integer :: ipiv, jpiv, ipib, jpib,ipsum,ilon,jlat,i,j
   real    :: t1,t2,a1,a2,a3,a4,lon_n, lat_n, asum
   logical :: gmsk(idm,jdm),ass
   real :: tmp1,tmp2,e,x,y
   integer :: indx, indy
   real    :: t,u, ripiv, rjpiv

   undef_regu=undef


   ! Check that the "native" projection option stays within model bounds
   ! (We couldn't check this until now)
   if (trim(cprojection)=='native') then
      if (firstlon>idm .or. lastlon>idm .or. &
          firstlat>jdm .or. lastlat>jdm ) then
          print *,'For native, bounds must be within grid'
          stop '(hyc2proj)'
      end if
   end if
      

   ! Error check
   if (dlon<0. .or. dlat<0.) then
      print *,'negative dlon or dlat!'
      stop '(regugrid_ini)'
   end if

   if (firstlon>lastlon .or. firstlat>lastlat) then
      print *,'last lon/lat less than first lon/lat!'
      stop '(regugrid_ini)'
   end if

   nlons=floor((lastlon-firstlon)/dlon) + 1
   nlats=floor((lastlat-firstlat)/dlat) + 1
   ini_called=.true.

   allocate(ongrid(nlons,nlats))
   allocate(ipivs (nlons,nlats))
   allocate(jpivs (nlons,nlats))
   allocate(a1s   (nlons,nlats))
   allocate(a2s   (nlons,nlats))
   allocate(a3s   (nlons,nlats))
   allocate(a4s   (nlons,nlats))
   allocate(lons  (nlons,nlats))
   allocate(lats  (nlons,nlats))
   allocate(xproj (idm,jdm))
   allocate(yproj (idm,jdm))

   if (trim(cprojection)=='regular') then
      do ilon=1,nlons
         lons(ilon,:)=firstlon+(ilon-1)*dlon
      end do
      do jlat=1,nlats
         lats(:,jlat)=firstlat+(jlat-1)*dlat
      end do
!  else if (trim(cprojection)=='ps_amsr') then
!     do jlat=1,nlats
!     do ilon=1,nlons
!        tmp1=firstlon+(ilon-1)*dlon
!        tmp2=firstlat+(jlat-1)*dlat
!        call ps_amsr(lons(ilon,jlat),lats(ilon,jlat),tmp2,tmp1)
!
!        !print *,tmp1,tmp2,lons(ilon,jlat),lats(ilon,jlat)
!     end do
!     end do
   else if (trim(cprojection)=='polar_stereographic') then
      do jlat=1,nlats
      do ilon=1,nlons
         tmp1=firstlon+(ilon-1)*dlon
         tmp2=firstlat+(jlat-1)*dlat
         call polar_stereographic(lons(ilon,jlat),lats(ilon,jlat), &
                                  tmp1,tmp2,'xy2ll')
      end do
      end do

      ! Also fetch grid
      do j=1,jdm
      do i=1,idm
         
         ! Because of intent(in)
         t=modlon(i,j)
         u=modlat(i,j)

         call polar_stereographic(t,u, &
                                  xproj(i,j),yproj(i,j),'ll2xy')
      end do
      end do
   end if


   if (.not. ini_called) then
      print *,'Initialize regugrid first!'
      stop '(to_regugrid)'
   end if

   ! Initialize conformal mapping
#if defined (INITCONFMAP)
   call initconfmap
#endif


   if (trim(cprojection)/='native') then

      do jlat=1,nlats
      do ilon=1,nlons

         ! Positions to interpolate to
         t1=lons(ilon,jlat)
         t2=lats(ilon,jlat)

         ! corresponding pivot points
#if defined (INITCONFMAP)
         call oldtonew(t2,t1,lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)
#elif defined(NORSEXGRID)
         call norsex_locate(2,1,2,ripiv,rjpiv,t1,t2)
         ipiv=floor(ripiv)
         jpiv=floor(rjpiv)
         !print *,ipiv,jpiv
#else
#error -
#endif

         ! Interpolation to regugrid
         if ((ipiv < 1).or.(ipiv > idm).or.(jpiv < 1).or.(jpiv > jdm)) then
            ongrid(ilon,jlat)=.false.
         else
            ipib=min(ipiv+1,idm)
            JPib=min(jpiv+1,jdm)
            ipsum=sum(ip(ipiv:ipiv+1,jpiv:jpiv+1))
            !if (ipsum >= 2) then
            if (ipsum >= 3) then
#if defined (INITCONFMAP)
               call bilincoeff(modlon,modlat,idm,jdm,t1,t2,ipiv,jpiv,&
                                a1,a2,a3,a4)
#elif defined(NORSEXGRID)
               t=(ripiv-ipiv)
               u=(rjpiv-jpiv)
               a1=(1-t)*(1-u)
               a2=t*(1-u)
               a3=t*u
               a4=(1-t)*u
#else
#error -
#endif
               
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
      do jlat=1,nlats
      do ilon=1,nlons
         indx=nint(firstlon)+(ilon-1)*nint(dlon)
         indy=nint(firstlat)+(jlat-1)*nint(dlat)
         lons(ilon,jlat)=indx
         lats(ilon,jlat)=indy
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
         !print *,'proj_ini ',indx,indy,ip(indx,indy)
      end do
      end do
   end if

end subroutine proj_ini



subroutine to_proj2d(varin,varout,idm,jdm)
   !use mod_dimensions
   implicit none
   integer, intent(in)  :: idm,jdm
   real,    intent(in)  :: varin (idm,jdm)
   real,    intent(out) :: varout(nlons,nlats)

   integer :: ilon,jlat,ipiv,jpiv,jpib,ipib,indx,indy
   real a1,a2,a3,a4
   if (trim(cprojection)/='native') then
      do jlat=1,nlats
      do ilon=1,nlons
      if (ongrid(ilon,jlat)) then
         a1=a1s(ilon,jlat)
         a2=a2s(ilon,jlat)
         a3=a3s(ilon,jlat)
         a4=a4s(ilon,jlat)
         ipiv=ipivs(ilon,jlat)
         jpiv=jpivs(ilon,jlat)
         ipib=min(ipiv+1,idm)
         JPib=min(jpiv+1,jdm)


         varout(ilon,jlat) = a1*varin(ipiv,jpiv) &
                           + a2*varin(ipib,jpiv) &
                           + a3*varin(ipib,jpib) &
                           + a4*varin(ipiv,jpib)

         !if (testflag) print *,varout(ilon,jlat),varin(ipiv,jpiv)
      else
         varout(ilon,jlat)=undef_regu
      end if
      end do
      end do
   else
      do jlat=1,nlats
      do ilon=1,nlons
      if (ongrid(ilon,jlat)) then
         indx = nint(firstlon)+(ilon-1)*nint(dlon)
         indy = nint(firstlat)+(jlat-1)*nint(dlat)
         !print *,ilon,jlat,indx,indy,varin(indx,indy)
         varout(ilon,jlat) = varin(indx,indy)
      else
         varout(ilon,jlat) = undef_regu
      end if
      end do
      end do
   end if
      
end subroutine to_proj2d


subroutine to_proj3d(varin,varout,idm,jdm,kdm)
   implicit none
   integer, intent(in)  :: idm,jdm,kdm
   real,    intent(in)  :: varin (idm,jdm,kdm)
   real,    intent(out) :: varout(nlons,nlats,kdm)


   integer :: k

   do k=1,kdm
      call to_proj2d(varin(:,:,k),varout(:,:,k),idm,jdm)
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
   if (trim(cprojection)/='regular' .and. trim(cprojection)/='ps_amsr' .and. &
       trim(cprojection)/='polar_stereographic' .and. trim(cprojection)/='native') then
      print *,'unknown projection "'//trim(cprojection)//'"'
      stop 
   end if
   print '(a)','Projection:'//trim(cprojection)
   read(10,*,iostat=ios) firstlon, tmpchar
   print '(a,f10.3)','SW Lon   :',firstlon
   read(10,*,iostat=ios) lastlon, tmpchar
   print '(a,f10.3)','NE lon   :',lastlon
   read(10,*,iostat=ios) dlon   , tmpchar
   print '(a,f10.3)','Delta lon:', dlon
   read(10,*,iostat=ios) firstlat, tmpchar
   print '(a,f10.3)','SW lat   :',firstlat
   read(10,*,iostat=ios) lastlat, tmpchar
   print '(a,f10.3)','NE lat   :',lastlat
   read(10,*,iostat=ios) dlat   , tmpchar
   print '(a,f10.3)','Delta lat:', dlat
   if (trim(cprojection)=='polar_stereographic') then
      read(10,*,iostat=ios) polar_stereographic_central_lon,tmpchar
      print '(a,f10.3)','polar_stereographic central longitude:', &
         polar_stereographic_central_lon
      read(10,*,iostat=ios) polar_stereographic_central_lat,tmpchar
      print '(a,f10.3)','polar_stereographic central latitude:', &
         polar_stereographic_central_lat
      read(10,*,iostat=ios) proj_conv_fac,tmpchar
      print '(a,f10.3)','polar_stereographic conv factor:', proj_conv_fac

      firstlon=firstlon*proj_conv_fac
      lastlon =lastlon *proj_conv_fac
      firstlat=firstlat*proj_conv_fac
      lastlat =lastlat *proj_conv_fac
      dlon    =dlon    *proj_conv_fac
      dlat    =dlat    *proj_conv_fac

      read(10,*,iostat=ios) gridrotate,tmpchar
      print '(a,l5)','polar_stereographic rotate to output grid:', gridrotate


   end if
   close(10)

   ! Safet check for the "native" option
   if (trim(cprojection)=='native') then
      if (abs(dlat-nint(dlat))>1e-6 .or.  abs(dlon-nint(dlon))>1e-6) then
         print *,'For native - dellat and dellon should be integers'
         stop '(get_proj)'
      else if (dlat<1. .or. dlon<1. .or. firstlon<1. .or. firstlat<1. .or. &
               lastlon<1. .or. lastlat<1.) then
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
