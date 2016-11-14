module mod_levitus
use mod_parameters

   private

   logical, parameter ::  silent=.true.
   real, dimension(:)  , allocatable :: levlon,levlat,levt,levz

   integer :: londimlen, latdimlen,zdimlenmonth,zdimlenannual,tdimlen
   logical :: woa2001=.false.,woa2005=.false.

   !character(len=80) :: save

   character(len=200), save :: pathwoa2001
   character(len=200), save :: pathwoa2005


   character(len=*), parameter ::  &
      woa2001sannual='WOA01.sal.annual.cdf',&
      woa2001smonth ='WOA01.sal.monthly.cdf',&
      woa2001tannual='WOA01.tem.annual.cdf', &
      woa2001tmonth ='WOA01.tem.monthly.cdf'
   character(len=*), parameter ::  &
      woa2005sannual='s00an1.nc',&
      woa2005tmonth ='t0112an1.nc',&
      woa2005smonth ='s0112an1.nc',&
      woa2005tannual='t00an1.nc'

   !character(len=250), save ::  &
   !   woa2001sannual,&
   !   woa2001smonth ,&
   !   woa2001tannual,&
   !   woa2001tmonth 
   !character(len=250), save ::  &
   !   woa2005sannual,&
   !   woa2005tmonth ,&
   !   woa2005smonth ,&
   !   woa2005tannual


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   public :: levitus_setup, stationsInterpLevitus, levitus_interp3d

contains



! Set up arrays for reading levitus data
   subroutine levitus_setup(cwoa)
      use netcdf
      implicit none
      character(len=*), intent(in) :: cwoa
      character(len=80):: filenamemonth,filenameannual
      integer :: ncresult, dimid,ncidmonth,ncidannual
      integer :: lonvarid, latvarid, zvarid, tvarid
      integer, dimension(nf90_max_dims) :: vdims ! Variable dimensions - not kept
      integer                           :: vndim ! number of variable dims - not kept
      real, allocatable :: levz2(:)
      integer :: k
      character(len=200) :: cenv


      if (trim(cwoa)=='WOA2001') then
         woa2001=.true.
         woa2005=.false.
         call getenv('WOA2001_PATH',cenv)
         if (trim(cenv)=='') then
            print *,'environment variable WOA2001_PATH not set '
            stop
         end if
         filenameannual=trim(cenv)//'/'//woa2001tannual
         filenamemonth =trim(cenv)//'/'//woa2001tmonth
      elseif (trim(cwoa)=='WOA2005') then
         woa2005=.true.
         woa2001=.false.
         call getenv('WOA2005_PATH',cenv)
         if (trim(cenv)=='') then
            print *,'environment variable WOA2005_PATH not set '
            stop
         end if
         filenameannual=trim(cenv)//'/'//woa2005tannual
         filenamemonth =trim(cenv)//'/'//woa2005tmonth
      else
         print *,'cwoa is ',cwoa
         stop '(levitus_setup -- invalid woa value)'
      end if


      if (.not. silent)print *,'reading Levitus fields from '//trim(filenamemonth)
      ncresult = nf90_open(trim(filenamemonth),nf90_nowrite,ncidmonth) 
      if (ncresult/= nf90_noerr) then
         call ncerr(ncresult)
      end if

      if (.not. silent)print *,'reading deep Levitus fields from '//trim(filenameannual)
      ncresult = nf90_open(trim(filenameannual),nf90_nowrite,ncidannual) 
      if (ncresult/= nf90_noerr) then
         call ncerr(ncresult)
      end if

      ! Dimension IDs and lengths (assuming the names are right)
      if (woa2001) then
         call getdimms(ncidmonth,'X',dimid,londimlen )
         call getdimms(ncidmonth,'Y',dimid,latdimlen )
         call getdimms(ncidannual,'Z',dimid,zdimlenannual ) ! NB - annual !
         call getdimms(ncidmonth,'Z',dimid,zdimlenmonth ) ! NB - month !
         call getdimms(ncidmonth,'T',dimid,tdimlen )
      elseif (woa2005) then
         call getdimms(ncidmonth,'lon',dimid,londimlen )
         call getdimms(ncidmonth,'lat',dimid,latdimlen )
         call getdimms(ncidannual,'depth',dimid,zdimlenannual ) ! NB - annual !
         call getdimms(ncidmonth,'depth',dimid,zdimlenmonth ) ! NB - month !
         call getdimms(ncidmonth,'time',dimid,tdimlen )
      else
         print *,'(levitus_setup: woa flag)'
      end if
      !print *,londimlen,latdimlen,zdimlen,tdimlen

      ! Allocate fields
      if (.not. allocated(levlon )) allocate(levlon(londimlen))
      if (.not. allocated(levlat )) allocate(levlat(latdimlen))
      if (.not. allocated(levz   )) allocate(levz  (  zdimlenannual))
      if (.not. allocated(levz2  )) allocate(levz2 (  zdimlenmonth))
      if (.not. allocated(levt   )) allocate(levt  (  tdimlen))


      ! Read z, t, and lon and lat levels
      if (woa2001) then
         call getvarstat(ncidmonth,'X',lonvarid,vndim  ,vdims)
         call getvarstat(ncidmonth,'Y',latvarid,vndim  ,vdims)
         call getvarstat(ncidannual,'Z',  zvarid,vndim  ,vdims) ! Annual
         call getvarstat(ncidmonth,'T',  tvarid,vndim  ,vdims)
      elseif (woa2005) then
         call getvarstat(ncidmonth,'lon',lonvarid,vndim  ,vdims)
         call getvarstat(ncidmonth,'lat',latvarid,vndim  ,vdims)
         call getvarstat(ncidannual,'depth',  zvarid,vndim  ,vdims) ! Annual
         call getvarstat(ncidmonth,'time',  tvarid,vndim  ,vdims)
      else
         print *,'(levitus_setup: woa flag)'
      end if
      call ncerr( nf90_get_var(ncidmonth,   tvarid, levt  ) )
      call ncerr( nf90_get_var(ncidannual,   zvarid, levz  ) ) ! Annual
      call ncerr( nf90_get_var(ncidannual,   zvarid, levz2 ) ) ! month
      call ncerr( nf90_get_var(ncidmonth, lonvarid, levlon) )
      call ncerr( nf90_get_var(ncidmonth, latvarid, levlat) )
      if (woa2005) then
         do k=1,12
            levt(k)=k-0.5
         end do
      end if
      ncresult = nf90_close(ncidmonth)
      ncresult = nf90_close(ncidannual)
      !print *,levz
      !print *,levt
      do k=1,zdimlenmonth
         !print *,levz(k),levz2(k)
         if (abs(levz(k)-levz2(k))>1e-4) then
            print *,'annual and monthly levitus difference..'
            print *,'index ',k
            print *,levz(k),levz2(k)
            stop
         end if
      end do
   end subroutine


   subroutine stationsInterpLevitus(varname, &
                             csjuld,cslon,cslat,cspres,clevvar,nprof,nzlevel,&
                             month)
   use netcdf
   use mod_year_info
   implicit none
   character(len=*), intent(in) :: varname
   character(len=80) :: filenamemonth, filenameannual
   character(len=20) :: varannual, varmonth
   integer, intent(in) :: nprof,nzlevel
   real, intent(in), dimension(nprof) :: csjuld, cslon, cslat
   real, intent(in), dimension(nzlevel,nprof) :: cspres
   real, dimension(nzlevel,nprof) :: clevvar
   integer, intent(in),optional :: month

   real, allocatable,dimension(:,:) :: levprof, levfld, tmp,a
   integer, dimension(:), allocatable :: ipiv,jpiv
   real ripiv,rjpiv,t,u,lon1,lat1
   integer :: ip, ipb, jp, jpb
   real :: fillvalue_var


   integer :: i,lwind,upind,k,k2
   integer :: iyear, iday, imonth
   real    :: rmonth,wlw,wup
   integer :: varidmonth, varidannual, ncidmonth, ncidannual
   integer :: ncresult
   integer :: vndim,vdims(nf90_max_dims)

   call checkinput(varname,varannual,varmonth,filenameannual,filenamemonth)
   if (present(month)) then
      call getMonthIndex(csjuld(1),lwind,upind,wup,wlw,month)
   else
      call getMonthIndex(csjuld(1),lwind,upind,wup,wlw)
   end if

   ! Profile positions on lev grid + bilin coeff
   allocate(ipiv(nprof))
   allocate(jpiv(nprof))
   allocate(a   (nprof,4))
   do k=1,nprof
      call getPivots(cslon(k),cslat(k),ipiv(k),jpiv(k),a(k,:))
   end do


   ! Get levitus data for profiles - on levitus vert grid
   allocate(levprof(zdimlenannual,nprof))
   allocate(tmp(londimlen,latdimlen))
   allocate(levfld(londimlen,latdimlen))

   !print *,'reading Levitus fields from '//trim(filenamemonth)
   call openLevitusNC(trim(filenamemonth),trim(varmonth),varidmonth,ncidmonth,fillvalue_var)
   call openLevitusNC(trim(filenameannual),trim(varannual),varidannual,ncidannual,fillvalue_var)

   do k=1,zdimlenannual
      ! time interpolated saln for this level interpolation is done in subroutines
      levfld=0.
      call getVarLevitusNC(ncidmonth,varidmonth,ncidannual,varidannual,k,lwind,wlw,levfld)
      call getVarLevitusNC(ncidmonth,varidmonth,ncidannual,varidannual,k,upind,wup,levfld)
      do k2=1,nprof
         ! Horizontal interpolation - quick fix version
         call interpHorizontal(ipiv(k2),jpiv(k2),a(k2,:),fillvalue_var,levfld,levprof(k,k2))
      end do
   end do
   ncresult = nf90_close(ncidannual)
   ncresult = nf90_close(ncidmonth)

   ! Final step - interpolate onto coriolis vertical levels
   do i=1,nprof
      do k=1,nzlevel
         if (cspres(k,i)/=undef) then

            !call findLayers(cspres(k,i),upind,lwind)
            call findLayers(cspres(k,i),upind,lwind,wup,wlw)
 
            !print *,upind,lwind,wup,wlw
            if (abs(wup+wlw-1.)<1e-4) then
               clevvar(k,i)=levprof(lwind,i)*wlw+levprof(upind,i)*wup
            else
               clevvar(k,i)=undef
            end if
         else
            clevvar(k,i)=undef
         end if
      end do
   end do
   end subroutine

   subroutine levitus_interp3d(varname, juld,plon,plat,pres,levint, &
                               idm,jdm,kdm, month)
   use netcdf
   use mod_year_info
   implicit none
   character(len=*), intent(in) :: varname
   character(len=80) :: filenamemonth, filenameannual
   character(len=20) :: varannual, varmonth
   integer, intent(in) :: idm,jdm,kdm
   real, intent(in) :: juld
   real, intent(in), dimension(idm,jdm)      :: plon, plat
   real, intent(in), dimension(kdm)          :: pres
   real, intent(out), dimension(idm,jdm,kdm) :: levint
   integer, intent(in),optional :: month

   real   , dimension(londimlen,latdimlen)   :: levfldlw, levfldup
   real   , dimension(idm,jdm,4) :: a
   integer, dimension(idm,jdm)   :: ipiv,jpiv
   real   , dimension(idm,jdm)   :: int_levfldlw, int_levfldup
   real :: fillvalue_var

   integer :: i,j,k,lwind,upind
   integer :: z_upind, z_lwind
   real    :: z_wlw, z_wup
   integer :: iyear, iday, imonth
   real    :: rmonth,wlw,wup, a2(4)
   integer :: varidmonth, varidannual, ncidmonth, ncidannual
   integer :: ncresult, ios


   call checkinput(varname,varannual,varmonth,filenameannual,filenamemonth)
   if (present(month)) then
      call getMonthIndex(juld,lwind,upind,wup,wlw,month)
   else
      call getMonthIndex(juld,lwind,upind,wup,wlw)
   end if

   ! Profile positions on lev grid + bilin coeff
   do j=1,jdm
   do i=1,idm
      call getPivots(plon(i,j),plat(i,j),ipiv(i,j),jpiv(i,j),a(i,j,1:4))
   end do
   end do


   !print *,'reading Levitus fields from '//trim(filenamemonth)
   call openLevitusNC(trim(filenamemonth),trim(varmonth),varidmonth,ncidmonth,fillvalue_var)
   call openLevitusNC(trim(filenameannual),trim(varannual),varidannual,ncidannual,fillvalue_var)

   ! Final step - interpolate onto coriolis vertical levels
   do k=1,kdm
      call findLayers(pres(k),z_upind,z_lwind,z_wup,z_wlw)
      !print *,'z layers ',k,z_upind, z_lwind, z_wup,z_wlw

      ! Time interpolation of lev level z_lwind
      levfldlw=0.
      call getVarLevitusNC(ncidmonth,varidmonth,ncidannual,varidannual,z_lwind,lwind,wlw,levfldlw)
      call getVarLevitusNC(ncidmonth,varidmonth,ncidannual,varidannual,z_lwind,upind,wup,levfldlw)
      !print *,'lowlevel :',z_lwind,maxval(levfldlw),minval(levfldlw)

      ! Time interpolation of lev level z_upind
      levfldup=0.
      call getVarLevitusNC(ncidmonth,varidmonth,ncidannual,varidannual,z_upind,lwind,wlw,levfldup)
      call getVarLevitusNC(ncidmonth,varidmonth,ncidannual,varidannual,z_upind,upind,wup,levfldup)
      !print *,'hilevel :',z_upind,maxval(levfldup),minval(levfldup)

      if (abs(z_wup+z_wlw-1.)<1e-4) then

         do j=1,jdm
         do i=1,idm

            ! Horizontal interpolation - not optimal without compiler help :-)
            call interpHorizontal(ipiv(i,j),jpiv(i,j),a(i,j,:),fillvalue_var,levfldlw,int_levfldlw(i,j))
            call interpHorizontal(ipiv(i,j),jpiv(i,j),a(i,j,:),fillvalue_var,levfldup,int_levfldup(i,j))

            ! Final step - vertical interpolation 
            levint(i,j,k)=int_levfldlw(i,j)*z_wlw+int_levfldup(i,j)*z_wup

         end do
         end do
         !print *,'levint lw',maxval(int_levfldlw),minval(int_levfldlw)
         !print *,'levint up',maxval(int_levfldup),minval(int_levfldup)
         !print *,'final',maxval(levint(:,:,k)),minval(levint(:,:,k))

      else
         levint(:,:,k)=undef
      end if
   end do
   end subroutine





   subroutine checkinput(varname,varannual,varmonth,filenameannual,filenamemonth)
   implicit none
   character(len=*), intent(in)  :: varname
   character(len=*), intent(out) :: varannual,varmonth
   character(len=*), intent(out) :: filenameannual,filenamemonth
   character(len=250) :: cenv
   if (woa2001) then
      call getenv('WOA2001_PATH',cenv)
      if (trim(cenv)=='') then
         print *,'environment variable WOA2001_PATH not set '
         stop
      end if
      if (varname=='salinity') then
         filenameannual=trim(cenv)//'/'//woa2001sannual
         filenamemonth =trim(cenv)//'/'//woa2001smonth
         varannual='salinity'
         varmonth ='salinity'
      elseif (varname=='temperature') then
         !filenameannual=woa2001tannual
         !filenamemonth =woa2001tmonth
         filenameannual=trim(cenv)//'/'//woa2001tannual
         filenamemonth =trim(cenv)//'/'//woa2001tmonth
         varannual='temperature'
         varmonth ='temperature'
      else
         print *,'Unknown field '//trim(varname)
      end if
   elseif (woa2005) then
      call getenv('WOA2005_PATH',cenv)
      if (trim(cenv)=='') then
         print *,'environment variable WOA2005_PATH not set '
         stop
      end if
      if (varname=='salinity') then
         !filenameannual=woa2005sannual
         !filenamemonth =woa2005smonth
         filenameannual=trim(cenv)//'/'//woa2005sannual
         filenamemonth =trim(cenv)//'/'//woa2005smonth
         varannual='s00an1'
         varmonth ='s0112an1'
      elseif (varname=='temperature') then
         !filenameannual=woa2005tannual
         !filenamemonth =woa2005tmonth
         filenameannual=trim(cenv)//'/'//woa2005tannual
         filenamemonth =trim(cenv)//'/'//woa2005tmonth
         varannual='t00an1'
         varmonth ='t0112an1'
      else
         print *,'Unknown field '//trim(varname)
      end if
   else
      print *,'clim flag not set '
      stop '(levitus_interp)'
   end if
   end subroutine
      

   subroutine getMonthIndex(csjuld,lwind,upind,wup,wlw,month)
      use mod_year_info
      implicit none
      real,intent(in) :: csjuld
      real,    intent(out) :: wup, wlw
      integer, intent(out) :: upind,lwind
      integer, intent(in), optional :: month


      integer :: iyear,imonth,iday
      real    :: rmonth
      integer :: i

      ! Get month,   day and year
      call juliantodate(floor(csjuld),iyear,imonth,iday,1950,1,1)
      rmonth=(imonth-1)+min(iday/30.,1.) ! Close enough

      ! find closest months in lev data
      lwind=tdimlen !"lower index"
      do i=1,tdimlen
         if (levt(i)<=rmonth) then
            lwind=i
         end if
      end do
         
      upind=1 ! Upper index
      do i=tdimlen,1,-1
         if (levt(i)>rmonth) then
            upind=i
         end if
      end do


      ! Weights for upper and lower time index
      if (levt(upind)-levt(lwind)>0.) then
         wup=(rmonth-levt(lwind))/(levt(upind)-levt(lwind))
      else
         wup=(rmonth+12.-levt(lwind))/(levt(upind)+12.-levt(lwind))
      end if
      wlw=1.-wup
      !print *,lwind,upind,rmonth,wlw,wup
      !print *,levt(lwind),levt(upind)

      if (present(month)) then
         upind=month ; wup=1.
         lwind=1     ; wlw=0.
      end if
   end subroutine getMonthIndex
        
   subroutine getPivots(lon,lat,ipiv,jpiv,a)
      implicit none
      real   , intent(in) :: lon,lat
      integer, intent(out) :: ipiv,jpiv
      real   , intent(out) :: a(4)

      real :: lon1,lat1, ripiv, t, rjpiv,u

      lon1=lon
      lat1=lat

      if(lon1<levlon(1)) lon1=lon1+360.

      ! index into lon/lat matrix 
      ripiv=(lon1-levlon(1))/(levlon(2)-levlon(1))+1.
      ipiv=floor(ripiv);
      t=ripiv-ipiv

      rjpiv=(lat1-levlat(1))/(levlat(2)-levlat(1))+1.
      jpiv=floor(rjpiv);
      u=rjpiv-jpiv

      a(1)=(1-t)*(1-u)
      a(2)=    t*(1-u)
      a(3)=    t*u
      a(4)=(1-t)*u     
      !if (ipiv<0) print *,lon1,lat1,ipiv
   end subroutine


   subroutine openLevitusNC(fname,varname,varid,ncid,fillvalue)
   use netcdf
   implicit none
      character(len=*), intent(in)  :: fname,varname
      integer         , intent(out) :: ncid,varid
      real            , intent(out) :: fillvalue

      integer :: ncresult
      integer :: vndim,vdims(nf90_max_dims)

      ncresult = nf90_open(trim(fname),nf90_nowrite,ncid) 
      if (ncresult/= nf90_noerr) then
         call ncerr(ncresult)
      end if
      call getvarstat(ncid,trim(varname),varid,vndim  ,vdims)
      ncresult=nf90_get_att(ncid, varid   , 'missing_value', fillvalue)
   end subroutine


      ! time interpolated saln for this level
   subroutine getVarLevitusNC(ncidmonth,varidmonth,ncidannual,varidannual,k,ind,w,levfld)
   use netcdf
   implicit none
      integer, intent(in) :: ncidmonth,varidmonth,ncidannual,varidannual,k,ind
      real,    intent(in) :: w
      real, intent(inout) :: levfld(londimlen,latdimlen)

      real :: tmp(londimlen,latdimlen)

      if (k<=zdimlenmonth) then
         call ncerr( nf90_get_var(ncidmonth,varidmonth, tmp    , start=(/1,1,k,ind/)))
         levfld=levfld+tmp*w
      else
         call ncerr( nf90_get_var(ncidannual,varidannual, tmp    , start=(/1,1,k/)))
         levfld=levfld+tmp*w
      end if
   end subroutine
         
   subroutine interpHorizontal(ipin,jpin,ain,fillvalue,levfld,levprof)
      implicit none
      integer, intent(in)  :: ipin,jpin
      real,    intent(in)  :: ain(4),fillvalue,levfld(londimlen,latdimlen)
      real,    intent(out) :: levprof

      integer ipb,jpb,ip,jp
      real :: a(4)

      ip=ipin
      jp=jpin
      a=ain

      ! Horizontal interpolation
      ipb=mod(ip,londimlen)+1
      jp =min(max(1,jp  ),latdimlen)
      jpb=min(max(1,jp+1),latdimlen)
      if ((abs(levfld(ip  ,jp  )-fillvalue))/fillvalue <1e-4) a(1)=0.
      if ((abs(levfld(ipb ,jp  )-fillvalue))/fillvalue <1e-4) a(2)=0.
      if ((abs(levfld(ipb ,jpb )-fillvalue))/fillvalue <1e-4) a(3)=0.
      if ((abs(levfld(ip  ,jpb )-fillvalue))/fillvalue <1e-4) a(4)=0.
      if (sum(a) > .1) then
         a= a/sum(a)
         levprof = levfld(ip  ,jp  )*a(1) + &
                   levfld(ipb ,jp  )*a(2) + &
                   levfld(ipb ,jpb )*a(3) + &
                   levfld(ip  ,jpb )*a(4)   
      else
         levprof = undef
      end if
   end subroutine


   subroutine findLayers(cspres,upind,lwind,wup,wlw)
      implicit none
      real, intent(in) :: cspres !, levz(zdimlenannual)
      integer, intent(out) :: lwind, upind
      real,    intent(out) :: wlw, wup

      integer :: k2

      ! find levitus depth index right above level cspres
      k2=1
      do while (cspres>=levz(k2) .and.k2<=zdimlenannual)
         k2=k2+1
      end do
      k2=k2-1
      upind=k2
      lwind=min(k2+1,zdimlenannual)

      if (upind/=lwind) then
         wup=(levz(lwind)-cspres)/(levz(lwind)-levz(upind))
         wlw=1.-wup
      else
         wup=0.
         wlw=0.
      end if
   end subroutine

   

   ! Give netcdf error message
   subroutine ncerr(error)
      use netcdf
      implicit none
      integer, intent(in) :: error
      if (error/=nf90_noerr) then
         print *,'mod_levitus: '//nf90_strerror(error)
         stop
      end if
   end subroutine ncerr




   subroutine getdimms(ncid,dname,dmid,dmlen)
      use netcdf
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: dname
      integer,          intent(out) :: dmid,dmlen
      call ncerr( nf90_inq_dimid(ncid,dname,dmid) )
      call ncerr( nf90_inquire_dimension(ncid,dmid,len=dmlen) )
   end subroutine

   subroutine getvarstat(ncid,varname,varid,var_ndim,var_dimids)
      use netcdf
      implicit none
      integer,                          intent(in)  :: ncid
      character(len=*),                 intent(in)  :: varname
      integer,                          intent(out) :: varid,var_ndim
      integer, dimension(nf90_max_dims),intent(out) :: var_dimids
      call ncerr( nf90_inq_varid(ncid,varname,varid) )
      call ncerr( nf90_inquire_variable(ncid,varid,ndims=var_ndim,dimids=var_dimids) )
   end subroutine getvarstat

end module mod_levitus
