module mod_levitus

   private

   real, dimension(:)  , allocatable :: levlon,levlat,levt,levz

   integer :: londimlen, latdimlen,zdimlenmonth,zdimlenannual,tdimlen
   logical :: woa2001=.false.,woa2005=.false.

   character(len=80) :: save
! TODO - make this more flexible
! TODO - Replace with new mod_levitus from hyc2proj
!   character(len=*), parameter ::  &
!      woa2001sannual='/work2/knutali/WOA2001/WOA01.sal.annual.cdf',&
!      woa2001smonth ='/work2/knutali/WOA2001/WOA01.sal.monthly.cdf',&
!!      woa2001tannual='/work2/knutali/WOA2001/WOA01.tem.annual.cdf', &
!      woa2001tmonth ='/work2/knutali/WOA2001/WOA01.tem.monthly.cdf'
!   character(len=*), parameter ::  &
!      woa2005sannual='/work2/knutali/WOA2005/s00an1.nc',&
!      woa2005tmonth ='/work2/knutali/WOA2005/t0112an1.nc',&
!      woa2005smonth ='/work2/knutali/WOA2005/s0112an1.nc',&
!      woa2005tannual='/work2/knutali/WOA2005/t00an1.nc'

   character(len=*), parameter ::  &
      woa2001sannual='/work/shared/nersc/msc/WOA2001/WOA01.sal.annual.cdf',&
      woa2001smonth ='/work/shared/nersc/msc/WOA2001/WOA01.sal.monthly.cdf',&
      woa2001tannual='/work/shared/nersc/msc/WOA2001/WOA01.tem.annual.cdf', &
      woa2001tmonth ='/work/shared/nersc/msc/WOA2001/WOA01.tem.monthly.cdf'
   character(len=*), parameter ::  &
      woa2005sannual='/work/shared/nersc/msc/WOA2005/s00an1.nc',&
      woa2005tmonth ='/work/shared/nersc/msc/WOA2005/t0112an1.nc',&
      woa2005smonth ='/work/shared/nersc/msc/WOA2005/s0112an1.nc',&
      woa2005tannual='/work/shared/nersc/msc/WOA2005/t00an1.nc'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   public :: levitus_setup, levitus_interp

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


      if (trim(cwoa)=='WOA2001') then
         woa2001=.true.
         woa2005=.false.
         filenameannual=woa2001tannual
         filenamemonth =woa2001tmonth
      elseif (trim(cwoa)=='WOA2005') then
         woa2005=.true.
         woa2001=.false.
         filenameannual=woa2005tannual
         filenamemonth =woa2005tmonth
      else
         print *,'cwoa is ',cwoa
         stop '(levitus_setup -- invalid woa value)'
      end if


      print *,'reading Levitus fields from '//trim(filenamemonth)
      ncresult = nf90_open(trim(filenamemonth),nf90_nowrite,ncidmonth) 
      if (ncresult/= nf90_noerr) then
         call ncerr(ncresult)
      end if

      print *,'reading deep Levitus fields from '//trim(filenameannual)
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


   subroutine levitus_interp(varname, &
                             csjuld,cslon,cslat,cspres,clevvar,nprof,nzlevel,undef)
   use netcdf
   use mod_year_info
   implicit none
   character(len=*), intent(in) :: varname
   character(len=80) :: filenamemonth, filenameannual
   character(len=20) :: varannual, varmonth
   integer, intent(in) :: nprof,nzlevel
   real, intent(in), dimension(nprof) :: csjuld, cslon, cslat
   real, intent(in), dimension(nzlevel,nprof) :: cspres
   real, pointer,   dimension(:,:) :: clevvar
   real, intent(in) :: undef

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

   if (woa2001) then
      if (varname=='salinity') then
         filenameannual=woa2001sannual
         filenamemonth =woa2001smonth
         varannual='salinity'
         varmonth ='salinity'
      elseif (varname=='temperature') then
         filenameannual=woa2001tannual
         filenamemonth =woa2001tmonth
         varannual='temperature'
         varmonth ='temperature'
      else
         print *,'Unknown field '//trim(varname)
      end if
   elseif (woa2005) then
      if (varname=='salinity') then
         filenameannual=woa2005sannual
         filenamemonth =woa2005smonth
         varannual='s00an1'
         varmonth ='s0112an1'
      elseif (varname=='temperature') then
         filenameannual=woa2005tannual
         filenamemonth =woa2005tmonth
         varannual='t00an1'
         varmonth ='t0112an1'
      else
         print *,'Unknown field '//trim(varname)
      end if
   else
      print *,'clim flag not set '
      stop '(levitus_interp)'
   end if

   ! Get month,   day and year
   do i=1,nprof
      call juliantodate(floor(csjuld(i)),iyear,imonth,iday,1950,1,1)
      rmonth=(imonth-1)+min(iday/30.,1.) ! Close enough
      !print *,iyear, imonth, iday,rmonth
   end do

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


   ! Profile positions on lev grid + bilin coeff
   allocate(ipiv(nprof))
   allocate(jpiv(nprof))
   allocate(a   (nprof,4))
   do k=1,nprof

      lon1=cslon(k)
      lat1=cslat(k)

      if(lon1<0) lon1=lon1+360.

      ! index into lon/lat matrix 
      ripiv=(lon1-levlon(1))/(levlon(2)-levlon(1))+1.
      ipiv(k)=floor(ripiv);
      t=ripiv-ipiv(k)

      rjpiv=(lat1-levlat(1))/(levlat(2)-levlat(1))+1.
      jpiv(k)=floor(rjpiv);
      u=rjpiv-jpiv(k)

      a(k,1)=(1-t)*(1-u)
      a(k,2)=    t*(1-u)
      a(k,3)=    t*u
      a(k,4)=(1-t)*u     

      !print *,lon1,lat1,ripiv,rjpiv,t,u
   end do




   ! Get levitus data for profiles - on levitus vert grid
   allocate(levprof(zdimlenannual,nprof))
   allocate(tmp(londimlen,latdimlen))
   allocate(levfld(londimlen,latdimlen))

   print *,'reading Levitus fields from '//trim(filenamemonth)
   ncresult = nf90_open(trim(filenamemonth),nf90_nowrite,ncidmonth) 
   if (ncresult/= nf90_noerr) then
      call ncerr(ncresult)
   end if
   call getvarstat(ncidmonth,trim(varmonth),varidmonth,vndim  ,vdims)
   ncresult=nf90_get_att(ncidmonth, varidmonth    , 'missing_value', fillvalue_var)

   print *,'reading deep Levitus fields from '//trim(filenameannual)
   ncresult = nf90_open(trim(filenameannual),nf90_nowrite,ncidannual) 
   if (ncresult/= nf90_noerr) then
      call ncerr(ncresult)
   end if
   call getvarstat(ncidannual,trim(varannual),varidannual,vndim  ,vdims)

   do k=1,zdimlenannual

      !print *,k,zdimlenmonth,zdimlenannual

      ! time interpolated saln for this level
      if (k<=zdimlenmonth) then
         call ncerr( nf90_get_var(ncidmonth,varidmonth, tmp    , start=(/1,1,k,lwind/)))
         levfld=tmp*wlw
      else
         call ncerr( nf90_get_var(ncidannual,varidannual, tmp    , start=(/1,1,k/)))
         levfld=tmp*wlw
      end if

      if (k<=zdimlenmonth) then
         call ncerr( nf90_get_var(ncidmonth,varidmonth, tmp    , start=(/1,1,k,upind/)))
         levfld=tmp*wup+levfld
      else
         call ncerr( nf90_get_var(ncidannual,varidannual, tmp    , start=(/1,1,k/)))
         levfld=tmp*wup+levfld
      end if

      do k2=1,nprof
         
         ! Horizontal interpolation
         ip=ipiv(k2)
         jp=jpiv(k2)
         ipb=mod(ip,size(levlon,1))+1
         jp =min(max(1,jp  ),size(levlat,1))
         jpb=min(max(1,jp+1),size(levlat,1))

         if (abs(levfld(ip  ,jp  )-fillvalue_var) <1e-4) a(k2,1)=0.
         if (abs(levfld(ip+1,jp  )-fillvalue_var) <1e-4) a(k2,2)=0.
         if (abs(levfld(ip+1,jp+1)-fillvalue_var) <1e-4) a(k2,3)=0.
         if (abs(levfld(ip  ,jp+1)-fillvalue_var) <1e-4) a(k2,4)=0.

         if (sum(a(k2,:)) > .1) then
            a(k2,:)= a(k2,:)/sum(a(k2,:))

            levprof(k,k2) = levfld(ip  ,jp  )*a(k2,1) + &
                            levfld(ip+1,jp  )*a(k2,2) + &
                            levfld(ip+1,jp+1)*a(k2,3) + &
                            levfld(ip  ,jp+1)*a(k2,4)   
         else
            levprof(k,k2) = undef
         end if

      end do
   end do
   ncresult = nf90_close(ncidannual)
   ncresult = nf90_close(ncidmonth)


   ! Final step - interpolate onto coriolis vertical levels
   allocate(clevvar(nzlevel,nprof))
   do i=1,nprof
      do k=1,nzlevel
         if (cspres(k,i)/=undef) then

            ! find levitus depth index right above level cspres
            k2=1
            do while (cspres(k,i)>=levz(k2) .and.k2<=nzlevel)
               k2=k2+1
            end do
            k2=k2-1

            upind=k2
            lwind=min(k2+1,zdimlenannual)

            if (upind/=lwind) then
               wup=(levz(lwind)-cspres(k,i))/(levz(lwind)-levz(upind))
               wlw=1.-wup
               !print *,levz(upind),cspres(k,i),levz(lwind),wup,wlw
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
