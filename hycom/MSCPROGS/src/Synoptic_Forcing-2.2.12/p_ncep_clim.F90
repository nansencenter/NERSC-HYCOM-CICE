program ncep_clim
   ! Program to create a climatology from the NCEP files
   use mod_year_info22
   use mod_atm_func
   use m_ncvar_dims
   use m_ncvar_read
   use netcdf
   implicit none

   ! Field names
   character(len=200), save :: ncep_path='./NCEP/'

   ! Variable names
   character(len=*), parameter :: vicec='icec'
   character(len=*), parameter :: vprec='prate'
   character(len=*), parameter :: vtair='air'
   character(len=*), parameter :: vtcdc='tcdc'
   character(len=*), parameter :: vuwnd='uwnd'
   character(len=*), parameter :: vvwnd='vwnd'
   character(len=*), parameter :: vlmsk='land'
   character(len=*), parameter :: vhgt ='hgt'
   character(len=*), parameter :: vpres='pres'
   character(len=*), parameter :: vshum='shum'
   character(len=*), parameter :: vssrd='dswrf'

   ! Variable names of derived quantities
   character(len=*), parameter :: vwspd='windspeed'
   character(len=*), parameter :: vtaux='tauewd'
   character(len=*), parameter :: vtauy='taunwd'
   character(len=*), parameter :: vrhum='relhum'
   character(len=*), parameter :: vvpmx='vapmix'

   ! vapor pressure parameters
   real, parameter :: & 
      aice     =9.5  ,&    ! --
      bice     =7.66 ,&    ! k                    ..
      awater   =7.5  ,&    ! --                   ..
      bwater   =35.86,&    ! k                    ..
      cd       = 0.0012, &
      airdns   =    1.2


   character(len=4) :: clyy, cyy

   ! file names
   character(len=200) :: cenv
   character(len=80) ::  &
      ficec,fprec,ftair,ftcdc,fuwnd, fvwnd, flmsk, fpres,  &
      fshum, fhgt, fssrd, fssrd2 

   real, dimension(:,:,:), allocatable :: tmplon, tmplat
   real, allocatable, dimension(:,:,:) ::  &
      ncep_lmsk  , ncep_hgt   , ncep_shum  , ncep_dswrf , &
      nceptmp    , ncep_icec  , ncep_precip, ncep_tair  , &
      ncep_ccov  , ncep_uwind , ncep_vwind , ncep_wspd  , &
      ncep_taux  , ncep_tauy  , ncep_pres  , ncep_rhum  , &
      ncep_vpmx
   real, allocatable, dimension(:,:,:) ::  &
      ncep_clim_shum  , ncep_clim_dswrf , &
      ncep_clim_icec  , ncep_clim_precip, ncep_clim_tair  , &
      ncep_clim_ccov  , ncep_clim_uwind , ncep_clim_vwind , ncep_clim_wspd  , &
      ncep_clim_taux  , ncep_clim_tauy  , ncep_clim_pres  , ncep_clim_rhum  , &
      ncep_clim_vpmx

   integer :: ncep_clim_cnt(12)

   real, allocatable :: lon(:),lat(:),tmpdswrf(:,:)

   real :: rfac,wndfac,cd_new,w4,wfact,dlon,dlat
   real :: satvapp, ax, bx,sathumid
   integer :: nlon,nlat,ndims,recdim 
   integer , dimension(NF90_MAX_VAR_DIMS) :: dimsizes
   integer  :: i,j,k,ix,jy,irec, maxrec, maxrec2
   real     :: fday
   logical :: allocerr, ldump, ldump2
   logical, save :: lfirst=.true.
   integer :: iyear, imonth, iday, iyear2
   integer, parameter :: lp=6
   integer*4, external :: system
   integer*4 :: ret


   ! Check NCEPR path on first pass
   if (lfirst) then
      call getenv('NCEP_PATH',cenv)
      if (trim(cenv)/='') then ! prefer this path if present
         ncep_path=trim(cenv)//'/'
      end if
      ret=system('[ -d '//trim(ncep_path)//' ]')
      if (ret /=0 ) then
         print *
         print *,'The directory  '//trim(ncep_path)//' does not exist.'
         print *,'Make sure that '
         print *,' - You have linked the NCEP data into this catalogue, and '
         print *,'   that the variable NCEP_PATH is empty'
         print *,' - OR, the variable NCEP_PATH is set to the location of '
         print *,'   the NCEP data'
         call exit(1)
      end if
   end if

   flmsk=trim(ncep_path)//'land.sfc.gauss.nc'
   fhgt =trim(ncep_path)//'hgt.sfc.gauss.nc' ! New - on gauss grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Get longitude  & latitude for NCEP data -- Two different grids 
   call ncvar_dims(flmsk,'lat',dimsizes,ndims,recdim)
   nlat=dimsizes(1) ; allocate(lat(nlat)) ; allocate(tmplat(nlat,1,1))
   call ncvar_read(flmsk,'lat', tmplat, nlat,1,1,1,1)
   lat=tmplat(:,1,1)
   deallocate(tmplat)

   call ncvar_dims(flmsk,'lon',dimsizes,ndims,recdim)
   nlon=dimsizes(1) ; allocate(lon(nlon)) ; allocate(tmplon(nlon,1,1))
   call ncvar_read(flmsk,'lon',tmplon,nlon,1,1,1,1)
   lon=tmplon(:,1,1)
   deallocate(tmplon)

   dlon=lon(2)-lon(1)
   do i=2,nlon
      if ( lon(i)-lon(i-1) /= dlon) then
         write(lp,'(a)') 'Longitude is not uniform !!'
         stop '(forfun_ncep)'
      end if
   end do

   ! Allocate snapshot fields
   allocate(ncep_lmsk  (nlon,nlat,1))
   allocate(ncep_hgt   (nlon,nlat,1)) 
   allocate(nceptmp    (nlon,nlat,1))
   allocate(ncep_icec  (nlon,nlat,1))
   allocate(ncep_precip(nlon,nlat,1))
   allocate(ncep_tair  (nlon,nlat,1))
   allocate(ncep_ccov  (nlon,nlat,1))
   allocate(ncep_uwind (nlon,nlat,1))
   allocate(ncep_vwind (nlon,nlat,1))
   allocate(ncep_wspd  (nlon,nlat,1))
   allocate(ncep_taux  (nlon,nlat,1))
   allocate(ncep_tauy  (nlon,nlat,1))
   allocate(ncep_pres  (nlon,nlat,1))
   allocate(ncep_shum  (nlon,nlat,1))
   allocate(ncep_rhum  (nlon,nlat,1))
   allocate(ncep_vpmx  (nlon,nlat,1))
   allocate(ncep_dswrf (nlon,nlat,1))

   ! Allocate climatology fields
   allocate(ncep_clim_icec  (nlon,nlat,12))
   allocate(ncep_clim_precip(nlon,nlat,12))
   allocate(ncep_clim_tair  (nlon,nlat,12))
   allocate(ncep_clim_ccov  (nlon,nlat,12))
   allocate(ncep_clim_uwind (nlon,nlat,12))
   allocate(ncep_clim_vwind (nlon,nlat,12))
   allocate(ncep_clim_wspd  (nlon,nlat,12))
   allocate(ncep_clim_taux  (nlon,nlat,12))
   allocate(ncep_clim_tauy  (nlon,nlat,12))
   allocate(ncep_clim_pres  (nlon,nlat,12))
   allocate(ncep_clim_shum  (nlon,nlat,12))
   allocate(ncep_clim_rhum  (nlon,nlat,12))
   allocate(ncep_clim_vpmx  (nlon,nlat,12))
   allocate(ncep_clim_dswrf (nlon,nlat,12))

   ncep_clim_icec  =0.
   ncep_clim_precip=0.
   ncep_clim_tair  =0.
   ncep_clim_ccov  =0.
   ncep_clim_uwind =0.
   ncep_clim_vwind =0.
   ncep_clim_wspd  =0.
   ncep_clim_taux  =0.
   ncep_clim_tauy  =0.
   ncep_clim_pres  =0.
   ncep_clim_shum  =0.
   ncep_clim_rhum  =0.
   ncep_clim_vpmx  =0.
   ncep_clim_dswrf =0.
   ncep_clim_cnt =0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP land mask
   write(lp, *) 'Using NCEP land mask'
   call ncvar_read(flmsk,vlmsk,ncep_lmsk   ,nlon,nlat, 1,1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP geopot height
   call ncvar_read(fhgt,vhgt,ncep_hgt   ,nlon,nlat, 1,1,1)
   write(lp, *) 'Using NCEP geopot height (orography pressure correction'

   ! Loop over files in NCEP directory
   !do iyear=1949,1949 ! test
   do iyear=1949,2007 

      print *,'New year',iyear

      ! Which files to open
      write(cyy,'(i4.4)') iyear
      write(clyy,'(i4.4)') iyear - 1
      ficec=trim(ncep_path)//'icec.sfc.gauss.'//cyy//'.nc'
      fprec=trim(ncep_path)//'prate.sfc.gauss.'//cyy//'.nc'
      ftair=trim(ncep_path)//'air.2m.gauss.'//cyy//'.nc'
      ftcdc=trim(ncep_path)//'tcdc.eatm.gauss.'//cyy//'.nc'
      fuwnd=trim(ncep_path)//'uwnd.10m.gauss.'//cyy//'.nc'
      fvwnd=trim(ncep_path)//'vwnd.10m.gauss.'//cyy//'.nc'
      fpres=trim(ncep_path)//'pres.sfc.gauss.'//cyy//'.nc'
      fshum=trim(ncep_path)//'shum.2m.gauss.'//cyy//'.nc'
      fssrd=trim(ncep_path)//'dswrf.sfc.gauss.'//cyy//'.nc'
      fssrd2=trim(ncep_path)//'dswrf.sfc.gauss.'//clyy//'.nc'

      ! Get number of records in file
      call ncvar_dims(ficec,vicec,dimsizes,ndims,recdim)
      maxrec=dimsizes(recdim)

      lfirst=.true.

      ! Now loop over file records
      do irec=1,maxrec,50 ! test
      !do irec=1,maxrec

         ! Time relative to 1-1-iyear 00:00:00
         fday = (irec-1)*0.25 ! 6hr time step -> floating point day

         ! Actual date
         call juliantodate(floor(fday),iyear2,imonth,iday,iyear,1,1)

         ldump=lfirst
         ldump2=lfirst.or.ldump.or.mod(irec,80)==0

         if (ldump) then
            write(lp,'(a)')'####################################################'
            write(lp,'(a)')'Reading NCEP synoptic forcing fields '
         end if
         if (ldump2) then
            write(lp,'(a,i5,a,i4,a,i3,a,i2,a,i2,a,i3)')'Record=',irec,' Year=',iyear2, &
               ' Day=',floor(fday),' hour=',nint(24*(fday-floor(fday))), &
               ' month=',imonth,' dayinmonth=',iday
         end if



         !! NCEP ice concentration
         if (ldump) write(lp, *) 'Using NCEP ice conc (involved in relhum)'
         call ncvar_read(ficec,vicec,ncep_icec,nlon,nlat,1,irec,irec)
         ncep_icec=max(0.,min(1.,ncep_icec))

         ! NCEP precipitation   kg/(m^2 s) -> m/month (rho=1000)
         ! convert prcipitation  from  kg/(m2*s) to m^3 m^-2 s^-1
         if (ldump) write(lp, *) 'Using NCEP synoptic precipitation '
         call ncvar_read(fprec,vprec,ncep_precip,nlon,nlat,1,irec,irec)
         ncep_precip=max(0.,ncep_precip*1e-3)

         ! NCEP air temperature   K 
         if (ldump) write(lp,'(a)')' Using NCEP synoptic temperature'
         call ncvar_read(ftair,vtair,ncep_tair,nlon,nlat,1,irec,irec)

         ! Correct temperature for dry adiabatic lapse rate 
         ncep_tair=ncep_tair + ncep_hgt * 0.0098

         if (ldump) then
            print *,'max/min geopot hgt adj (tair) :',maxval(ncep_hgt)*0.0098, &
                                              minval(ncep_hgt)*0.0098
            write(lp,'(a)') ' Using NCEP synoptic specific humidity'
            call flush(lp)
         end if

         ! NCEP clouds   % -> []  (0-1)
         if (ldump) write(lp,'(a)') ' Using NCEP synoptic clouds '
         call ncvar_read(ftcdc,vtcdc,ncep_ccov,nlon,nlat,1,irec,irec)
         ncep_ccov=0.01*ncep_ccov
         ncep_ccov=min(1.0,max(ncep_ccov,0.))

         ! NCEP winds
         if (ldump) write(lp,'(a)') ' Using NCEP synoptic winds '
         call ncvar_read(fuwnd,vuwnd,ncep_uwind,nlon,nlat,1,irec,irec)
         call ncvar_read(fvwnd,vvwnd,ncep_vwind,nlon,nlat,1,irec,irec)

         ! Derived  wind stress and absolute windspeed
         !$OMP PARALLEL DO PRIVATE(j,i,wndfac,w4,wfact,cd_new,im1,jm1,ip1,jp1) &
         !$OMP SCHEDULE(STATIC,jblk)
         do j=1,nlat
         do i=1,nlon

            ncep_wspd(i,j,1)=sqrt(ncep_uwind(i,j,1)*ncep_uwind(i,j,1)+ &
                                  ncep_vwind(i,j,1)*ncep_vwind(i,j,1))

            wndfac=(1.+sign(1.,ncep_wspd(i,j,1)-11.))*.5
            cd_new=(0.49+0.065*ncep_wspd(i,j,1))*1.0e-3*wndfac+cd*(1.-wndfac)

            wfact=ncep_wspd(i,j,1)*airdns*cd_new
            ncep_taux(i,j,1)=ncep_uwind(i,j,1)*wfact
            ncep_tauy(i,j,1)=ncep_vwind(i,j,1)*wfact
         end do
         end do

         !NCEP pressure
         if (ldump) write(lp,'(a)') ' Using NCEP synoptic SLP '
         call ncvar_read(fpres,vpres,ncep_pres,nlon,nlat,1,irec,irec)
         ncep_pres=0.01*ncep_pres                  ! Convert from Pa to mbar 
         ncep_pres=ncep_pres+ncep_hgt*9.8*1.3*0.01 ! Simple correction for geopot. height

         if (ldump) then
            print *,'max/min geopot hgt adj (pressure) :',maxval(ncep_hgt)*9.8*1.3*0.01, &
                                              minval(ncep_hgt)*9.8*1.3*0.01
            write(lp,'(a)') ' Using NCEP synoptic specific humidity'
            call flush(lp)
         end if

         ! Correct shum with hgt too ? (For saturation)
         call ncvar_read(fshum,vshum,ncep_shum,nlon,nlat,1,irec,irec)
         do jy=1,nlat
         do ix=1,nlon

            ! Saturation Vapour pressure (sea level)
            satvapp=satvap(ncep_tair(ix,jy,1))

            ! Saturation specific humidity at sea level pressure
            sathumid=humid(ncep_pres(ix,jy,1)*100.,satvapp) ! ncep_pres in mBar

            ! Correction of shum for saturation 
            ncep_shum(ix,jy,1)  =min(ncep_shum(ix,jy,1),sathumid)

            ! Vapor mixing ratio
            ncep_vpmx(ix,jy,1)=ncep_shum(ix,jy,1)/(1.-ncep_shum(ix,jy,1))

            ! Relative humidty
            ncep_rhum(ix,jy,1)= spechum_to_relhum(ncep_shum(ix,jy,1),sathumid)
            ncep_rhum(ix,jy,1)= max(0.,min(1., ncep_rhum(ix,jy,1)))

         end do
         end do


         !NCEP downwelling shortwave radiation
         if (ldump) then
            write(lp,'(a)') ' Using NCEP synoptic DSWRF '
            print *,'NB: DSWRF not testet thoroughly'
            write(lp,'(a)')'####################################################'
         end if
         call ncvar_read(fssrd,vssrd,ncep_dswrf,nlon,nlat,1,irec,irec)

         ! ssrd is averaged over 6 next hrs, read last ssrd field as well
         if (irec==1)  then
            call ncvar_dims(fssrd2,vssrd,dimsizes,ndims,recdim)
            maxrec2=dimsizes(recdim)
            call ncvar_read(fssrd2,vssrd,nceptmp,nlon,nlat,1,maxrec2,maxrec2)
         else 
            call ncvar_read(fssrd ,vssrd,nceptmp,nlon,nlat,1,irec-1,irec-1)
         end if
         ncep_dswrf=(nceptmp+ncep_dswrf)*0.5
         ncep_dswrf=ncep_dswrf*.5



         ! Add to climatology fields
         ncep_clim_icec  (:,:,imonth) = ncep_clim_icec  (:,:,imonth) + ncep_icec  (:,:,1)
         ncep_clim_precip(:,:,imonth) = ncep_clim_precip(:,:,imonth) + ncep_precip(:,:,1)
         ncep_clim_tair  (:,:,imonth) = ncep_clim_tair  (:,:,imonth) + ncep_tair  (:,:,1)
         ncep_clim_ccov  (:,:,imonth) = ncep_clim_ccov  (:,:,imonth) + ncep_ccov  (:,:,1)
         ncep_clim_uwind (:,:,imonth) = ncep_clim_uwind (:,:,imonth) + ncep_uwind (:,:,1)
         ncep_clim_vwind (:,:,imonth) = ncep_clim_vwind (:,:,imonth) + ncep_vwind (:,:,1)
         ncep_clim_wspd  (:,:,imonth) = ncep_clim_wspd  (:,:,imonth) + ncep_wspd  (:,:,1)
         ncep_clim_taux  (:,:,imonth) = ncep_clim_taux  (:,:,imonth) + ncep_taux  (:,:,1)
         ncep_clim_tauy  (:,:,imonth) = ncep_clim_tauy  (:,:,imonth) + ncep_tauy  (:,:,1)
         ncep_clim_pres  (:,:,imonth) = ncep_clim_pres  (:,:,imonth) + ncep_pres  (:,:,1)
         ncep_clim_shum  (:,:,imonth) = ncep_clim_shum  (:,:,imonth) + ncep_shum  (:,:,1)
         ncep_clim_rhum  (:,:,imonth) = ncep_clim_rhum  (:,:,imonth) + ncep_rhum  (:,:,1)
         ncep_clim_vpmx  (:,:,imonth) = ncep_clim_vpmx  (:,:,imonth) + ncep_vpmx  (:,:,1)
         ncep_clim_dswrf (:,:,imonth) = ncep_clim_dswrf (:,:,imonth) + ncep_dswrf (:,:,1)
         ncep_clim_cnt(imonth)= ncep_clim_cnt(imonth) + 1

         lfirst=.false.

      end do ! Loop over file records
      print *
      print *
   end do ! Loop over files

   do imonth=1,12

      ! Add to climatology fields
      ncep_clim_icec  (:,:,imonth) = ncep_clim_icec  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_precip(:,:,imonth) = ncep_clim_precip(:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_tair  (:,:,imonth) = ncep_clim_tair  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_ccov  (:,:,imonth) = ncep_clim_ccov  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_uwind (:,:,imonth) = ncep_clim_uwind (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_vwind (:,:,imonth) = ncep_clim_vwind (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_wspd  (:,:,imonth) = ncep_clim_wspd  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_taux  (:,:,imonth) = ncep_clim_taux  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_tauy  (:,:,imonth) = ncep_clim_tauy  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_pres  (:,:,imonth) = ncep_clim_pres  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_shum  (:,:,imonth) = ncep_clim_shum  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_rhum  (:,:,imonth) = ncep_clim_rhum  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_vpmx  (:,:,imonth) = ncep_clim_vpmx  (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_dswrf (:,:,imonth) = ncep_clim_dswrf (:,:,imonth) / ncep_clim_cnt(imonth)
      ncep_clim_cnt(imonth)= ncep_clim_cnt(imonth) 

   end do


   ! Create climatology files
   print *,minval(ncep_clim_icec),maxval(ncep_clim_icec)
   call clmtonc('ncep_clim_icec.nc','icec',ncep_clim_icec,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_prcp.nc','prcp',ncep_clim_precip,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_tair.nc','tair',ncep_clim_tair ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_ccov.nc','ccov',ncep_clim_ccov ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_uwnd.nc','uwnd',ncep_clim_uwind,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_vwnd.nc','vwnd',ncep_clim_vwind,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_wspd.nc','wspd',ncep_clim_wspd ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_taux.nc','taux',ncep_clim_taux ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_tauy.nc','tauy',ncep_clim_tauy ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_pres.nc','pres',ncep_clim_pres ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_shum.nc','shum',ncep_clim_shum ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_rhum.nc','rhum',ncep_clim_rhum ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_vpmx.nc','vpmx',ncep_clim_vpmx ,lon,lat,nlon,nlat,12)
   call clmtonc('ncep_clim_dswr.nc','dswr',ncep_clim_dswrf,lon,lat,nlon,nlat,12)

   contains

      subroutine clmtonc(fname,vname,fld,lon,lat,nlon,nlat,nrec)
      use m_handle_err
      implicit none
      integer, intent(in) :: nlon, nlat, nrec
      real, intent(in) :: fld(nlon,nlat,nrec) 
      real, intent(in) :: lon(nlon)
      real, intent(in) :: lat(nlat)
      character(len=*), intent(in) :: fname, vname

      integer :: xdimid, ydimid, tdimid, ncid, varid, varidlon, varidlat

      if (NF90_CREATE(fname,NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(norsexclim)'
      end if

      call handle_err(NF90_DEF_DIM(ncid,'longitude',nlon,xdimid))
      call handle_err(NF90_DEF_DIM(ncid,'latitude' ,nlat,ydimid))
      call handle_err(NF90_DEF_DIM(ncid,'nrec'     ,nrec,tdimid))
      call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_Float, &
         (/xdimid/),varidlon))
      call handle_err(NF90_DEF_VAR(ncid,'latitude' ,NF90_Float, &
         (/ydimid/),varidlat))
      call handle_err(NF90_DEF_VAR(ncid,vname,NF90_Float, &
         (/xdimid,ydimid,tdimid/),varid))
      call handle_err(NF90_PUT_ATT(ncid,varid,'comment','NCEP Climatology'))
      call handle_err(NF90_enddef(ncid))
      call handle_err(NF90_put_var(ncid,varidlon,lon))
      call handle_err(NF90_put_var(ncid,varidlat,lat))
      call handle_err(NF90_put_var(ncid,varid   ,fld))
      call handle_err(NF90_close(ncid))
      end subroutine




!      real function vapp(ax,bx,tx)
!      implicit none
!        real, intent(in) :: ax,bx,tx
!         vapp=611.*10.**(ax*(tx-273.16)/(tx-bx))
!      end function vapp
!
!      real function humid(px,ex)
!      implicit none
!        real, intent(in) :: px,ex
!         humid=.622*ex/(px-.378*ex)
!      end function humid

   END program ncep_clim
