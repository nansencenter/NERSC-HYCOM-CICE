      module mod_forcing_nersc
      use mod_xc
c --- TODO: taux/tauy/uwind,vwind should perhaps be on p-grid
      implicit none
c
      character(len=5), save :: rforce   ! Forcing (ecmwf/ncep ..)
      character(len=5), save :: clmflag  ! Climate forcing to use (era40 or "old")
c
c --- Forcing parameters -- particular to nersc 
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::   
     &  uwind  ,    ! Winds            [m/s]
     &  vwind  ,    ! Winds            [m/s]
     &  clouds ,    ! Cloud Cover      [0-1]
     &  slp    ,    ! Cloud Cover      [0-1]
     &  relhum      ! Relative humidity[0-1]
c
c --- Units to connect to for reading fields above
      integer, parameter :: 
     &   unit_uwind  =895, 
     &   unit_vwind  =896, 
     &   unit_clouds =897, 
     &   unit_relhum =898, 
     &   unit_slp    =899   
c
C --- Time info 
      real*8,save ::  
     & dtime_uwind   ,
     & dtime_vwind   ,
     & dtime_relhum  ,
     & dtime_clouds  ,
     & dtime_slp
c
C --- Parameters
      real,parameter :: cd=0.0012
      real,parameter :: airdns  =    1.2
      real,parameter :: slp0=1012.5
c
c
C --- Relaxation time scales 
      real, save :: time_srelax
      real, save :: time_trelax
c
c
c --- Climatology index used in hf forcing
      real   , save,private :: hf_cwgt
      real   , save,private :: hf_wc0,hf_wc1
      integer, save,private :: hf_lc0,hf_lc1
      integer, save,private :: hf_mc0,hf_mc1, hf_lastmc1
      logical, save,private :: hf_cinit, hf_cread
c
c --- Last time a forcing field was read
      real*8 ,save :: dtimefrc_last
c
c --- time step of hf forcing
      real*8 ,save :: hf_rdtime
c
c --- Lon and lat parameters for uniform grids
c --- Kept for active hf forcing
      real, save,private    :: flon, flat
      real, save, private   :: dlon, dlat
c
c --- Keeps gaussian latitudes for gaussian grids
      real, save, private, allocatable   :: gausslat(:)
c
      private :: rdncrec, rdncrec_1
#if defined (FORCING_INLINE)
     &  , humid
#endif
      contains
c
#if defined (FORCING_INLINE)
c --- Reads high-frequency forcing fields from various input files.
c --- Replacement for forfunh in forfun.F
c --- TODO: Add various options for seatmp, surtmp etc
      subroutine nersc_forcing_inline(dtime,fmonth)
      implicit none
      real*8,  intent(in) :: dtime,fmonth
      real*8 :: dtimefrc,dtimefrc2
      integer :: i,j
      logical :: init
      logical, parameter :: force_diag=.true.
c
      include 'common_blocks.h'
c
      init=.false.
      if (w0<-1.) then
c ---   linear interpolation in time, so slots 3 and 4 are zero.
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
              taux(i,j,3) = 0.0
              taux(i,j,4) = 0.0
              tauy(i,j,3) = 0.0
              tauy(i,j,4) = 0.0
            wndspd(i,j,3) = 0.0
            wndspd(i,j,4) = 0.0
            airtmp(i,j,3) = 0.0
            airtmp(i,j,4) = 0.0
            vapmix(i,j,3) = 0.0
            vapmix(i,j,4) = 0.0
            precip(i,j,3) = 0.0
            precip(i,j,4) = 0.0
            radflx(i,j,3) = 0.0
            radflx(i,j,4) = 0.0
             swflx(i,j,3) = 0.0
             swflx(i,j,4) = 0.0
            surtmp(i,j,3) = 0.0
            surtmp(i,j,4) = 0.0
            seatmp(i,j,3) = 0.0
            seatmp(i,j,4) = 0.0
#if defined (NERSC_VERSION)
            uwind (i,j,3) = 0.0
            uwind (i,j,4) = 0.0
            vwind (i,j,3) = 0.0
            vwind (i,j,4) = 0.0
            slp   (i,j,3) = 0.0
            slp   (i,j,4) = 0.0
            clouds(i,j,3) = 0.0
            clouds(i,j,4) = 0.0
            relhum(i,j,3) = 0.0
            relhum(i,j,4) = 0.0
#endif /* NERSC_VERSION */
          enddo
        enddo
        init=.true.
      end if
c
c --- Set weights for inear interpolation of climatology
      call set_hfclmwgt(fmonth,init)
c
      hf_rdtime=6.d0/24.d0
      dtimefrc =floor(dtime/hf_rdtime)*hf_rdtime
      dtimefrc2=dtimefrc+hf_rdtime
      if (mnproc==1) 
     & print '(a,4f16.4)','dtimefrc, dtime, dtimefrc2,dtimefrclast',
     &         dtimefrc,dtime,dtimefrc2,dtimefrc_last

      if (trim(rforce)=='era40') then
         if (w0<-1) then
            call read_era40(dtimefrc ,1)
            call read_era40(dtimefrc2,2)
         else if (dtimefrc_last < dtimefrc2) then
            call shifthf()
            call read_era40(dtimefrc2,2)
         end if
       elseif (trim(rforce)=='era-i') then
         if (w0<-1) then
            call read_erai(dtimefrc ,1)
            call read_erai(dtimefrc2,2)
         else if (dtimefrc_last < dtimefrc2) then
            call shifthf()
            call read_erai(dtimefrc2,2)
         end if
      elseif (trim(rforce)=='ecnc') then
         if (w0<-1) then
            call read_ecmwf_nc(dtimefrc ,1)
            call read_ecmwf_nc(dtimefrc2,2)
         else if (dtimefrc_last < dtimefrc2) then
            call shifthf()
            call read_ecmwf_nc(dtimefrc2,2)
         end if
      elseif (trim(rforce)=='ecnc2') then
         if (w0<-1) then
            call read_ecmwf_nc2(dtimefrc ,1)
            call read_ecmwf_nc2(dtimefrc2,2)
         else if (dtimefrc_last < dtimefrc2) then
            call shifthf()
            call read_ecmwf_nc2(dtimefrc2,2)
         end if
      elseif (trim(rforce)=='ecmwf') then
         if (w0<-1) then
            call read_ecmwf(dtimefrc ,1)
            call read_ecmwf(dtimefrc2,2)
         else if (dtimefrc_last < dtimefrc2) then
            call shifthf()
            call read_ecmwf(dtimefrc2,2)
         end if
      elseif (trim(rforce)=='metno') then
         if (w0<-1) then
            call read_ecmwf(dtimefrc ,1)
            call read_ecmwf(dtimefrc2,2)
         else if (dtimefrc_last < dtimefrc2) then
            call shifthf()
            call read_ecmwf(dtimefrc2,2)
         end if
      elseif (trim(rforce)=='ncepr') then
         if (w0<-1) then
            call read_ncepr(dtimefrc ,1)
            call read_ncepr(dtimefrc2,2)
         else if (dtimefrc_last < dtimefrc2) then
            call shifthf()
            call read_ncepr(dtimefrc2,2)
         end if
      else
         write(lp,'(a)') 'So far only era40 supports the '//
     &                   'inline forcing option'
         call xcstop('(xcstop)')
      end if
c
c --- These options are not yet supported, and we should stop if they 
c --- are active
      if     (sstflg.eq.3) then
         if (mnproc==1) write(lp,'(a)') 'sstflg=3 not supported'//
     &      'for inline forcing'
         call xcstop('(mod_forcing_nersc:nersc_forcing_inline)')
      elseif (flxoff) then
         if (mnproc==1) write(lp,'(a)') 'flux offset not supported'//
     &      'for inline forcing'
         call xcstop('(mod_forcing_nersc:nersc_forcing_inline)')
      elseif (lwflag.eq.2 .or. sstflg.eq.2   .or.
     &         icmflg.eq.2 .or. ticegr.eq.0.0     ) then
         if (mnproc==1) write(lp,'(a)') 'lwflg=2, sstflg=2 iclmflg=2 '//
     &      'and ticegr=0.0 not supported for inline forcing'
         call xcstop('(mod_forcing_nersc:nersc_forcing_inline)')
      end if
c
c --- Diagnostics    
      if (force_diag) then
      if (w0<-1) then
         call inline_diag('hforig',dtimefrc ,1)
         call inline_diag('hforig',dtimefrc2,2)
      else if (dtimefrc_last < dtimefrc2) then
         call inline_diag('hforig',dtimefrc2,2)
      end if
      end if
c --- weights (linear)
      w0=(dtimefrc2 - dtime    )/(dtimefrc2-dtimefrc)
      w1=(dtime     - dtimefrc )/(dtimefrc2-dtimefrc)
      w2=0.
      w3=0.
c --- time last read
      dtimefrc_last=dtimefrc2
c
      end subroutine
c
      subroutine shifthf()
      implicit none
      integer :: i,j
      include 'common_blocks.h'
c --- Move old forcing to slot 1
      do j= 1-nbdy,jj+nbdy
      do i= 1-nbdy,ii+nbdy
         precip(i,j,1) = precip(i,j,2)
         airtmp(i,j,1) = airtmp(i,j,2)
         clouds(i,j,1) = clouds(i,j,2)
         uwind (i,j,1) = uwind (i,j,2)
         vwind (i,j,1) = vwind (i,j,2)
         taux  (i,j,1) = taux  (i,j,2)
         tauy  (i,j,1) = tauy  (i,j,2)
         wndspd(i,j,1) = wndspd(i,j,2)
         slp   (i,j,1) = slp   (i,j,2)
         relhum(i,j,1) = relhum(i,j,2)
         vapmix(i,j,1) = vapmix(i,j,2)
         swflx (i,j,1) = swflx (i,j,2)
      enddo
      enddo
      end subroutine
c
c --- Sets climate weights used for linear interpolation
c --- of climatology when running with hf forcing
      subroutine set_hfclmwgt(x,init)
      implicit none
      real, intent   (in) :: x
      logical, intent(in) :: init
      hf_cwgt=x
      hf_cinit=init

      hf_mc0=mod(11+floor(hf_cwgt),12)+1
      hf_mc1=mod(hf_mc0,12)+1
      !hf_wc0=ceiling(hf_cwgt)-hf_cwgt
      !hf_wc0=hf_mc0+1.d0-hf_cwgt
      hf_wc1=hf_cwgt-floor(hf_cwgt)
      hf_wc0=1.d0-hf_wc1
      hf_cread=.false.
c --- On init - read two cons. months into slot hf_lc0 and hf_lc1
      if (hf_cinit) then
         hf_lc0=1
         hf_lc1=2
         hf_lastmc1=hf_mc1
c --- On new month - read data into slothf_lc1
      else if (hf_mc1 /= hf_lastmc1 ) then
         if (hf_mc0 == hf_lastmc1) then
            hf_cread=.true.
            hf_lc0=hf_lc1
            hf_lc1=mod(hf_lc0,2)+1
            hf_lastmc1=hf_mc1
c ---    This shuld not happen
         else
           if (mnproc==1)  then
              print *,'hf_cwgt',hf_cwgt
              print *,'months ',hf_mc0,hf_mc1
              print *,'weights',hf_wc0,hf_wc1
              print *,'indices',hf_lc0,hf_lc1
              print *,'lastrd ',hf_lastmc1
              print *,'read ? ',hf_cread
            end if
            call xcstop('(set_hfclmwgt)')
            stop '(set_hfclmwgt)'
          end if
      end if
      if (mnproc==1)  then
        print *,'hf_cwgt',hf_cwgt
        print *,'months ',hf_mc0,hf_mc1
        print *,'weights',hf_wc0,hf_wc1
        print *,'indices',hf_lc0,hf_lc1
        print *,'lastrd ',hf_lastmc1
        print *,'read ? ',hf_cread
        print *,'floor  (hf_cwgt)',floor(hf_cwgt)
        print *,'ceiling(hf_cwgt)',ceiling(hf_cwgt)
      end if
      end subroutine
c
c
c
c --- Routine for reading era40 data from netcdf files, inline version
c --- It is assumed that era40 files are organized in years and variable
c --- names. Each netcdf file must have one record every 6 hour for 
c --- that year. These records must come in chronological order,
c --- starting Jan 1 at 00:00:00
c --- TODO : Now all threads read the same netcdf file, 
      subroutine read_era40(dtime,slot)
      use mod_xc
      use mod_common_ice
      use mod_hycom_nersc
      use netcdf
      implicit none
      real*8 ,  intent(in) :: dtime
      integer,  intent(in) :: slot
c
c --- Data catalogue name
      character(len=*), parameter :: era40_path='./ERA40/'
c
c --- Variable names - ERA40 1.125 by 1.125 
      character(len=*), parameter :: vuwnd='U10M_sfc'
      character(len=*), parameter :: vvwnd='V10M_sfc'
      character(len=*), parameter :: vtair='T2M_sfc'
      character(len=*), parameter :: vdewp='D2M_sfc'
      character(len=*), parameter :: vprec='TP'
      character(len=*), parameter :: vtcc ='TCC_sfc'
      character(len=*), parameter :: vblh ='BLH_sfc'
      character(len=*), parameter :: vlmsk='LSM_sfc'
      character(len=*), parameter :: vpres='MSL_sfc'
      character(len=*), parameter :: vssrd='SSRD_sfc'
c
c --- Hardcoded grid sizes 
      integer, parameter :: nlon=320, nlat=161
c
      integer  :: i,j,k,ix,jy,irec,maxrec, nlon2, nlat2
      real, dimension(nlon, nlat) :: 
     &   tmp, tmpu, tmpv, tmpw, tmp2, tmptx, tmpty
      real               :: e40_lon(nlon), e40_lat(nlat)
      real               :: rfac,wndfac,cd_new,w4,wfact,svpair, svpdew
      real               :: dlon,dlat,flon,flat,llon,llat,rtmp
      integer            :: dimid, varid, ncid
      integer            :: iyear, ihour, iday
      logical            :: lfirst=.true.
      character(len=80)  :: ncfil, ncfilu, ncfilv
      character(len=4)   :: cyy, cyyp1
      integer, parameter :: itype=1 ! Sets interpolation type bicubic
      !integer, parameter :: itype=0 ! Sets interpolation type bilinear
#ifdef KARA2002
      real               :: tmpw2,CD0,CD1,airdns2
#endif
c
      include 'common_blocks.h'
c
      call forday(dtime,yrflag,iyear,iday,ihour)
      irec  =  (iday-1)*4 + ihour/6 +1
      if (mnproc==1) 
     &  write(lp,'(a,i6,a,i4,x,i3,x,i2)')'ERA40 data from record: ',
     &  irec, '  time= ',iyear,iday,ihour
c
c --- Set file names
      write(cyy  ,'(i4.4)') iyear
      write(cyyp1,'(i4.4)') iyear+1
c
c --- Check dimensions against hardcoded ones
c --- Get longitude  & latitude for data  - use land mask file
      ncfil=era40_path//'ans.LSM.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_dimid(ncid, "lon", dimid))
      call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlon2))
      call ncerr(nf90_inq_dimid(ncid, "lat", dimid))
      call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlat2))
      call ncerr(nf90_inq_dimid(ncid, "record", dimid))
      call ncerr(nf90_inquire_dimension(ncid, dimid, len=maxrec))
      if (nlon/=nlon2.or.nlat/=nlat2) then
         if (mnproc==1) then
            write(lp,'(a)') 'Dimension mismatch in ERA40 data'
            write(lp,'(a,2i5)') 'nlon :',nlon,nlon2
            write(lp,'(a,2i5)') 'nlat :',nlat,nlat2
         end if
         call xcstop('(read_era40_inline)')
         stop '(read_era40_inline)'
      end if
C 
C --- New approach 
      call ncerr(nf90_inq_varid(ncid, "Lo1", varid))
      call ncerr(nf90_get_var  (ncid, varid, flon))
      call ncerr(nf90_inq_varid(ncid, "La1", varid))
      call ncerr(nf90_get_var  (ncid, varid, flat))
      call ncerr(nf90_inq_varid(ncid, "Di", varid))
      call ncerr(nf90_get_var  (ncid, varid, dlon))
      call ncerr(nf90_inq_varid(ncid, "Dj", varid))
      call ncerr(nf90_get_var  (ncid, varid, dlat))
      call ncerr(nf90_close    (ncid))
c
c --- ERA40 precipitation   m/ 6hr -> m/month (rho=1000)
      if (mnproc==1)print *,'ERA40 precip'
      ncfil=era40_path//'fcs.6h.'//cyy//'.TP.nc'   ! Forecast field !
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vprec, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      !call era40_fix(vprec,tmp,nlon,nlat,flon,flat,dlon,dlat)
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat,  
     &             precip(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      precip(:,:,slot)=precip(:,:,slot)/(6.*3600.) ! Precipitation is accumulated field...
c
c --- ERA40 air temperature   K 
      if (mnproc==1) print *,'ERA40 airtmp'
      ncfil=era40_path//'ans.6h.'//cyy//'.2T.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vtair, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &                  airtmp(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      do j= 1-nbdy,jj+nbdy
      do i= 1-nbdy,ii+nbdy
         airtmp(i,j,slot) = airtmp(i,j,slot)-t0deg
      end do
      end do
c
c --- ERA40 clouds   []  (0-1)
      if (mnproc==1)print *,'ERA40 tcc'
      ncfil =era40_path//'ans.6h.'//cyy//'.TCC.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vtcc, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      call interpug(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
     &             clouds(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
c ---  ERA40 winds
      if (mnproc==1)print *,'ERA40 u/v'
      ncfilu=era40_path//'ans.6h.'//cyy//'.10U.nc'
      call ncerr(NF90_open(trim(ncfilu),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vuwnd, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmpu,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      ncfilv=era40_path//'ans.6h.'//cyy//'.10V.nc'
      call ncerr(NF90_open(trim(ncfilv),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vvwnd, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmpv,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
c
!#ifdef ERA40_WIND_POLE_FIX
!         ! Rotatate ERA40 winds to x/y components on plane normal to NP
!         xwnd=cos(elon*radian)*(-era40_vwnd) - sin(elon*radian)*era40_uwnd
!         ywnd=sin(elon*radian)*(-era40_vwnd) + cos(elon*radian)*era40_uwnd
!
!         ! Average northernmost wind components - averaging strategy?
!         !xwnd(:,nlat,1)=sum(xwnd(:,nlat-1,1))/nlon
!         !ywnd(:,nlat,1)=sum(ywnd(:,nlat-1,1))/nlon
!         xwnd(:,nlat,1)=xwnd(:,nlat-1,1)
!         ywnd(:,nlat,1)=ywnd(:,nlat-1,1)
!         era40_vwnd(:,nlat,1)=-xwnd(:,nlat,1)*cos(elon(:,nlat,1)*radian) - &
!                               ywnd(:,nlat,1)*sin(elon(:,nlat,1)*radian)
!         era40_uwnd(:,nlat,1)=-xwnd(:,nlat,1)*sin(elon(:,nlat,1)*radian) + &
!                               ywnd(:,nlat,1)*cos(elon(:,nlat,1)*radian)
!#endif
c
#ifdef KARA2002
c --- Interpolate to p-points
      call interpug(tmpu(:,1:nlat-1) ,nlon,nlat-1,flon,flat,dlon,dlat, 
     &                 uwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpv(:,1:nlat-1) ,nlon,nlat-1,flon,flat,dlon,dlat,
     &                 vwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call xctilr(airtmp(1-nbdy,1-nbdy,  slot),1,1 ,nbdy,nbdy, halo_ps)
      call xctilr(temp  (1-nbdy,1-nbdy,1,slot),1,kk,nbdy,nbdy, halo_ps)
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         wndspd(i,j,slot)=sqrt(uwind(i,j,slot)*uwind(i,j,slot)
     &   +vwind(i,j,slot)*vwind(i,j,slot))
         tmpw2=MAX(2.5,MIN(32.5,wndspd(i,j,slot)))
         CD0 = 1.0E-3*(.692 + .0710*tmpw2 - .000700*tmpw2**2)
         if (temp(i,j,1,slot)<100.) then
            CD1 = 1.0E-3*(.083 - .0054*tmpw2 + .000093*tmpw2**2)
            cd_new= CD0 + CD1*(temp(i,j,1,slot) - airtmp(i,j,slot))
         else
            cd_new= CD0 
         endif
         airdns2=(1013.0 * 100.0 ) / (287.1 * (airtmp(i,j,slot)+t0deg))
         wfact=wndspd(i,j,slot)*airdns2*cd_new
         taux(i,j,slot)=uwind(i,j,slot)*wfact
         tauy(i,j,slot)=vwind(i,j,slot)*wfact
!         if (i.eq.itest.and.j.eq.jtest) then
!            print *,'TST TA-SST',airtmp(i,j,slot),temp(i,j,1,slot)
!         endif
      end do
      end do
#else

c --- Wind speed on original (colocated) grid
      do j=1,nlat
      do i=1,nlon
         tmpw(i,j)=sqrt(tmpu(i,j)*tmpu(i,j)+ tmpv(i,j)*tmpv(i,j))
         wndfac=(1.+sign(1.,tmpw(i,j)-11.))*.5
         cd_new=(0.49+0.065*tmpw(i,j))*1.0e-3*wndfac+cd*(1.-wndfac)
         wfact=tmpw(i,j)*airdns*cd_new
         tmptx(i,j)=tmpu(i,j)*wfact
         tmpty(i,j)=tmpv(i,j)*wfact
      end do
      end do
c --- Interpolate to p-points
      call interpug(tmpu(:,1:nlat-1) ,nlon,nlat-1,flon,flat,dlon,dlat, 
     &                 uwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpv(:,1:nlat-1) ,nlon,nlat-1,flon,flat,dlon,dlat,
     &                 vwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmptx(:,1:nlat-1),nlon,nlat-1,flon,flat,dlon,dlat, 
     &                  taux(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpty(:,1:nlat-1),nlon,nlat-1,flon,flat,dlon,dlat, 
     &                  tauy(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpw(:,1:nlat-1),nlon,nlat-1,flon,flat,dlon,dlat, 
     &                  wndspd(1-nbdy,1-nbdy,slot),plon,plat,itype)   
#endif
c --- Rotate winds/stress to model u/v points
      if (wndflg==1) then
        call rotate2(uwind(1-nbdy,1-nbdy,slot),
     &               vwind(1-nbdy,1-nbdy,slot),
     &               plat,plon,'l2m',.false.)
        call rotate2(taux(1-nbdy,1-nbdy,slot),tauy(1-nbdy,1-nbdy,slot), 
     &               plat(1-nbdy,1-nbdy  ), plon(1-nbdy,1-nbdy  ),
     &               'l2m',.false.)
        call xctilr(uwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_uv)
        call xctilr(vwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_vv)
        call xctilr(taux (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_uv)
        call xctilr(tauy (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_vv)
c --- Rotate winds/stress to model p points
      else if (wndflg==2 .or. wndflg==3) then
        call rotate2(uwind(1-nbdy,1-nbdy,slot),
     &               vwind(1-nbdy,1-nbdy,slot),
     &               plat,plon,'l2m',.true.)
        call rotate2(taux(1-nbdy,1-nbdy,slot),tauy(1-nbdy,1-nbdy,slot), 
     &               plat, plon,  'l2m',.true.)
        call xctilr(uwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(vwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(taux (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(tauy (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
      end if
! --  Pressure
      if (mnproc==1)print *,'ERA40 slp'
      ncfil=era40_path//'ans.6h.'//cyy//'.MSL.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vpres, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      call interpug(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
     &                 slp(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
! --- Relative humidity + vapor mixing ratio
      if (mnproc==1)print *,'ERA40 relhum'
      ncfil=era40_path//'ans.6h.'//cyy//'.2D.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vdewp, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      call interpug(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
     &                  relhum(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         svpair=satvap(airtmp(i,j,slot)+t0deg)
         svpdew=satvap(relhum(i,j,slot))
         rtmp= relhumid(svpair,svpdew,slp(i,j,slot))*0.01
         relhum(i,j,slot)=rtmp
c ---    These may be out of range for non-monotonic interpolation
         relhum(i,j,slot)  =max(0.0,min(relhum(i,j,slot),1.0))
         precip(i,j,slot)  =max(0.0,precip(i,j,slot))
         clouds(i,j,slot)  =max(0.0,min(clouds(i,j,slot),1.0))
         vapmix(i,j,slot) = rhtovpmix(relhum(i,j,slot),
     &                             airtmp(i,j,slot)+t0deg,slp(i,j,slot))
         slp(i,j,slot)=slp(i,j,slot)*0.01 ! Pa -> mBar (HPa)
      enddo
      enddo
c
C --- SW radiation - this one is  off for now
      if (mnproc==1)print *,'ERA40 ssr'
      if (.false.) then
         ncfil =era40_path//'fcs.6h.'//cyy//'.SSRD.nc' ! Forecast field !
         call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
         call ncerr(nf90_inq_varid(ncid, vssrd, varid))
         call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
         call ncerr(nf90_close    (ncid))
c
c ---   ssrd is accumulated over 6hrs, midpoint is 6hr-3hr
         if (irec==maxrec)  then
            ncfil=era40_path//'fcs.6h.'//cyyp1//'.SSRD.nc' ! Forecast field !
            call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
            call ncerr(nf90_inq_varid(ncid, vssrd, varid))
            call ncerr(nf90_get_var  (ncid, varid, tmp2,
     &                                start=(/1,1,1/)))
            call ncerr(nf90_close    (ncid))
         else if (irec<maxrec) then
            ncfil =era40_path//'fcs.6h.'//cyy//'.SSRD.nc' ! Forecast field !
            call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
            call ncerr(nf90_inq_varid(ncid, vssrd, varid))
            call ncerr(nf90_get_var  (ncid, varid, tmp2,
     &                                start=(/1,1,irec+1/)))
            call ncerr(nf90_close    (ncid))
         else
            print *,'You shouldnt be here ...'
            print *,irec,maxrec
            call xcstop('(read_era40)')
            stop '(read_era40)'
         end if
         tmp=(tmp+tmp2)*0.5d0
c
         call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &                 swflx(1-nbdy,1-nbdy,slot),plon,plat,itype)   
         swflx(:,:,slot)=swflx(:,:,slot)/(6*3600) ! Wm^-2 * s - > Wm^-2
         swflx(:,:,slot)=max(0.,swflx(:,:,slot))
      else
         if (mnproc==1) then
          write(lp,'(a)') 'TODO - treat swfl in read_era40_inline'
         end if
      end if

C --- These variables are also available
C --- fblh =era40_path//'fcs.6h.'//cyy//'.BLH.nc'  ! Forecast field !
C --- ficec=era40_path//'ans.6h.'//rt%cyy//'.CI.nc'
C --- fZ   =era40_path//'ans.Z.nc'
c
      lfirst=.false.
c
      END subroutine read_era40
c
c --- Routine for reading era40 data from netcdf files, inline version
c --- It is assumed that era40 files are organized in years and variable
c --- names. Each netcdf file must have one record every 6 hour for 
c --- that year. These records must come in chronological order,
c --- starting Jan 1 at 00:00:00
c --- TODO : Now all threads read the same netcdf file, 
c
c
c --- Routine for reading era-i data from netcdf files, inline version
c --- It is assumed that era-i files are organized in years and variable
c --- names. Each netcdf file must have one record every 6 hour for 
c --- that year. These records must come in chronological order,
c --- starting Jan 1 at 00:00:00
c --- TODO : Now all threads read the same netcdf file, 
      subroutine read_erai(dtime,slot)
      use mod_xc
      use mod_common_ice
      use mod_hycom_nersc
      use netcdf
      implicit none
      real*8 ,  intent(in) :: dtime
      integer,  intent(in) :: slot
c
c --- Data catalogue name
      character(len=*), parameter :: erai_path='./ERA-I/'
c
c --- Variable names - ERA-I 0.5 by 0.5   
      character(len=*), parameter :: vuwnd='10U'
      character(len=*), parameter :: vvwnd='10V'
      character(len=*), parameter :: vtair='2T'
      character(len=*), parameter :: vdewp='2D'
      character(len=*), parameter :: vprec='TP'
      character(len=*), parameter :: vtcc ='TCC'
      character(len=*), parameter :: vblh ='BLH'
      character(len=*), parameter :: vlmsk='LSM_sfc'
      character(len=*), parameter :: vpres='MSL'
      character(len=*), parameter :: vssrd='SSRD'
c
c --- Hardcoded grid sizes 
      integer, parameter :: nlon=720, nlat=361
c
      integer  :: i,j,k,ix,jy,irec,maxrec, nlon2, nlat2
      real, dimension(nlon, nlat) :: 
     &   tmp, tmpu, tmpv, tmpw, tmp2, tmptx, tmpty
      real               :: e40_lon(nlon), e40_lat(nlat)
      real               :: rfac,wndfac,cd_new,w4,wfact,svpair, svpdew
      real               :: dlon,dlat,flon,flat,llon,llat,rtmp
      integer            :: dimid, varid, ncid
      integer            :: iyear, ihour, iday
      logical            :: lfirst=.true.
      character(len=80)  :: ncfil, ncfilu, ncfilv
      character(len=4)   :: cyy, cyyp1
      integer, parameter :: itype=1 ! Sets interpolation type bicubic
#ifdef KARA2002
      real               :: tmpw2,CD0,CD1,airdns2
#endif
c      integer, parameter :: itype=0 ! Sets interpolation type bilinear
c
      include 'common_blocks.h'
c
      call forday(dtime,yrflag, iyear,iday,ihour)
      irec  =  (iday-1)*4 + ihour/6 +1
      if (mnproc==1) 
     &  write(lp,'(a,i6,a,i4,x,i3,x,i2)')'ERA-I data from record: ',
     &  irec, '  time= ',iyear,iday,ihour
c
c --- Set file names
      write(cyy  ,'(i4.4)') iyear
      write(cyyp1,'(i4.4)') iyear+1

c --- Check dimensions against hardcoded ones
c --- Get longitude  & latitude for data  - use land mask file
c      ncfil=erai_path//'ans.LSM.nc'
c      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
c      print *, 'toto2'
!      call ncerr(nf90_inq_dimid(ncid, "lon", dimid))
!      call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlon2))
!      call ncerr(nf90_inq_dimid(ncid, "lat", dimid))
!      call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlat2))
!      call ncerr(nf90_inq_dimid(ncid, "record", dimid))
!      print *, 'toto3'
!      call ncerr(nf90_inquire_dimension(ncid, dimid, len=maxrec))
!      print *, 'toto4'
!
!      if (nlon/=nlon2.or.nlat/=nlat2) then
!         if (mnproc==1) then
!            write(lp,'(a)') 'Dimension mismatch in ERA40 data'
!            write(lp,'(a,2i5)') 'nlon :',nlon,nlon2
!            write(lp,'(a,2i5)') 'nlat :',nlat,nlat2
!         end if
!         call xcstop('(read_erai_inline)')
!         stop '(read_erai_inline)'
!      end if
C 
C --- New approach 
C      call ncerr(nf90_inq_varid(ncid, "Lo1", varid))
C      call ncerr(nf90_get_var  (ncid, varid, flon))
      flon=-180
      flat=90
      dlon=0.5
      dlat=-0.5
C      call ncerr(nf90_inq_varid(ncid, "La1", varid))
C      call ncerr(nf90_get_var  (ncid, varid, flat))
C      call ncerr(nf90_inq_varid(ncid, "Di", varid))
C      call ncerr(nf90_get_var  (ncid, varid, dlon))
C      call ncerr(nf90_inq_varid(ncid, "Dj", varid))
C      call ncerr(nf90_get_var  (ncid, varid, dlat))
c      call ncerr(nf90_close    (ncid))
c
c --- ERA-I precipitation   m/ 6hr -> m/month (rho=1000)
      if (mnproc==1)print *,'ERA-I precip'
      ncfil=erai_path//'fcs.6h.'//cyy//'.TP.nc'   ! Forecast field !
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vprec, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
!      call erai_fix(vprec,tmp,nlon,nlat,flon,flat,dlon,dlat)
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat,  
     &             precip(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      precip(:,:,slot)=precip(:,:,slot)/(6.*3600.) ! Precipitation is accumulated field...
c
c --- ERA-I air temperature   K 
      if (mnproc==1) print *,'ERA-I airtmp'
      ncfil=erai_path//'ans.6h.'//cyy//'.2T.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vtair, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &                  airtmp(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      do j= 1-nbdy,jj+nbdy
      do i= 1-nbdy,ii+nbdy
         airtmp(i,j,slot) = airtmp(i,j,slot)-t0deg
      end do
      end do
c
c --- ERA-I clouds   []  (0-1)
      if (mnproc==1)print *,'ERA-I tcc'
      ncfil =erai_path//'ans.6h.'//cyy//'.TCC.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vtcc, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      call interpug(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
     &             clouds(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
c ---  ERA-I winds
      if (mnproc==1)print *,'ERA-I u/v'
      ncfilu=erai_path//'ans.6h.'//cyy//'.10U.nc'
      call ncerr(NF90_open(trim(ncfilu),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vuwnd, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmpu,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      ncfilv=erai_path//'ans.6h.'//cyy//'.10V.nc'
      call ncerr(NF90_open(trim(ncfilv),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vvwnd, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmpv,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
c
c --- Wind speed on original (colocated) grid
#ifdef KARA2002
      call interpug(tmpu(:,1:nlat-1) ,nlon,nlat-1,flon,flat,dlon,dlat, 
     &                 uwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpv(:,1:nlat-1) ,nlon,nlat-1,flon,flat,dlon,dlat,
     &                 vwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call xctilr(airtmp(1-nbdy,1-nbdy,  slot),1,1 ,nbdy,nbdy, halo_ps)
      call xctilr(temp  (1-nbdy,1-nbdy,1,slot),1,kk,nbdy,nbdy, halo_ps)
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         wndspd(i,j,slot)=sqrt(uwind(i,j,slot)*uwind(i,j,slot)
     &   +vwind(i,j,slot)*vwind(i,j,slot))
         tmpw2=MAX(2.5,MIN(32.5,wndspd(i,j,slot)))
         CD0 = 1.0E-3*(.692 + .0710*tmpw2 - .000700*tmpw2**2)
         if (temp(i,j,1,slot)<100.) then
            CD1 = 1.0E-3*(.083 - .0054*tmpw2 + .000093*tmpw2**2)
            cd_new= CD0 + CD1*(temp(i,j,1,slot) - airtmp(i,j,slot))
         else
            cd_new= CD0 
         endif
         airdns2=(1013.0 * 100.0 ) / (287.1 * (airtmp(i,j,slot)+t0deg))
         wfact=wndspd(i,j,slot)*airdns2*cd_new
         taux(i,j,slot)=uwind(i,j,slot)*wfact
         tauy(i,j,slot)=vwind(i,j,slot)*wfact
!         if (i.eq.itest.and.j.eq.jtest) then
!            print *,'TST TA-SST',airtmp(i,j,slot),temp(i,j,1,slot)
!         endif
      end do
      end do
#else
      do j=1,nlat
      do i=1,nlon
         tmpw(i,j)=sqrt(tmpu(i,j)*tmpu(i,j)+ tmpv(i,j)*tmpv(i,j))
         wndfac=(1.+sign(1.,tmpw(i,j)-11.))*.5
         cd_new=(0.49+0.065*tmpw(i,j))*1.0e-3*wndfac+cd*(1.-wndfac)
         wfact=tmpw(i,j)*airdns*cd_new
         tmptx(i,j)=tmpu(i,j)*wfact
         tmpty(i,j)=tmpv(i,j)*wfact
      end do
      end do
c --- Interpolate to p-points
      call interpug(tmpu(:,1:nlat-1) ,nlon,nlat-1,flon,flat,dlon,dlat, 
     &                 uwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpv(:,1:nlat-1) ,nlon,nlat-1,flon,flat,dlon,dlat,
     &                 vwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmptx(:,1:nlat-1),nlon,nlat-1,flon,flat,dlon,dlat, 
     &                  taux(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpty(:,1:nlat-1),nlon,nlat-1,flon,flat,dlon,dlat, 
     &                  tauy(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpw(:,1:nlat-1),nlon,nlat-1,flon,flat,dlon,dlat, 
     &                  wndspd(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
#endif
c --- Rotate winds/stress to model u/v points
      if (wndflg==1) then
        call rotate2(uwind(1-nbdy,1-nbdy,slot),
     &               vwind(1-nbdy,1-nbdy,slot),
     &               plat,plon,'l2m',.false.)
        call rotate2(taux(1-nbdy,1-nbdy,slot),tauy(1-nbdy,1-nbdy,slot), 
     &               plat(1-nbdy,1-nbdy  ), plon(1-nbdy,1-nbdy  ),
     &               'l2m',.false.)
        call xctilr(uwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_uv)
        call xctilr(vwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_vv)
        call xctilr(taux (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_uv)
        call xctilr(tauy (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_vv)
c --- Rotate winds/stress to model p points
      else if (wndflg==2 .or. wndflg==3) then
        call rotate2(uwind(1-nbdy,1-nbdy,slot),
     &               vwind(1-nbdy,1-nbdy,slot),
     &               plat,plon,'l2m',.true.)
        call rotate2(taux(1-nbdy,1-nbdy,slot),tauy(1-nbdy,1-nbdy,slot), 
     &               plat, plon,  'l2m',.true.)
        call xctilr(uwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(vwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(taux (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(tauy (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
      end if
! --  Pressure
      if (mnproc==1)print *,'ERA-I slp'
      ncfil=erai_path//'ans.6h.'//cyy//'.MSL.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vpres, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      call interpug(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
     &                 slp(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
! --- Relative humidity + vapor mixing ratio
      if (mnproc==1)print *,'ERA-I relhum'
      ncfil=erai_path//'ans.6h.'//cyy//'.2D.nc'
      call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_varid(ncid, vdewp, varid))
      call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
      call ncerr(nf90_close    (ncid))
      call interpug(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
     &                  relhum(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         svpair=satvap(airtmp(i,j,slot)+t0deg)
         svpdew=satvap(relhum(i,j,slot))
         rtmp= relhumid(svpair,svpdew,slp(i,j,slot))*0.01
         relhum(i,j,slot)=rtmp
c ---    These may be out of range for non-monotonic interpolation
         relhum(i,j,slot)  =max(0.0,min(relhum(i,j,slot),1.0))
         precip(i,j,slot)  =max(0.0,precip(i,j,slot))
         clouds(i,j,slot)  =max(0.0,min(clouds(i,j,slot),1.0))
         vapmix(i,j,slot) = rhtovpmix(relhum(i,j,slot),
     &                             airtmp(i,j,slot)+t0deg,slp(i,j,slot))
         slp(i,j,slot)=slp(i,j,slot)*0.01 ! Pa -> mBar (HPa)
      enddo
      enddo
c
C --- SW radiation - this one is  off for now
      if (mnproc==1)print *,'ERA-I ssr'
      if (.false.) then
         ncfil =erai_path//'fcs.6h.'//cyy//'.SSRD.nc' ! Forecast field !
         call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
         call ncerr(nf90_inq_varid(ncid, vssrd, varid))
         call ncerr(nf90_get_var  (ncid, varid, tmp,start=(/1,1,irec/)))
         call ncerr(nf90_close    (ncid))
c
c ---   ssrd is accumulated over 6hrs, midpoint is 6hr-3hr
         if (irec==maxrec)  then
            ncfil=erai_path//'fcs.6h.'//cyyp1//'.SSRD.nc' ! Forecast field !
            call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
            call ncerr(nf90_inq_varid(ncid, vssrd, varid))
            call ncerr(nf90_get_var  (ncid, varid, tmp2,
     &                                start=(/1,1,1/)))
            call ncerr(nf90_close    (ncid))
         else if (irec<maxrec) then
            ncfil =erai_path//'fcs.6h.'//cyy//'.SSRD.nc' ! Forecast field !
            call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
            call ncerr(nf90_inq_varid(ncid, vssrd, varid))
            call ncerr(nf90_get_var  (ncid, varid, tmp2,
     &                                start=(/1,1,irec+1/)))
            call ncerr(nf90_close    (ncid))
         else
            print *,'You shouldnt be here ...'
            print *,irec,maxrec
            call xcstop('(read_erai)')
            stop '(read_erai)'
         end if
         tmp=(tmp+tmp2)*0.5d0
c
         call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &                 swflx(1-nbdy,1-nbdy,slot),plon,plat,itype)   
         swflx(:,:,slot)=swflx(:,:,slot)/(6*3600) ! Wm^-2 * s - > Wm^-2
         swflx(:,:,slot)=max(0.,swflx(:,:,slot))
      else
         if (mnproc==1) then
          write(lp,'(a)') 'TODO - treat swfl in read_erai_inline'
         end if
      end if

C --- These variables are also available
C --- fblh =erai_path//'fcs.6h.'//cyy//'.BLH.nc'  ! Forecast field !
C --- ficec=erai_path//'ans.6h.'//rt%cyy//'.CI.nc'
C --- fZ   =erai_path//'ans.Z.nc'
c
      lfirst=.false.
c
      END subroutine read_erai
c
c --- Routine for reading era40 data from netcdf files, inline version
c --- It is assumed that era40 files are organized in years and variable
c --- names. Each netcdf file must have one record every 6 hour for 
c --- that year. These records must come in chronological order,
c --- starting Jan 1 at 00:00:00
c --- TODO : Now all threads read the same netcdf file, 

      subroutine read_ecmwf_nc(dtime,slot)
      use mod_xc
      use netcdf
      use mod_common_ice
      use mod_hycom_nersc
      implicit none
      real*8, intent(in) :: dtime
      integer, intent(in) :: slot
c
      ! Data catalogue name
      character(len=*), parameter :: ecmwf_path ='./Ecmwf.nc/'
      character(len=*), parameter :: era40_cpath='./ERA40-Clim/'
      !character(len=*), parameter :: ecmwf_path='./ECMWFR_T799/'
c
      ! Variable names - Ecmwf 0.25 x 0.25
      character(len=*), parameter :: vuwnd='U10M'
      character(len=*), parameter :: vvwnd='V10M'
      character(len=*), parameter :: vtair='T2M'
      character(len=*), parameter :: vdewp='D2M'
      character(len=*), parameter :: vprec='TP'
      character(len=*), parameter :: vtcc ='TCC'
      character(len=*), parameter :: vpres='MSL'
c
c --- Hardcoded grid sizes 
      integer, parameter :: nlon=1440, nlat=721
c
      real, dimension(nlon, nlat) :: 
     &   tmp, tmpu, tmpv, tmpw, tmp2, tmptx, tmpty, lmsk
      real               :: eclon(nlon), eclat(nlat), rtmp
c
      character(len= 4)  :: cyy, cyyp1
      character(len=80)  :: ncfil, ncfilu, ncfilv
      integer  :: i,j,k,irec,maxrec
      real :: rfac,wndfac,cd_new,w4,wfact,svpair, svpdew,wspd
      integer :: nlon2,nlat2
      integer :: im1,ip1,jm1,jp1,iyear, iday, ihour
      integer :: ncid, dimid, varid, ierr
      logical :: first=.true.
c --- Files are kept open, so we store file and var ids
      integer, save :: ncidu=-9999, ncidv=-9999,nciddp=-9999,
     &   ncidat=-9999,ncidpr=-9999, ncidcc=-9999
      integer, save :: varidu,varidv,variddp,varidat,varidpr,
     &   varidcc
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     &   tu, tv
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: 
     &   cprecip
c     integer, parameter :: itype=1 ! Sets interpolation type bicubic
      integer, parameter :: itype=0 ! Sets interpolation type bilinear
#ifdef KARA2002
      real               :: tmpw2,CD0,CD1,airdns2
#endif
c
      include 'common_blocks.h'
c
      call forday(dtime,yrflag, iyear,iday,ihour)
      irec  =  (iday-1)*4 + ihour/6 +1
      if (mnproc==1) 
     &  write(lp,'(a,i6,a,i4,x,i3,x,i2)')'Ecnc data from record: ',
     &  irec, '  time= ',iyear,iday,ihour
c
c --- Set file names
      write(cyy  ,'(i4.4)') iyear
      write(cyyp1,'(i4.4)') iyear+1
c
c
c --- Check dimensions against hardcoded ones
c --- Get longitude  & latitude for data  - use land mask file
c --- Only on first pass
      if (first) then
         ncfil=ecmwf_path//'ec_atmo_geo_mad_T2M_'//cyy//'.nc'
         call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
         call ncerr(nf90_inq_dimid(ncid, "lon", dimid))
         call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlon2))
         call ncerr(nf90_inq_dimid(ncid, "lat", dimid))
         call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlat2))
         call ncerr(nf90_inq_dimid(ncid, "time", dimid))
         call ncerr(nf90_inquire_dimension(ncid, dimid, len=maxrec))
         if (nlon/=nlon2.or.nlat/=nlat2) then
            if (mnproc==1) then
               write(lp,'(a)') 'Dimension mismatch in Ecnc data'
               write(lp,'(a,2i5)') 'nlon :',nlon,nlon2
               write(lp,'(a,2i5)') 'nlat :',nlat,nlat2
            end if
            call xcstop('(read_ecmwf_nc)')
            stop '(read_ecmwf_nc)'
         end if
         call ncerr(nf90_inq_varid(ncid, "lon", varid))
         call ncerr(nf90_get_var  (ncid, varid, eclon))
         call ncerr(nf90_inq_varid(ncid, "lat", varid))
         call ncerr(nf90_get_var  (ncid, varid, eclat))
         call ncerr(nf90_close    (ncid))
         flon=eclon(1) ; flat=eclat(1) ;
         dlon=eclon(2)-eclon(1) ; dlat=eclat(2)-eclat(1)
      end if
c
c
c --- Ecmwf precipitation   m/ 6hr -> m/month (rho=1000)
c --- TODO: ust be handled by climatology 
c
c ---    Read one precipitation field from climatology. cprecip is kept
c ---    locally.
      if (mnproc==1)print *,'Ecnc precip (climatology)'
      if (hf_cinit) then
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc0),'precip',hf_mc0)
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc1),'precip',hf_mc1)
      else if (hf_cread) then
       cprecip(:,:,hf_lc0)=cprecip(:,:,hf_lc1)
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc1),'precip',hf_mc1)
      end if
c ---    Interpolate climatology onto hf forcing field times
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         precip(i,j,slot) = hf_wc0*cprecip(i,j,hf_lc0) +
     &                   hf_wc1*cprecip(i,j,hf_lc1) 
      end do
      end do
c
c --- Ecnc air temperature   K 
      if (mnproc==1) print *,'Ecnc airtmp'
      ncfil=ecmwf_path//'ec_atmo_geo_mad_T2M_'//cyy//'.nc'
      call rdncrec(ncfil,vtair,tmp,nlon,nlat,irec)
c     call bilin_ecmwf2(tmp,nlon,nlat,flon,flat,dlon,dlat, 
c    &                  airtmp(1-nbdy,1-nbdy,slot),plon,plat)   
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              airtmp(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      do j= 1-nbdy,jj+nbdy
      do i= 1-nbdy,ii+nbdy
         airtmp(i,j,slot) = airtmp(i,j,slot)-t0deg
      end do
      end do
c
c --- Ecnc clouds   []  (0-1)
      if (mnproc==1)print *,'Ecnc tcc'
      ncfil =ecmwf_path//'ec_atmo_geo_mad_TCC_'//cyy//'.nc'
      call rdncrec(ncfil,vtcc,tmp,nlon,nlat,irec)
c     call bilin_ecmwf2(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
c    &                 clouds(1-nbdy,1-nbdy,slot),plon,plat)   
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              clouds(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
c ---  Ecnc winds
      if (mnproc==1)print *,'Ecnc u/v'
      ncfilu=ecmwf_path//'ec_atmo_geo_mad_U10M_'//cyy//'.nc'
      ncfilv=ecmwf_path//'ec_atmo_geo_mad_V10M_'//cyy//'.nc'
      call rdncrec(ncfilu,vuwnd,tmpu,nlon,nlat,irec)
      call rdncrec(ncfilv,vvwnd,tmpv,nlon,nlat,irec)
c
c
c --- New approach. The original grid is quite large - so
c --- for small MPI tiles this should be more effective. 
c
c --- Interpolate u/v wind to p-points
      call interpug(tmpu,nlon,nlat,flon,flat,dlon,dlat, 
     &              uwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpv,nlon,nlat,flon,flat,dlon,dlat, 
     &              vwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
#ifdef KARA2002
      call xctilr(airtmp(1-nbdy,1-nbdy,  slot),1,1 ,nbdy,nbdy, halo_ps)
      call xctilr(temp  (1-nbdy,1-nbdy,1,slot),1,kk,nbdy,nbdy, halo_ps)
#endif
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         wndspd(i,j,slot)=sqrt(uwind(i,j,slot)*uwind(i,j,slot)
     &   +vwind(i,j,slot)*vwind(i,j,slot))
#ifdef KARA2002
         tmpw2=MAX(2.5,MIN(32.5,wndspd(i,j,slot)))
         CD0 = 1.0E-3*(.692 + .0710*tmpw2 - .000700*tmpw2**2)
         if (temp(i,j,1,slot)<100.) then
            CD1 = 1.0E-3*(.083 - .0054*tmpw2 + .000093*tmpw2**2)
            cd_new= CD0 + CD1*(temp(i,j,1,slot) - airtmp(i,j,slot))
         else
            cd_new= CD0 
         endif
         airdns2=(1013.0 * 100.0 ) / (287.1 * (airtmp(i,j,slot)+t0deg))
         wfact=wndspd(i,j,slot)*airdns2*cd_new
#else
         wndfac=(1.+sign(1.,wndspd(i,j,slot)-11.))*.5
         cd_new=(0.49+0.065*wndspd(i,j,slot))*1.0e-3*wndfac
     &   +cd*(1.-wndfac)
         wfact=wndspd(i,j,slot)*airdns*cd_new
#endif
         taux(i,j,slot)=uwind(i,j,slot)*wfact
         tauy(i,j,slot)=vwind(i,j,slot)*wfact
      end do
      end do

c --- Rotate winds/stress to model u/v points
      if (wndflg==1) then
        call rotate2(uwind(1-nbdy,1-nbdy,slot),
     &               vwind(1-nbdy,1-nbdy,slot),
     &               plat,plon,'l2m',.false.)
        call rotate2(taux(1-nbdy,1-nbdy,slot),tauy(1-nbdy,1-nbdy,slot), 
     &               plat(1-nbdy,1-nbdy  ), plon(1-nbdy,1-nbdy  ),
     &               'l2m',.false.)
        call xctilr(uwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_uv)
        call xctilr(vwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_vv)
        call xctilr(taux (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_uv)
        call xctilr(tauy (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_vv)
c --- Rotate winds/stress to model p points
      else if (wndflg==2 .or. wndflg==3) then
        call rotate2(uwind(1-nbdy,1-nbdy,slot),
     &               vwind(1-nbdy,1-nbdy,slot),
     &               plat,plon,'l2m',.true.)
        call rotate2(taux(1-nbdy,1-nbdy,slot),tauy(1-nbdy,1-nbdy,slot), 
     &               plat, plon,  'l2m',.true.)
        call xctilr(uwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(vwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(taux (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(tauy (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
      end if

c --  Pressure
      if (mnproc==1)print *,'Ecnc slp'
      ncfil=ecmwf_path//'ec_atmo_geo_mad_MSL_'//cyy//'.nc'
      call rdncrec(ncfil,vpres,tmp,nlon,nlat,irec)
C     call bilin_ecmwf2(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
C    &                 slp(1-nbdy,1-nbdy,slot),plon,plat)   
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              slp(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
c --- Relative humidity ( + vapor mixing ratio )
      if (mnproc==1)print *,'Ecnc relhum'
      ncfil=ecmwf_path//'ec_atmo_geo_mad_D2M_'//cyy//'.nc'
      call rdncrec(ncfil,vdewp,tmp,nlon,nlat,irec)
C     call bilin_ecmwf2(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
C    &                  relhum(1-nbdy,1-nbdy,slot),plon,plat)   
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              relhum(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
c --- Downward shortwave radiation
      if (mnproc==1)print *,'Ecnc ssr'
      if (.false.) then
      else
         if (mnproc==1) then
          write(lp,'(a)') 'TODO - treat swfl in read_era40_inline'
         end if
      end if
c
c --- Final processing
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         svpair=satvap(airtmp(i,j,slot)+t0deg)
         svpdew=satvap(relhum(i,j,slot))
         rtmp= relhumid(svpair,svpdew,slp(i,j,slot))*0.01
         relhum(i,j,slot)=rtmp
c ---    These may be out of range for non-monotonic interpolation
         relhum(i,j,slot)  =max(0.0,min(relhum(i,j,slot),1.0))
         precip(i,j,slot)  =max(0.0,precip(i,j,slot))
         clouds(i,j,slot)  =max(0.0,min(clouds(i,j,slot),1.0))
         vapmix(i,j,slot) = rhtovpmix(relhum(i,j,slot),
     &                             airtmp(i,j,slot)+t0deg,slp(i,j,slot))
         slp   (i,j,slot)=     slp(i,j,slot)*0.01  ! Pa -> mBar (HPa)
      enddo
      enddo
      first=.false.
      END subroutine read_ecmwf_nc


c --- Routine for reading era40 data from netcdf files, inline version
c --- It is assumed that era40 files are organized in years and variable
c --- names. Each netcdf file must have one record every 6 hour for 
c --- that year. These records must come in chronological order,
c --- starting Jan 1 at 00:00:00
c --- TODO : Now all threads read the same netcdf file, 

      subroutine read_ecmwf_nc2(dtime,slot)
      use mod_xc
      use netcdf
      use mod_common_ice
      use mod_hycom_nersc
      implicit none
      real*8, intent(in) :: dtime
      integer, intent(in) :: slot
c
      ! Data catalogue name
      character(len=*), parameter :: ecmwf_path ='./Ecmwf.nc/'
      character(len=*), parameter :: era40_cpath='./ERA40-Clim/'
      !character(len=*), parameter :: ecmwf_path='./ECMWFR_T799/'
c
      ! Variable names - Ecmwf 1/8 x 1/8 
      character(len=*), parameter :: vuwnd='U10M'
      character(len=*), parameter :: vvwnd='V10M'
      character(len=*), parameter :: vtair='T2M'
      character(len=*), parameter :: vdewp='D2M'
      character(len=*), parameter :: vprec='TP'
      character(len=*), parameter :: vtcc ='TCC'
      character(len=*), parameter :: vpres='MSL'
c
c --- Hardcoded grid sizes 
      integer, parameter :: nlon=2880, nlat=1441
c
      real, dimension(nlon, nlat) :: 
     &   tmp, tmpu, tmpv, tmpw, tmp2, tmptx, tmpty, lmsk
      real               :: eclon(nlon), eclat(nlat), rtmp
c
      character(len= 4)  :: cyy, cyyp1
      character(len=80)  :: ncfil, ncfilu, ncfilv
      integer  :: i,j,k,irec,maxrec
      real :: rfac,wndfac,cd_new,w4,wfact,svpair, svpdew,wspd
      integer :: nlon2,nlat2
      integer :: im1,ip1,jm1,jp1,iyear, iday, ihour
      integer :: ncid, dimid, varid, ierr
      logical :: first=.true.
c --- Files are kept open, so we store file and var ids
      integer, save :: ncidu=-9999, ncidv=-9999,nciddp=-9999,
     &   ncidat=-9999,ncidpr=-9999, ncidcc=-9999
      integer, save :: varidu,varidv,variddp,varidat,varidpr,
     &   varidcc
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     &   tu, tv
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: 
     &   cprecip
c     integer, parameter :: itype=1 ! Sets interpolation type bicubic
      integer, parameter :: itype=0 ! Sets interpolation type bilinear
#ifdef KARA2002
      real               :: tmpw2,CD0,CD1,airdns2
#endif
c
      include 'common_blocks.h'
c
      call forday(dtime,yrflag, iyear,iday,ihour)
      irec  =  (iday-1)*4 + ihour/6 +1
      if (mnproc==1) 
     &  write(lp,'(a,i6,a,i4,x,i3,x,i2)')'Ecnc2 data from record: ',
     &  irec, '  time= ',iyear,iday,ihour
c
c --- Set file names
      write(cyy  ,'(i4.4)') iyear
      write(cyyp1,'(i4.4)') iyear+1
c
c --- Check dimensions against hardcoded ones
c --- Get longitude  & latitude for data  - use land mask file
c --- Only on first pass
      if (first) then
         ncfil=ecmwf_path//'ec_atmo_geo_la_T2M_'//cyy//'.nc'
         call ncerr(NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid))
         call ncerr(nf90_inq_dimid(ncid, "lon", dimid))
         call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlon2))
         call ncerr(nf90_inq_dimid(ncid, "lat", dimid))
         call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlat2))
         call ncerr(nf90_inq_dimid(ncid, "time", dimid))
         call ncerr(nf90_inquire_dimension(ncid, dimid, len=maxrec))
         if (nlon/=nlon2.or.nlat/=nlat2) then
            if (mnproc==1) then
               write(lp,'(a)') 'Dimension mismatch in Ecnc2 data'
               write(lp,'(a,2i5)') 'nlon :',nlon,nlon2
               write(lp,'(a,2i5)') 'nlat :',nlat,nlat2
            end if
            call xcstop('(read_ecmwf_nc)')
            stop '(read_ecmwf_nc)'
         end if
         call ncerr(nf90_inq_varid(ncid, "lon", varid))
         call ncerr(nf90_get_var  (ncid, varid, eclon))
         call ncerr(nf90_inq_varid(ncid, "lat", varid))
         call ncerr(nf90_get_var  (ncid, varid, eclat))
         call ncerr(nf90_close    (ncid))
         flon=eclon(1) ; flat=eclat(1) ;
         dlon=eclon(2)-eclon(1) ; dlat=eclat(2)-eclat(1)
      end if
c
c
c --- Ecmwf precipitation   m/ 6hr -> m/month (rho=1000)
c --- TODO: ust be handled by climatology 
c
c ---    Read one precipitation field from climatology. cprecip is kept
c ---    locally.
      if (mnproc==1)print *,'Ecnc2 precip (climatology)'
      if (hf_cinit) then
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc0),'precip',hf_mc0)
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc1),'precip',hf_mc1)
      else if (hf_cread) then
       cprecip(:,:,hf_lc0)=cprecip(:,:,hf_lc1)
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc1),'precip',hf_mc1)
      end if
c ---    Interpolate climatology onto hf forcing field times
      do j=1-nbdy,jdm+nbdy
      do i=1-nbdy,idm+nbdy
         precip(i,j,slot) = hf_wc0*cprecip(i,j,hf_lc0) +
     &                   hf_wc1*cprecip(i,j,hf_lc1) 
      end do
      end do
c
c --- Ecnc air temperature   K 
      if (mnproc==1) print *,'Ecnc2 airtmp'
      ncfil=ecmwf_path//'ec_atmo_geo_la_T2M_'//cyy//'.nc'
      call rdncrec(ncfil,vtair,tmp,nlon,nlat,irec)
c     call bilin_ecmwf2(tmp,nlon,nlat,flon,flat,dlon,dlat, 
c    &                  airtmp(1-nbdy,1-nbdy,slot),plon,plat)   
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              airtmp(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      do j= 1-nbdy,jj+nbdy
      do i= 1-nbdy,ii+nbdy
         airtmp(i,j,slot) = airtmp(i,j,slot)-t0deg
      end do
      end do
c
c --- Ecnc clouds   []  (0-1)
      if (mnproc==1)print *,'Ecnc2 tcc'
      ncfil =ecmwf_path//'ec_atmo_geo_la_TCC_'//cyy//'.nc'
      call rdncrec(ncfil,vtcc,tmp,nlon,nlat,irec)
c     call bilin_ecmwf2(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
c    &                 clouds(1-nbdy,1-nbdy,slot),plon,plat)   
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              clouds(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
c ---  Ecnc winds
      if (mnproc==1)print *,'Ecnc2 u/v'
      ncfilu=ecmwf_path//'ec_atmo_geo_la_U10M_'//cyy//'.nc'
      ncfilv=ecmwf_path//'ec_atmo_geo_la_V10M_'//cyy//'.nc'
      call rdncrec(ncfilu,vuwnd,tmpu,nlon,nlat,irec)
      call rdncrec(ncfilv,vvwnd,tmpv,nlon,nlat,irec)
c
c
c --- New approach. The original grid is quite large - so
c --- for small MPI tiles this should be more effective. 
c
c --- Interpolate u/v wind to p-points
      call interpug(tmpu,nlon,nlat,flon,flat,dlon,dlat, 
     &              uwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
      call interpug(tmpv,nlon,nlat,flon,flat,dlon,dlat, 
     &              vwind(1-nbdy,1-nbdy,slot),plon,plat,itype)   
#ifdef KARA2002
      call xctilr(airtmp(1-nbdy,1-nbdy,  slot),1,1 ,nbdy,nbdy, halo_ps)
      call xctilr(temp  (1-nbdy,1-nbdy,1,slot),1,kk,nbdy,nbdy, halo_ps)
#endif
      do j=1-nbdy,jdm+nbdy
      do i=1-nbdy,idm+nbdy
         wndspd(i,j,slot)=sqrt(uwind(i,j,slot)*uwind(i,j,slot)
     &   +vwind(i,j,slot)*vwind(i,j,slot))
#ifdef KARA2002
         tmpw2=MAX(2.5,MIN(32.5,wndspd(i,j,slot)))
         CD0 = 1.0E-3*(.692 + .0710*tmpw2 - .000700*tmpw2**2)
         if (temp(i,j,1,slot)<100.) then
            CD1 = 1.0E-3*(.083 - .0054*tmpw2 + .000093*tmpw2**2)
            cd_new= CD0 + CD1*(temp(i,j,1,slot) - airtmp(i,j,slot))
         else
            cd_new= CD0 
         endif
         airdns2=(1013.0 * 100.0 ) / (287.1 * (airtmp(i,j,slot)+t0deg))
         wfact=wndspd(i,j,slot)*airdns2*cd_new
#else
         wndfac=(1.+sign(1.,wndspd(i,j,slot)-11.))*.5
         cd_new=(0.49+0.065*wndspd(i,j,slot))*1.0e-3*wndfac
     &   +cd*(1.-wndfac)
         wfact=wndspd(i,j,slot)*airdns*cd_new
#endif
         taux(i,j,slot)=uwind(i,j,slot)*wfact
         tauy(i,j,slot)=vwind(i,j,slot)*wfact
      end do
      end do

c --- Rotate winds/stress to model u/v points
      if (wndflg==1) then
        call rotate2(uwind(1-nbdy,1-nbdy,slot),
     &               vwind(1-nbdy,1-nbdy,slot),
     &               plat,plon,'l2m',.false.)
        call rotate2(taux(1-nbdy,1-nbdy,slot),tauy(1-nbdy,1-nbdy,slot), 
     &               plat(1-nbdy,1-nbdy  ), plon(1-nbdy,1-nbdy  ),
     &               'l2m',.false.)
        call xctilr(uwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_uv)
        call xctilr(vwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_vv)
        call xctilr(taux (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_uv)
        call xctilr(tauy (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_vv)
c --- Rotate winds/stress to model p points
      else if (wndflg==2 .or. wndflg==3) then
        call rotate2(uwind(1-nbdy,1-nbdy,slot),
     &               vwind(1-nbdy,1-nbdy,slot),
     &               plat,plon,'l2m',.true.)
        call rotate2(taux(1-nbdy,1-nbdy,slot),tauy(1-nbdy,1-nbdy,slot), 
     &               plat, plon,  'l2m',.true.)
        call xctilr(uwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(vwind(1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(taux (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
        call xctilr(tauy (1-nbdy,1-nbdy,slot),1,1,nbdy,nbdy,halo_ps)
      end if

c --  Pressure
      if (mnproc==1)print *,'Ecnc2 slp'
      ncfil=ecmwf_path//'ec_atmo_geo_la_MSL_'//cyy//'.nc'
      call rdncrec(ncfil,vpres,tmp,nlon,nlat,irec)
C     call bilin_ecmwf2(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
C    &                 slp(1-nbdy,1-nbdy,slot),plon,plat)   
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              slp(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
c --- Relative humidity ( + vapor mixing ratio )
      if (mnproc==1)print *,'Ecnc2 relhum'
      ncfil=ecmwf_path//'ec_atmo_geo_la_D2M_'//cyy//'.nc'
      call rdncrec(ncfil,vdewp,tmp,nlon,nlat,irec)
C     call bilin_ecmwf2(tmp ,nlon,nlat,flon,flat,dlon,dlat, 
C    &                  relhum(1-nbdy,1-nbdy,slot),plon,plat)   
      call interpug(tmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              relhum(1-nbdy,1-nbdy,slot),plon,plat,itype)   
c
c --- Downward shortwave radiation
      if (mnproc==1)print *,'Ecnc2 ssr'
      if (.false.) then
      else
         if (mnproc==1) then
          write(lp,'(a)') 'TODO - treat swfl in read_era40_inline'
         end if
      end if
c
c --- Final processing
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         svpair=satvap(airtmp(i,j,slot)+t0deg)
         svpdew=satvap(relhum(i,j,slot))
         rtmp= relhumid(svpair,svpdew,slp(i,j,slot))*0.01
         relhum(i,j,slot)=rtmp
c ---    These may be out of range for non-monotonic interpolation
         relhum(i,j,slot)  =max(0.0,min(relhum(i,j,slot),1.0))
         precip(i,j,slot)  =max(0.0,precip(i,j,slot))
         clouds(i,j,slot)  =max(0.0,min(clouds(i,j,slot),1.0))
         vapmix(i,j,slot) = rhtovpmix(relhum(i,j,slot),
     &                             airtmp(i,j,slot)+t0deg,slp(i,j,slot))
         slp   (i,j,slot)=     slp(i,j,slot)*0.01  ! Pa -> mBar (HPa)
      enddo
      enddo
      first=.false.
      END subroutine read_ecmwf_nc2


c --- Reads NCEP reanalyzed fields
      subroutine read_ncepr(dtime,slot)
      use mod_xc
      use netcdf
      use mod_common_ice
      use mod_hycom_nersc
      implicit none
      real*8, intent(in) :: dtime ! HYCOM time
      integer, intent(in) :: slot ! forcing slot to fill
c --- NCEP path
      character(len=200), save :: ncep_path='./NCEP/'
c
c --- Variable names
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
c 
c --- Hardcoded grid size
      integer, parameter :: nlon=192, nlat=94
c
c --- file names
      character(len=80) ::  ficec,fprec,ftair,ftcdc,fuwnd, fvwnd, 
     &   flmsk, fpres, fshum, fhgt, fssrd, fssrd2 
c
c --- Temporary fields on NCEP grid
      real, dimension(nlon,nlat) :: ncepu, ncepv, nceptaux, nceptauy,
     &   nceptmp, ncep_icec, ncep_lmsk, ncep_hgt,tmpdswrf,ncepwspd
      real :: nceplon(nlon)
c
      character(len=4) :: clyy, cyy
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::tmp_icec
      real :: satvapp, sathumid,wfact
      integer :: ndims,recdim,dimid,varid, ncid
      integer , dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      integer :: im1,ip1,jm1,jp1, maxrec, i,j,k,irec,nlon2,nlat2
      integer :: iyear, iday, ihour
      logical ::  ldump
      logical, save :: lfirst=.true.
      !integer, parameter :: itype=0 ! itype=1 means bicubic, 0=bilinear
      integer, parameter :: itype=1 ! itype=1 means bicubic, 0=bilinear
      logical, parameter :: lncepdswr=.false. ! use downwelling sw radiation
      include 'common_blocks.h'

c --- Calculate record number from year, ordinal day and hour
      call forday(dtime,yrflag, iyear,iday,ihour)
      irec  =  (iday-1)*4 + ihour/6 +1
      if (mnproc==1) 
     &  write(lp,'(a,i6,a,i4,x,i3,x,i2)')'NCEP data from record: ',
     &  irec, '  time= ',iyear,iday,ihour
C
      ldump=lfirst.or.mod(irec,40)==0
      if (mnproc==1.and.ldump) then
         print *
         write(lp,'(a)')'##############################################'
         write(lp,'(a)')'Reading NCEP synoptic forcing fields '
         write(lp,'(a,i5,a,i4,a,i3,a,i2)')'Record=',irec,' Year=',
     &             iyear,' Day=',iday,' hour=',ihour
      end if
C
c --- Which file to open
      write(clyy,'(i4.4)') iyear-1
      write(cyy ,'(i4.4)') iyear
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
      flmsk=trim(ncep_path)//'land.sfc.gauss.nc'
      fhgt =trim(ncep_path)//'hgt.sfc.gauss.nc' ! New - on gauss grid
c 
c --- Get dimensions, use air temperature file
      call ncerr(NF90_open(trim(ftair),NF90_NOCLOBBER,ncid))
      call ncerr(nf90_inq_dimid(ncid, "time", dimid))
      call ncerr(nf90_inquire_dimension(ncid, dimid, len=maxrec))
      if (lfirst) then
         call ncerr(nf90_inq_dimid(ncid, "lon", dimid))
         call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlon2))
         call ncerr(nf90_inq_dimid(ncid, "lat", dimid))
         call ncerr(nf90_inquire_dimension(ncid, dimid, len=nlat2))
         if (nlon/=nlon2.or.nlat/=nlat2) then
            if (mnproc==1) then
               write(lp,'(a)') 'Dimension mismatch in ERA40 data'
               write(lp,'(a,2i5)') 'nlon :',nlon,nlon2
               write(lp,'(a,2i5)') 'nlat :',nlat,nlat2
            end if
            call xcstop('(read_era40_inline)')
            stop '(read_era40_inline)'
         end if
c
c ---    NCEP Grid info
         allocate(gausslat(nlat))
         call ncerr(nf90_inq_varid(ncid, "lon", varid))
         call ncerr(nf90_get_var  (ncid, varid, nceplon))
         call ncerr(nf90_inq_varid(ncid, "lat", varid))
         call ncerr(nf90_get_var  (ncid, varid, gausslat))
C
C ---    NB: flat, dlat are just place holders. gausslat will be used
         flon=nceplon(1) ; flat=gausslat(1) ;
         dlon=nceplon(2)-nceplon(1) ; dlat=gausslat(2)-gausslat(1)
      end if
      call ncerr(nf90_close    (ncid))
c
c --- NCEP land mask. Not used yet
      if (mnproc==1.and.ldump) write(lp, *) 'Reading NCEP land mask'
      call rdncrec(flmsk,vlmsk,ncep_lmsk,nlon,nlat,1)
c
c --- NCEP geopotential height
      call rdncrec(fhgt,vhgt,ncep_hgt,nlon,nlat,1)
      if (mnproc==1.and.ldump) write(lp, *) 
     &   'Using NCEP geopot height (orography pressure correction'
c
c --- NCEP ice concentration
      if (mnproc==1.and.ldump) write(lp, *) 
     &   'Using NCEP ice conc (involved in relhum)'
      call rdncrec(ficec,vicec,ncep_icec,nlon,nlat,irec)
      call interpug(ncep_icec,nlon,nlat,flon,flat,dlon,dlat, 
     &              tmp_icec(1-nbdy,1-nbdy),plon,plat,itype,gausslat)   
c
c --- NCEP precipitation   kg/(m^2 s) -> m/month (rho=1000)
      if (mnproc==1.and.ldump) write(lp, *)  
     &   'Using NCEP synoptic precipitation '
      call rdncrec(fprec,vprec,nceptmp,nlon,nlat,irec)
      call interpug(nceptmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              precip(1-nbdy,1-nbdy,slot),plon,plat,itype,gausslat)
c --- convert prcipitation  from  kg/(m2*s) to m^3/m^2s
      precip(:,:,slot)=  precip(:,:,slot)*1e-3
c
c --- NCEP air temperature   K -> K
      if (mnproc==1.and.ldump) write(lp,'(a)')
     &   ' Using NCEP synoptic temperature'
      call rdncrec(ftair,vtair,nceptmp,nlon,nlat,irec)
      call interpug(nceptmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              airtmp(1-nbdy,1-nbdy,slot),plon,plat,itype,gausslat)
c
c --- NCEP clouds   % -> []  (0-1)
      if (mnproc==1.and.ldump) write(lp,'(a)')
     &      ' Using NCEP synoptic clouds '
      call rdncrec(ftcdc,vtcdc,nceptmp,nlon,nlat,irec)
      if (mnproc==1)print *,'max cc',maxval(nceptmp)
      call interpug(nceptmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              clouds(1-nbdy,1-nbdy,slot),plon,plat,itype,gausslat)
      clouds(:,:,slot)=0.01*clouds(:,:,slot)
      clouds(:,:,slot)=min(1.0,max(clouds(:,:,slot),0.))
c
c --- NCEP winds 
      if (mnproc==1.and.ldump) write(lp,'(a)') 
     &      ' Using NCEP synoptic winds '
      call rdncrec(fuwnd,vuwnd,ncepu,nlon,nlat,irec)
      call rdncrec(fvwnd,vvwnd,ncepv,nlon,nlat,irec)
c
c --- Calculate derived  wind stress and absolute windspeed on NCEP grid
c$OMP PARALLEL DO PRIVATE(j,i,wfact) &
c$OMP SCHEDULE(STATIC,jblk)
      do j=1,nlat
      do i=1,nlon
c
c ---    Wind speed
         ncepwspd(i,j)=sqrt(ncepu(i,j)*ncepu(i,j)+ 
     &                      ncepv(i,j)*ncepv(i,j))

c ---    Wind drag factor from Large&Pond 1981
C        cd_new=cd_lp81(wndspd(i,j,slot))
C        wndfac=(1.+sign(1.,wndspd(i,j,slot)-11.))*.5
C        cd_new=(0.49+0.065*wndspd(i,j,slot))*1.0e-3*wndfac+
C    &      cd*(1.-wndfac)
C
C        wfact=wndspd(i,j,slot)*airdns*cd_new
c ---    Drag coeff moved to function
         wfact=ncepwspd(i,j)*airdns*cd_lp81(ncepwspd(i,j))
         nceptaux(i,j)=ncepu(i,j)*wfact
         nceptauy(i,j)=ncepv(i,j)*wfact
      end do
      end do
c$OMP END PARALLEL DO
c
c
c --- Interpolate wind, speed, ++ to p-points
      call interpug(ncepu,nlon,nlat,flon,flat,dlon,dlat, 
     &           uwind(1-nbdy,1-nbdy,slot),PLON,PLAT,itype,gausslat)
      call interpug(ncepv,nlon,nlat,flon,flat,dlon,dlat, 
     &           vwind(1-nbdy,1-nbdy,slot),PLON,PLAT,itype,gausslat) 
      call interpug(nceptaux,nlon,nlat,flon,flat,dlon,dlat, 
     &           taux(1-nbdy,1-nbdy,slot) ,PLON,PLAT,itype,gausslat)
      call interpug(nceptauy,nlon,nlat,flon,flat,dlon,dlat, 
     &           tauy(1-nbdy,1-nbdy,slot) ,PLON,PLAT,itype,gausslat) 
      call interpug(ncepwspd,nlon,nlat,flon,flat,dlon,dlat, 
     &            wndspd(1-nbdy,1-nbdy,slot) ,PLON,PLAT,itype,gausslat) 
c
c --- Rotate winds and stresses to model grid 
c --- wndflg=1 : Rotate winds/stress to model u/v points
      if (wndflg==1) then
      call rotate2(uwind(1-nbdy,1-nbdy,slot),vwind(1-nbdy,1-nbdy,slot),
     &             plat (1-nbdy,1-nbdy)     ,plon (1-nbdy,1-nbdy  ),
     &             'l2m',keeppoint=.false.)
      call rotate2(taux (1-nbdy,1-nbdy,slot),tauy (1-nbdy,1-nbdy,slot),
     &             plat (1-nbdy,1-nbdy)     ,plon (1-nbdy,1-nbdy  ),
     &             'l2m',keeppoint=.false.)
c --- wndflg=2 : Rotate winds/stress to model p points
      else if (wndflg==2) then
      call rotate2(uwind(1-nbdy,1-nbdy,slot),vwind(1-nbdy,1-nbdy,slot),
     &             plat (1-nbdy,1-nbdy)     ,plon (1-nbdy,1-nbdy  ),
     &             'l2m',keeppoint=.true.)
      call rotate2(taux (1-nbdy,1-nbdy,slot),tauy (1-nbdy,1-nbdy,slot),
     &             plat (1-nbdy,1-nbdy)     ,plon (1-nbdy,1-nbdy  ),
     &             'l2m',keeppoint=.true.)
      end if

c --- NCEP relative humidity and pressure
      if (mnproc==1.and.ldump) write(lp,'(a)') 
     &      ' Using NCEP synoptic SLP + SHUM'
      call rdncrec(fpres,vpres,nceptmp,nlon,nlat,irec)
      nceptmp=0.01*nceptmp                  ! Convert from Pa to mbar 
      nceptmp=nceptmp+ncep_hgt*9.8*1.3*0.01 ! Simple correction for geopot. height
      call interpug(nceptmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              slp(1-nbdy,1-nbdy,slot),plon,plat,itype,gausslat)   

      if (ldump.and.mnproc==1) print *,'max/min geopot hgt adj:',
     &   maxval(ncep_hgt)*9.8*1.3*0.01,minval(ncep_hgt)*9.8*1.3*0.01

      if (mnproc==1.and.ldump) then
         write(lp,'(a)') ' Using NCEP synoptic specific humidity'
         call flush(lp)
      end if

c --- Read specific humidity and interpolate to grid
      call rdncrec(fshum,vshum,nceptmp,nlon,nlat,irec)
      call interpug(nceptmp,nlon,nlat,flon,flat,dlon,dlat, 
     &              relhum(1-nbdy,1-nbdy,slot),plon,plat,itype,gausslat)
c
c --- Convert to relative humidity and vapor mixing ration
      do j=2-nbdy,jj+nbdy-1
      do i=2-nbdy,ii+nbdy-1
c
c ---    !! Saturation Vapour pressure (sea level)
         satvapp=satvap(airtmp(i,j,slot))
c
c ---    ! Vapor mixing ratio (synrelhum is shum at this point)
         vapmix(i,j,slot)=relhum(i,j,slot) / (1.-relhum(i,j,slot))
c
c ---    ! Saturation humidity
         sathumid=humid(slp(i,j,slot)*100.,satvapp)
c
c ---    ! Relative humidty
C ---    TODO: Check this...
         relhum(i,j,slot)=spechum_to_relhum(relhum(i,j,slot),sathumid)
         relhum(i,j,slot)=min(max(0.0,relhum(i,j,slot)),1.)
c
c ---    Air temperature must be in celsius on exit
         airtmp(i,j,slot)=airtmp(i,j,slot)-t0deg
      end do
      end do

C --- NCEP downwelling shortwave radiation
      if (lncepdswr) then
         if (mnproc==1.and.ldump)
     &      write(lp,'(a)')' Using NCEP synoptic DSWRF '
         call rdncrec(fssrd,vssrd,tmpdswrf,nlon,nlat,irec)
c
c ---    ssrd is averaged over 6 next hrs, read last ssrd field as well
         if (irec==1)  then
            call ncerr(NF90_open(trim(fssrd2),NF90_NOCLOBBER,ncid))
            call ncerr(nf90_inq_dimid(ncid, "time", dimid))
            call ncerr(nf90_inquire_dimension(ncid, dimid, len=maxrec))
            call ncerr(NF90_close(ncid))
            call rdncrec(fssrd2,vssrd,nceptmp,nlon,nlat,maxrec)
         else 
            call rdncrec(fssrd2,vssrd,nceptmp,nlon,nlat,irec-1)
         end if
         nceptmp(:,:)=nceptmp(:,:)+tmpdswrf
         nceptmp=nceptmp*.5
         call interpug(nceptmp,nlon,nlat,flon,flat,dlon,dlat, 
     &                 swflx(1-nbdy,1-nbdy,slot),
     &                 plon,plat,itype,gausslat)   
      end if
      if (mnproc==1.and.ldump)  then
         print '(a)',
     &      '####################################################'
      end if
      lfirst=.false.
      END subroutine read_ncepr



c --- #######################################################################
c --- Specific humidity from vapor pressure
c --- #######################################################################
c ... TODO: reconsider name
      real function humid(px,es)
c --- Function calculates the specific humidity based on the total presure and 
c --- the vapor ! pressure (ie partial pressure of vapor). Nondimensional. Both
c --- Input is in Pa (doesnt really matter as long as they have the same unit)
c --- Based on formulas in Gill (1982).
      implicit none
      real,intent(in) :: px,es
      humid=.622*es/(px-.378*es)
      end function humid

      real function spechum_to_relhum(q,qs)
c --- Function converts from specific humidity to relative humidity. 
      implicit none
      real, intent(in) :: q   ! humidity
      real, intent(in) :: qs  ! saturation humidity
      real :: aa,bb
      !aa=q *(1-qs)
      !bb=qs*(1-q )
      !spechum_to_relhum=aa/bb
      spechum_to_relhum=q/qs
      end function




C!   real function vapp(ax,bx,tx)
C!   implicit none
C!     real, intent(in) :: ax,bx,tx
C!      vapp=611.*10.**(ax*(tx-273.16)/(tx-bx))
C!   end function vapp
C!
C!   real function humid(px,ex)
C!   implicit none
C!     real, intent(in) :: px,ex
C!      humid=.622*ex/(px-.378*ex)
C!   end function humid



c
      subroutine read_ecmwf(dtime,slot)
      use mod_xc
      use mod_common_ice
      use mod_hycom_nersc
      implicit none
      real*8, intent(in) :: dtime
      integer, intent(in) :: slot
c --- Input Data dimensions (to be decided)
      integer :: nxd, nyd
c
      integer mo,im1,ip1,jm1,jp1
      real fac,w4,wfact,rtmp
      logical ex
      real lonref,dlon
      real latref,dlat
      real, dimension(:, :), allocatable :: 
     &   tmp, tmpu, tmpv, tmpw, tmptx, tmpty
      integer i,j,irec,k, ihour, iday, iyear
      integer imax(2),imin(2)
      real maxtem,maxdew,maxmsl,maxuuu,maxvvv
      real mintem,mindew,minmsl,minuuu,minvvv
      character(len=9) path0
      character(len=4) cyy
      real wndfac,cd_new, svpair, svpdew
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: 
     &   cprecip,cclouds
      include 'common_blocks.h'
c
      call forday(dtime,yrflag, iyear,iday,ihour)
      irec=(iday-1)*4+ihour/6+1
      write(cyy  ,'(i4.4)') iyear
c --- Set up for forcing
      if (rforce=='ecmwf') then
         path0(1:9)='./Ecmwfr/'
         nxd=320
         nyd=160
      else if (rforce=='metno') then
         path0(1:9)='./Met.no/'
         nxd=720
         nyd=361
      else
         write(lp,*) 'read_ecmwf only understands "metno" or "ecmwf" '
         write(lp,*) 'forcing ID'
         call flush(lp)
         call xcstop('(read_ecmwf)')
         stop '(read_ecmwf)'
      end if
c
      allocate( tmp  (nxd,nyd) )
      allocate( tmpu (nxd,nyd) )
      allocate( tmpv (nxd,nyd) )
      allocate( tmpw (nxd,nyd) )
      allocate( tmptx(nxd,nyd) )
      allocate( tmpty(nxd,nyd) )
c
c --- read ecmwf headers
      inquire(file=trim(path0)//'ecmwf_'//cyy//'.hdr',exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(lp,*),
     &     'read_ecmwf:  file does not exist :',trim(path0),'ecmwf_',
     &     cyy,'.hdr'
         call xcstop ('(read_ecmwf)')
         stop '(read_ecmwf)'
      endif
      open(10,file=trim(path0)//'ecmwf_'//cyy//'.hdr',status='old')
         read(10,'(t14,i5)')i
         read(10,'(t14,i5)')j
         read(10,'(t14,2f9.4)')lonref,dlon
         read(10,'(t14,2f9.4)')latref,dlat
      close(10)

c --- Read input files
c
c --- Read one precipitation field from climatology. cprecip is kept
c --- locally.
      if (mnproc==1)print *,'Ecmwf precip (climatology)'
      if (hf_cinit) then
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc0),'precip',hf_mc0)
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc1),'precip',hf_mc1)
      else if (hf_cread) then
       cprecip(:,:,hf_lc0)=cprecip(:,:,hf_lc1)
       call hf_readclm(cprecip(1-nbdy,1-nbdy,hf_lc1),'precip',hf_mc1)
      end if
c --- Interpolate climatology onto hf forcing field times
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         precip(i,j,slot) = hf_wc0*cprecip(i,j,hf_lc0) +
     &                   hf_wc1*cprecip(i,j,hf_lc1) 
      end do
      end do
c
c --- Read one cloud field from climatology. ccloud is kept
c --- locally.
      if (mnproc==1)print *,'Ecmwf clouds (climatology)'
      if (hf_cinit) then
       call hf_readclm(cclouds(1-nbdy,1-nbdy,hf_lc0),'clouds',hf_mc0)
       call hf_readclm(cclouds(1-nbdy,1-nbdy,hf_lc1),'clouds',hf_mc1)
      else if (hf_cread) then
       cclouds(:,:,hf_lc0)=cclouds(:,:,hf_lc1)
       call hf_readclm(cclouds(1-nbdy,1-nbdy,hf_lc1),'clouds',hf_mc1)
      end if
c --- Interpolate climatology onto hf forcing field times
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         clouds(i,j,slot) = hf_wc0*cclouds(i,j,hf_lc0) +
     &                   hf_wc1*cclouds(i,j,hf_lc1) 
      end do
      end do
c
c --- air temperature
      call getuf(trim(path0)//'ecmwf_'//cyy//'_tem.uf',
     &   irec,100.,400., tmp,nxd,nyd)
      if (mnproc==1)print *,'back again'
      call bilin_ecmwf2(tmp,nxd,nyd,lonref,latref,dlon,dlat,
     &                  airtmp(1-nbdy,1-nbdy,slot),plon,plat)
c 
c --- dewpoint temperature
      call getuf(trim(path0)//'ecmwf_'//cyy//'_dew.uf',
     &   irec,100.,400., tmp,nxd,nyd)
      call bilin_ecmwf2(tmp,nxd,nyd,lonref,latref,dlon,dlat,
     &                  relhum(1-nbdy,1-nbdy,slot),plon,plat)
c
c --- SLP             
      call getuf(trim(path0)//'ecmwf_'//cyy//'_msl.uf',
     &   irec,100.,1400., tmp,nxd,nyd)
      call bilin_ecmwf2(tmp,nxd,nyd,lonref,latref,dlon,dlat,
     &                  slp(1-nbdy,1-nbdy,slot),plon,plat)
c
c --- Winds
      call getuf(trim(path0)//'ecmwf_'//cyy//'_uuu.uf',
     &   irec,-100.,100., tmpu,nxd,nyd)
      call getuf(trim(path0)//'ecmwf_'//cyy//'_vvv.uf',
     &   irec,-100.,100., tmpv,nxd,nyd)
c
c --- Wind speed on original (colocated) grid
      do j=1,nyd
      do i=1,nxd
         tmpw(i,j)=sqrt(tmpu(i,j)*tmpu(i,j)+ tmpv(i,j)*tmpv(i,j))
         wndfac=(1.+sign(1.,tmpw(i,j)-11.))*.5
         cd_new=(0.49+0.065*tmpw(i,j))*1.0e-3*wndfac+cd*(1.-wndfac)
         wfact=tmpw(i,j)*airdns*cd_new
         tmptx(i,j)=tmpu(i,j)*wfact
         tmpty(i,j)=tmpv(i,j)*wfact
      end do
      end do
c
c --- Interpolate to u/v points
      call bilin_ecmwf2(tmpu ,nxd,nyd,lonref,latref,dlon,dlat, 
     &                 uwind(1-nbdy,1-nbdy,slot),ulon,ulat)   
      call bilin_ecmwf2(tmpv ,nxd,nyd,lonref,latref,dlon,dlat,
     &                 vwind(1-nbdy,1-nbdy,slot),vlon,vlat)   
      call bilin_ecmwf2(tmptx,nxd,nyd,lonref,latref,dlon,dlat, 
     &                  taux(1-nbdy,1-nbdy,slot),ulon,ulat)   
      call bilin_ecmwf2(tmpty,nxd,nyd,lonref,latref,dlon,dlat, 
     &                  tauy(1-nbdy,1-nbdy,slot),vlon,vlat)   
      call bilin_ecmwf2(tmpw,nxd,nyd,lonref,latref,dlon,dlat, 
     &                  wndspd(1-nbdy,1-nbdy,slot),plon,plat)   
c
c --- Rotate winds/stress to model grid 
      call rotate(uwind(1-nbdy,1-nbdy,slot), vwind(1-nbdy,1-nbdy,slot), 
     &            plat (1-nbdy,1-nbdy)  , plon (1-nbdy,1-nbdy  ),
     &            idm+2*nbdy,jdm+2*nbdy,'l2m')
      call rotate(taux(1-nbdy,1-nbdy,slot), tauy(1-nbdy,1-nbdy,slot),  
     &            plat(1-nbdy,1-nbdy  ), plon(1-nbdy,1-nbdy  ),
     &            idm+2*nbdy,jdm+2*nbdy,'l2m')
c
c --- Final processing
      do j=1-nbdy,jj+nbdy
      do i=1-nbdy,ii+nbdy
         svpair=satvap(airtmp(i,j,slot))
         svpdew=satvap(relhum(i,j,slot)) ! relhum contains dewp temp here
         !relhum(i,j,slot)= relhumid(svpair,svpdew,slp(i,j,slot)*100.)*0.01
         rtmp= relhumid(svpair,svpdew,slp(i,j,slot)*100.)*0.01
         relhum(i,j,slot)=rtmp
         relhum(i,j,slot)  =max(0.0,min(relhum(i,j,slot),1.0))
         vapmix(i,j,slot) = rhtovpmix(relhum(i,j,slot),
     &                             airtmp(i,j,slot),slp(i,j,slot)*100.)
         airtmp(i,j,slot) = airtmp(i,j,slot)-t0deg ! K  -> C
      enddo
      enddo
c
c --- Diagnostic output
      write(lp,'(a,i6,a,i4,x,i3,x,i2,a,i3,i3)')
     &   'Ecmwfr data from record: ',irec,'  rt= ',iyear,iday-1,
     &   ihour
      call flush(lp)
c --- 
      deallocate( tmp  )
      deallocate( tmpu )
      deallocate( tmpv )
      deallocate( tmpw )
      deallocate( tmptx)
      deallocate( tmpty)
      end subroutine read_ecmwf
c
c
c --- "Helper" routine used in read_ecmwf
      subroutine getuf(fname,irec,vmin,vmax,fld,nx,ny)
      implicit none
      character(len=*), intent(in)  :: fname
      integer         , intent(in)  :: irec,nx,ny
      real            , intent(in)  :: vmax,vmin
      real            , intent(out) :: fld(nx,ny)
      real*4, dimension(nx,ny) ::  tmpin
      logical :: ex
      integer :: imin(2),imax(2),j 
      real :: minfld, maxfld
c
c --- Check for existence
      inquire(file=fname,exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(lp,*), 
     &      'read_ecmwf:  file does not exist :',fname
         call xcstop ('(read_ecmwf)')
         stop '(read_ecmwf)'
      endif
c
c --- Open and read (all nodes!). To improve read and scatter to all nodes
      inquire(iolength=j) tmpin
      open(21,file=trim(fname),form='unformatted',access='direct',
     &        recl=j)
      read(21,rec=irec)tmpin; 
      close(21)
      fld=tmpin
c
c --- Consistency check
      imax(1:2)=maxloc(fld); 
      if (mnproc==1) print *,'getuf - maxloc',imax
      maxfld=fld(imax(1),imax(2))
      if (maxfld > vmax)  then
         if (mnproc==1) print '(a,g13.5,2i5)',
     &      'read_ecmwf: Max value and location : ',maxfld,imax(:)
         call xcstop('(read_ecmwf: max fld too large)')
         stop 'Read_ecmwf: max fld too large'
      end if
      imin(1:2)=minloc(fld); 
      if (mnproc==1) print *,'getuf - minloc',imin
      minfld=fld(imin(1),imin(2))
      if (minfld < vmin)  then
         if (mnproc==1) print '(a,g13.5,2i5)',
     &      'read_ecmwf: Min value and location : ',minfld,imin(:)
         call xcstop('(read_ecmwf: min fld too small)')
         stop 'Read_ecmwf: min fld too small'
      end if
      end subroutine
c
c
      subroutine hf_readclm(fld,cfld,mnth)
      use mod_xc
      use mod_za
      implicit none
c --- Read one Field from climatology and put it into 
c --- fld.
      character(len=*), intent(in) :: cfld
      integer         , intent(in) :: mnth
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
      character preambl(5)*79
      character (len=80) :: cline
      integer :: i,ios,k,imnth
      real :: hmina, hmaxa, hminb, hmaxb
      include 'common_blocks.h'

c
      if (mnproc==1) write(lp,'(a)') 'Reading '//cfld
     &   //' climatology(hf)'
      call zaiopf(trim(flnmfor)//'forcing.'//cfld//'.a', 'old', 333)
      if (mnproc==1) then
         open (unit=uoff+333,
     &      file=trim(flnmfor)//'forcing.'//cfld//'.b',
     &      status='old', action='read')
         read (uoff+333,'(a79)') preambl
      end if
      !call prembl_print(preambl)
      do k=1,mnth-1
         call zaiosk(333)
         call zagetc(cline,ios,uoff+333)
      end do
      call zagetc(cline,ios,uoff+333)
      i = index(cline,'=')
      read (cline(i+1:),*) imnth,hminb,hmaxb
      if  (imnth/=mnth) then
         write(lp,'(a)') 'Month mismatch is this montlhy data?'
         call xcstop('(hf_readclm)')
         stop '(hf_readclm)'
      end if
      call zaiord(fld,ip,.false., hmina,hmaxa, 333)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,i3 / a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    cfld//' error - .a and .b files not consistent:',
     &    'iunit = ',333,
     &    cline,
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif !1st tile
        call xcstop('(hf_readclm)')
               stop '(hf_readclm)'
      endif
      call zaiocl(333)
      if (mnproc==1) close(333+uoff)
      end subroutine
c
#endif /* FORCING_INLINE */
c
c --- Dumps diagnostics. 
      subroutine inline_diag(prefix,dtime,slot)
      use mod_hycom_nersc
      use netcdf
      implicit none
      characteR(len=*), intent(in) :: prefix
      real*8 , intent(in) :: dtime
      integer, intent(in) :: slot
      character(len=80) :: ncfil
      integer :: iyear, iday, ihour
      real    ,dimension(itdm,jtdm) :: gfld
      include 'common_blocks.h'
c
      call forday(dtime, yrflag, iyear,iday,ihour)
      write(ncfil,'(a,i4.4,"_",i3.3,"_",i2.2,".nc")') 
     &      prefix,iyear,iday,ihour
c   
      call xcaget(gfld,depths,0)
      call ncdraft(ncfil,gfld,'depths','clobber')
c   
      call xcaget(gfld,plon,0)
      call ncdraft(ncfil,gfld,'plon','')
c   
      call xcaget(gfld,plat,0)
      call ncdraft(ncfil,gfld,'plat','')
c   
      call xcaget(gfld,airtmp(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'airtmp','')
c
      call xcaget(gfld,clouds(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'clouds','')
c
      call xcaget(gfld,precip(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'precip','')
c
      call xcaget(gfld,relhum(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'relhum','')
c
      call xcaget(gfld,slp(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'slp','')
c
      call xcaget(gfld,vapmix(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'vapmix','')
c
      call xcaget(gfld,vwind(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'vwind','')
c
      call xcaget(gfld,uwind(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'uwind','')
c
      call xcaget(gfld,wndspd(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'wndspd','')
c
      call xcaget(gfld,taux(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'taux','')
c
      call xcaget(gfld,tauy(1-nbdy,1-nbdy,slot),0)
      call ncdraft(ncfil,gfld,'tauy','')
      end subroutine
c
c
c --- convert relative humidity to vapor mixing ratio
      real function rhtovpmix(relhum,t,slp)
      implicit none
      real, intent(in) :: relhum ! relative humidity (range 0-1)
      real, intent(in) :: t      ! air temperature [K]
      real, intent(in) :: slp    ! sea level pressure [Pa]
      rhtovpmix = relhum*satvap(t)/slp
      end function
c
c
c --- This function calculates the saturation vapour pressure
c --- [Pa] from the temperature [deg K].
c --- Modified: Anita Jacob, June '97
c --- Input: t: temperature [deg K]
c --- Output: satvap: saturation vapour pressure at temp. t
c --- es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual
      real function satvap(t)
      implicit none
      real, intent(in):: t
      real :: aa,bb,c1,c3,c4,t00,cc
      data c1/610.78/,t00/273.16/
c
      if (t < t00) then
         c3 = 21.875
         c4 = 7.66
      else
         c3 = 17.269
         c4 = 35.86
      endif
c
      aa = c3 * (t - t00)
      bb = t - c4
      cc=aa/bb
      if (cc < -20.0) then
         satvap=0.0
      else
         satvap = c1 * exp(aa/bb)
      endif
      end function satvap
c
c
      real function relhumid(sva,svd,msl)
c --- This routine calculates the relative humidity by the 
c --- dew point temperature and the mean sea level pressure.
c --- Modified: Anita Jacob, June '97
c
c --- Input: sva: saturatn vapour press at air temp [K]
c --- svd: saturatn vapour press at dew pt temp [K]
c --- msl: pressure at mean sea level [Pa]
c --- Output: relhumid: Relative Humidity
c
c --- We use the Tetens formula:
c --- es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual
c ---              es(Tdew)        p - es(Tair)
c --- RH = 100 *  -----------   *  ------------
c ---             p - es(tdew)       es(Tair)
      implicit none
      real, intent(in):: sva
      real, intent(in):: svd
      real, intent(in):: msl
      real :: aaa,bbb
      aaa=msl - svd
      aaa = svd/aaa
      bbb = (msl - sva)/sva
      relhumid = 100. * aaa * bbb
      end function relhumid
c
c
c
c --- Routines checks error value of netcdf operation. Exits
c --- if an error occured
      subroutine ncerr(errcode)
      use mod_xc 
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      integer :: errcode2
      real a(1)
      a(1)=errcode
      call xcastr(a,1)
      errcode2=int(a(1))
      if (errcode2/=NF90_NOERR) then
         if (mnproc==1) then
            write(lp,'(a)') NF90_STRERROR(errcode2)
            print *,errcode2
         end if
         stop '(ncerr)'
         call xcstop('(ncerr)')
      end if
      end subroutine

Cc --- -----------------------------------------------------------------
Cc --- "interpug" - interpolation routine for "uniformly spaced" (in lon
Cc --- lat space) and global data sets which are periodic in longitude. 
Cc --- Interpolates data from data grid (old) to input grid (new) with 
Cc --- positions specified by newlon, newlat
Cc ---
Cc --- Interpolation methods:
Cc --- itype == 0 : bilinear interpolation. Continous, but discontinuous
Cc --- 1st derivatives. Always monotonic
Cc ---
Cc --- itype == 1 : bicubic interpolation. Continous 0-th and 1st order 
Cc --- derivatives. 2nd order cross derivative continous at corners. Not 
Cc --- Monotonic
Cc --- -----------------------------------------------------------------
C      subroutine interpug(old,onx,ony,olonref,olatref,odlon,odlat,
C     &                   new,newlon,newlat,itype)   
C      use mod_xc
C      implicit none
Cc --- Dims of old (data) grid, as well as reference lon/lat and grid increment
C      integer, intent(in) :: onx, ony,itype 
C      real, intent(in)    :: olonref, olatref    
C      real, intent(in)    :: odlon  , odlat
C      real, intent(in)    :: old(onx,ony)! old grid
Cc --- Longitude/Latitude of new grid, as well as new (interpolated) values
C      real, intent(in), dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
C     &   newlon, newlat
C      real, intent(out)   :: new   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)  
Cc
C      integer i,j,ipos,jpos,ia,ib,ja,jb,ib2,jb2,maxlati,minlati,numerr
C      real    aa,bb,maxlat,lon,minlat,rnumerr, newlon2
C      real :: f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,
C     &        fxy00,fxy10,fxy01,fxy11, a1,a2,a3,a4
C      real :: rhs(16), coeffs(16)
Cc
Cc
Cc --- Set "minimumlat" and "maximumlat" index..
C      maxlat=olatref+(ony-1)*odlat
C      maxlat=max(olatref,maxlat) ; ! For odlat<0
C      minlat=olatref+(ony-1)*odlat
C      minlat=min(minlat,olatref)   ! For odlat<0
C      if (olatref<olatref+(ony-1)*odlat) then
C         minlati=1
C         maxlati=ony
C      else
C         minlati=ony
C         maxlati=1
C      end if
Cc
Cc --- Start interpolation
C      numerr=0
CCC$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,jb,lon,ib, newlon2,
CCc&OMP& ia,ja, ib2, jb2, f00,f10,f01,f11,fx00,fx10,fx01,fx11,
CCc&OMP& fy00,fy10,fy01,fy11,fxy00,fxy10,fxy01,fxy11,
CCc&OMP& a1,a2,a3,a4,rhs,coeffs)
CCC$OMP&SHARED (lerr) SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
C      do j=1-nbdy,jdm+nbdy
C      do i=1-nbdy,idm+nbdy
Cc ---    New index in old data. New point is between
Cc ---    ipos, ib and jpos, jb
C         newlon2 = mod(newlon(i,j)+360.d0,360.d0)
C         ipos =int((newlon2-olonref)/odlon+1.0d0)
C         ia =mod(onx+ipos-2,onx)+1    ! Periodic
C         ib =mod(ipos,onx)+1          ! Periodic
C         ib2=mod(ib  ,onx)+1          ! Periodic
C         lon=olonref+(ipos-1)*odlon
Cc
C         jpos =int((newlat(i,j)-olatref)/odlat+1.0d0)
C         if (jpos>ony .or. jpos < 1 ) then
C            write(lp,'(a,2i5,2f10.2)') 'jpos error ',
C     &         jpos,ony,newlat(i,j),olatref+odlat*(ony-1)
C            numerr=numerr+1
C         endif
C         jpos =min(max(1,jpos  ),ony)
C         ja   =min(max(1,jpos-1),ony)
C         jb   =min(max(1,jpos+1),ony)
C         jb2  =min(max(1,jpos+2),ony)
CC        if (j+j0==250) print '(a,5i4)','i',i,ia,ipos,ib,ib2
Cc
Cc ---    Grid distance new point -> ipos, jpos
C         aa=(newlon2 - olonref-float(ipos-1)*odlon)/odlon
C         bb=(newlat(i,j) - olatref-float(jpos-1)*odlat)/odlat
Cc
Cc ---    Catch errors - but dont stop until after loop
C         if ((aa > 1.0).or.(aa < 0.0)) then
C            !write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon(i,j),lon+odlon
C            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2,lon+odlon
C            print *,'bilin_ecmwf2: invalid aa',aa
C            numerr=numerr+1
C         endif
C         if ((bb > 1.0).or.(bb < 0.0)) then
C            print *,'bilin_ecmwf2: invalid bb',bb
C            numerr=numerr+1
C         endif
Cc
Cc ---   Set up bilinear weights
C        if (itype==0) then
Cc ---      Bilinear weights
C           a1=(1.0-aa)*(1.0-bb)
C           a2=aa*(1.0-bb)
C           a3=aa*bb
C           a4=(1.0-aa)*bb
Cc ---      New data value
C           new(i,j) = a1*old(ipos,jpos)+a2*old(ib  ,jpos)+
C     &                a3*old(ib  ,jb  )+a4*old(ipos,jb  ) 
Cc ---   TODO: most of this can be re-used if (ipos,jpos) hasnt changed
Cc ---   Set up function and derivatives at ecmwf nodes
C        elseif (itype==1) then
C          f00  =old(ipos,jpos)
C          f10  =old(ib  ,jpos)
C          f01  =old(ipos,jb  )
C          f11  =old(ib  ,jb  )
Cc ---     X derivative with gridspacing 1 - TODO: modify with latitude
C          fx00 = 0.5*(old(ib  ,jpos) - old(ia  ,jpos))/odlon
C          fx10 = 0.5*(old(ib2 ,jpos) - old(ipos,jpos))/odlon
C          fx01 = 0.5*(old(ib  ,jb  ) - old(ia  ,jb  ))/odlon
C          fx11 = 0.5*(old(ib2 ,jb  ) - old(ipos,jb  ))/odlon
Cc ---     Y derivative with gridspacing 1 
C          fy00 = 0.5*(old(ipos,jb  ) - old(ipos,ja  ))/odlat
C          fy10 = 0.5*(old(ib  ,jb  ) - old(ib  ,ja  ))/odlat
C          fy01 = 0.5*(old(ipos,jb2 ) - old(ipos,jpos))/odlat
C          fy11 = 0.5*(old(ib  ,jb2 ) - old(ib  ,jpos))/odlat
Cc ---     Cross derivative with gridspacing 1  - TODO: modify with latitude
C          fxy00=0.25*( old(ib ,jb )-old(ib ,ja  )-
C     &                (old(ia  ,jb )-old(ia  ,ja  )))/(odlon*odlat)
C          fxy10=0.25*( old(ib2,jb )-old(ib2,ja  )-
C     &                (old(ipos,jb )-old(ipos,ja  )))/(odlon*odlat)
C          fxy01=0.25*( old(ib ,jb2)-old(ib ,jpos)-
C     &                (old(ia  ,jb2)-old(ia  ,jpos)))/(odlon*odlat)
C          fxy11=0.25*( old(ib2,jb2)-old(ib2,jpos)-
C     &                (old(ipos,jb2)-old(ipos,jpos)))/(odlon*odlat)
C          rhs=(/f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,
C     &        fxy00,fxy10,fxy01,fxy11/)
Cc ---     Solve matrix for cubic coeffs
C          coeffs=cubiccoeff(rhs)
Cc ---     Calculate solution
C          new(i,j)=cubicsol(coeffs,aa,bb)
C        end if
C      end do
C      end do
C      rnumerr=numerr ;call xcmaxr(rnumerr); numerr=int(rnumerr)
C      if (numerr>0) then
C         if (mnproc==1)write(lp,'(a)')'Error(s) occured in interpglob..'
C         call xcstop('(interpglob)')
C         stop '(interpglob)'
C      end if
C      end subroutine interpg
C
C
C
C
C
Cc
Cc --- Routine calculates bilinear weights for interpolation of field
Cc --- "old" to "new" locations. Also does the interpolation. A Periodic
Cc --- "old" grid is assumed
Cc
C      subroutine bilin_ecmwf2(old,onx,ony,olonref,olatref,odlon,odlat,
C     &                       new,newlon,newlat)   
C      use mod_xc
C      implicit none
C      integer, intent(in) :: onx         ! x-dimension of old field
C      integer, intent(in) :: ony         ! y-dimension of old field
C      real, intent(in)    :: old(onx,ony)! old grid
C      real, intent(in)    :: olonref     ! Lon - reference point
C      real, intent(in)    :: olatref     ! Lat - reference point
C      real, intent(in)    :: odlon       ! Lon - grid spacing in old grid
C      real, intent(in)    :: odlat       ! Lat - grid spacing in old grid
Cc
C      real, intent(out)   :: new   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)  ! New interpolated field
C      real, intent(in)    :: newlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! Longitudes for new grid
C      real, intent(in)    :: newlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! Latitudes for new grid
Cc
C      integer i,j,ia,ib,ja,jb,ifalse,maxlati,minlati, numerr
C      integer ipos    !  index of i-pivot grid point in old grid
C      integer jpos    !  index of j-pivot grid point in old grid
C      real    aa,bb,a1,a2,a3,a4,maxlat,lon,minlat,rnumerr, newlon2
Cc
C      maxlat=olatref+(ony-1)*odlat
C      maxlat=max(olatref,maxlat) ; ! For odlat<0
C      minlat=olatref+(ony-1)*odlat
C      minlat=min(minlat,olatref)   ! For odlat<0
Cc
Cc --- Set "minimumlat" and "maximumlat" index..
C      if (olatref<olatref+(ony-1)*odlat) then
C         minlati=1
C         maxlati=ony
C      else
C         minlati=ony
C         maxlati=1
C      end if
Cc
Cc
Cc --- Start interpolation
C      numerr=0
CC$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,jb,lon,ib, newlon2,
Cc&OMP&                    aa,bb,a1,a2,a3,a4) 
CC$OMP&SHARED (lerr) SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
C      do j=1-nbdy,jdm+nbdy
C      do i=1-nbdy,idm+nbdy
Cc ---    New index in old data. New point is between
Cc ---    ipos, ipos+1 and jpos, jpos+1
C         newlon2 = mod(newlon(i,j)+360.d0,360.d0)
C         !ipos =int((newlon(i,j)-olonref)/odlon+1.0)
C         ipos =int((newlon2-olonref)/odlon+1.0d0)
C         ib=mod(ipos,onx)+1          ! Periodic
C         lon=olonref+(ipos-1)*odlon
C         !print *,newlon(i,j),newlon2
Cc
C         jpos =int((newlat(i,j)-olatref)/odlat+1.0d0)
C         if (jpos>ony .or. jpos < 1 ) then
C            write(lp,'(a,2i5,2f10.2)') 'jpos error ',
C     &         jpos,ony,newlat(i,j),olatref+odlat*(ony-1)
C            numerr=numerr+1
C         endif
C         jpos =min(max(1,jpos  ),ony)
C         jb   =min(max(1,jpos+1),ony)
C         !if (j+j0==250) print *,i,ipos,ib,jpos,newlon2
Cc
Cc ---    Grid distance new point -> ipos, jpos
C         aa=(newlon2 - olonref-float(ipos-1)*odlon)/odlon
C         bb=(newlat(i,j) - olatref-float(jpos-1)*odlat)/odlat
Cc
Cc ---    Catch errors - but dont stop until after loop
C         if ((aa > 1.0).or.(aa < 0.0)) then
C            !write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon(i,j),lon+odlon
C            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2,lon+odlon
C            print *,'bilin_ecmwf2: invalid aa',aa
C            numerr=numerr+1
C         endif
C         if ((bb > 1.0).or.(bb < 0.0)) then
C            print *,'bilin_ecmwf2: invalid bb',bb
C            numerr=numerr+1
C         endif
C
Cc ---    Bilinear weights
C         a1=(1.0-aa)*(1.0-bb)
C         a2=aa*(1.0-bb)
C         a3=aa*bb
C         a4=(1.0-aa)*bb
Cc ---    New data value
C         new(i,j) = a1*old(ipos,jpos)+a2*old(ib  ,jpos)+
C     &              a3*old(ib  ,jb  )+a4*old(ipos,jb  ) 
C      enddo
C      enddo
CC$OMP END PARALLEL DO
C      rnumerr=numerr
C      call xcmaxr_0o(rnumerr)
C      numerr=int(rnumerr)
C
C      if (numerr>0) then
C         if (mnproc==1) then
C            write(lp,'(a)') 'Error(s) occured in bilin_ecmwf2..'
C            call flush(lp)
C         end if
C         call xcstop('(bilin_ecmwf2)')
C         stop '(bilin_ecmwf2)'
C      end if
C      end subroutine bilin_ecmwf2
Cc
C      subroutine bicubic(old,onx,ony,olonref,olatref,odlon,odlat,
C     &                       new,newlon,newlat)   
C      use mod_xc
C      implicit none
C      integer, intent(in) :: onx         ! x-dimension of old field
C      integer, intent(in) :: ony         ! y-dimension of old field
C      real, intent(in)    :: old(onx,ony)! old grid
C      real, intent(in)    :: olonref     ! Lon - reference point
C      real, intent(in)    :: olatref     ! Lat - reference point
C      real, intent(in)    :: odlon       ! Lon - grid spacing in old grid
C      real, intent(in)    :: odlat       ! Lat - grid spacing in old grid
Cc
C      real, intent(out)   :: new   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)  ! New interpolated field
C      real, intent(in)    :: newlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! Longitudes for new grid
C      real, intent(in)    :: newlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! Latitudes for new grid
Cc
C      integer i,j,ia,ib,ja,jb, ib2, jb2
C      integer :: maxlati,minlati, numerr
C      integer ipos    !  index of i-pivot grid point in old grid
C      integer jpos    !  index of j-pivot grid point in old grid
C      real    aa,bb,maxlat,lon,minlat,rnumerr, newlon2
C      real :: f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,
C     &        fxy00,fxy10,fxy01,fxy11
C      real :: rhs(16), coeffs(16)
Cc
C      maxlat=olatref+(ony-1)*odlat
C      maxlat=max(olatref,maxlat) ; ! For odlat<0
C      minlat=olatref+(ony-1)*odlat
C      minlat=min(minlat,olatref)   ! For odlat<0
C
C      !call xcstop ('(bicubic needs more testing )')     
C
Cc
Cc --- Set "minimumlat" and "maximumlat" index..
C      if (olatref<olatref+(ony-1)*odlat) then
C         minlati=1
C         maxlati=ony
C      else
C         minlati=ony
C         maxlati=1
C      end if
Cc
Cc
Cc --- Start interpolation
C      numerr=0
CCC$OMP PARALLEL DO PRIVATE(i,j,ipos,jpos,jb,lon,ib, newlon2,
CCc&OMP& ia,ja, ib2, jb2, f00,f10,f01,f11,fx00,fx10,fx01,fx11,
CCc&OMP& fy00,fy10,fy01,fy11,fxy00,fxy10,fxy01,fxy11,
CCc&OMP& rhs,coeffs)
CCC$OMP&SHARED (lerr) SCHEDULE(STATIC,jblk) REDUCTION(+:numerr)
C      do j=1-nbdy,jdm+nbdy
C      do i=1-nbdy,idm+nbdy
Cc ---    New index in old data. New point is between
Cc ---    ipos, ipos+1 and jpos, jpos+1
C         newlon2 = mod(newlon(i,j)+360.d0,360.d0)
C         ipos =int((newlon2-olonref)/odlon+1.0d0)
C         ia =mod(onx+ipos-2,onx)+1          ! Periodic
C         ib =mod(ipos,onx)+1          ! Periodic
C         ib2=mod(ib  ,onx)+1          ! Periodic
C         lon=olonref+(ipos-1)*odlon
Cc
C         jpos =int((newlat(i,j)-olatref)/odlat+1.0d0)
C         if (jpos>ony .or. jpos < 1 ) then
C            write(lp,'(a,2i5,2f10.2)') 'jpos error ',
C     &         jpos,ony,newlat(i,j),olatref+odlat*(ony-1)
C            numerr=numerr+1
C         endif
C         jpos =min(max(1,jpos  ),ony)
C         ja   =min(max(1,jpos-1),ony)
C         jb   =min(max(1,jpos+1),ony)
C         jb2  =min(max(1,jpos+2),ony)
C         if (j+j0==250) print '(a,5i4)','i',i,ia,ipos,ib,ib2
C
Cc
Cc ---    Grid distance new point -> ipos, jpos
C         aa=(newlon2 - olonref-float(ipos-1)*odlon)/odlon
C         bb=(newlat(i,j) - olatref-float(jpos-1)*odlat)/odlat
Cc
Cc ---    Catch errors - but dont stop until after loop
C         if ((aa > 1.0).or.(aa < 0.0)) then
C            !write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon(i,j),lon+odlon
C            write(*,'(3i5,3f10.2)')i,j,ipos,lon,newlon2,lon+odlon
C            print *,'bilin_ecmwf2: invalid aa',aa
C            numerr=numerr+1
C         endif
C         if ((bb > 1.0).or.(bb < 0.0)) then
C            print *,'bilin_ecmwf2: invalid bb',bb
C            numerr=numerr+1
C         endif
Cc ---   TODO: most of this can be re-used if (ipos,jpos) hasnt changed
Cc ---   Set up function and derivatives at ecmwf nodes
C        f00  =old(ipos,jpos)
C        f10  =old(ib  ,jpos)
C        f01  =old(ipos,jb  )
C        f11  =old(ib  ,jb  )
Cc ---   X derivative with gridspacing 1 - TODO: modify with latitude
C        fx00 = 0.5*(old(ib  ,jpos) - old(ia  ,jpos))/odlon
C        fx10 = 0.5*(old(ib2 ,jpos) - old(ipos,jpos))/odlon
C        fx01 = 0.5*(old(ib  ,jb  ) - old(ia  ,jb  ))/odlon
C        fx11 = 0.5*(old(ib2 ,jb  ) - old(ipos,jb  ))/odlon
Cc ---   Y derivative with gridspacing 1 
C        fy00 = 0.5*(old(ipos,jb  ) - old(ipos,ja  ))/odlat
C        fy10 = 0.5*(old(ib  ,jb  ) - old(ib  ,ja  ))/odlat
C        fy01 = 0.5*(old(ipos,jb2 ) - old(ipos,jpos))/odlat
C        fy11 = 0.5*(old(ib  ,jb2 ) - old(ib  ,jpos))/odlat
Cc ---   Cross derivative with gridspacing 1  - TODO: modify with latitude
C        fxy00=0.25*(
C     &        old(ib ,jb )-old(ib ,ja  )-(old(ia  ,jb )-old(ia  ,ja  ))
C     &        )/(odlon*odlat)
C        fxy10=0.25*(
C     &        old(ib2,jb )-old(ib2,ja  )-(old(ipos,jb )-old(ipos,ja  ))
C     &        )/(odlon*odlat)
C        fxy01=0.25*(
C     &        old(ib ,jb2)-old(ib ,jpos)-(old(ia  ,jb2)-old(ia  ,jpos))
C     &        )/(odlon*odlat)
C        fxy11=0.25*(
C     &        old(ib2,jb2)-old(ib2,jpos)-(old(ipos,jb2)-old(ipos,jpos))
C     &        )/(odlon*odlat)
Cc ---   Testing
Cc       fx00 =0.
Cc       fx10 =0.
Cc       fx01 =0.
Cc       fx11 =0.
Cc       fy00 =0.
Cc       fy10 =0.
Cc       fy01 =0.
Cc       fy11 =0.
Cc       fxy00=0.
Cc       fxy10=0.
Cc       fxy01=0.
Cc       fxy11=0.
Cc ---   Form rhs
C        rhs=(/f00,f10,f01,f11,fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11,
C     &        fxy00,fxy10,fxy01,fxy11/)
Cc ---   Solve matrix for cubic coeffs
C        coeffs=cubiccoeff(rhs)
Cc ---   Calculate solution
C        new(i,j)=cubicsol(coeffs,aa,bb)
C      enddo
C      enddo
CC$OMP END PARALLEL DO
C      rnumerr=numerr
C      call xcmaxr_0o(rnumerr)
C      numerr=int(rnumerr)
C      !call xcstop('bicubic')
C
C      if (numerr>0) then
C         if (mnproc==1) then
C            write(lp,'(a)') 'Error(s) occured in bicubic..'
C            call flush(lp)
C         end if
C         call xcstop('(bicubic)')
C         stop '(bicubic)'
C      end if
C      end subroutine bicubic
C
C      function cubiccoeff(rhs)
C      implicit  none
C      real, dimension(16), intent(in) :: rhs
C      real, dimension(16) :: cubiccoeff
CC     real, dimension(16*16), parameter ::
CC    & invcb=(/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
CC    &         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
CC    &        -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
CC    &         2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
CC    &         0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
CC    &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
CC    &         0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
CC    &         0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
CC    &        -3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0,
CC    &         0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0,
CC    &         9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1,
CC    &        -6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1,
CC    &         2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
CC    &         0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0,
CC    &        -6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1,
CC    &         4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1/)
CC     real, dimension(16,16),parameter:: invcb2=transpose(reshape(invcb,(/16,16/)))
Cc --- invcb2 = reshape(invcb,(/16,16/)) TODO - fix this
C      !invcb2 = transpose(reshape(invcb,(/16,16/)))
Cc 
CC --- Matrix is sparse - so there is room for reducing the number
Cc --- of operations (i.e avoid matmul)
C      real, dimension(16*16), parameter ::
C     &invcb=(/ 1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4,
C     &         0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4,
C     &         0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4,
C     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4,
C     &         0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2,
C     &         0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2,
C     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2,
C     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2,
C     &         0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2,
C     &         0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2,
C     &         0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2,
C     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2,
C     &         0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1,
C     &         0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1,
C     &         0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1,
C     &         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1/)
C      real,parameter,dimension(16,16):: invcb2=reshape(invcb,(/16,16/))  
C      cubiccoeff=matmul(invcb2,rhs)
C      end function
C
C      real function cubicsol(coeffs,aa,bb)
C      implicit none
C      real, intent(in) :: aa,bb,coeffs(16)
C      real :: coeffs2(4,4)
C      integer :: i,j
C
C      coeffs2=reshape(coeffs,(/4,4/))
C      cubicsol=0.
C      do j=1,4
C      do i=1,4
Cc ---    Double check indices
C         cubicsol=cubicsol+coeffs2(i,j)*aa**(i-1)*bb**(j-1)
C      end do
C      end do
C      end function cubicsol
C
C
c
c
      subroutine era40_fix(vname,fld,nlon,nlat,flon,flat,dlon,dlat)
      use mod_xc
      implicit none
c
      character(len=*), intent(in)    :: vname
      real            , intent(in)    :: flon,flat,dlon,dlat
      integer         , intent(in)    :: nlon,nlat
      real            , intent(inout) :: fld(nlon,nlat)
c
c --- Variable names - ERA40 1.125 by 1.125 
      character(len=*), parameter :: vprec='TP'
c
      real :: latera40, weight,wrad,wrad2
      integer :: i,j
c
      if (vname==vprec) then
c ---    ERA40 precipitation is biased in the tropics (See Biasotti et al., J.Clim)
c ---    bias is approx 50 % so real precip is ~ 2/3 of precip fields
c ---    This puts a latitudinal weighting on precipitation, centered on the
c ---    mean ITCZ-position (~5 degrees N). TODO: Adjust to seasonal ITCZ
c ---    position
         do j=1,nlat
            latera40=flat+(j-1)*dlat
c
            ! Weight=1 at 5N, 0 at 30N, 25S
            weight=abs(latera40-5.)
            wrad=20. ! 0-weight distance from 5N
            wrad2=10. ! region around 5N where weight=1
            if (weight<wrad2) then
               weight=1.
            else if (weight<wrad) then
               weight=.5+.5*cos((weight-wrad2)*3.14156/(wrad-wrad2))
            else 
               weight=0.
            end if
c
            do i=1,nlon
               fld(i,j)=fld(i,j)*(1. - weight*0.33)
            end do
         end do
      else
         if (mnproc==1) print *,'No ERA40 fix for variable '//vname
         call xcstop('(era40_fix)')
         stop '(era40_fix)'
      end if
c
      END subroutine era40_fix
c
      real function cd_lp81(wspd)
c --- define WINDDRAG results in the use of a winddependent winddrag.
c --- compute air-water wind stresses  tau_x  and  tau_y  from the surface wind
c --- components. surface wind components have dim [m/s], and are defined on a
c --- c-grid;  airdns  has dim [kg/m^3], and  cd  is dimensionless.
c --- tau_x  and  tau_y  have dim [kg/(m s^2) = N/m^2]
c ---
c --- include wind dependent drag coefficient according to
c --- large & pond, jpo, 11, 324-336, 1981
c ---
c ---   10^3 cd = 1.2             for winds < 11 m/s
c ---           = .49+0.065*wnd   for winds > 11 m/s
      implicit none
      real, intent(in) :: wspd
      real :: wndfac
      wndfac=(1.+sign(1.,wspd-11.))*.5
      cd_lp81=(0.49+0.065*wspd)*1.0e-3*wndfac+cd*(1.-wndfac)
      end function
c
c --- Helper routine for reading fields from netcdf forcing
c --- all processes call this
      subroutine rdncrec(ncfil,vname,fld,nx,ny,irec)
      use netcdf
      integer         ,intent(in) :: nx,ny,irec
      character(len=*),intent(in) :: ncfil, vname
      real  , dimension(nx,ny), intent(out) :: fld
      real*4, dimension(nx*ny) :: bcvec
      integer :: ierr
c
      ierr=0
c --- one cpu only
      if (mnproc==1) ierr=rdncrec_1(ncfil,vname,fld, nx,ny,irec)
      call ncerr(ierr)
c --- TODO: broadcast results to all cpus. This operation is
c --- probably costly for large grids. Not-so-naive options might 
c --- be considered
      if (mnproc==1) bcvec=reshape(fld,(/nx*ny/))
      call xcastr4(bcvec,1)  ! New - r*4 version of hycom xcastr
      fld=reshape(bcvec,(/nx,ny/))
      end subroutine
c
c
c --- Helper routine for reading fields from netcdf forcing
c --- one process only
      integer function rdncrec_1(ncfil,vname,fld,nx,ny,irec)
      use netcdf
      implicit none
      integer         ,intent(in) :: nx,ny,irec
      character(len=*),intent(in) :: ncfil, vname
      real, dimension(nx,ny), intent(out) :: fld
      integer :: varid,ncid,ierr
      real    :: scfac, offset
c
      rdncrec_1=NF90_open(trim(ncfil),NF90_NOCLOBBER,ncid)
      if (rdncrec_1/=NF90_NOERR)  return
      rdncrec_1=nf90_inq_varid(ncid, vname, varid)
      if (rdncrec_1/=NF90_NOERR)  return
      rdncrec_1=nf90_get_var  (ncid, varid, fld,start=(/1,1,irec/))
      if (rdncrec_1/=NF90_NOERR)  return
c
c --- Also need to check for these attributes. But do not return error 
c --- if we can't locate them
      ierr=nf90_get_att  (ncid, varid, 'scale_factor',scfac)
      if (ierr/=NF90_NOERR) scfac=1.
      ierr=nf90_get_att  (ncid, varid, 'add_offset',offset)
      if (ierr/=NF90_NOERR) offset=0.
c
c --- Scale final field
      fld=fld*scfac+offset
c 
      rdncrec_1=nf90_close    (ncid)
      if (rdncrec_1/=NF90_NOERR)  return
      end function


      end module mod_forcing_nersc
