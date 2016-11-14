module mod_clim_atm
use mod_clim_ocn
! KAL - a small attempt to clean up messy code - moved climatology-related 
!     - routines in here. 
!
! -- read_clim_all  - reads all climatologies in one swoop and places them 
!                     in the "clm" variables of mod_forcing_nersc.
!                     which climatology to read is conrolled by "clmflag" in
!                     mod_forcing_nersc
! -- read_clim      - Reads a single climatology field. Called repeatedly
!                     by read_clim_all. Also controlled by "clmflag" in 
!                     mod_forcing_nersc
! -- nersc_clm_calc - Processes the climatology read by read_clim_all when
!                     clmflag is "old"
! -- era40_clm_calc - Processes the climatology read by read_clim_all when
!                     clmflag is "era40"
! -- clim_diag      - tecplot diagnostic of climatology - independent of clmflg
! -- syn_diag       - tecplot diagnostic of syn. forcing - independent of frcflg
!                   - This doesnt really belong here, but dunno where to put it
!
! NB - All these routines require that mod_grid and mod_forcing is initialized 


character(len=200), save :: path_era40_clim='./ERA40-Clim/'
character(len=200), save :: path_ncep_clim='./NCEP-Clim/'

! Temporal variables - used when reading climatology
real, allocatable, dimension(:,:,:)  :: &
   clmuwind, clmvwind, clmwndspd, clmairtmp, clmrelhum, clmprecip, &
   clmclouds, clmsst,  clmsss, clmtaux, clmtauy, clmvapmix, &
   clmradflx, clmshwflx, clmslp

contains


   subroutine init_clmvars_nersc()
   use mod_xc
   implicit none
   allocate(clmuwind (idm,jdm,2))
   allocate(clmvwind (idm,jdm,2))
   allocate(clmwndspd(idm,jdm,2))
   allocate(clmtaux  (idm,jdm,2))
   allocate(clmtauy  (idm,jdm,2))
   allocate(clmvapmix(idm,jdm,2))
   allocate(clmairtmp(idm,jdm,2))
   allocate(clmrelhum(idm,jdm,2))
   allocate(clmprecip(idm,jdm,2))
   allocate(clmclouds(idm,jdm,2))
   allocate(clmradflx(idm,jdm,2))
   allocate(clmshwflx(idm,jdm,2))
   allocate(clmslp   (idm,jdm,2))
   allocate(clmsss   (idm,jdm,2))
   allocate(clmsst   (idm,jdm,2))
   clmuwind=0.
   clmvwind=0.
   clmwndspd=0.
   clmtaux=0.
   clmtauy=0.
   clmvapmix=0.
   clmairtmp=0.
   clmrelhum=0.
   clmprecip=0.
   clmclouds=0.
   clmradflx=0.
   clmshwflx=0.
   clmsss=0.
   clmsst=0.
   end subroutine

   subroutine read_clim_all(iplace,mo)
   use mod_grid
   use mod_forcing_nersc
   implicit none
   integer, intent(in) :: iplace ! Slot to place data in climate arrays
   integer, intent(in) :: mo     ! Month to retrieve data from

   logical :: ex
   real*8 :: tmr0,tmr1
   real*8, external :: wtime
   logical, parameter :: timer=.false.
   logical, save :: first=.true.
   character(len=200) :: cenv
   integer*4, external :: system
   integer*4 :: ret


   if (clmflag=='old') then
      if (first) then
         print '(a)','Will look for "old" climatology under ./Data/'
      end if

      CALL read_clim(mo,clmuwind (:,:,iplace), &
             './Data/uwind.forc','',clmflag)       !m/s

      CALL read_clim(mo,clmvwind (:,:,iplace), &
             './Data/vwind.forc','',clmflag)       !m/s

      CALL read_clim(mo,clmwndspd(:,:,iplace), &
             './Data/wndab.forc','',clmflag)       !m/s

      CALL read_clim(mo,clmairtmp(:,:,iplace), &
             './Data/airtp.forc','',clmflag)       !kelvin

      CALL read_clim(mo,clmrelhum(:,:,iplace), &
             './Data/humid.forc','',clmflag)       !rh fraction 0-1

      CALL read_clim(mo,clmprecip(:,:,iplace), &
             './Data/precp.forc','',clmflag)       !mm/month precepitation

      CALL read_clim(mo,clmclouds(:,:,iplace), &
             './Data/cloud.forc','',clmflag)       !mm/month precepitation
 
      clmslp(:,:,iplace)=slp0

      call sss_sst_nodc(clmsst(:,:,iplace), clmsss(:,:,iplace),mo)


      call nersc_clm_calc(iplace,mo)
      print *


   else if (clmflag=='era40') then
      ! Check ERA40 path on first pass
      if (first) then
         call getenv('ERA40_CLIM_PATH',cenv)
         if (trim(cenv)/='') then ! prefer this path if present
            path_era40_clim=trim(cenv)//'/'
         end if
         ret=system('[ -d '//trim(path_era40_clim)//' ]')
         if (ret /=0 ) then
            print *
            print *,'The directory  '//trim(path_era40_clim)//' does not exist.'
            print *,'Make sure that '
            print *,' - You have linked the ERA40-Clim data into this catalogue, and '
            print *,'   that the variable ERA40_CLIM_PATH is empty'
            print *,' - OR, the variable ERA40_CLIM_PATH is set to the location of '
            print *,'   the ERA40-Clim data'
            call exit(1)
         end if
      end if

      !tmr0=wtime()
      CALL read_clim(mo,clmuwind (:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_U10M_sfc.nc','U10M_sfc',clmflag)       !m/s
      !tmr1=wtime()
      !print *,'One clm field read took ',tmr1-tmr0

      CALL read_clim(mo,clmvwind (:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_V10M_sfc.nc','V10M_sfc',clmflag)       !m/s

      CALL read_clim(mo,clmwndspd(:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_wndspd.nc'  ,'wndspd',clmflag)         !m/s

      CALL read_clim(mo,clmairtmp(:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_T2M_sfc.nc' ,'T2M_sfc',clmflag)       !kelvin

      CALL read_clim(mo,clmrelhum(:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_relhum.nc'  ,'relhum',clmflag)        !rh fraction 0-1

      CALL read_clim(mo,clmprecip(:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_TP.nc'      ,'TP',clmflag)            !mm/month precepitation

      CALL read_clim(mo,clmclouds(:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_TCC_sfc.nc' ,'TCC_sfc',clmflag)       !cc fraction 0-1

      ! Climatologies were calculated for these fields as well 
      CALL read_clim(mo,clmtaux (:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_taux.nc'    ,'taux',clmflag)          !Nm-2

      CALL read_clim(mo,clmtauy (:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_tauy.nc'    ,'tauy',clmflag)          !Nm-2

      CALL read_clim(mo,clmslp (:,:,iplace), trim(path_era40_clim)//&
             'era40_climatology_MSL_sfc.nc' ,'MSL_sfc',clmflag)       !N m-2

      !call sss_sst_nodc(clmsst(:,:,iplace),clmsss(:,:,iplace),mo)
      call sss_sst_WOA2005(clmsst(:,:,iplace),clmsss(:,:,iplace),mo)

      if (timer)tmr0=wtime()
      call era40_clm_calc(iplace,mo)
      if (timer)tmr1=wtime()
      if (timer)print *,'One era40 clm calc took ',tmr1-tmr0

      !write (lp,'(a,i3)') 'Read ERA40 clim for month ',mo


   else if (clmflag=='ncepr') then ! NEW - Climatology from NCEP fields

      ! Check NCEP path on first pass
      if (first) then
         call getenv('NCEP_CLIM_PATH',cenv)
         if (trim(cenv)/='') then ! prefer this path if present
            path_ncep_clim=trim(cenv)//'/'
         end if
         ret=system('[ -d '//trim(path_ncep_clim)//' ]')
         if (ret /=0 ) then
            print *
            print *,'The directory  '//trim(path_ncep_clim)//' does not exist.'
            print *,'Make sure that '
            print *,' - You have linked the NCEP-Clim data into this catalogue, and '
            print *,'   that the variable NCEP_CLIM_PATH is empty'
            print *,' - OR, the variable NCEP_CLIM_PATH is set to the location of '
            print *,'   the NCEP-Clim data'
            call exit(1)
         end if
      end if

      CALL read_clim(mo,clmuwind (:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_uwnd.nc','uwnd',clmflag)       !m/s

      CALL read_clim(mo,clmvwind (:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_vwnd.nc','vwnd',clmflag)       !m/s

      CALL read_clim(mo,clmwndspd(:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_wspd.nc'  ,'wspd',clmflag)         !m/s

      CALL read_clim(mo,clmairtmp(:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_tair.nc' ,'tair',clmflag)       !kelvin

      CALL read_clim(mo,clmrelhum(:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_rhum.nc'  ,'rhum',clmflag)        !rh fraction 0-1

      CALL read_clim(mo,clmprecip(:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_prcp.nc'  ,'prcp',clmflag)            !mm/month precepitation

      CALL read_clim(mo,clmclouds(:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_ccov.nc' ,'ccov',clmflag)       !cc fraction 0-1

      ! Climatologies were calculated for these fields as well 
      CALL read_clim(mo,clmtaux (:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_taux.nc'    ,'taux',clmflag)          !Nm-2

      CALL read_clim(mo,clmtauy (:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_tauy.nc'    ,'tauy',clmflag)          !Nm-2

      CALL read_clim(mo,clmslp  (:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_pres.nc' ,'pres',clmflag)       !mBar

      CALL read_clim(mo,clmvapmix(:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_vpmx.nc' ,'vpmx',clmflag)       ![]

      CALL read_clim(mo,clmshwflx(:,:,iplace), trim(path_ncep_clim)//&
             'ncep_clim_dswr.nc' ,'dswr',clmflag)       !W m-2

      call sss_sst_WOA2005(clmsst(:,:,iplace),clmsss(:,:,iplace),mo)
      call ncepr_clm_calc(iplace,mo)

   else 
      write(lp,*) 'Unknown Clim flag '//clmflag
      print *, '(m_forfun_nersc)'
      call exit(1)
   end if
   first=.false.
   end subroutine read_clim_all



subroutine read_clim(nmo,field,filen,vname,clmflag)
! Load data from file 'filen' into array 'field'
! and dumps it in integer format if 'diag'=1.
!KAL -- changed to read one month at a time
!KAL -- Based on old NERSC redclm routine, which only supported "old"
!KAL -- climatology type. Extended the routine to read era 40 nc climatology
!KAL -- TODO: Some cleanup...
   use mod_xc
   use mod_grid
   use m_bilin_ecmwf2
   use m_bilin_ncep_gauss
   use m_ncvar_read
   use m_era40_fix
   use netcdf
   implicit none

   integer, intent(in) :: nmo             ! Month number
   real, intent(out) :: field (idm,jdm)   ! Climatology field
   character(len=*), intent(in) :: filen  ! File name
   character(len=*), intent(in) :: vname  ! Variable name
   character(len=*), intent(in) :: clmflag! Climate flag

   integer nxl,nyl,nzl,fact
   integer :: iutil (idm,jdm)
   real grdn,ypn,xpn,ypo,rfact,vmax,vmin,dd,dw,offs,fl
   integer :: i,j,k,iind
   logical :: ex,exa, exb
   real    :: hmin,hmax,hmin2,hmax2
   character*40 :: filena,filenb
   real :: dummy(1,1,1)
   real, dimension(:,:,:), allocatable :: tmpdata
   real :: mlon (idm,jdm)
   integer :: nlon,nlat
   real    :: dlon,dlat,flon,flat,llon,llat
   integer :: index1, index2
   real*8 :: tmr0,tmr1
   real*8, external :: wtime
   integer, parameter :: diag=0
   real, allocatable :: tmplon(:), tmplat(:)
   integer :: ncid, varid, varidlon, varidlat, dimids(NF90_MAX_DIMS), ndims



if (clmflag=='era40') then

   ! Get lon/lat info:
   call ncvar_read(filen,'Ni' ,dummy(1,1,1)  , 1,1,1,1,1) ; nlon=dummy(1,1,1)
   call ncvar_read(filen,'Di' ,dummy(1,1,1)  , 1,1,1,1,1) ; dlon=dummy(1,1,1)
   call ncvar_read(filen,'Lo1',dummy(1,1,1)  , 1,1,1,1,1) ; flon=dummy(1,1,1)
   call ncvar_read(filen,'Lo2',dummy(1,1,1)  , 1,1,1,1,1) ; llon=dummy(1,1,1)
   call ncvar_read(filen,'Nj' ,dummy(1,1,1)  , 1,1,1,1,1) ; nlat=dummy(1,1,1)
   call ncvar_read(filen,'Dj' ,dummy(1,1,1)  , 1,1,1,1,1) ; dlat=dummy(1,1,1)
   call ncvar_read(filen,'La1',dummy(1,1,1)  , 1,1,1,1,1) ; flat=dummy(1,1,1)
   call ncvar_read(filen,'La2',dummy(1,1,1)  , 1,1,1,1,1) ; llat=dummy(1,1,1)

   allocate(tmpdata(nlon,nlat,1))
   call ncvar_read(filen,trim(vname),tmpdata,nlon,nlat,1,nmo,nmo)

   ! fix for era40 climatology precipitation
   if (trim(vname)=='TP') then
      call era40_fix(vname,tmpdata,nlon,nlat,flon,flat,dlon,dlat)
   end if

   ! Do bilinear interpolation of data
   mlon = plon ; where (mlon < 0.) mlon=mlon+360.
   call bilin_ecmwf2(tmpdata,nlon,nlat,flon,flat,dlon,dlat, &
                          field,mlon,plat,depths)   

   deallocate(tmpdata)

elseif (clmflag=='ncepr') then
   call nf90_handle_err(NF90_OPEN(filen,NF90_NOCLOBBER,ncid),'Can not open '//trim(filen))
   call nf90_handle_err(NF90_INQ_VARID(ncid,vname,varid),'Can not find '//trim(vname))
   call nf90_handle_err(NF90_INQ_VARID(ncid,'longitude',varidlon),'Can not find longitude')
   call nf90_handle_err(NF90_INQ_VARID(ncid,'latitude',varidlat),'Can not find latitude')

   call nf90_handle_err(NF90_INQUIRE_VARIABLE(ncid,varidlon,ndims=ndims,dimids=dimids),'')
   if (ndims>1 ) then
      print *,'mod_clim_atm:read_clim:Dimension error - longitude'
      call exit(1)
   else
      call nf90_handle_err(NF90_INQUIRE_DIMENSION(ncid,dimids(1),len=nlon),'')
   end if

   call nf90_handle_err(NF90_INQUIRE_VARIABLE(ncid,varidlat,ndims=ndims,dimids=dimids),'')
   if (ndims>1 ) then
      print *,'mod_clim_atm:read_clim:Dimension error - latitude'
      call exit(1)
   else
      call nf90_handle_err(NF90_INQUIRE_DIMENSION(ncid,dimids(1),len=nlat),'')
   end if

   allocate(tmplon(nlon))
   allocate(tmplat(nlat))
   allocate(tmpdata(nlon,nlat,1))
   call nf90_handle_err(NF90_GET_VAR(ncid,varidlon,tmplon, start=(/1/),count=(/nlon/)),'')
   call nf90_handle_err(NF90_GET_VAR(ncid,varidlat,tmplat, start=(/1/),count=(/nlat/)),'')
   call nf90_handle_err(NF90_GET_VAR(ncid,varid,tmpdata(:,:,1), start=(/1,1,nmo/),count=(/nlon, nlat,1/)),'')
   call nf90_handle_err(NF90_CLOSE(ncid),'')

   mlon = plon ; where (mlon < 0.) mlon=mlon+360.
   call bilin_ncep_gauss(tmpdata,nlon,nlat,tmplon,tmplat,field,mlon,plat)
   deallocate(tmplon,tmplat,tmpdata)

elseif (clmflag=='old') then

   vmax=-9999.
   vmin= 9999.

   ! Check for presence of original ".forc" files
   inquire(file=filen,exist=ex) 

   if (.not.ex) then
      write(lp,'(a)') 'Can not find file '//filen
      call flush(lp)
      call xcstop('(read_clim)')
      stop '(read_clim)'
   end if

   if (mnproc==1) then
      write(*,'(''I read: '',a,i5)')trim(filen)//' for month ',nmo
      call flush(lp)
   end if
   !close(10)
   open(10,FILE=filen,FORM='formatted',STATUS= 'UNKNOWN')
   read(10,200)nxl,nyl,nzl,fact
   if(idm.NE.nxl.OR.jdm.NE.nyl.OR.nzl.NE.12) then
      CLOSE(10)
      if (mnproc==1) then
         write(lp,*)'Wrong dimension in file ',filen    
         write(lp,*)'File dimension: (nx,ny,nmo)',nxl,nyl,nzl
         write(lp,*)'Model dimension:(nx,ny,nmo)',idm,jdm,12
         call flush(lp)
      end if
      call xcSTOP('read_clim')
      STOP 'read_clim'
   end if
   rfact=1./real(fact)


   do k=1,nmo ! Cycle all months
      read(10,201) iutil
      ! Calculate field
      if (k==nmo) then
      do i=1,idm
         do j=1,jdm
            field(i,j)=real(iutil(i,j))*rfact
            fl=.5+SIGN(.5,field(i,j)+998.)        !=0 if no dat
            vmin=MIN(vmin,field(i,j)*fl+(1.-fl)*vmin)
            vmax=MAX(vmax,field(i,j)*fl+(1.-fl)*vmax)
         enddo
      enddo
      end if
   end do
   !print *,vmin,vmax
   close(10)

! --- ---------------------------- Diagnostics -------------
   dd=vmax-vmin
   offs=0.
   if(vmin.LT.0.)offs=-vmin
   if(vmax.LT.1.AND.dd.LT.1.)dw=100.
   if(vmax.GE.1.AND.dd.LT.100.)dw=1.  
   if(dd.GT.100.AND.dd.LT.1000.)dw=.1
   if(dd.GT.1000.AND.dd.LT.10000.)dw=.01
   if(dd.GT.10000.AND.dd.LT.100000.)dw=.001
   if(vmin.GT.200.AND.vmax.LT.400) THEN
     offs=-223.15
     dw=1.
   endif
   !write(*,*)'File: ',filen(:)   ! Already stated
   !write(*,202)vmin,vmax,offs,dw ! Who cares ...

   if(diag.EQ.1)THEN
    do k=1,nzl
     write(*,203)k,filen
     do i=idm,1,-1
      fl=-999.*(.5-SIGN(.5,offs+field(i,j)-.1))
      write(*,'(67I2)')i,(INT(-999.*(.5-SIGN(.5,offs+field(i,j)-.1))&
                +(offs+field(i,j))*dw) ,j=jdm-3,3,-1)
     enddo
    enddo
   endif
else
   print *,'unknown forcing '//clmflag
   call exit(1)
end if

 200  format(4I8)
 201  format(15I8)
 202  format('Min: ',f9.2,' Max: ',f9.2,' Offs: ',f9.2,' Fact ',f9.2)
 203  format('Month: ',i2,'File: ',a10)
      continue

! KAL -- Based on old "reddat"  -- Changes below
! KAL -- Returns one month. Produces binary file versions of .forc files -- 22062005
! KAL -- Allows for use of  Produces binary file versions of .forc files -- 22062005
! KAL -- Removed .ab file creation - not needed in offline version       -- 02092008
end subroutine read_clim




      subroutine nersc_clm_calc(islot,month)
      use mod_xc
      use mod_forcing_nersc
      use mod_atm_func
      use mod_parameters
      use mod_grid
      use mod_year_info22
      use m_qsw0
      implicit none
      integer, intent(in) :: islot,month

      integer :: imargin,jm1,jp1,im1,ip1,i,j
      real :: wfact,vapp,shum,w4, fac
      real*8 :: dtime
      type(year_info) :: rttmp



! ---    ----------------------------------------------
! ---    Calculate properties dependent on what we read
! ---    ----------------------------------------------
         ! Calc sw rad forcing
         call year_day(real(0.,kind=8),0,rttmp,yrflag,.true.)
         dtime=real(rttmp%daysinyear)*(month/12.-1./24.)
         call year_day(dtime,refyear,rttmp,yrflag,.true.)
         !dtime=rttmp%idd+rttmp%ihh/24.+delt1/86400.
         dtime=rttmp%idd+rttmp%ihh/24.
         call qsw0(clmshwflx(:,:,islot),cawdir,clmclouds(:,:,islot), &
            plat/radian,dtime,rttmp%daysinyear)

! --- ------------------------------------------------------------------------
! --- compute air-water wind stresses  tau_x  and  tau_y  from the surface wind 
! --- components. surface wind components have dim [m/s], and are defined on a 
! --- c-grid;  airdns  has dim [kg/m^3], and  cd  is dimensionless. 
! ---  tau_x  and  tau_y  have dim [kg/(m s^2) = N/m^2]
! --- ------------------------------------------------------------------------

         fac=airdns*cd
         imargin=0
         ! KAL: This calculates monthly average stress from monthly average
         !      winds, which is .. kinda wrong..
!$OMP PARALLEL DO PRIVATE(j,l,w4,im1,ip1,jm1,jp1,wfact,vapp,shum)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm


            im1=max(1,i-1)
            jm1=max(1,j-1)
            ip1=min(i+1,idm)
            jp1=min(j+1,jdm)

            !print *,i,j,idm,jdm


            !windstress in u-dir
            w4=.25*(clmvwind(im1,jp1,islot)+clmvwind(i,jp1,islot)+  &
                    clmvwind(im1,j  ,islot)+clmvwind(i,j  ,islot))     
            wfact=sqrt(clmuwind(i,j,islot)*clmuwind(i,j,islot)+w4*w4) &
                   *fac 
            clmtaux(i,j,islot)=clmuwind(i,j,islot)*wfact               

            !windstress in v-dir  
            w4=.25*(clmuwind(i,jm1,islot)+clmuwind(ip1,jm1,islot)  &
                   +clmuwind(ip1,j,islot)+clmuwind(i,j,islot    ))
            !call xcsync(flush_lp)
            wfact=sqrt(clmvwind(i,j,islot)*clmvwind(i,j,islot)+w4*w4) &
                  *fac
            clmtauy(i,j,islot)=clmvwind(i,j,islot)*wfact               

            !from mm/month to m/month 
            clmprecip(i,j,islot)=clmprecip(i,j,islot)*1.e-3            

            !from m/month to m/s
            clmprecip(i,j,islot)=clmprecip(i,j,islot)/(30*86400)

            !! From Relative humidity to mixing ratio
            !vapp = satvappw(clmairtmp(i,j,islot))
            !!print *,mnproc,i,j,islot,vapp,clmairtmp(i,j,islot)
            !shum = clmrelhum(i,j,islot)*humid(1e5,vapp)
            !clmvapmix(i,j,islot) = shum/(1.-shum)
            clmvapmix(i,j,islot) = rhtovpmix(clmrelhum(i,j,islot), &
                                   clmairtmp(i,j,islot),slp0*100)

            ! Calculate Radiative flux - simple but not used
            clmradflx(i,j,islot)= (clmairtmp(i,j,islot))**4    - &
                                  (clmsst   (i,j,islot)+t0deg)**4
            clmradflx(i,j,islot)=clmradflx(i,j,islot)*stefanb*0.95
            clmradflx(i,j,islot)=clmradflx(i,j,islot) &
                                +clmshwflx(i,j,islot)
         enddo
         enddo
!$OMP END PARALLEL DO


      end subroutine nersc_clm_calc



   ! Final postprocessing of era40 fields. Called after all era40 fields
   ! are read.
   ! ------------------------------------------------------------------
   subroutine era40_clm_calc(islot,month)
   use mod_xc
   use mod_parameters
   use mod_forcing_nersc
   use mod_atm_func
   use mod_year_info22
   use mod_grid
   use m_qsw0
   use m_rotate
   implicit none
   integer, intent(in) :: islot,month
   !
   integer :: imargin,jm1,jp1,im1,ip1,i,j
   real :: wfact,vapp,shum,w4, fac
   real*8 :: dtime
   type(year_info) :: rttmp
   real*8 :: tmr0,tmr1
   real*8, external :: wtime
   logical, parameter :: timer=.false.



   ! --- Calc sw rad forcing
   ! --- NB: This operation is quite costly
   if (timer) tmr0=wtime()
   call year_day(real(0.,kind=8),0,rttmp,yrflag,.true.)
   dtime=real(rttmp%daysinyear)*(month/12.-1./24.)
   call year_day(dtime,refyear,rttmp,yrflag,.true.)
   !dtime=rttmp%idd+rttmp%ihh/24.+delt1/86400.
   dtime=rttmp%idd+rttmp%ihh/24.
   call qsw0(clmshwflx(:,:,islot),cawdir,clmclouds(:,:,islot), &
             plat/radian,dtime,rttmp%daysinyear)
    if (timer) tmr1=wtime()
    if (timer) print *,'One shwfl field calc took ',tmr1-tmr0


   fac=airdns*cd
   ! KAL: This calculates monthly average stress from monthly average
   !      winds, which is .. kinda wrong..
   !$OMP PARALLEL DO PRIVATE(j,l,w4,im1,ip1,jm1,jp1,wfact,vapp,shum)
   !$OMP&         SCHEDULE(STATIC,jblk)
   do j=1,jdm
   do i=1,idm


      ! From Relative humidity to mixing ratio
      clmvapmix(i,j,islot) = rhtovpmix(clmrelhum(i,j,islot), &
                             clmairtmp(i,j,islot),slp0*100)

      ! ERA40 Precipitation is accumulated over 6hrs
      clmprecip(i,j,islot)= clmprecip(i,j,islot)/(6*3600.)

      ! Avoid neg precip (very small, scaling errors?)
      clmprecip(i,j,islot)=max(0.,clmprecip(i,j,islot))

      ! ERA40 SLP is in Pa
      clmslp(i,j,islot)= clmslp(i,j,islot)*0.01

      ! Calculate Radiative flux
      clmradflx(i,j,islot)= (clmairtmp(i,j,islot))**4    - &
                   (clmsst   (i,j,islot)+t0deg)**4
      clmradflx(i,j,islot)=clmradflx(i,j,islot)*stefanb*0.95
      clmradflx(i,j,islot)=clmradflx(i,j,islot) &
                          +clmshwflx(i,j,islot)
   enddo
   enddo
   !$OMP END PARALLEL DO

   ! Rotate velocity model grid
   call rotate(clmuwind(:,:,islot), clmvwind(:,:,islot), plat, plon, idm,jdm,'l2m')

   ! Rotate stresses to model grid
   call rotate(clmtaux(:,:,islot),  clmtauy(:,:,islot) , plat, plon, idm,jdm,'l2m')
   end subroutine era40_clm_calc



   ! Final postprocessing of era40 fields. Called after all era40 fields
   ! are read.
   ! ------------------------------------------------------------------
   subroutine ncepr_clm_calc(islot,month)
   use mod_xc
   use mod_parameters
   use mod_forcing_nersc
   use mod_atm_func
   use mod_year_info22
   use mod_grid
   use m_qsw0
   use m_rotate
   implicit none
   integer, intent(in) :: islot,month
   !
   integer :: imargin,jm1,jp1,im1,ip1,i,j
   real :: wfact,vapp,shum,w4, fac

   !$OMP PARALLEL DO PRIVATE(j,l,w4,im1,ip1,jm1,jp1,wfact,vapp,shum)
   !$OMP&         SCHEDULE(STATIC,jblk)
   do j=1,jdm
   do i=1,idm

      ! Calculate Radiative flux
      clmradflx(i,j,islot)= (clmairtmp(i,j,islot))**4    - &
                            (clmsst   (i,j,islot)+t0deg)**4
      clmradflx(i,j,islot)=clmradflx(i,j,islot)*stefanb*0.95
      clmradflx(i,j,islot)=clmradflx(i,j,islot) &
                          +clmshwflx(i,j,islot)
   enddo
   enddo
   !$OMP END PARALLEL DO

   ! Rotate velocity model grid
   call rotate(clmuwind(:,:,islot), clmvwind(:,:,islot), plat, plon, idm,jdm,'l2m')

   ! Rotate stresses to model grid
   call rotate(clmtaux(:,:,islot),  clmtauy(:,:,islot) , plat, plon, idm,jdm,'l2m')
   end subroutine ncepr_clm_calc







! #######################################################################
! Diagnostic for climatology forcing (tecplot file)
! #######################################################################
   subroutine clim_diag(mo,fdiag,ifdiag)
   use mod_xc
   use mod_forcing_nersc
   use mod_grid
   use netcdf
   implicit none
   integer, intent(in) :: mo,ifdiag
   logical, intent(in) ::  fdiag

   integer :: i,j
   integer :: idmid,jdmid, ncid, tdmid, ierr, varid

   ! ---    --------------------------------------
   ! ---    used to be diagnostic ouput to tecplot ascii files
   ! ---    Diagnostic to netcdf files in stead - ascii 
   ! ---    output is very slow for large grids
   ! ---    --------------------------------------
   if (fdiag)THEN
      if (mo==1) then
         ierr=nf90_create('clmforce.nc',NF90_CLOBBER,ncid)
         ierr=nf90_def_dim(ncid,'idm',   idm,idmid)
         ierr=nf90_def_dim(ncid,'jdm',   jdm,jdmid)
         ierr=nf90_def_dim(ncid,'month', NF90_UNLIMITED,tdmid)
         ierr=nf90_def_var(ncid,'longitude', NF90_FLOAT, (/idmid,jdmid/),varid)
         ierr=nf90_def_var(ncid,'latitude ', NF90_FLOAT, (/idmid,jdmid/),varid)
         ierr=nf90_def_var(ncid,'depths'   , NF90_FLOAT, (/idmid,jdmid/),varid)
         ierr=nf90_def_var(ncid,'wndspd'   , NF90_FLOAT, (/idmid,jdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , 'm s-1')
         ierr=nf90_def_var(ncid,'airtmp'   , NF90_FLOAT, (/idmid,jdmid,tdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , 'C')
         ierr=nf90_def_var(ncid,'clouds'   , NF90_FLOAT, (/idmid,jdmid,tdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , '')
         ierr=nf90_def_var(ncid,'relhum'   , NF90_FLOAT, (/idmid,jdmid,tdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , '')
         ierr=nf90_def_var(ncid,'precip'   , NF90_FLOAT, (/idmid,jdmid,tdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , 'mm d-1')
         ierr=nf90_def_var(ncid,'taux'     , NF90_FLOAT, (/idmid,jdmid,tdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , 'N m-2')
         ierr=nf90_def_var(ncid,'tauy'     , NF90_FLOAT, (/idmid,jdmid,tdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , 'N m-2')
         ierr=nf90_def_var(ncid,'sss'      , NF90_FLOAT, (/idmid,jdmid,tdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , 'psu')
         ierr=nf90_def_var(ncid,'sst'      , NF90_FLOAT, (/idmid,jdmid,tdmid/),varid)
         ierr=nf90_put_att(ncid,varid, 'units'   , 'C')
         ierr=nf90_enddef(ncid)

         ! NB - var output must be same order as defined
         ierr=nf90_put_var(ncid,1,plon)
         ierr=nf90_put_var(ncid,2,plon)
         ierr=nf90_put_var(ncid,3,depths)
         ierr=nf90_close(ncid)
      end if
      ierr=nf90_open('clmforce.nc',NF90_WRITE,ncid)

      ! NB - var output must be same order as defined
      ierr=nf90_put_var(ncid, 4,clmwndspd(:,:,1),start=(/1,1,mo/))
      ierr=nf90_put_var(ncid, 5,clmairtmp(:,:,1),start=(/1,1,mo/))
      ierr=nf90_put_var(ncid, 6,clmclouds(:,:,1),start=(/1,1,mo/))
      ierr=nf90_put_var(ncid, 7,clmrelhum(:,:,1),start=(/1,1,mo/))
      ierr=nf90_put_var(ncid, 8,clmprecip(:,:,1)*86400.,start=(/1,1,mo/))
      ierr=nf90_put_var(ncid, 9,clmtaux  (:,:,1),start=(/1,1,mo/))
      ierr=nf90_put_var(ncid,10,clmtauy  (:,:,1),start=(/1,1,mo/))
      ierr=nf90_put_var(ncid,11,clmsss   (:,:,1),start=(/1,1,mo/))
      ierr=nf90_put_var(ncid,12,clmsst   (:,:,1),start=(/1,1,mo/))
      ierr=nf90_close(ncid)
   end if ! diag flag
   ! ---    --------------------------------------
   ! ---    End Diagnostic ouput to tecforce.dat
   ! ---    --------------------------------------
  99  FORMAT(30I4) 
 100  FORMAT(10(1x,e10.4)) 
   end subroutine clim_diag








! #######################################################################
! Diagnostic for synoptic forcing (tecplot file)
! #######################################################################
      subroutine syn_diag(lsyntst,rttmp)
      use mod_xc
      use mod_forcing_nersc
      use mod_year_info22
      use mod_grid
      implicit none
      logical, intent(in) :: lsyntst
      type(year_info) :: rttmp

      character(len=80) :: syntstf
      integer :: i,j


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnose synoptic forcing
      IF(lsyntst)THEN
         syntstf='synforce_y'//rttmp%cyy//'d'//rttmp%cdd//'h' //rttmp%chh//'.tec'
         print *,'Diagnostics to '//trim(syntstf)
         OPEN(10,FILE=trim(syntstf),STATUS='UNKNOWN')
         WRITE(10,'(''TITLE= "Forcing fields"'')')
         write(10,'(a)') &
         'VARIABLES="i" "j" "lon" "lat" "depths[m]" "wndspd[m/s]"  &
         "airtmp[C]" "clouds[0-1]" "relhum[0-1]" "vapmix[0-1]"  &
         "precip[mm/d]" "taux[N/m^2]" "tauy[N/m^2]" &
         "shwflx[W/m^2]" "radfl[W/m^2]x" "slp[mBar]"'
         WRITE(10,'(''ZONE I='',I3,'', J='',I3,'', F=BLOCK'')')  idm,jdm
         WRITE(10,99)((i,i=1,idm),j=1,jdm)
         WRITE(10,99)((j,i=1,idm),j=1,jdm)
         WRITE(10,100)((plon     (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((plat     (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((depths   (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synwndspd(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synairtmp(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synclouds(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synrelhum(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synvapmix(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synprecip(i,j)*86400.*1000.,i=1,idm),j=1,jdm)
         WRITE(10,100)((syntaux  (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((syntauy  (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synshwflx(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synradflx(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synslp   (i,j),i=1,idm),j=1,jdm)
         CLOSE(10)
      end if
  99  FORMAT(30I4) 
 100  FORMAT(10(1x,e10.4)) 
      end subroutine syn_diag




! ------- ------------ ----------- -------------- ------------ ---
! ------  Change notes:
! ------- ------------ ----------- -------------- ------------ ---
! --- Cleaned up main routine - created many subroutines            KAL -- 06022005
! --- Got rid of global arrays for climatology (reduces mem usage)  KAL -- 06022005
! --- Introduced preprocessing of clim. - placed in ".a" files      KAL -- 06022005
! --- Clim files now placed in .ab-files -- more robust             KAL -- 22022005
! ------- ------------ ----------- -------------- ------------ ---


end module mod_clim_atm
