module mod_read_erai

contains

   subroutine read_erai(rt,plon,plat,depths)
      use mod_xc
      use mod_year_info22
      use mod_forcing_nersc
      use mod_atm_func
      use m_bilin_ecmwf2
      use m_ncvar_dims
      use m_ncvar_read
      use m_rotate
      use netcdf
      use mod_za
      implicit none

      ! Data catalogue name
      character(len=200), save :: erai_path='./ERA-I/'

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


      type (year_info), intent(in) :: rt
      real, dimension(idm,jdm), intent(in) :: plon,plat,depths
      character(len=5) :: cirec


      ! File names
      character(len=200) :: fprec, ftair, ftcc, fuwnd, fvwnd, fblh, &
         flmsk, fpres, fdewp, fssrd, fssrd2, ficec , fZ
      character(len=100), save  :: cline
      character(len=200)  :: cenv
      logical :: lfirst=.true.,ex

      integer  :: i,j,k,ix,jy,irec,maxrec,margin_loc

      real, dimension(:,:,:), allocatable :: &
         erai_prec, erai_pres, erai_dewp, erai_tair, &
         erai_tcc,  erai_uwnd, erai_vwnd, erai_ssrd
      real, allocatable ::  erai_lmsk(:,:,:), erai_rhum(:,:), tmpfld(:,:), &
         xwnd(:,:,:),ywnd(:,:,:), elon(:,:,:),elat(:,:,:)

      real, dimension(idm,jdm) :: mlon,tstspd
      integer, dimension(idm,jdm) :: tstip
      real :: rfac,wndfac,cd_new,w4,wfact,svpair, svpdew
      real :: dlon,dlat,flon,flat,llon,llat,dummy(1,1,1)
      real :: vpair_x,vpair_w,vpair_i,humid_i,humid_w,sphum,vpair_1
      integer :: nlon,nlat,ndims,recdim
      integer , dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      real :: satvapp, ax, bx,eraicd,sathumid
      integer :: im1,ip1,jm1,jp1
      character(len=4)  :: cyy
      real :: radian
      real :: xmin,xmax
      integer*4, external :: system
      integer*4 :: ret


      ! Check ERAi path on first pass
      if (lfirst) then
         call getenv('ERAI_PATH',cenv)
         if (trim(cenv)/='') then ! prefer this path if present
            erai_path=trim(cenv)//'/'
         end if
         ret=system('[ -d '//trim(erai_path)//' ]')
         if (ret /=0 ) then
            print *
            print *,'The directory  '//trim(erai_path)//' does not exist.'
            print *,'Make sure that '
            print *,' - You have linked the ERAI data into this catalogue, and '
            print *,'   that the variable ERAI_PATH is empty'
            print *,' - OR, the variable ERAI_PATH is set to the location of '
            print *,'   the ERAI data'
            call exit(1)
         end if
      end if

      ! We are reading record number [erai starts on 1.Jan 00:00]
      ! Should have a more robust checking here ...
      irec  = rt%ihh/6 + rt%idd*4 +1

      ! Diagnostic output
      if (mod(irec,1)==0.or. (lfirst)) then
            write(lp,'(a,i6,a,i4,x,i3,x,i2,a,a)',advance='yes')'ERAi data from record: ',irec,&
            '  rt= ',rt%iyy,rt%idd,rt%ihh,' -- Using ',trim(cline)
      end if

      radian=asin(1.)/90.

      mlon = plon ; where (mlon < 0.) mlon=mlon+360.

      ! Set file names
      write(cyy,'(i4.4)') rt%iyy+1
      fuwnd=trim(erai_path)//'ans.6h.'//rt%cyy//'.10U.nc'
      fvwnd=trim(erai_path)//'ans.6h.'//rt%cyy//'.10V.nc'
      fdewp=trim(erai_path)//'ans.6h.'//rt%cyy//'.2D.nc'
      ftair=trim(erai_path)//'ans.6h.'//rt%cyy//'.2T.nc'
      fblh =trim(erai_path)//'fcs.6h.'//rt%cyy//'.BLH.nc'  ! Forecast field !
      ficec=trim(erai_path)//'ans.6h.'//rt%cyy//'.CI.nc'
      fpres=trim(erai_path)//'ans.6h.'//rt%cyy//'.MSL.nc'
      fssrd =trim(erai_path)//'fcs.6h.'//rt%cyy//'.SSRD.nc' ! Forecast field !
      fssrd2=trim(erai_path)//'fcs.6h.'//cyy//'.SSRD.nc' ! Forecast field !
      fprec=trim(erai_path)//'fcs.6h.'//rt%cyy//'.TP.nc'   ! Forecast field !
      ftcc =trim(erai_path)//'ans.6h.'//rt%cyy//'.TCC.nc'
      flmsk=trim(erai_path)//'ans.LSM.nc'
      fZ   =trim(erai_path)//'ans.Z.nc'


      ! Check the allocation status of arrays
      if ((.not.allocated(synslp   ).and. lsynslp   ) .or. &
          (.not.allocated(synrelhum).and. lsynrelhum) .or. &
          (.not.allocated(synuwind ).and. lsynwnd   ) .or. &
          (.not.allocated(synvwind ).and. lsynwnd   ) .or. &
          (.not.allocated(syntaux  ).and. lsynwnd   ) .or. &
          (.not.allocated(syntauy  ).and. lsynwnd   ) .or. &
          (.not.allocated(synwndspd).and. lsynwnd   ) .or. &
          (.not.allocated(synairtmp).and. lsynairtmp) .or. &
          (.not.allocated(synprecip).and. lsynprecip) .or. &
          (.not.allocated(synshwflx).and. lsynshwflx)) then
         if (mnproc==1) write(lp,'(a)') 'ERAi Allocation error'
         call xcstop('(read_erai)')
         stop '(read_erai)'
      end if

      if (lfirst) cline=''




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Get longitude  & latitude for data  - use land mask file
      call ncvar_dims(fssrd,vssrd,dimsizes,ndims,recdim) ; maxrec=dimsizes(recdim)
      !print *,irec,maxrec!,dimsizes,recdim
     ! call ncvar_read(flmsk,'Ni' ,dummy(1,1,1)  , 1,1,1,1,1) ; nlon=dummy(1,1,1)
     ! call ncvar_read(flmsk,'Di' ,dummy(1,1,1)  , 1,1,1,1,1) ; dlon=dummy(1,1,1)
     ! call ncvar_read(flmsk,'Lo1',dummy(1,1,1)  , 1,1,1,1,1) ; flon=dummy(1,1,1)
     ! call ncvar_read(flmsk,'Lo2',dummy(1,1,1)  , 1,1,1,1,1) ; llon=dummy(1,1,1)
     ! call ncvar_read(flmsk,'Nj' ,dummy(1,1,1)  , 1,1,1,1,1) ; nlat=dummy(1,1,1)
     ! call ncvar_read(flmsk,'Dj' ,dummy(1,1,1)  , 1,1,1,1,1) ; dlat=dummy(1,1,1)
     ! call ncvar_read(flmsk,'La1',dummy(1,1,1)  , 1,1,1,1,1) ; flat=dummy(1,1,1)
     ! call ncvar_read(flmsk,'La2',dummy(1,1,1)  , 1,1,1,1,1) ; llat=dummy(1,1,1)
     nlon=720
     nlat=361
     flon=-180
     flat=90
     dlon=0.5
     dlat=-0.5


      !print *,nlon,dlon,flon,llon,flon+(nlon-1)*dlon
      !print *,nlat,dlat,flat,llat,flat+(nlat-1)*dlat
      !stop


      ! Allocate temp fields
      allocate(erai_prec(nlon,nlat,1))
      allocate(erai_tair(nlon,nlat,1))
      allocate(erai_tcc (nlon,nlat,1))
      allocate(erai_uwnd(nlon,nlat,1))
      allocate(erai_vwnd(nlon,nlat,1))
      allocate(erai_pres(nlon,nlat,1))
      allocate(erai_ssrd(nlon,nlat,1))
      allocate(erai_dewp(nlon,nlat,1))
      allocate(erai_lmsk(nlon,nlat,1))
      allocate(erai_rhum(nlon,nlat))
      allocate(xwnd(nlon,nlat,1))
      allocate(ywnd(nlon,nlat,1))
      allocate(elon(nlon,nlat,1))
      allocate(elat(nlon,nlat,1))
      allocate(tmpfld(nlon,nlat))

      do j=1,nlat
      do i=1,nlon
         elon(i,j,1)=flon+dlon*(i-1)
         elat(i,j,1)=flat+dlat*(j-1)
      end do
      end do
      !print *,flat,flat+dlat*(nlat-1)
      !print *,flon,flon+dlon*(nlon-1)
      !call exit(1)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ERAi land mask
      call ncvar_read(flmsk,vlmsk,erai_lmsk,nlon,nlat, 1,1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ERAi precipitation   m/ 6hr -> m/month (rho=1000)
      if (lsynprecip) then
         if (lfirst) cline= trim(cline)//' precipitation '
         call ncvar_read(fprec,vprec,erai_prec,nlon,nlat,1,irec,irec)



         call bilin_ecmwf2(erai_prec,nlon,nlat,flon,flat,dlon,dlat, &
                          synprecip,mlon,plat,depths)   

         ! Precipitation is accumulated field...
         synprecip=synprecip/(6.*3600.)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ERAi air temperature   K 
      if (lsynairtmp) then
         if (lfirst) cline= trim(cline)//' temperature '
         call ncvar_read(ftair,vtair,erai_tair,nlon,nlat,1,irec,irec)
         call bilin_ecmwf2(erai_tair,nlon,nlat,flon,flat,dlon,dlat, &
                          synairtmp,mlon,plat,depths)   
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP clouds   []  (0-1)
      if (lsynclouds) then
         if (lfirst) cline= trim(cline)//' clouds '
         call ncvar_read(ftcc,vtcc,erai_tcc,nlon,nlat,1,irec,irec)
         call bilin_ecmwf2(erai_tcc ,nlon,nlat,flon,flat,dlon,dlat, &
                          synclouds,mlon,plat,depths)   
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP winds
      if (lsynwnd) then
         if (lfirst) cline= trim(cline)//' winds '

         call ncvar_read(fuwnd,vuwnd,erai_uwnd,nlon,nlat,1,irec,irec)
         call ncvar_read(fvwnd,vvwnd,erai_vwnd,nlon,nlat,1,irec,irec)

!#ifdef ERAi_WIND_POLE_FIX
!         ! Rotate ERAi winds to x/y components on plane normal to NP
!         xwnd=cos(elon*radian)*(-erai_vwnd) - sin(elon*radian)*erai_uwnd
!         ywnd=sin(elon*radian)*(-erai_vwnd) + cos(elon*radian)*erai_uwnd
!
!         ! Average northernmost wind components - averaging strategy?
!         !xwnd(:,nlat,1)=sum(xwnd(:,nlat-1,1))/nlon
!         !ywnd(:,nlat,1)=sum(ywnd(:,nlat-1,1))/nlon
!         xwnd(:,nlat,1)=xwnd(:,nlat-1,1)
!         ywnd(:,nlat,1)=ywnd(:,nlat-1,1)
!         erai_vwnd(:,nlat,1)=-xwnd(:,nlat,1)*cos(elon(:,nlat,1)*radian) - &
!                               ywnd(:,nlat,1)*sin(elon(:,nlat,1)*radian)
!         erai_uwnd(:,nlat,1)=-xwnd(:,nlat,1)*sin(elon(:,nlat,1)*radian) + &
!                               ywnd(:,nlat,1)*cos(elon(:,nlat,1)*radian)
!         !erai_uwnd(:,nlat,1)=erai_uwnd(:,nlat-1,1)
!         !erai_vwnd(:,nlat,1)=erai_vwnd(:,nlat-1,1)
!#endif

!KAL     ERAi wind data goes all to 90N - this neglects the values there
!        call bilin_ecmwf2(erai_uwnd(:,:,1) ,nlon,nlat,flon,flat,dlon,dlat, &
!                         synuwind,mlon,plat,depths)   
!        call bilin_ecmwf2(erai_vwnd(:,:,1) ,nlon,nlat,flon,flat,dlon,dlat, &
!                         synvwind,mlon,plat,depths)   
         call bilin_ecmwf2(erai_uwnd(:,1:nlat-1,1) ,nlon,nlat-1,flon,flat,dlon,dlat, &
                          synuwind,mlon,plat,depths)   
         call bilin_ecmwf2(erai_vwnd(:,1:nlat-1,1) ,nlon,nlat-1,flon,flat,dlon,dlat, &
                          synvwind,mlon,plat,depths)   

         ! Rotate winds to model grid 
         if (rotateflag==1) then
            call rotate(synuwind, &
                        synvwind, &
                        plat    , &
                        mlon    ,idm,jdm,'l2m')
                        !print *,'plain rotate'
         elseif (rotateflag==2) then
            call rotate_sphere(synuwind, &
                        synvwind, &
                        plat    , &
                        mlon    ,idm,jdm,'l2m')
                        !print *,'sphere rotate'
         end if

         ! Derived  wind stress and absolute windspeed
         eraicd=0.0012

         !$OMP PARALLEL DO PRIVATE(j,i,wndfac,w4,wfact,cd_new,im1,jm1,ip1,jp1) &
         !$OMP SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm

            im1=max(1,i-1)
            jm1=max(1,j-1)
            ip1=min(idm,i+1)
            jp1=min(jdm,j+1)

            synwndspd(i,j)=sqrt(synuwind(i,j)*synuwind(i,j)+ &
                             synvwind(i,j)*synvwind(i,j))

            wndfac=(1.+sign(1.,synwndspd(i,j)-11.))*.5
            cd_new=(0.49+0.065*synwndspd(i,j))*1.0e-3*wndfac+cd*(1.-wndfac)


            w4   =.25*(synvwind(im1,jp1)+synvwind(i,jp1)+ &
                       synvwind(im1,j  )+synvwind(i,j  ))
            wfact=sqrt(synuwind(i,j)*synuwind(i,j)+w4*w4)*airdns*cd_new
            syntaux(i,j)=synuwind(i,j)*wfact

            w4   =.25*(synuwind(i  ,jm1)+synuwind(ip1,jm1)+&
                       synuwind(ip1,j  )+synuwind(i  ,j  ))
            wfact=sqrt(synvwind(i,j)*synvwind(i,j)+w4*w4)*airdns*cd_new
            syntauy(i,j)=synvwind(i,j)*wfact

         end do
         end do
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pressure
      if (lsynslp) then
         if (lfirst) cline= trim(cline)//' slp '
         call ncvar_read(fpres,vpres,erai_pres,nlon,nlat,1,irec,irec)
         call bilin_ecmwf2(erai_pres ,nlon,nlat,flon,flat,dlon,dlat, &
                          synslp,mlon,plat,depths)   
         synslp=synslp*0.01 ! Pa -> mBar (HPa)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Relative humidity
      if (lsynrelhum) then

         if (lfirst) cline= trim(cline)//' relhum '
         call ncvar_read(fdewp,vdewp,erai_dewp,nlon,nlat,1,irec,irec)

         ! SLP is needed to calc relative humidity (alternative - use slp0)
         if (.not.lsynslp) then
            write(lp,'(a)') ' ERROR: ERAi slp is needed to calculate relhum ...  '
            call xcstop ('(read_erai)')
            stop '(read_erai)'
         end if

         ! relative humidity
         !$OMP PARALLEL DO PRIVATE (i,j,svpair,svpdew)
         do j=1,nlat
         do i=1,nlon
            svpair=satvap(erai_tair(i,j,1))
            svpdew=satvap(erai_dewp(i,j,1))
            erai_rhum(i,j)=relhumid(svpair,svpdew,erai_pres(i,j,1))*0.01

            ! Ahem ...
            erai_rhum(i,j)  =max(0.0,min(erai_rhum(i,j),1.0))
         enddo
         enddo
         !$OMP END PARALLEL DO
         call bilin_ecmwf2(erai_rhum ,nlon,nlat,flon,flat,dlon,dlat, &
                          synrelhum,mlon,plat,depths)   
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SW radiation

      if (lsynshwflx) then
         if (lfirst) cline= trim(cline)//' ssrd '

         call ncvar_read(fssrd,vssrd,erai_ssrd,nlon,nlat,1,irec,irec)
         tmpfld=erai_ssrd(:,:,1)

         ! ssrd is accumulated over 6hrs, midpoint is 6hr-3hr
         if (irec==maxrec)  then
            call ncvar_read(fssrd2,vssrd,erai_ssrd,nlon,nlat,1,1,1)
         else if (irec<maxrec) then
            call ncvar_read(fssrd ,vssrd,erai_ssrd,nlon,nlat,1,irec+1,irec+1)
         else
            print *,'You shouldnt be here ...'
            call xcstop('(read_erai)')
            stop '(read_erai)'
         end if
         erai_ssrd(:,:,1)=erai_ssrd(:,:,1)+tmpfld
         erai_ssrd=erai_ssrd/2.


         call bilin_ecmwf2(erai_ssrd,nlon,nlat,flon,flat,dlon,dlat, &
                          synshwflx,mlon,plat,depths)   
         synshwflx=synshwflx/(6*3600) ! Wm^-2 * s - > Wm^-2
      end if


      ! Clean up workspace
      deallocate(erai_prec)
      deallocate(erai_tair)
      deallocate(erai_tcc )
      deallocate(erai_uwnd)
      deallocate(erai_vwnd)
      deallocate(erai_pres)
      deallocate(erai_ssrd)
      deallocate(erai_dewp)
      deallocate(erai_lmsk)
      deallocate(erai_rhum)
      deallocate(tmpfld)
      lfirst=.false.
   END subroutine read_erai
end module mod_read_erai
