module m_read_era40

contains

   subroutine read_era40(rt,plon,plat,depths)
      use mod_xc
      use mod_year_info 
      use mod_forcing_nersc
      use m_bilin_ecmwf2
      use m_ncvar_dims
      use m_ncvar_read
      use m_rotate
      use m_rotate_sphere
      use m_satvap
      use m_relhumid
      use m_era40_fix
      use netcdf
      use mod_za
      implicit none

      ! Data catalogue name
      character(len=*), parameter :: era40_path='./ERA40/'

      ! Variable names - ERA40 1.125 by 1.125 
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

      type (year_info), intent(in) :: rt
      real, dimension(idm,jdm), intent(in) :: plon,plat,depths
      character(len=5) :: cirec


      ! File names
      character(len=80) :: fprec, ftair, ftcc, fuwnd, fvwnd, fblh, &
         flmsk, fpres, fdewp, fssrd, fssrd2, ficec,fz

      integer  :: i,j,k,ix,jy,irec,maxrec,margin_loc

      real, dimension(:,:,:), allocatable :: &
         era40_prec, era40_pres, era40_dewp, era40_tair, &
         era40_tcc,  era40_uwnd, era40_vwnd, era40_ssrd
      real, allocatable ::  era40_lmsk(:,:,:), era40_rhum(:,:), tmpfld(:,:), &
         xwnd(:,:,:),ywnd(:,:,:), elon(:,:,:),elat(:,:,:)

      real, dimension(idm,jdm) :: mlon,tstspd
      integer, dimension(idm,jdm) :: tstip
      real :: rfac,wndfac,cd_new,w4,wfact,svpair, svpdew
      real :: dlon,dlat,flon,flat,llon,llat,dummy(1,1,1)
      real :: vpair_x,vpair_w,vpair_i,humid_i,humid_w,sphum,vpair_1
      integer :: nlon,nlat,ndims,recdim
      integer , dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      real :: satvapp, ax, bx,era40cd,sathumid
      logical :: lfirst=.true.
      integer :: im1,ip1,jm1,jp1


      character(len=4)  :: cyy

      real :: radian
      real :: xmin,xmax

      radian=asin(1.)/90.

      !print *,1/radian


      !if (mnproc==1) print *,'ERA40 -- add time stamp check!'

      !call xcstop('read_rea40 not changed to using tiles')
      !stop 'read_rea40 not changed to using tiles'

      ! We are reading record number [era40 starts on 1.Jan 00:00]
      ! Should have a more robust checking here ...
      irec  = rt%ihh/6 + rt%idd*4 +1

      ! Diagnostic output
      if (mod(irec,40)==0.or.lfirst) then
         if (mnproc==1) then
            WRITE(lp,*)
            if (lfirst) then
               write(lp,'(a,i6,a,i4,x,i3,x,i2)',advance='yes')'ERA40 data from record: ',irec,&
               '  rt= ',rt%iyy,rt%idd,rt%ihh
            else
               write(lp,'(a,i6,a,i4,x,i3,x,i2)',advance='no') 'ERA40 data from record: ',irec,&
               '  rt= ',rt%iyy,rt%idd,rt%ihh
            end if
         end if
      else
         if (mnproc==1) write(lp,'(a1)',advance='no')'.'
         call flush(lp)
      end if

      mlon = plon ; where (mlon < 0.) mlon=mlon+360.

      ! Set file names
      write(cyy,'(i4.4)') rt%iyy+1
      fuwnd=era40_path//'ans.6h.'//rt%cyy//'.10U.nc'
      fvwnd=era40_path//'ans.6h.'//rt%cyy//'.10V.nc'
      fdewp=era40_path//'ans.6h.'//rt%cyy//'.2D.nc'
      ftair=era40_path//'ans.6h.'//rt%cyy//'.2T.nc'
      fblh =era40_path//'fcs.6h.'//rt%cyy//'.BLH.nc'  ! Forecast field !
      ficec=era40_path//'ans.6h.'//rt%cyy//'.CI.nc'
      fpres=era40_path//'ans.6h.'//rt%cyy//'.MSL.nc'
      fssrd =era40_path//'fcs.6h.'//rt%cyy//'.SSRD.nc' ! Forecast field !
      fssrd2=era40_path//'fcs.6h.'//cyy//'.SSRD.nc' ! Forecast field !
      fprec=era40_path//'fcs.6h.'//rt%cyy//'.TP.nc'   ! Forecast field !
      ftcc =era40_path//'ans.6h.'//rt%cyy//'.TCC.nc'
      flmsk=era40_path//'ans.LSM.nc'
      fZ   =era40_path//'ans.Z.nc'


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
         if (mnproc==1) write(lp,'(a)') 'ERA40 Allocation error'
         call xcstop('(read_era40)')
         stop '(read_era40)'
      end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Get longitude  & latitude for data  - use land mask file
      call ncvar_dims(fssrd,vssrd,dimsizes,ndims,recdim) ; maxrec=dimsizes(recdim)
      !print *,irec,maxrec!,dimsizes,recdim
      call ncvar_read(flmsk,'Ni' ,dummy(1,1,1)  , 1,1,1,1,1) ; nlon=dummy(1,1,1)
      call ncvar_read(flmsk,'Di' ,dummy(1,1,1)  , 1,1,1,1,1) ; dlon=dummy(1,1,1)
      call ncvar_read(flmsk,'Lo1',dummy(1,1,1)  , 1,1,1,1,1) ; flon=dummy(1,1,1)
      call ncvar_read(flmsk,'Lo2',dummy(1,1,1)  , 1,1,1,1,1) ; llon=dummy(1,1,1)
      call ncvar_read(flmsk,'Nj' ,dummy(1,1,1)  , 1,1,1,1,1) ; nlat=dummy(1,1,1)
      call ncvar_read(flmsk,'Dj' ,dummy(1,1,1)  , 1,1,1,1,1) ; dlat=dummy(1,1,1)
      call ncvar_read(flmsk,'La1',dummy(1,1,1)  , 1,1,1,1,1) ; flat=dummy(1,1,1)
      call ncvar_read(flmsk,'La2',dummy(1,1,1)  , 1,1,1,1,1) ; llat=dummy(1,1,1)

      !print *,nlon,dlon,flon,llon,flon+(nlon-1)*dlon
      !print *,nlat,dlat,flat,llat,flat+(nlat-1)*dlat
      !stop


      ! Allocate temp fields
      allocate(era40_prec(nlon,nlat,1))
      allocate(era40_tair(nlon,nlat,1))
      allocate(era40_tcc (nlon,nlat,1))
      allocate(era40_uwnd(nlon,nlat,1))
      allocate(era40_vwnd(nlon,nlat,1))
      allocate(era40_pres(nlon,nlat,1))
      allocate(era40_ssrd(nlon,nlat,1))
      allocate(era40_dewp(nlon,nlat,1))
      allocate(era40_lmsk(nlon,nlat,1))
      allocate(era40_rhum(nlon,nlat))
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
! ERA40 land mask
      call ncvar_read(flmsk,vlmsk,era40_lmsk,nlon,nlat, 1,1,1)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ERA40 precipitation   m/ 6hr -> m/month (rho=1000)
      if (lsynprecip) then
         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using ERA40 precipitation '
         call ncvar_read(fprec,vprec,era40_prec,nlon,nlat,1,irec,irec)

         call era40_fix(vprec,era40_prec(:,:,1),nlon,nlat,flon,flat,dlon,dlat)


         call bilin_ecmwf2(era40_prec,nlon,nlat,flon,flat,dlon,dlat, &
                          synprecip,mlon,plat,depths)   

         ! Precipitation is accumulated field...
         synprecip=synprecip/(6.*3600.)

         ! 

         ! convert prcipitation  from  m/6hr to m/(month)
         !rfac=1./(4*30)
         !synprecip=rfac*synprecip
      end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ERA40 air temperature   K 
      if (lsynairtmp) then
         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using ERA40 temperature'
         call ncvar_read(ftair,vtair,era40_tair,nlon,nlat,1,irec,irec)
         call bilin_ecmwf2(era40_tair,nlon,nlat,flon,flat,dlon,dlat, &
                          synairtmp,mlon,plat,depths)   
      end if





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP clouds   []  (0-1)
      if (lsynclouds) then
         if (mnproc==1.and.lfirst) write(lp,'(a)')' Using ERA40 clouds '
         call ncvar_read(ftcc,vtcc,era40_tcc,nlon,nlat,1,irec,irec)
         call bilin_ecmwf2(era40_tcc ,nlon,nlat,flon,flat,dlon,dlat, &
                          synclouds,mlon,plat,depths)   
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP winds
      if (lsynwnd) then
         if (mnproc==1.and.lfirst) write(lp,'(a)')' Using ERA40 winds '

         call ncvar_read(fuwnd,vuwnd,era40_uwnd,nlon,nlat,1,irec,irec)
         call ncvar_read(fvwnd,vvwnd,era40_vwnd,nlon,nlat,1,irec,irec)

!#ifdef ERA40_WIND_POLE_FIX
!         ! Rotate ERA40 winds to x/y components on plane normal to NP
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
!         !era40_uwnd(:,nlat,1)=era40_uwnd(:,nlat-1,1)
!         !era40_vwnd(:,nlat,1)=era40_vwnd(:,nlat-1,1)
!#endif


!        call bilin_ecmwf2(era40_uwnd(:,:,1) ,nlon,nlat,flon,flat,dlon,dlat, &
!                         synuwind,mlon,plat,depths)   
!        call bilin_ecmwf2(era40_vwnd(:,:,1) ,nlon,nlat,flon,flat,dlon,dlat, &
!                         synvwind,mlon,plat,depths)   
!KAL     ERA40 wind data goes all to 90N - this neglects the values there
         call bilin_ecmwf2(era40_uwnd(:,1:nlat-1,1) ,nlon,nlat-1,flon,flat,dlon,dlat, &
                          synuwind,mlon,plat,depths)   
         call bilin_ecmwf2(era40_vwnd(:,1:nlat-1,1) ,nlon,nlat-1,flon,flat,dlon,dlat, &
                          synvwind,mlon,plat,depths)   

         !call zaiopf('uicetst.a','unknown',321)
         !tstspd=synuwind**2+synvwind**2
         !call zaiowr(tstspd,tstip,.false.,xmin,xmax, 321, .false.)

         ! Rotate winds to model grid 
         if (.true.) then
            call rotate(synuwind, &
                        synvwind, &
                        plat    , &
                        mlon    ,idm,jdm,'l2m')
            !print *,'plain rotate'
         else
            ! needs more tests
            call rotate_sphere(synuwind, &
                        synvwind, &
                        plat    , &
                        mlon    ,idm,jdm,'l2m')
            !print *,'sphere rotate'
         end if

         !tstspd=synuwind**2+synvwind**2
         !call zaiowr(tstspd,tstip,.false.,xmin,xmax, 321, .false.)
         !call zaiocl(321)
         !call exit(1)

         ! Derived  wind stress and absolute windspeed
         era40cd=0.0012

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
         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using ERA40 SLP '
         call ncvar_read(fpres,vpres,era40_pres,nlon,nlat,1,irec,irec)
         call bilin_ecmwf2(era40_pres ,nlon,nlat,flon,flat,dlon,dlat, &
                          synslp,mlon,plat,depths)   
         synslp=synslp*0.01 ! Pa -> mBar (HPa)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Relative humidity
      if (lsynrelhum) then

         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using ERA40-derived relative humidity'
         call ncvar_read(fdewp,vdewp,era40_dewp,nlon,nlat,1,irec,irec)

         ! SLP is needed to calc relative humidity (alternative - use slp0)
         if (.not.lsynslp) then
            write(lp,'(a)') ' ERROR: ERA40 slp is needed to calculate relhum ...  '
            call xcstop ('(read_era40)')
            stop '(read_era40)'
         end if

         ! relative humidity
         !$OMP PARALLEL DO PRIVATE (i,j,svpair,svpdew)
         do j=1,nlat
         do i=1,nlon
            svpair=satvap(era40_tair(i,j,1))
            svpdew=satvap(era40_dewp(i,j,1))
            era40_rhum(i,j)=relhumid(svpair,svpdew,era40_pres(i,j,1))*0.01

            ! Ahem ...
            era40_rhum(i,j)  =max(0.0,min(era40_rhum(i,j),1.0))
         enddo
         enddo
         !$OMP END PARALLEL DO
         call bilin_ecmwf2(era40_rhum ,nlon,nlat,flon,flat,dlon,dlat, &
                          synrelhum,mlon,plat,depths)   
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SW radiation

      if (lsynssr) then
         if (mnproc==1.and.lfirst) write(lp,'(a)') ' Using ERA40 shortwave radiation '

         call ncvar_read(fssrd,vssrd,era40_ssrd,nlon,nlat,1,irec,irec)
         tmpfld=era40_ssrd(:,:,1)

         ! ssrd is accumulated over 6hrs, midpoint is 6hr-3hr
         if (irec==maxrec)  then
            call ncvar_read(fssrd2,vssrd,era40_ssrd,nlon,nlat,1,1,1)
         else if (irec<maxrec) then
            call ncvar_read(fssrd ,vssrd,era40_ssrd,nlon,nlat,1,irec+1,irec+1)
         else
            print *,'You shouldnt be here ...'
            call xcstop('(read_era40)')
            stop '(read_era40)'
         end if
         era40_ssrd(:,:,1)=era40_ssrd(:,:,1)+tmpfld
         era40_ssrd=era40_ssrd/2.


         call bilin_ecmwf2(era40_ssrd,nlon,nlat,flon,flat,dlon,dlat, &
                          synshwflx,mlon,plat,depths)   
         synshwflx=synshwflx/(6*3600) ! Wm^-2 * s - > Wm^-2
      end if


      ! Clean up workspace
      deallocate(era40_prec)
      deallocate(era40_tair)
      deallocate(era40_tcc )
      deallocate(era40_uwnd)
      deallocate(era40_vwnd)
      deallocate(era40_pres)
      deallocate(era40_ssrd)
      deallocate(era40_dewp)
      deallocate(era40_lmsk)
      deallocate(era40_rhum)
      deallocate(tmpfld)

      lfirst=.false.


   END subroutine read_era40
end module m_read_era40
