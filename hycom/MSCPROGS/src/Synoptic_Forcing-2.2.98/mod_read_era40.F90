module mod_read_era40

contains

   subroutine read_era40(rt,plon,plat,depths)
      use mod_xc
      use mod_year_info22
      use mod_forcing_nersc
      use mod_atm_func
      use m_interpug
      use m_ncvar_dims
      use m_ncvar_read
      use m_era40_fix
      use netcdf
      use mod_za
      implicit none

      ! Data catalogue name
      character(len=200), save :: era40_path='./ERA40/'

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
      character(len=200) :: fprec, ftair, ftcc, fuwnd, fvwnd, fblh, &
         flmsk, fpres, fdewp, fssrd, fssrd2, ficec , fZ
      character(len=100), save  :: cline
      character(len=200)  :: cenv
      logical :: lfirst=.true.,ex

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
      integer :: im1,ip1,jm1,jp1
      character(len=4)  :: cyy
      real :: radian
      real :: xmin,xmax
      integer*4, external :: system
      integer*4 :: ret
      integer, parameter :: itype=0 ! Sets interpolation type bilinear


      ! Check ERA40 path on first pass
      if (lfirst) then
         call getenv('ERA40_PATH',cenv)
         if (trim(cenv)/='') then ! prefer this path if present
            era40_path=trim(cenv)//'/'
         end if
         ret=system('[ -d '//trim(era40_path)//' ]')
         if (ret /=0 ) then
            print *
            print *,'The directory  '//trim(era40_path)//' does not exist.'
            print *,'Make sure that '
            print *,' - You have linked the ERA40 data into this catalogue, and '
            print *,'   that the variable ERA40_PATH is empty'
            print *,' - OR, the variable ERA40_PATH is set to the location of '
            print *,'   the ERA40 data'
            call exit(1)
         end if
      end if

      ! We are reading record number [era40 starts on 1.Jan 00:00]
      ! Should have a more robust checking here ...
      irec  = rt%ihh/6 + rt%idd*4 +1

      ! Diagnostic output
      if (mod(irec,1)==0.or. (lfirst)) then
            write(lp,'(a,i6,a,i4,x,i3,x,i2,a,a)',advance='yes')'ERA40 data from record: ',irec,&
            '  rt= ',rt%iyy,rt%idd,rt%ihh,' -- Using ',trim(cline)
      end if

      radian=asin(1.)/90.

      mlon = plon ; where (mlon < 0.) mlon=mlon+360.

      ! Set file names
      write(cyy,'(i4.4)') rt%iyy+1
      fuwnd=trim(era40_path)//'ans.6h.'//rt%cyy//'.10U.nc'
      fvwnd=trim(era40_path)//'ans.6h.'//rt%cyy//'.10V.nc'
      fdewp=trim(era40_path)//'ans.6h.'//rt%cyy//'.2D.nc'
      ftair=trim(era40_path)//'ans.6h.'//rt%cyy//'.2T.nc'
      fblh =trim(era40_path)//'fcs.6h.'//rt%cyy//'.BLH.nc'  ! Forecast field !
      ficec=trim(era40_path)//'ans.6h.'//rt%cyy//'.CI.nc'
      fpres=trim(era40_path)//'ans.6h.'//rt%cyy//'.MSL.nc'
      fssrd =trim(era40_path)//'fcs.6h.'//rt%cyy//'.SSRD.nc' ! Forecast field !
      fssrd2=trim(era40_path)//'fcs.6h.'//cyy//'.SSRD.nc' ! Forecast field !
      fprec=trim(era40_path)//'fcs.6h.'//rt%cyy//'.TP.nc'   ! Forecast field !
      ftcc =trim(era40_path)//'ans.6h.'//rt%cyy//'.TCC.nc'
      flmsk=trim(era40_path)//'ans.LSM.nc'
      fZ   =trim(era40_path)//'ans.Z.nc'


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

      if (lfirst) cline=''




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
         if (lfirst) cline= trim(cline)//' precipitation '
         call ncvar_read(fprec,vprec,era40_prec,nlon,nlat,1,irec,irec)

         call era40_fix(vprec,era40_prec(:,:,1),nlon,nlat,flon,flat,dlon,dlat)


         call interpug(era40_prec,nlon,nlat,flon,flat,dlon,dlat, &
                          synprecip,mlon,plat,itype)   

         ! Precipitation is accumulated field...
         synprecip=synprecip/(6.*3600.)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ERA40 air temperature   K 
      if (lsynairtmp) then
         if (lfirst) cline= trim(cline)//' temperature '
         call ncvar_read(ftair,vtair,era40_tair,nlon,nlat,1,irec,irec)
         call interpug(era40_tair,nlon,nlat,flon,flat,dlon,dlat, &
                          synairtmp,mlon,plat,itype)   
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP clouds   []  (0-1)
      if (lsynclouds) then
         if (lfirst) cline= trim(cline)//' clouds '
         call ncvar_read(ftcc,vtcc,era40_tcc,nlon,nlat,1,irec,irec)
         call interpug(era40_tcc ,nlon,nlat,flon,flat,dlon,dlat, &
                          synclouds,mlon,plat,itype)   
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCEP winds
      if (lsynwnd) then
         if (lfirst) cline= trim(cline)//' winds '

         call ncvar_read(fuwnd,vuwnd,era40_uwnd,nlon,nlat,1,irec,irec)
         call ncvar_read(fvwnd,vvwnd,era40_vwnd,nlon,nlat,1,irec,irec)

         call interpug(era40_uwnd(:,1:nlat-1,1) ,nlon,nlat-1,flon,flat,dlon,dlat, &
                          synuwind,mlon,plat,itype)   
         call interpug(era40_vwnd(:,1:nlat-1,1) ,nlon,nlat-1,flon,flat,dlon,dlat, &
                          synvwind,mlon,plat,itype)   

! --- Rotate winds/stress to model u/v points
         if (wndflag==1) then
            call rotate2(synuwind, synvwind, plat, mlon,'l2m',.false.)
         else
! --- Rotate winds/stress to model p-points
            call rotate2(synuwind, synvwind, plat, mlon,'l2m',.true.)
         end if

         ! Derived  wind stress and absolute windspeed
         do j=1,jdm
         do i=1,idm
            synwndspd(i,j)=sqrt(synuwind(i,j)*synuwind(i,j)+ &
                             synvwind(i,j)*synvwind(i,j))
#ifdef KARALIGHT
            wndfac=MAX(2.5,MIN(32.5,synwndspd(i,j)))
            cd_new = 1.0E-3*(.692 + .0710*wndfac - .000700*wndfac**2)
            wfact=synwndspd(i,j)*airdns*cd_new
            syntaux(i,j)=synuwind(i,j)*wfact
            syntauy(i,j)=synvwind(i,j)*wfact
#else
            wndfac=(1.+sign(1.,synwndspd(i,j)-11.))*.5
            cd_new=(0.49+0.065*synwndspd(i,j))*1.0e-3*wndfac+cd*(1.-wndfac)
            wfact=synwndspd(i,j)*airdns*cd_new
            syntaux(i,j)=synuwind(i,j)*wfact
            syntauy(i,j)=synvwind(i,j)*wfact
#endif
! ustar - new. Not quite right if in u/v points
            synustar(i,j) = sqrt((syntaux(i,j)**2+syntauy(i,j)**2)*1e-3)

         end do
         end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pressure
      if (lsynslp) then
         if (lfirst) cline= trim(cline)//' slp '
         call ncvar_read(fpres,vpres,era40_pres,nlon,nlat,1,irec,irec)
         call interpug(era40_pres ,nlon,nlat,flon,flat,dlon,dlat, &
                          synslp,mlon,plat,itype)   
         synslp=synslp*0.01 ! Pa -> mBar (HPa)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Relative humidity
      if (lsynrelhum) then

         if (lfirst) cline= trim(cline)//' relhum '
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
         call interpug(era40_rhum ,nlon,nlat,flon,flat,dlon,dlat, &
                          synrelhum,mlon,plat,itype)   
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SW radiation

      if (lsynshwflx) then
         if (lfirst) cline= trim(cline)//' ssrd '

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


         call interpug(era40_ssrd,nlon,nlat,flon,flat,dlon,dlat, &
                          synshwflx,mlon,plat,itype)   
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
end module mod_read_era40
