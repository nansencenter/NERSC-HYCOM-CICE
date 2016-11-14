module m_read_ecmwf_nc


! Time info is kept -- speeds up reading
real, allocatable, save, dimension(:), private :: ect_all
integer, save, private :: maxrec_all

contains

   subroutine read_ecmwf_nc(rt,plon,plat,depths)
      use mod_xc
      use mod_za
      use mod_year_info22
      use mod_forcing_nersc
      use mod_atm_func
      use m_bilin_ecmwf2
      use m_ncvar_dims
      use m_ncvar_read
      use m_rotate
      use netcdf
      implicit none

      ! Data catalogue name
      character(len=200), save :: ecmwf_path='./Ecmwf.nc/'
      character(len=200) :: cenv

      ! Variable names - Ecmwf 0.25 x 0.25
      character(len=*), parameter :: vuwnd='U10M'
      character(len=*), parameter :: vvwnd='V10M'
      character(len=*), parameter :: vtair='T2M'
      character(len=*), parameter :: vdewp='D2M'
      character(len=*), parameter :: vprec='TP'
      character(len=*), parameter :: vtcc ='TCC'
      character(len=*), parameter :: vpres='MSL'

      type (year_info), intent(in) :: rt
      real, dimension(idm,jdm), intent(in) :: plon,plat,depths

      ! File names
      character(len=80) :: fprec, ftair, ftcc, fuwnd, fvwnd, fblh, &
         flmsk, fpres, fdewp, fssrd, fssrd2, ficec,fz

      integer  :: i,j,k,irec,maxrec,margin_loc,irec2

      real, dimension(:,:,:), allocatable :: &
         ecmwf_prec, ecmwf_pres, ecmwf_dewp, ecmwf_tair, &
         ecmwf_tcc,  ecmwf_uwnd, ecmwf_vwnd, ecmwf_ssrd
      real, allocatable ::  ecmwf_lmsk(:,:,:), ecmwf_rhum(:,:), tmpfld(:,:), &
         xwnd(:,:,:),ywnd(:,:,:), elon(:,:,:),elat(:,:,:)
      real, dimension(:), pointer :: eclon, eclat,eclon2,eclat2

      real, dimension(idm,jdm) :: mlon
      integer, dimension(idm,jdm) :: ip
      real :: hmin,hmax
      real :: rfac,wndfac,cd_new,w4,wfact,svpair, svpdew, ecmwfcd
      real :: dlon,dlat,flon,flat,llon,llat,dummy(1,1,1),flon2,flat2,dlon2,dlat2
      integer :: nlon,nlat,nlon2,nlat2
      integer :: im1,ip1,jm1,jp1

      logical,save :: lfirst=.true.
      logical      :: ldisplay
      integer,save :: lastyear=-1
      integer*4, external :: system
      integer*4 :: ret

      ldisplay=.false.

      ! Check ECMWF path on first pass
      if (lfirst) then
         call getenv('ECNC_PATH',cenv)
         if (trim(cenv)/='') then ! prefer this path if present
            ecmwf_path=trim(cenv)//'/'
         end if
         ret=system('[ -d '//trim(ecmwf_path)//' ]')
         if (ret /=0 ) then
            print *
            print *,'The directory  '//trim(ecmwf_path)//' does not exist.'
            print *,'Make sure that '
            print *,' - You have linked the ECNC data into this catalogue (as'
            print *,'   Ecmwf.nc), and  that the variable ECNC_PATH is empty'
            print *,' - OR, the variable ECNC_PATH is set to the location of '
            print *,'   the ECNC data'
            call exit(1)
         end if
      end if

      ! De-allocate time variable - new info is available ....
      if (.not. lfirst .and. (rt%iyy/=lastyear) ) then
         print *,'New year -- removing old time info'
         deallocate(ect_all)
      !else
      !   print *,'First time or last year = this year'
      end if



      ! mlon must conform to the lon range of the input data 
      mlon = plon ; 
      !where (mlon >  180.) mlon=mlon-360.
      !where (mlon < -180.) mlon=mlon+360.

      ! Set file names
      fuwnd=trim(ecmwf_path)//'ec_atmo_geo_mad_U10M_'//rt%cyy//'.nc'
      fvwnd=trim(ecmwf_path)//'ec_atmo_geo_mad_V10M_'//rt%cyy//'.nc'
      ftair=trim(ecmwf_path)//'ec_atmo_geo_mad_T2M_'//rt%cyy//'.nc'
      fdewp=trim(ecmwf_path)//'ec_atmo_geo_mad_D2M_'//rt%cyy//'.nc'
      fpres=trim(ecmwf_path)//'ec_atmo_geo_mad_MSL_'//rt%cyy//'.nc'
      ftcc =trim(ecmwf_path)//'ec_atmo_geo_mad_TCC_'//rt%cyy//'.nc'


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
         write(lp,'(a)') 'Ecmwf Allocation error'
         call exit(1)
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ecmwf air temperature   K 
      if (lsynairtmp) then
         call ecnc_get_info(rt,trim(ftair),vtair,eclon,eclat,irec,nlon,nlat,flon,flat,dlon,dlat)

         ! Diagnostic output
         if (mod(irec,4)==0.or.lfirst) then
            ldisplay=.true.
            WRITE(lp,*)
            write(lp,'(a,i6,a,i4,x,i3,x,i2)',advance='yes')'Ecmwf data from nc record: ',irec,&
            '  rt= ',rt%iyy,rt%idd,rt%ihh
         else
            if (mnproc==1) write(lp,'(a1)',advance='no')'.'
            ldisplay=.false.
         end if

         if (ldisplay) write(lp,'(a)') ' Using Ecmwf temperature'

         if  (irec/=-1) then
            allocate(ecmwf_tair(nlon,nlat,1))

            call ncvar_read(ftair,vtair,ecmwf_tair,nlon,nlat,1,irec,irec)
            call bilin_ecmwf2(ecmwf_tair,nlon,nlat,flon,flat,dlon,dlat, &
                             synairtmp,mlon,plat,depths)   
            deallocate(eclon,eclat) 

            ! tair is needed later
            if (.not. lsynrelhum) deallocate(ecmwf_tair)

         else
            print *,'Could not find correct time record '//trim(ftair)
            call exit(1)
         end if
      end if





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ECMWF clouds   []  (0-1)
      if (lsynclouds) then

         call ecnc_get_info(rt,ftcc,vtcc,eclon,eclat,irec,nlon,nlat,flon,flat,dlon,dlat)

         if (irec/=-1) then
            allocate(ecmwf_tcc(nlon,nlat,1))
            if (ldisplay) write(lp,'(a)')' Using Ecmwf clouds '
            call ncvar_read(ftcc,vtcc,ecmwf_tcc,nlon,nlat,1,irec,irec)
            call bilin_ecmwf2(ecmwf_tcc ,nlon,nlat,flon,flat,dlon,dlat, &
                             synclouds,mlon,plat,depths)   
            ip=1
            deallocate(eclon,eclat,ecmwf_tcc)

         else
            print *,'Could not find correct time record '//trim(ftcc)
            call exit(1)
         end if
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ECMWF winds
      if (lsynwnd) then
         if (ldisplay) write(lp,'(a)')' Using Ecmwf winds '
         call ecnc_get_info(rt,fuwnd,vuwnd,eclon,eclat,irec,nlon,nlat,flon,flat,dlon,dlat)
         call ecnc_get_info(rt,fvwnd,vvwnd,eclon2,eclat2,irec2,nlon2,nlat2,flon2,flat2,dlon2,dlat2)

         ! Safety check -- irec and irec2 might be different
         if (nlon/=nlon2 .or.nlat/=nlat2 .or. flon/=flon2 .or. dlon/=dlon2) then
            print *,'Winds from different grid sizes ...'
            call exit(1)
         end if

         if (irec2/=-1.and.irec/=-1) then

            allocate(ecmwf_uwnd(nlon,nlat,1))
            allocate(ecmwf_vwnd(nlon,nlat,1))

            call ncvar_read(fuwnd,vuwnd,ecmwf_uwnd,nlon,nlat,1,irec,irec)
            call ncvar_read(fvwnd,vvwnd,ecmwf_vwnd,nlon,nlat,1,irec2,irec2)


            call bilin_ecmwf2(ecmwf_uwnd ,nlon,nlat,flon,flat,dlon,dlat, &
                             synuwind,mlon,plat,depths)   
            call bilin_ecmwf2(ecmwf_vwnd ,nlon,nlat,flon,flat,dlon,dlat, &
                             synvwind,mlon,plat,depths)   

            ! Rotate winds to model grid 
            call rotate(synuwind, &
                        synvwind, &
                        plat    , &
                        mlon    ,idm,jdm,'l2m')

            ! Derived  wind stress and absolute windspeed
            ecmwfcd=0.0012

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
            ip=1
            deallocate(eclon,eclat,ecmwf_uwnd,ecmwf_vwnd,eclon2,eclat2)
         else
            print *,'Could not find correct time record '//trim(fuwnd)//' '// &
                    trim(fvwnd)
            call exit(1)
         end if
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pressure

         
      if (lsynslp) then
         if (ldisplay) write(lp,'(a)')' Using Ecmwf SLP '

         call ecnc_get_info(rt,fpres,vpres,eclon,eclat,irec,nlon,nlat,flon,flat,dlon,dlat)

         if (irec/=-1) then
            allocate(ecmwf_pres(nlon,nlat,1))
            call ncvar_read(fpres,vpres,ecmwf_pres,nlon,nlat,1,irec,irec)
            call bilin_ecmwf2(ecmwf_pres ,nlon,nlat,flon,flat,dlon,dlat, &
                             synslp,mlon,plat,depths)   
            synslp=synslp*0.01 ! Pa -> mBar (HPa)
            ip=1
            deallocate(eclon,eclat)

            if (.not. lsynrelhum) deallocate(ecmwf_pres)

         else
            print *,'Could not find correct time record '//trim(fpres)
            call exit(1)
         end if
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Relative humidity
      if (lsynrelhum) then
         if (ldisplay) write(lp,'(a)') ' Using Ecmwf-derived relative humidity'


         ! SLP is needed to calc relative humidity (alternative - use slp0)
         if (.not.lsynslp) then
            write(lp,'(a)') ' ERROR: Ecmwf slp is needed to calculate relhum ...  '
            call exit(1)
         end if

         call ecnc_get_info(rt,fdewp,vdewp,eclon,eclat,irec,nlon,nlat,flon,flat,dlon,dlat)

         if (irec/=-1) then
            allocate(ecmwf_dewp(nlon,nlat,1))
            allocate(ecmwf_rhum(nlon,nlat))

            call ncvar_read(fdewp,vdewp,ecmwf_dewp,nlon,nlat,1,irec,irec)

            ! relative humidity
            !$OMP PARALLEL DO PRIVATE (i,j,svpair,svpdew)
            do j=1,nlat
            do i=1,nlon
               svpair=satvap(ecmwf_tair(i,j,1))
               svpdew=satvap(ecmwf_dewp(i,j,1))
               !print *,i,j,svpair,svpdew,ecmwf_pres(i,j,1)
               ecmwf_rhum(i,j)=relhumid(svpair,svpdew,ecmwf_pres(i,j,1))*0.01

               ! Ahem ...
               ecmwf_rhum(i,j)  =max(0.0,min(ecmwf_rhum(i,j),1.0))
            enddo
            enddo
            !$OMP END PARALLEL DO

            call bilin_ecmwf2(ecmwf_rhum ,nlon,nlat,flon,flat,dlon,dlat, &
                             synrelhum,mlon,plat,depths)   

            !call zaiopf('tst7.a','replace',999)
            !call zaiowr(synrelhum,ip,.false.,hmin,hmax,999,.true.)
            !call zaiocl(999)
            deallocate(eclon,eclat,ecmwf_dewp,ecmwf_tair,ecmwf_pres)

         else
            print *,'Could not find correct time record '//trim(fdewp)
            call exit(1)
         end if

      end if
 
       lfirst=.false.
       lastyear=rt%iyy


   END subroutine read_ecmwf_nc

   subroutine ecnc_get_info(rt,filein,varin,eclon,eclat,irec,nlon,nlat,flon,flat,dlon,dlat)
   use mod_year_info22
   use m_ncvar_dims
   use m_ncvar_read
   use netcdf
   implicit none

   type(year_info),  intent(in) :: rt
   character(len=*), intent(in) :: filein, varin
   real, pointer, dimension(:) :: eclon,eclat
   integer,intent(out) :: nlon,nlat,irec
   real,intent(out) :: flon,flat,dlon,dlat

   integer, dimension(NF90_MAX_DIMS) :: dimsizes
   integer :: ndims, recdim,maxrec

   real, dimension(:,:,:), allocatable :: ec_time, ec_lon, ec_lat
   logical :: lmatch,ex
   integer :: i, ndays, nhours, iyear,imonth,iday

   inquire(exist=ex,file=trim(filein))
   !print *,ex,trim(filein)
   if (.not.ex) then
      print *,'ecnc_get_info: Could not find file '//trim(filein)
      call exit(1)
   end if


   ! Get variable dimensions
   call ncvar_dims(filein,varin,dimsizes,ndims,recdim) ; maxrec=dimsizes(recdim)
   !print *,dimsizes(1:ndims),maxrec
   allocate(ec_time(1,1,maxrec))

   ! This assumes all time records are the same in all files!!!
   ! -- done because reading the time variable is slow (for some reason)
   if (.not.allocated(ect_all)) then
      call ncvar_read(filein,'time',ec_time,1,1,maxrec,1,maxrec)
      allocate(ect_all(maxrec))
      ect_all=ec_time(1,1,:)
      maxrec_all=maxrec
      !print *,'hei1'
   else
      ec_time(1,1,:)=ect_all(:)
      !print *,'hei2'
   end if

   ! Security check
   if (maxrec_all/=maxrec) then
      print *,'ecnc_get_info: Maxrec error -- ',maxrec,maxrec_all
      call exit(1)
   end if


   ! Get time, lon, lat
   nlon=dimsizes(1)
   nlat=dimsizes(2)
   allocate(ec_lon (1,1,nlon  ))
   allocate(ec_lat (1,1,nlat  ))
   allocate(eclon(nlon))
   allocate(eclat(nlat))
   call ncvar_read(filein,'lon' ,ec_lon ,nlon,1,1,1,1)
   !print *,'lon done'
   call ncvar_read(filein,'lat' ,ec_lat ,nlat,1,1,1,1)
   !print *,'lat done'
   eclon=ec_lon(1,1,:)
   eclat=ec_lat(1,1,:)
   

   flon=eclon(1)
   flat=eclat(1)
   dlon=eclon(2)-eclon(1)
   dlat=eclat(2)-eclat(1)

   ! ec_time is in hours since 1950-1-1 - this is days
   irec=-1
   lmatch=.false.
   do i=1,maxrec
      ndays  = floor(ec_time(1,1,i)/24.)
      nhours = (ec_time(1,1,i)-ndays*24)
      call juliantodate(ndays,iyear,imonth,iday,1950,1,1)

      !print '(a,5i5,f10.2)','nc file ',i,iyear,imonth,iday,nhours,ec_time(1,1,i)
      !print '(a,5i5)','rt      ',i,rt%iyy,rt%imm,rt%idm,rt%ihh

      ! KAL --- Remember that "rt%idm" starts from zero
      if (rt%iyy  ==iyear .and. rt%imm==imonth .and.  &
          rt%idm+1==iday  .and. rt%ihh==nhours ) then
         lmatch=.true.
         irec=i
      end if

      !print *,rt%iyy,rt%imm,rt%idm,rt%ihh
      !print *,iyear,imonth,iday,nhours,lmatch,irec
   end do

   deallocate(ec_time,ec_lon,ec_lat)

   end subroutine

end module m_read_ecmwf_nc
