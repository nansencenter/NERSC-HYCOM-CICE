module m_read_ecmwf
contains
subroutine read_ecmwf(rt,plon,plat,depths)
   use mod_xc
   use mod_forcing_nersc
   use mod_year_info22
   use mod_atm_func
   use m_interpug
   implicit none
   !character(len=5), intent(in) :: rforce ! used from mod_forcing_nersc
   type(year_info) , intent(in) :: rt
   real, dimension(idm,jdm), intent(in) :: plon,plat,depths

   !integer, parameter :: nxd=320     ! data i-dim DO NOT CHANGE
   !integer, parameter :: nyd=160     ! data j-dim DO NOT CHANGE
   integer :: nxd, nyd

   integer mo,im1,ip1,jm1,jp1
   real fac,w4,wfact
   logical ex


   real*4, allocatable :: tmp(:,:) 
   real, allocatable :: tem(:,:) 
   real, allocatable :: dew(:,:) 
   real, allocatable :: msl(:,:) 
   real, allocatable :: uuu(:,:) 
   real, allocatable :: vvv(:,:) 

   real, allocatable :: svpair(:,:) 
   real, allocatable :: svphum(:,:) 
   real, allocatable :: svpdew(:,:) 
   real, allocatable :: rhum(:,:) 

   real :: mlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

   real lonref,dlon
   real latref,dlat
   logical, save :: lfirst=.true.

   integer i,j,irec,k
   integer imax(2),imin(2)
   real maxtem,maxdew,maxmsl,maxuuu,maxvvv
   real mintem,mindew,minmsl,minuuu,minvvv

   character(len=3) month(12)
   character(len=2) year
   character(len=3) desc(4)
   character(len=7) tag7
   character(len=2) tag2
   character(len=40) string
   character(len= 3) cproc

   real wndfac,cd_new
   !real,parameter :: airdns  =    1.2
   !real,parameter :: cd      =    0.0012

!  define WINDDRAG results in the use of a winddependent winddrag.
!  compute air-water wind stresses  tau_x  and  tau_y  from the surface wind
!  components. surface wind components have dim [m/s], and are defined on a
!  c-grid;  airdns  has dim [kg/m^3], and  cd  is dimensionless.
!   tau_x  and  tau_y  have dim [kg/(m s^2) = N/m^2]
!
!  include wind dependent drag coefficient according to
!   large & pond, jpo, 11, 324-336, 1981
!
!   10^3 cd = 1.2             for winds < 11 m/s
!           = .49+0.065*wnd   for winds > 11 m/s
! 
   character(len=200) :: cenv
   character(len=200), save :: path0
   integer*4, external :: system
   integer*4 :: ret
   integer, parameter :: itype=0 ! Sets interpolation type bilinear

   ! Set up for forcing
   if (rforce=='ecmwf') then
      nxd=320
      nyd=160
      ! Check ECMWF path on first pass
      if (lfirst) then
         path0='./Ecmwfr/'
         call getenv('ECMWF_PATH',cenv)
         if (trim(cenv)/='') then ! prefer this path if present
            path0=trim(cenv)//'/'
         end if
         ret=system('[ -d '//trim(path0)//' ]')
         if (ret/=0 ) then
            print *
            print *,'The directory  '//trim(path0)//' does not exist.'
            print *,'Make sure that '
            print *,' - You have linked the Ecmwfr data into this catalogue, and '
            print *,'   that the variable ECMWF_PATH is empty'
            print *,' - OR, the variable ECMWF_PATH is set to the location of '
            print *,'   the Ecmwf data'
            call exit(1)
         end if
      end if
   else if (rforce=='metno') then
      nxd=720
      nyd=361
      ! Check ECMWF path on first pass
      if (lfirst) then
         path0='./Met.no/'
         call getenv('METNO_PATH',cenv)
         if (trim(cenv)/='') then ! prefer this path if present
            path0=trim(cenv)//'/'
         end if
         ret=system('[ -d '//trim(path0)//' ]')
         if (ret /=0 ) then
            print *
            print *,'The directory  '//trim(path0)//' does not exist.'
            print *,'Make sure that '
            print *,' - You have linked the Met.no data into this catalogue, and '
            print *,'   that the variable METNO_PATH is empty'
            print *,' - OR, the variable METNO_PATH is set to the location of '
            print *,'   the Met.no data'
            call exit(1)
         end if
      end if
   else
      write(lp,*) 'read_ecmwf only understands "metno" or "ecmwf" '
      write(lp,*) 'forcing ID'
      call flush(lp)
      call xcstop('(read_ecmwf)')
      stop '(read_ecmwf)'
   end if





   ! Check that vars in "mod_synforc" are allocated
   if (.not.allocated(synslp)    .or. .not. allocated(syntaux)   .or. &
       .not.allocated(syntauy)   .or. .not. allocated(synuwind)  .or. &
       .not.allocated(synvwind)  .or. .not. allocated(synwndspd) .or. &
       .not.allocated(synairtmp) .or. .not. allocated(synrelhum) ) then
       write(lp,*) 'Unallocated variables in mod_synforc'
       call flush(lp)
       call xcstop ('(read_ecmwf)')
       stop '(read_ecmwf)'
   end if

   !write(lp,*)'read_ecmw: remember that read_ecvmw works on tiles now....'
   !call flush(lp)


   allocate( tmp(nxd,nyd) )
   allocate( tem(nxd,nyd) )
   allocate( dew(nxd,nyd) )
   allocate( msl(nxd,nyd) )
   allocate( uuu(nxd,nyd) )
   allocate( vvv(nxd,nyd) )
   allocate( svpair(nxd,nyd) )
   allocate( svphum(nxd,nyd) )
   allocate( svpdew(nxd,nyd) )
   allocate( rhum(nxd,nyd) )


   mlon=plon
   where (plon < 0.0) mlon=plon+360.0

   inquire(iolength=j)tmp
   !print *,trim(path0),'ecmwf_',rt%cyy,'_xxx.uf'


! --- Open input files

   inquire(file=trim(path0)//'ecmwf_'//rt%cyy//'_tem.uf',exist=ex)
   if (.not.ex) then
      print '(a)','read_ecmwf:  file does not exist :'//trim(path0)//'ecmwf_'//rt%cyy//'_tem.uf'
      call xcstop ('(read_ecmwf)')
      stop
   endif
   open(21,file=trim(path0)//'ecmwf_'//rt%cyy//'_tem.uf',form='unformatted',access='direct',recl=j)

   inquire(file=trim(path0)//'ecmwf_'//rt%cyy//'_dew.uf',exist=ex)
   if (.not.ex) then
      print '(a)','read_ecmwf:  file does not exist :'//trim(path0)//'ecmwf_'//rt%cyy//'_dew.uf'
      call xcstop ('(read_ecmwf)')
      stop
   endif
   open(22,file=trim(path0)//'ecmwf_'//rt%cyy//'_dew.uf',form='unformatted',access='direct',recl=j)

   inquire(file=trim(path0)//'ecmwf_'//rt%cyy//'_msl.uf',exist=ex)
   if (.not.ex) then
      print '(a)','read_ecmwf:  file does not exist :'//trim(path0)//'ecmwf_',rt%cyy//'_msl.uf'
      call xcstop ('(read_ecmwf)')
      stop
   endif
   open(23,file=trim(path0)//'ecmwf_'//rt%cyy//'_msl.uf',form='unformatted',access='direct',recl=j)

   inquire(file=trim(path0)//'ecmwf_'//rt%cyy//'_uuu.uf',exist=ex)
   if (.not.ex) then
      print '(a)','read_ecmwf:  file does not exist :'//trim(path0)//'ecmwf_'//rt%cyy//'_uuu.uf'
      call xcstop ('(read_ecmwf)')
      stop
   endif
   open(24,file=trim(path0)//'ecmwf_'//rt%cyy//'_uuu.uf',form='unformatted',access='direct',recl=j)

   inquire(file=trim(path0)//'ecmwf_'//rt%cyy//'_vvv.uf',exist=ex)
   if (.not.ex) then
      print '(a)','read_ecmwf:  file does not exist :'//trim(path0)//'ecmwf_'//rt%cyy//'_vvv.uf'
      call xcstop ('(read_ecmwf)')
      stop
   endif
   open(25,file=trim(path0)//'ecmwf_'//rt%cyy//'_vvv.uf',form='unformatted',access='direct',recl=j)


   inquire(file=trim(path0)//'ecmwf_'//rt%cyy//'.hdr',exist=ex)
   if (.not.ex) then
      print '(a)','read_ecmwf:  file does not exist :'//trim(path0)//'ecmwf_'//rt%cyy//'.hdr'
      call xcstop ('(read_ecmwf)')
      stop
   endif

! --- read ecmwf headers
   open(10,file=trim(path0)//'ecmwf_'//rt%cyy//'.hdr',status='old')
      read(10,'(t14,i5)')i
      read(10,'(t14,i5)')j
      read(10,'(t14,2f9.4)')lonref,dlon
      read(10,'(t14,2f9.4)')latref,dlat
   close(10)
         
   irec=rt%idd*4+rt%ihh/6+1

! --- Read ecmwf  data files
   read(21,rec=irec)tmp; tem=dble(tmp) 
   read(22,rec=irec)tmp; dew=dble(tmp)
   read(23,rec=irec)tmp; msl=dble(tmp)
   read(24,rec=irec)tmp; uuu=dble(tmp)
   read(25,rec=irec)tmp; vvv=dble(tmp)

      
! --- Diagnostic output
   write(lp,'(a,i6,a,i4,x,i3,x,i2,a,i3,i3)')'Ecmwfr data from record: ',irec,&
   '  rt= ',rt%iyy,rt%idd,rt%ihh,' month,dayinmonth=',rt%imm,rt%idm+1
   call flush(lp)

! Check fields for consistency
   imax(1:2)=maxloc(tem); maxtem=tem(imax(1),imax(2))
   if (maxtem > 400.0)  then
      if (mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Max value and location for tem: ',maxtem,imax(:)
      call xcstop('(read_ecmwf: max tem too large)')
      stop 'Read_ecmwf: max tem too large'
   end if
   imin(1:2)=minloc(tem); mintem=tem(imin(1),imin(2))
   if (mintem < 100.0)  then
      if (mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Min value and location for tem: ',mintem,imin(:)
      call xcstop('(read_ecmwf: min tem too large)')
      stop 'Read_ecmwf: Min tem too low'
   end if
   imax(1:2)=maxloc(dew); maxdew=dew(imax(1),imax(2))
   if (maxdew > 400.0) then
      if(mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Max value and location for dew: ',maxdew,imax(:)
      call xcstop('(read_ecmwf: max dew too large)')
      stop 'Read_ecmwf: max dew too large'
   end if
   imin(1:2)=minloc(dew); mindew=dew(imin(1),imin(2))
   if (mindew < 100.0) then
      if(mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Min value and location for dew: ',mindew,imin(:)
      call xcstop('(read_ecmwf: min dew too large)')
      stop 'Read_ecmwf: Min dew too low'
   end if
   imax(1:2)=maxloc(msl); maxmsl=msl(imax(1),imax(2))
   if (maxmsl > 1400.0) then
      if(mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Max value and location for msl: ',maxmsl,imax(:)
      call xcstop('(read_ecmwf: max msl too large)')
      stop 'Read_ecmwf: max msl too large'
   end if
   imin(1:2)=minloc(msl); minmsl=msl(imin(1),imin(2))
   if (minmsl < 600.0) then
      if (mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Min value and location for msl: ',minmsl,imin(:)
      call xcstop('(read_ecmwf: min msl too large)')
      stop 'Read_ecmwf: Min msl too low'
   end if
   imax(1:2)=maxloc(uuu); maxuuu=uuu(imax(1),imax(2))
   if (maxuuu > 200.0) then
      if (mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Max value and location for uuu: ',maxuuu,imax(:)
      call xcstop('(read_ecmwf: max uuu too large)')
      stop 'Read_ecmwf: max uuu too large'
   end if
   imin(1:2)=minloc(uuu); minuuu=uuu(imin(1),imin(2))
   if (minuuu < -200.0) then
      if (mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Min value and location for uuu: ',minuuu,imin(:)
      call xcstop('(read_ecmwf: min uuu too large)')
      stop 'Read_ecmwf: Min uuu too low'
   end if
   imax(1:2)=maxloc(vvv); maxvvv=vvv(imax(1),imax(2))
   if (maxvvv > 200.0) then
      if(mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Max value and location for vvv: ',maxvvv,imax(:)
      call xcstop('(read_ecmwf: max vvv too large)')
      stop 'Read_ecmwf: max vvv too large'
   end if
   imin(1:2)=minloc(vvv); minvvv=vvv(imin(1),imin(2))
   if (minvvv < -200.0) then
      if (mnproc==1) print '(a,g13.5,2i5)','read_ecmwf: Min value and location for vvv: ',minvvv,imin(:)
      call xcstop('(read_ecmwf: min vvv too large)')
      stop 'Read_ecmwf: Min vvv too low'
   end if


   !if (mnproc==1) then
   !   write(lp,*)' Skipping bilin calls in read_ecmwf'
   !   call flush(lp)
   !end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sea level pressure
   call interpug(msl,nxd,nyd,lonref,latref,dlon,dlat,synslp,mlon,plat,itype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! air temperature
   call interpug(tem,nxd,nyd,lonref,latref,dlon,dlat,synairtmp,mlon,plat,itype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! windstresses and windspeed
   call interpug(uuu,nxd,nyd,lonref,latref,dlon,dlat,synuwind,mlon,plat,itype)
   call interpug(vvv,nxd,nyd,lonref,latref,dlon,dlat,synvwind,mlon,plat,itype)
   ! --- Rotate winds/stress to model u/v points
   if (wndflag==1) then
       call rotate2(synuwind, synvwind, plat, mlon,'l2m',.false.)
   else
   ! --- Rotate winds/stress to model p-points
       call rotate2(synuwind, synvwind, plat, mlon,'l2m',.true.)
   end if
   do j=1,jdm
      jm1=max(1,j-1)
      jp1=min(j+1,jdm)
      do i=1,idm
         synwndspd(i,j)=sqrt(synuwind(i,j)*synuwind(i,j)+ &
                             synvwind(i,j)*synvwind(i,j))
#ifdef KARALIGHT
         wndfac=MAX(2.5,MIN(32.5,synwndspd(i,j)))
         cd_new = 1.0E-3*(.692 + .0710*wndfac - .000700*wndfac**2)
#else

         wndfac=(1.+sign(1.,synwndspd(i,j)-11.))*.5
         cd_new=(0.49+0.065*synwndspd(i,j))*1.0e-3*wndfac+cd*(1.-wndfac)
#endif
         syntaux(i,j)=synuwind(i,j)*synwndspd(i,j)*airdns*cd_new
         syntauy(i,j)=synvwind(i,j)*synwndspd(i,j)*airdns*cd_new
      enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! relative humidity
   do j=1,nyd
   do i=1,nxd
      msl(i,j)=100.0*msl(i,j)             ! hPa to Pa
      svpair(i,j)=satvap(tem(i,j))
      svpdew(i,j)=satvap(dew(i,j))
      rhum  (i,j)=relhumid(svpair(i,j),svpdew(i,j),msl(i,j))*0.01
      ! Ahem ...
      rhum  (i,j)=max(0.0,min(rhum(i,j),1.0))
   enddo
   enddo
   call interpug(rhum,nxd,nyd,lonref,latref,dlon,dlat,synrelhum(:,:),mlon,plat,itype)
   lfirst=.false.


end subroutine read_ecmwf
end module m_read_ecmwf
