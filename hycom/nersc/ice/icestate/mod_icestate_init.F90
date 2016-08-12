module mod_icestate_init
! This module contains procedures for initializing ICESTATE thermo-
! dynamics when not using ICESTATE restart files.
contains

! ======================================================================
! ===================== icestate_MICOM_init ============================
! ======================================================================
!
!  This routine is used when ICESTATE is started first time from 
!  MICOM thermodynamics. It sets reasonable values for thicknesses
!  temperature profiles etc. etc. Only one iceclass will exist after
!  this routine is done, then the redistribution begins ...
!
   subroutine icestate_MICOM_init(h,f,hs,ts,hml,sml,tml,rsn)
   use mod_icestate
   use mod_icestate_tools
   implicit none

   real,    intent(in), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
      h  , & ! ice thickness      [m]
      f  , &! ice concentr.       []
      hs , &! snow thickness     [m]
      ts , &! surf. temp.        [K]
      hml, &! ML thickness       [m]
      tml, &! ML temperature     [C]
      sml   ! ML salinity      [ppm]
   real,    intent(in)                     :: rsn! snow density  [kg/m^3]

   ! Local variables
   integer     :: i,j,l,hl,hk,p
   type(t_ice) :: icem(nthick)
   real    :: t_f, ficetot, omficetot, rksn 
#ifdef ICESTATE_TEST 
   character (len= 3) :: cnn,cll 
   character (len=30) :: frmt,frmt2

   write (cnn,'(i3.3)') nthick
   frmt= trim(cnn)//'f14.8'
   write (cll,'(i2.2)') nlaymax
   frmt2= trim(cll)//'f14.8'
#endif


   ! Set ice thicknesses etc...
   ! Note that MICOMs fraction (f(i,j)) is the fraction of thick ice.
   ! Here we set h_thick_ice ~ average_h/f, rest is open water
   ! Remember that values were initialized in icestate_init
   !do j=1-margin,jj+margin
   !do i=1-margin,ii+margin
   do j=1-margin,jj+margin
   do l=1,isp(j)
   do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))



#ifdef ICESTATE_TEST
if (i==itst.and.j==jtst) print '(a,2i3,a,'//trim(frmt)//')','pstadv :',i,j,' Redist fice             :',icem%fice
if (i==itst.and.j==jtst) print '(a,2i3,a,'//trim(frmt)//')','pstadv :',i,j,' Redist hice             :',icem%hice
#endif
#if defined(SSNOWD)
      icem(:)        = clear_ice(icem(:),nlay(:),cvsnw(:))
#else
      icem(:)        = clear_ice(icem(:),nlay(:))
#endif 
      ficetot        = f(i,j)
      omficetot      = 1.-f(i,j)
      icem(:)%tsrf   = ts(i,j) 
      icem(:)%qstore = 0.
      t_f = 273.216  - .057 * sml(i,j)
      !t_f = freeze_temp(sml(i,j))

      if (ficetot>epsil1.and.ip(i,j)>0.) then 
         t_f = 273.216  - .057 * sml(i,j)
         !t_f = freeze_temp(sml(i,j))


         hk = nthick
         do while (hk>=1)

            if ( h(i,j)>thickl(hk)) then
               !MICOM's thickness fits into this category...
               icem(hk)%fice   = f(i,j)
               icem(hk)%hice   = h(i,j)  
               icem(hk)%hsnw   = hs(i,j) 
               icem(hk)%rhosnw = rsn
               rksn            = rkice * (icem(hk)%rhosnw / rhow)**1.885
               icem(hk)%vtp(:) = vtp_lin(t_f,icem(hk)%fice,icem(hk)%hice,icem(hk)%hsnw, &
                                         icem(hk)%rhosnw,rksn,icem(hk)%tsrf,icem(hk)%vtp(:), &
                                         nlay(hk))
               icem(hk)%nlay   = nlay(hk) 
               icem(hk)%albs   = (albsmax+albsmin)/2.
#if defined (SSNOWD)
               !We assume the snow is 100% snow covered. 
               icem(hk)%hprcp  = hs(i,j)
               icem(hk)%hmelt  = 0.
               icem(hk)%cv     = cvsnw(hk)
#endif
#if defined (ICEAGE)
               icem(hk)%age    =0. 
#endif
               hk = 0  ! Step out of loop

            end if
            hk = hk -1

         end do
      endif 


      icestate(i,j)%ice(:) = icem(:)
      icestate(i,j)%hml    = hml(i,j)
      icestate(i,j)%sml    = sml(i,j)
      icestate(i,j)%tml    = tml(i,j) + t0deg

#ifdef ICESTATE_TEST
      if (i==itst.and.j==jtst) then
      print '(a,2i3,a,'//trim(frmt)//')','M_init :',i,j,' M_init fice             :',icem%fice
      print '(a,2i3,a,'//trim(frmt)//')','M_init :',i,j,' M_init hice             :',icem%hice
      print '(a,2i3,a,'//trim(frmt)//')','M_init :',i,j,' M_init hsnw             :',icem%hsnw
      print '(a,2i3,a,3f8.3)','M_init :',i,j,' sml,tml,hml :',icestate(i,j)%sml,&
      icestate(i,j)%tml,icestate(i,j)%hml
      do hk=1,nthick
      print '(a,2i3,a,i3,a,'//trim(frmt2)//')','M_init :',i,j,&
                                     'Pvtp class number ',hk,'    :',icem(hk)%vtp(1:nlay(hk))
      end do
      endif
#endif

   enddo
   enddo
   enddo
   end subroutine icestate_MICOM_init




! ======================================================================
! ========================= icestate_init ==============================
! ======================================================================
!  
!  Initialize all variables used in ICESTATE thermodynamics and re-
!  distribution. 
   subroutine icestate_init(tstep,depths,rlat,threfx)
   use mod_icestate
   use mod_icestate_tools
   use mod_icestate_fluxes
   use mod_icestate_redist
   use mod_icestate_srfbudget
   implicit none

   real,            intent(in) :: tstep             ! Timestep  [s]
   real,            intent(in) :: depths(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   ! Ocean depths
   real,            intent(in) :: rlat  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   ! Cell center latitude [rad]
   real,            intent(in) :: threfx            ! Specific volume of ref. density

   real, parameter :: densfac = 1.

   type(t_ice) :: icem(nthick)
   integer     :: i,j,k,start,hk

   ! Read numerical values of latent transfer coefficient
   call clat_turbm(clat)

   ! Set timestep 
   dtt = tstep

   ! Set latitude
   do j=1-nbdy,jj+nbdy
   do i=1-nbdy,ii+nbdy
      rlatm(i,j)=rlat(i,j)
   enddo
   enddo

   ! Set reference specific volume
   thref = threfx / densfac

   ! So the reference density of seawater is [PSI]
   rhoref = 1./thref 

   ! Initialize all fields 
   ! -- Redistribution
   repfunct(:)   = 0.
   stat(:)       = 0.

   ! -- Icefields and ML
#if defined(SSNOWD)
      icem(:)        = clear_ice(icem,nlay,cvsnw)
#else
      icem(:)        = clear_ice(icem,nlay)
#endif   
   do j=1-nbdy,jj+nbdy
   do i=1-nbdy,ii+nbdy
      icestate(i,j)%ice = icem
   enddo
   enddo
   icestate(:,:)%hml = 10.
   icestate(:,:)%sml = 30.
   icestate(:,:)%tml = t0deg

   ! Fluxes to/fro micom.
   Iustar  = 0.
   Ibuoyfl = 0.
   Isalflx = 0.
   Isalflx2= 0.
   Isurflx = 0.
   Ibuoyt  = 0.
   Ibuoysw = 0.
   Ibuoylw = 0.
   divu    = 0.


   ! Determine the coefficient repartition function (repfunct)
   call calc_rep(avkri,devkri,stat,repfunct,stepfunct)

   if (mnproc==1) print '(a,i3)','Number of ICESTATE ice categories :',nthick
   do hk = 1,nthick
      if (hk == nthick) then
         if (mnproc==1)  &
         print '(a,i3,a,i3,a,f10.2,a)',' Category ',hk,' has ',nlay(hk), &
               ' layer(s). Limits for ice category is: ',thickl(hk),'      infinity'
      else
         if (mnproc==1)  &
         print '(a,i3,a,i3,a,2f10.2)',' Category ',hk,' has ',nlay(hk), &
               ' layer(s). Limits for ice category is: ',thickl(hk), thickl(hk+1)
      end if
   end do


   ! Set the "thin ice" mask
   if (hklim>=1.and.hklim<=nthick) then
      thin          = .false.
      thin(1:hklim) = .true.
   else
      if (mnproc==1) then
         print *,'You have set parameter "hklim" to an unrealistic'
         print *,'value in mod_icestate_classes...'
      end if
      call xcstop ('(m_icestate_init)')
      stop 'icestate_init'
   end if
   end subroutine icestate_init


! ======================================================================
! ========================= icestate_read_infile =======================
! ======================================================================
!  
!  Initialize all variables used in ICESTATE thermodynamics and re-
!  distribution. This routine is independent of ICESTATE implemen-
!  tation in MICOM.
!
   subroutine icestate_read_infile
   use mod_icestate
   use mod_icestate_hpar,       ONLY : hice_min
   implicit none
   character(len=*), parameter :: fname_in  ='infile.icestate'
   real :: cversion = 1.1
   real :: fversion =0.0
   integer :: tmp1,tmp2,ios,hor_number
   logical :: ex, l_horizontal_melt

   hor_number=0
   inquire(file=fname_in,exist=ex)
   if (ex) then
      ios=0
      open(10,file=fname_in)
      read(10,*,err=200,end=200,iostat=ios) fversion     !; print *,fversion
      if (fversion/=cversion) then
         if (mnproc==1) then
            print '(a)','icestate_read_infile wrong version'
            print '(a,f4.1)','Version needed : ',cversion
            print '(a,f4.1)','Version of file: ',fversion
         end if
         call xcstop('icestate_read_infile')
         stop '(icestate_read_infile)'
      end if
      read(10,*,err=200,end=200,iostat=ios) albi_mlt     !; print *,albi_mlt
      read(10,*,err=200,end=200,iostat=ios) albi_max     !; print *,albi_max
      read(10,*,err=200,end=200,iostat=ios) albsmin      !; print *,albsmin
      read(10,*,err=200,end=200,iostat=ios) albsmax      !; print *,albsmax
      read(10,*,err=200,end=200,iostat=ios) hice_min     !; print *,hice_min
      read(10,*,err=200,end=200,iostat=ios) fice_max     !; print *,fice_max
      read(10,*,err=200,end=200,iostat=ios) snwlim       !; print *, dump
      read(10,*,err=200,end=200,iostat=ios) qst_frac     !; print *,evp_timers
      read(10,*,err=200,end=200,iostat=ios) hor_number      
      close(10)
200 end if

      
 
   ! Prints message if read failed
   if (ios/=0.or.(.not.ex)) then
      if (mnproc==1)  then
         print *
         print '(a)','Error while reading icestate diagnostics file "'//fname_in//'"'
         print '(a,f4.1)','It is probably wrongly formatted (Your versions is ',fversion
         print '(a,f4.1,a)','and newest version is ',cversion,') or it is nonexistent.'
         print '(a)','A starting point for a new input file can be pasted from '
         print '(a)', 'content between ------- markers below'
         print *
         print *
      end if

      fversion=cversion
      albi_mlt = .6
      albi_max = .71
      albsmin  = .71
      albsmax  = .85
      hice_min = .2
      fice_max = .999
      snwlim   = .5
      qst_frac = .3
      l_horizontal_melt = .true. ; hor_number=0
   end if
      
   if (mnproc==1) then
      write(*,'(a)') '---------------------------------------------------------------------------'
      write(*,'(f3.1,12xa)') fversion,               '# FILE    : Version number of this file'
      write(*,'(f5.3,10xa)')  albi_mlt,              '# ALBEDO  : Albedo value of melting ice []'
      write(*,'(f5.3,10xa)')  albi_max,              '# ALBEDO  : Albedo value of dry ice []'
      write(*,'(f5.3,10xa)')  albsmin ,              '# ALBEDO  : Minimum albedo value of snow []'
      write(*,'(f5.3,10xa)')  albsmax ,              '# ALBEDO  : Maximum albedo value of snow []'
      write(*,'(f5.3,10xa)')  hice_min,              '# FREEZE  : Initial thickness of frozen ice [m]'
      write(*,'(f5.3,10xa)')  fice_max,              '# LEAD    : Maximum value for ice concentration []'
      write(*,'(f5.3,10xa)')  snwlim  ,              '# SNWLIM  : Maximum allowed snow thickness [m]'
      write(*,'(f5.3,10xa)')  qst_frac,              '# QSTORE  : Max heat store in frac. of ice latent heat []'
      write(*,'(i2,13xa)') hor_number , &
         '# LATERAL MELT: 0=none, 1=standard (Drange94), 2=Hakkinen & Mellor  '
      write(*,'(a)') '---------------------------------------------------------------------------'
   end if

   ! Halt on error
   if (ios/=0.or.(.not.ex)) then
      if (mnproc==1) print *,'icestate_read_infile stopping after read error...'
      call xcstop('icestate_read_infile')
      stop '(icestate_read_infile)'
   end if


   hmelt_solar=.false.
   hmelt_hak  =.false.
   if (hor_number==0) then
      if (mnproc==1) print '(a)','ICESTATE:Lateral melt switched off'
   elseif (hor_number==1) then
      hmelt_solar=.true.
      if (mnproc==1) print '(a)','ICESTATE:Switched on lateral melt (Drange94)'
   else if (hor_number==2) then
      hmelt_hak=.true.
      if (mnproc==1) print '(a)','ICESTATE:Switched on lateral melt (Hakkinen & Mellor)'
   else
      if (mnproc==1) print '(a)','ICESTATE:Lateral melt parameterization = 0, 1 or 2'
      call xcstop('icestate_read_infile')
      stop '(icestate_read_infile)'
   end if
   if(mnproc==1) print *
   end subroutine icestate_read_infile



end module mod_icestate_init
