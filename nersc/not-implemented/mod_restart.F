      module mod_restart
! --- Removed old stuff used by the old file types

      private 
      public :: read_restart_NERSC, save_restart_NERSC
      contains


c --- Main routine for reading restart files. 
      subroutine read_restart_NERSC(day1,nstepx,dtimex)
      use mod_xc
      use mod_year_info, only : year_info, year_day, refyear
      use mod_random_forcing
      use mod_common_ice
#if defined (ICESTATE)
      !use mod_icestate_io
      use mod_icestate_init
#endif
#if defined (NOR05)
      use m_NOR05_sed_restart
#endif /* NOR05 */
#if defined (ECOSM)
      use m_ECOSM_sed_restart
#endif /* ECOSM */
      use mod_hycom_nersc ! Contains rungen, imem
      implicit none
      integer,          intent(out) :: nstepx ! nstep read from hycom rst
      real*8,           intent(out) :: dtimex ! time  read from hycom rst
      real*8,           intent(in ) :: day1   ! Specified start time
c
      type(year_info)  :: rttmp
      integer :: reslt
      logical :: ltmp
      include 'common_blocks.h'

      call year_day(day1,refyear,rttmp,yrflag)
      if (mnproc==1) print *,dtimex,day1,nstepx
      if (mnproc==1) print *,rttmp
c --- Read one hycom restart file
      call read_restart_mem(rttmp,nstepx,dtimex)
c
c --- First, initialize from model fields
      call iceinit
c
#if defined (ICE) || defined (ICESTATE)
c --- Read ice model restart - standard 1-category ice model
c --- We try to read this even if ICESTATE is defined 
      if (.not.read_restart_ice(rttmp)) then
         if (mnproc==1) write(lp,*)'Error reading ice restart file'
      end if
#endif  /* ICE || ICESTATE*/
c
#if defined (ICESTATE)
c --- First attempt to initialize ICESTATE. Always succeeds if iceinit
c --- was run or if standard ice file was read 
      call icestate_MICOM_init(hicem,ficem,hsnwm,tsrfm,
     &   dp(:,:,1,1)/onem,saln(:,:,1,1),temp(:,:,1,1),rhosnw)
c ---  Second attempt -- Only succeeds if ICESTATE restart file is present
       !ltmp = read_restart_icestate(rungen,rttmp,imem)
       reslt= icestate_read_restart(rttmp)
#endif /* ICESTATE*/
c
#if defined (EVP)
      reslt= evp_read_restart_new(rttmp)
      if (reslt/=0 .and. mnproc==1) then
         print *,'Can not read EVP restart file. I continue however...'
         write(lp,*) 'Result of evp_read_restart is ',reslt
      end if
#endif /* EVP */
c
#if defined (NOR05) || (ECOSM)
      if (trcrin) call read_restart_sed(rttmp)
#endif /* NOR05 or ECOSM */
c
c --- Read random restart
      if (randf) then
         call read_restart_rand_ab(rttmp)
      end if
      end subroutine




c --- Main routine for saving restart files
      subroutine save_restart_NERSC(nstepx,dtimex)
      use mod_xc
      use mod_year_info, only : year_info, year_day, refyear
      use mod_random_forcing
#if defined (NOR05)
      use m_NOR05_sed_restart
#endif /* NOR05 */
#if defined (ECOSM)
      use m_ECOSM_sed_restart
#endif /* ECOSM */
      implicit none
      integer,          intent(in) :: nstepx
      real*8,           intent(in) :: dtimex
c
      type(year_info)  :: rttmp
      integer :: itmp, ierr
      logical :: ltmp
      include 'common_blocks.h'
c
      call year_day(dtimex,refyear,rttmp,yrflag)
c --- Save hycom  restart file
      call save_restart_mem(rttmp,nstep,dtimex)
c
#if defined (ICE) 
c --- Save ice  restart file. Also done when running with ICESTATE
      ltmp=save_restart_ice(rttmp)
#endif
c
#if defined (ICESTATE)
c --- Save icestate  restart file. 
      itmp = icestate_save_restart(rttmp)
      ltmp =save_restart_ice_ab(rttmp)
#endif
c
#if defined (EVP)
c --- Save EVP  restart file. 
      itmp=evp_save_restart_new(rttmp)
#endif
#if defined (NOR05) || (ECOSM)
c --- Save sediment restart file, the 3D biological components are saved as tracers
      if (ntracr>0)  call save_restart_sed(rttmp)
#endif /* NOR05 or ECOSM*/
c
c --- Save random forcing 
      if (randf) call save_restart_rand_ab(rttmp)
c
#if defined (SINGLE_RESTART) /* We dump ssh when this flag is set  - its convenient*/
c --- Also put ssh field in restart file - useful for assimilation
      call restart_appendfield(trim(restartfile(rttmp)),
     &   'ssh     ',0,1,srfhgt/g  ,ierr)
c 
#endif
      end subroutine 

c
c
c --- ================================================================
c --- ================================================================
c ---  ----------- Routines for reading restart files ----------------
c --- ================================================================
c --- ================================================================
c
c

c --- ----------------------------------------------
c --- Read hycom model state from file.
c --- Now the only option is hycom-style .[ab] files
c --- ----------------------------------------------
      subroutine read_restart_mem(rt,nstep0,dtime0)
      use mod_xc
      use mod_year_info, only: year_info,dtime_start
      implicit none
      type(year_info), intent(in)  :: rt
      integer,          intent(out) :: nstep0
      real*8,           intent(out) :: dtime0
c
      character(len=80) :: filename
      logical ex
      integer iday, iyear, ihour
      include 'common_blocks.h'
c
      filename=restartfile(rt)
      call forday(dtime_start,yrflag,iyear,iday,ihour)
      if (mnproc==1) print *,'read_restart:',iyear,iday,ihour
c
c --- Check for existence. Stop 
      inquire(file=trim(filename)//'.a',exist=ex)
      if (.not.ex) then
         if(mnproc==1) 
     &      write(lp,*)' File '//trim(filename)//'.[ab]'//
     &                 ' does not exist'
         call xcstop('(m_read_restart:read_restart_mem)')
         stop '(m_read_restart:read_restart_mem)'
      end if 
c
c --- use hycom restart routine. Modify restart filename 
      flnmrsi=trim(filename)
      if (mnproc==1) write(lp,'(a)',advance='yes') 'reading: '//
     &   trim(filename)//'.[ab]'
      call restart_in(nstep0,dtime0,
     &   trim(flnmrsi)//'.a', trim(flnmrsi)//'.b')
      if (mnproc==1) then
         write(lp,'(a,f10.2)')'read_restart:Restart day is ',dtime0
         write(lp,'(a,i10)') 'read_restart:nstep is ',nstep0
         write(lp,'(a)')'read_restart:rt is '//rt%cyy//' '//rt%cmm
     &      //' '//rt%cdm
      end if
      end subroutine read_restart_mem



c --- Read of old ice format files. real*8, one record for each member.
c --- Should be made obsolete some time in the future. See
c --- SINGLE_RESTART flags which calls another ice restart routine. This
c --- routine appends ice fields to existing hycom files
      function read_restart_ice(rt)
      use mod_xc
      use mod_hycom_nersc
      use mod_year_info, only: year_info
      use mod_common_ice
      implicit none
      logical read_restart_ice
      type(year_info) , intent(in) :: rt
c
      real, dimension(itdm,jtdm) :: 
     &   lficem,lhicem,lhsnwm,lticem,ltsrfm,liceU,liceV
      real*8, dimension(itdm,jtdm)  ::  ! Hard-coded real*8 for file IO
     &   io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm,io_iceU,io_iceV
      character(len=27) icename
      logical ex
      integer i,j,iostat
      real x(10)

#if defined (SINGLE_RESTART) && defined (SINGLE_RESTART_ONLY)
c --- In this case this routine is entirely replaced by the new
c --- IO
      read_restart_ice=read_restart_ice_ab(rt)
      return
#endif
c
c --- Prepare to read old ice filetype
      read_restart_ice=.true.
      icename=fileice(rt)
      if(mnproc==1) write(lp,'(a,i4)',advance='yes') 
     &   'reading: '//trim(icename)//
     & '                     record=',imem
      inquire(file=icename,exist=ex)
      if (ex) then
c ---    TODO:Note that this is done by all tasks
         inquire(iolength=j)
     &      io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm  !,iceU,iceV
         open(10,file=icename,status='old',form='unformatted',
     &           access='direct',recl=j)
            iceU=0.0; iceV=0.0
            io_iceU=0.0; io_iceV=0.0
            read(10,rec=imem,err=200,iostat=iostat)io_ficem,io_hicem,
     &           io_hsnwm,io_ticem,io_tsrfm !,liceU,liceV
  200       continue
         close(10)
c ---    Bail on error
         if (iostat/=0) then
            if (mnproc==1) then
            write(lp,*) 'An error occured while reading '//trim(icename)
            write(lp,*) 'IOSTAT is :',iostat
            end if
            read_restart_ice=.false.
            call xcstop('(read_restart_ice)')
            stop '(read_restart_ice)'
c ---    Spread to all nodes
         else
            lficem=io_ficem
            lhicem=io_hicem
            lhsnwm=io_hsnwm
            lticem=io_ticem
            ltsrfm=io_tsrfm
            call xcaput(lficem,ficem,0) 
            call xcaput(lhicem,hicem,0)
            call xcaput(lhsnwm,hsnwm,0)
            call xcaput(lticem,ticem,0)
            call xcaput(ltsrfm,tsrfm,0)
         end if
c --- File not found - return with read flag== false
      else 
         if (mnproc==1) then
            write(lp,*) 'Restart file for ICE does not exist.'
            write(lp,*) 'I continue with fields from iniice'
         end if
         read_restart_ice=.false.
      endif
      end function read_restart_ice


c --- Read ice fields from hycom restart file. Called if the
c --- SINGLE_RESTART flag is set. These restart files are 
c --- more robust wrt error checking than the old ice files.
      logical function read_restart_ice_ab(rt)
      use mod_xc
      use mod_za
      use mod_common_ice
      use mod_year_info, only: year_info
      use mod_hycom_nersc
      implicit none
      type(year_info) , intent(in) :: rt

      character filebase*80
      integer :: ierr(9)
      real rand

      read_restart_ice_ab=.true.
      filebase=restartfile(rt)
      if(mnproc==1) write(lp,'(a)',advance='yes') 
     &   'read_restart_ice_ab: reading '//trim(filebase)//'.[ab]'
c
c --- Retrieve record indexes for the variables
      call restart_getfield(filebase,'ficem',0,1,ficem,ierr(1))
      call restart_getfield(filebase,'hicem',0,1,hicem,ierr(2))
      call restart_getfield(filebase,'hsnwm',0,1,hsnwm,ierr(3))
      call restart_getfield(filebase,'tsrfm',0,1,tsrfm,ierr(4))
      call restart_getfield(filebase,'ticem',0,1,ticem,ierr(5))
#if defined (TEST_ICE_AGE)
      call restart_getfield(filebase,'fy_frac ',0,1,fy_frac,ierr(6))
      call restart_getfield(filebase,'fy_age  ',0,1,fy_age ,ierr(6))
      call restart_getfield(filebase,'rdg_frac',0,1,rdg_frac,ierr(6))
#endif
#if defined (PARAM_EST)
      call restart_getfield(filebase,'msshb   ',0,1,msshb,ierr(9))
      call restart_getfield(filebase,'sstb    ',0,1,sstb,ierr(9))
#endif
#if defined (ALBSNW_EVOL)
      call restart_getfield(filebase,'albsnwm ',0,1,albsnwm,ierr(7))
#endif 
#if defined (SSNOWD_ICE)
      call restart_getfield(filebase,'hprcp   ',0,1,hprcp,ierr(8))
      call restart_getfield(filebase,'hmelt   ',0,1,hmelt,ierr(8))
#endif       
c
c --- These are essential fields so we halt if they are not found
      if (any(ierr(1:8)/=0)) then
         read_restart_ice_ab=.false.
         if (mnproc==1) write (lp,'(a)') 'An error occured when '// 
     &      ' reading the ice fields'
#if ! defined(TEST_ICE_AGE) && ! defined(SSNOWD_ICE) && ! defined (ALBSNW_EVOL)
         call xcstop('(read_restart_ice_ab')
#endif      
      end if
#if defined (PARAM_EST)
      if (ierr(9)/=0) then
         read_restart_ice_ab=.true.
             call RANDOM_SEED
             call random_number(rand)
             msshb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)=-0.2+rand*0.4
             call random_number(rand)
             sstb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)=-2+rand*4
      end if
#endif       
      end function read_restart_ice_ab


c --- Read random forcing files. The random forcing fields are
c --- located after hycom fields in the standard hycom restart file
      subroutine read_restart_rand_ab(rt)
      use mod_xc
      use mod_random_forcing
      use mod_year_info, only: year_info
      implicit none
      type(year_info) ,     intent(in)  :: rt
c
      character filebase*80
      integer :: ierr(10)
c 
c --- File name
      filebase=trim(restartfile(rt))
      if (mnproc==1) then
         write(lp,'(a)')'read_restart_rand_ab reading '
     &      //trim(filebase)//'.[ab]'
      end if
c
c --- Retrieve record indexes for the variables
      call restart_getfield(filebase,'rf_slp'  ,0,1,ran%slp   ,ierr( 1))
      call restart_getfield(filebase,'rf_taux' ,0,1,ran%taux  ,ierr( 2))
      call restart_getfield(filebase,'rf_tauy' ,0,1,ran%tauy  ,ierr( 3))
      call restart_getfield(filebase,'rf_wspd' ,0,1,ran%wndspd,ierr( 4))
      call restart_getfield(filebase,'rf_tair' ,0,1,ran%airtmp,ierr( 5))
      call restart_getfield(filebase,'rf_rlhum',0,1,ran%relhum,ierr( 6))
      call restart_getfield(filebase,'rf_cloud',0,1,ran%clouds,ierr( 7))
      call restart_getfield(filebase,'rf_prcp' ,0,1,ran%precip,ierr( 8))
      call restart_getfield(filebase,'rf_sss'  ,0,1,ran%sss,   ierr( 9))
      call restart_getfield(filebase,'rf_sst'  ,0,1,ran%sst,   ierr(10))
c
c --- These are not essential fields so we continue if they are not found
      if (any(ierr(1:10)/=0)) then
         if (mnproc==1) write (lp,'(a)') 'An error occured when '// 
     &      ' reading the random fields'
c ---    Init random forcing fields if an error occured. NB - also
c ---    resets random seed and fft dims
         call init_random_forcing()
      end if
      end subroutine read_restart_rand_ab



#if defined (EVP)
c --- EVP save/read routines moved in here because
c --- they are appended to HYCOM restart file
      function evp_read_restart_new(rt)
      use mod_xc
      use mod_evp
      use mod_common_ice, only : tauxice,tauyice
#if defined(TEST_ICE_AGE)
     &  ,fy_age,fy_frac,rdg_frac
#endif      
#if defined(PARAM_EST)
     &  ,msshb,sstb
#endif      
      use mod_year_info, only: year_info
      implicit none
      type(year_info) , intent(in) :: rt
c
      character filebase*80
      integer :: evp_read_restart_new
      integer :: ierr(7)
      real :: rand
      evp_read_restart_new=0
      filebase=trim(restartfile(rt))
      if (mnproc==1) then
         write(lp,'(a)')'evp_read_restart_new: reading '
     &      //trim(filebase)//'.[ab]'
      end if
c
C --- Retrieve record indexes for the variables
      call restart_getfield(filebase,'uice'   ,0,1,uvel   ,ierr(1))
      call restart_getfield(filebase,'vice'   ,0,1,vvel   ,ierr(2))
      call restart_getfield(filebase,'tauxice',0,1,tauxice,ierr(3))
      call restart_getfield(filebase,'tauyice',0,1,tauyice,ierr(4))
c 
#if defined (TEST_ICE_AGE)
      call restart_getfield(filebase,'fy_frac ',0,1,fy_frac ,ierr(5))
      call restart_getfield(filebase,'fy_age  ',0,1,fy_age  ,ierr(5))
      call restart_getfield(filebase,'rdg_frac',0,1,rdg_frac,ierr(5))
#endif 
#if defined (PARAM_EST)
      call restart_getfield(filebase,'msshb   ',0,1,msshb   ,ierr(6))
      call restart_getfield(filebase,'sstb    ',0,1,sstb    ,ierr(7))
#endif 
c --- These are not essential fields so we warn if they are not found
      if (any(ierr(1:4)/=0)) then
         if (mnproc==1) write (lp,'(a)') 'An error occured when '// 
     &      ' reading the evp fields'
         evp_read_restart_new=-1
      end if
#if defined (PARAM_EST)
      if (any(ierr(6:7)/=0)) then
             call RANDOM_SEED
             call random_number(rand)
             msshb(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy)=-0.2+rand*0.4
             call random_number(rand)
             sstb(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy)=-2+rand*4
      end if

#endif 
      end function evp_read_restart_new
#endif


#if defined (ICESTATE)
c --- ICESTATE save/read routines moved in here because
c --- they are appended to HYCOM restart file
      function icestate_read_restart(rt)
      use mod_xc
      use mod_year_info, only: year_info
      use mod_icestate
      implicit none
      type(year_info) , intent(in) :: rt ! Time and date info 
c
      character filebase*80
      character(len=2) :: ccat, clev
      integer :: icestate_read_restart, hk, hl
      integer :: e(nthick,11),e2(nthick,nlaymax)
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp
#if defined(ICEAGE)
      real, dimension(nthick) :: amax
#endif      
c
#if defined(ICEAGE)
      amax=0.
#endif
      icestate_read_restart=0
      e=-1
      e2=-1
      filebase=trim(restartfile(rt))
      if (mnproc==1) then
         write(lp,'(a)')'icestate_read_restart: reading '
     &      //trim(filebase)//'.[ab]'
      end if
c
C --- Retrieve record indexes for the variables
      do hk=1,nthick ! loop over thickness categories
c
c ---    Encode category
         write(ccat,'(i2.2)') hk
         call restart_getfield(filebase,'i_nlay'//ccat,0,1,tmp,e(hk,1))
         if (e(hk,1)==0) icestate%ice(hk)%nlay  = nint(tmp) ! Number of layers
c
         call restart_getfield(filebase,'i_qbrn'//ccat,0,1,tmp,e(hk,2))
         if (e(hk,2)==0) icestate%ice(hk)%qstore = tmp      ! Brine heat store
c
         call restart_getfield(filebase,'i_albs'//ccat,0,1,tmp,e(hk,3))
         if (e(hk,3)==0) icestate%ice(hk)%albs   = tmp      ! Surface albedo
c
         call restart_getfield(filebase,'i_fice'//ccat,0,1,tmp,e(hk,4))
         if (e(hk,4)==0) icestate%ice(hk)%fice   = tmp      ! Ice concentration
c
         call restart_getfield(filebase,'i_hice'//ccat,0,1,tmp,e(hk,5))
         if (e(hk,5)==0) icestate%ice(hk)%hice   = tmp      ! Ice thickness
c
         call restart_getfield(filebase,'i_hsnw'//ccat,0,1,tmp,e(hk,6))
         if (e(hk,6)==0) icestate%ice(hk)%hsnw   = tmp      ! Snow thickness
c
         call restart_getfield(filebase,'i_rsnw'//ccat,0,1,tmp,e(hk,7))
         if (e(hk,7)==0) icestate%ice(hk)%rhosnw   = tmp    ! Snow density
c
         call restart_getfield(filebase,'i_tsrf'//ccat,0,1,tmp,e(hk,8))
         if (e(hk,8)==0) icestate%ice(hk)%tsrf     = tmp    ! Surface temperature
#if defined(ICEAGE)   
          ! Ice age (s)
         call restart_getfield(filebase,'i_iage'//ccat,0,1,tmp,e(hk,9))
         if (e(hk,9)==0) icestate%ice(hk)%age = tmp
         amax(hk)=maxval(tmp)
         call xcmaxr(amax(hk))
#endif
#if defined(SSNOWD)
         call restart_getfield(filebase,'i_hpcp'//ccat,0,1,tmp,e(hk,10))
         if (e(hk,10)==0) icestate%ice(hk)%hprcp = tmp

         call restart_getfield(filebase,'i_hmlt'//ccat,0,1,tmp,e(hk,11))
         if (e(hk,11)==0) icestate%ice(hk)%hmelt = tmp
#endif          
c
         do hl=1,nlaymax   ! loop over vertical levels
            write(clev,'(i2.2)') hl
c ---       Read vertical temperature profile
            call restart_getfield(filebase,'i_t'//ccat//'_'//clev,
     &                            0,1,tmp,e2(hk,hl))
            if (e(hk,hl)==0) icestate%ice(hk)%vtp(hl)   = tmp    
         end do
      end do
#if defined(ICEAGE)   
c --- Maximum ice age in the restart file      
         age_restart=maxval(amax)
#endif      
c
c --- TODO:Can do error checking against arrays e and e2 here
c 
      end function icestate_read_restart
#endif


c
c
c --- ================================================================
c --- ================================================================
c ---  ----------- Routines for saving restart files ----------------
c --- ================================================================
c --- ================================================================
c
c


c --- Read of old ice format files. real*8, one record for each member.
c --- Should be made obsolete some time in the futore
      subroutine save_restart_mem(rt,nstep0,dtime0)
      use mod_year_info, only: year_info
      use mod_xc
      use mod_za
      implicit none
      type(year_info),  intent(in) :: rt
      integer,          intent(in) :: nstep0
      real*8,           intent(in) :: dtime0
c
      character(len=80) filename 
      include 'common_blocks.h'
c
c --- Use hycom restart routine. Modify restart filename
      filename=trim(restartfile(rt))
      flnmrso=trim(filename) 
      if (mnproc==1)  write(lp,*) 'Dumping to file '//
     &   trim(flnmrso)//'.[ab]'
      call restart_out(nstep0,dtime0, 
     &     trim(flnmrso),trim(flnmrso),.true.)
c
c --- Hmm .. not closed automatically ...
c     call zaiocl(12)
c     if (mnproc==1)  then
c        close(12)
c     end if
      end subroutine save_restart_mem
c
c
c
c --- Save of old ice format files. real*8, one record for each member.
c --- Should be made obsolete some time in the future. See
c --- SINGLE_RESTART flags which calls another ice restart routine. That
c --- routine appends ice fields to existing hycom files
      function save_restart_ice(rt)
      use mod_xc
      use mod_year_info, only: year_info
      use mod_common_ice
      use mod_hycom_nersc
      implicit none
      type(year_info),       intent(in)  :: rt
c
      real, dimension(itdm,jtdm)  :: 
     &   tficem,thicem,thsnwm,tticem,ttsrfm,ticeU,ticeV
      real*8, dimension(itdm,jtdm)  ::  ! Hard-coded real*8 for file IO
     &   io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm,io_iceU,io_iceV
      logical :: save_restart_ice
      character(len=27) icename
      integer i,j, iostat, ioproc
      real :: a(1)

#if defined (SINGLE_RESTART) /* We always dump when this flag is set */
      save_restart_ice=save_restart_ice_ab(rt)
#if defined (SINGLE_RESTART_ONLY) /* Short-circuit when this flag is set */
      return
#endif /* SINGLE_RESTART_ONLY */
#endif /* SINGLE_RESTART */

      ioproc=1

      ! NB get stop - signals valid for MPI runs as well !

      ! Put ice data into non-tiled arrays
      call xcaget(tficem,ficem,0) ; io_ficem=tficem
      call xcaget(thicem,hicem,0) ; io_hicem=thicem
      call xcaget(thsnwm,hsnwm,0) ; io_hsnwm=thsnwm
      call xcaget(tticem,ticem,0) ; io_ticem=tticem
      call xcaget(ttsrfm,tsrfm,0) ; io_tsrfm=ttsrfm
      call xcaget(ticeU ,iceU ,0) ; io_iceu =ticeu 
      call xcaget(ticeV ,iceV ,0) ; io_icev =ticev

      ! Only one node is alloved to access this file (only important for MPI)
      save_restart_ice=.true.
      if (mnproc==ioproc) then
         icename=fileice(rt)
         write(lp,'(a,i4)',advance='yes')'Saving : '//
     &      trim(icename)//'                     record=',imem
         !inquire(iolength=j)tficem,thicem,thsnwm,tticem,ttsrfm  !,ticeU,ticeV
         inquire(iolength=j)
     &      io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm  !,ticeU,ticeV
         !print *,'iolength=',j ! Not important
         open(10,file=icename,status='unknown',form='unformatted',
     &      access='direct',recl=j)
         write(10,rec=imem,err=200,iostat=iostat)io_ficem,io_hicem,
     &      io_hsnwm,io_ticem,io_tsrfm   !,io_iceU,io_iceV
200      continue
         close(10)
         if (iostat==0) then
            !write(lp,*) '  ICE done'
            save_restart_ice=.true.
         else
            write(lp,'(a,i4)')
     &         'An Error encountered while writing record ',imem
            save_restart_ice=.false.
            write(lp,'(a,i4)') 'IOSTAT is :',iostat
         end if
      end if

      ! KAL : Should communicate result to all processes in MPI...
      if (save_restart_ice) a(1)=1.
      call xcastr(a,1); save_restart_ice=nint(a(1))==1
      a(1)=iostat  
      call xcastr(a,1); iostat=int(a(1))

      ! KAL - new: bail out on save errors.. Prevent disaster in the
      ! making
      if (iostat/=0) then
         if (mnproc==1) then
            write(lp,*) 'Exiting on iostat save  errors...'
            write(lp,*) 'IOSTAT is ',iostat
         end if
         call xcstop('(save_restart_ice)')
         stop '(save_restart_ice)'
      end if

      end function save_restart_ice




      ! KAL - with the advent of ice types, fy fraction and ice age, 
      ! KAL - it is finally time to drop the old restart file setup for 
      ! KAL - ice variables. The .ab files are much more flexible for
      ! KAL - future changes. This version appends everything 
      ! KAL - to the HYCOM restart file
      logical function save_restart_ice_ab(rt)
      use mod_xc
      use mod_za
      use mod_common_ice
      use mod_year_info, only: year_info
      use mod_hycom_nersc
      implicit none
      type(year_info) , intent(in) :: rt
      character filebase*80
      integer :: ierr

      save_restart_ice_ab=.true.
      filebase=trim(restartfile(rt))
      if (mnproc==1) then
         write(lp,'(a)',advance='yes')'save_restart_ice_ab:  Saving  '//
     &      trim(filebase)//'.[ab]'
      end if
      call restart_appendfield(filebase,'ficem   ',0,1,ficem,ierr)
      call restart_appendfield(filebase,'hicem   ',0,1,hicem,ierr)
      call restart_appendfield(filebase,'hsnwm   ',0,1,hsnwm,ierr)
      call restart_appendfield(filebase,'tsrfm   ',0,1,tsrfm,ierr)
      call restart_appendfield(filebase,'ticem   ',0,1,ticem,ierr)
#if defined (TEST_ICE_AGE)
      call restart_appendfield(filebase,'fy_frac ',0,1,fy_frac,ierr)
      call restart_appendfield(filebase,'fy_age  ',0,1,fy_age ,ierr)
      call restart_appendfield(filebase,'rdg_frac',0,1,rdg_frac,ierr)
#endif
#if defined (PARAM_EST)
      call restart_appendfield(filebase,'msshb   ',0,1,msshb,ierr)
      call restart_appendfield(filebase,'sstb    ',0,1,sstb,ierr)
#endif
#if defined (ALBSNW_EVOL)   
      call restart_appendfield(filebase,'albsnwm ',0,1,albsnwm,ierr)
#endif  
#if defined (SSNOWD_ICE)   
      call restart_appendfield(filebase,'hprcp   ',0,1,hprcp,ierr)
      call restart_appendfield(filebase,'hmelt   ',0,1,hmelt,ierr)
#endif       
      if (ierr/=0) save_restart_ice_ab=.false.
      end function save_restart_ice_ab




      subroutine save_restart_rand_ab(rt)
      use mod_xc
      use mod_random_forcing
      use mod_year_info, only: year_info
      implicit none
      type(year_info)     , intent(in) :: rt
c
      character(len=80)   :: filebase
      integer             :: ierr(10)
c
      filebase=trim(restartfile(rt))
      if (mnproc==1) then 
         write(lp,'(a)',advance='yes') 
     &      'save_restart_rand_ab: saving '//
     &      trim(filebase)//'.[ab]'
      end if
      call restart_appendfield(filebase,'rf_slp  '
     &   ,0,1,ran%slp   ,ierr(1))
      call restart_appendfield(filebase,'rf_taux '
     &   ,0,1,ran%taux  ,ierr(2))
      call restart_appendfield(filebase,'rf_tauy '
     &   ,0,1,ran%tauy  ,ierr(3))
      call restart_appendfield(filebase,'rf_wspd '
     &   ,0,1,ran%wndspd,ierr(4))
      call restart_appendfield(filebase,'rf_tair '
     &   ,0,1,ran%airtmp,ierr(5))
      call restart_appendfield(filebase,'rf_rlhum'
     &   ,0,1,ran%relhum,ierr(6))
      call restart_appendfield(filebase,'rf_cloud'
     &   ,0,1,ran%clouds,ierr(7))
      call restart_appendfield(filebase,'rf_prcp '
     &   ,0,1,ran%precip,ierr(8))
      call restart_appendfield(filebase,'rf_sss  '
     &   ,0,1,ran%sss   ,ierr(9))
      call restart_appendfield(filebase,'rf_sst  '
     &   ,0,1,ran%sst   ,ierr(10))
      if (any(ierr/=0)) then
         if (mnproc==1) write(lp,'(a)') 
     &     'An error occured when saving random forcing'
         call xcstop ('(mod_restart:save_restart_rand_ab)')
      end if
      end subroutine save_restart_rand_ab


#if defined(EVP)
c --- EVP save/read routines moved in here because
c --- they are appended to HYCOM restart file
      integer function evp_save_restart_new(rt)
      use mod_xc
      use mod_za
      use mod_common_ice, only : tauxice,tauyice
#if defined(TEST_ICE_AGE)
     &  ,fy_age,fy_frac,rdg_frac
#endif      
#if defined(PARAM_EST)
     &  ,msshb,sstb
#endif      
      use mod_evp
      use mod_year_info, only: year_info
      implicit none
      type(year_info)     , intent(in) :: rt
c     
      integer :: ierr(9)
      character(len=80) :: filebase
c
      ierr=0
      evp_save_restart_new=0
      filebase=trim(restartfile(rt))
      if (mnproc==1) then
         write(lp,'(a)')'evp_save_restart_new: reading '
     &      //trim(filebase)//'.[ab]'
      end if
      call restart_appendfield(filebase,'uice    ',0,1,uvel    ,ierr(1))
      call restart_appendfield(filebase,'vice    ',0,1,vvel    ,ierr(2))
      call restart_appendfield(filebase,'tauxice ',0,1,tauxice ,ierr(3))
      call restart_appendfield(filebase,'tauyice ',0,1,tauyice ,ierr(4))
#if ! defined (EVP_REDUNDANT)
#if defined (TEST_ICE_AGE)
      call restart_appendfield(filebase,'fy_frac ',0,1,fy_frac,ierr(5))
      call restart_appendfield(filebase,'fy_age  ',0,1,fy_age ,ierr(6))
      call restart_appendfield(filebase,'rdg_frac',0,1,rdg_frac,ierr(7))
#endif
#if defined (PARAM_EST)
      call restart_appendfield(filebase,'msshb   ',0,1,msshb,ierr(8))
      call restart_appendfield(filebase,'sstb    ',0,1,sstb,ierr(9))
#endif
#endif
      evp_save_restart_new=maxval(abs(ierr))
      end function evp_save_restart_new
#endif

#if defined(ICESTATE)
c --- ICESTATE save routines moved in here because
c --- they are appended to HYCOM restart file
      integer function icestate_save_restart(rt)
      use mod_xc
      use mod_icestate
      use mod_year_info, only: year_info
      implicit none
      type(year_info)     , intent(in) :: rt
c     
      integer :: hk,hl,ierr
      character(len=80) :: filebase
      character(len= 2) :: ccat, clev
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp
c
      ierr=0
      icestate_save_restart=0
      filebase=trim(restartfile(rt))
      if (mnproc==1) then
         write(lp,'(a)')'icestate_save_restart: reading '
     &      //trim(filebase)//'.[ab]'
      end if

C --- Retrieve record indexes for the variables
      do hk=1,nthick ! loop over thickness categories
c
c ---    Encode category
         write(ccat,'(i2.2)') hk

         tmp=icestate%ice(hk)%nlay        ! number of vertical layers
         call restart_appendfield(filebase,'i_nlay'//ccat,0,1,tmp,ierr)
c
         tmp=icestate%ice(hk)%qstore      ! Brine heat store
         call restart_appendfield(filebase,'i_qbrn'//ccat,0,1,tmp,ierr)
c
         tmp=icestate%ice(hk)%albs        ! Surface albedo
         call restart_appendfield(filebase,'i_albs'//ccat,0,1,tmp,ierr)
c
         tmp=icestate%ice(hk)%fice        ! Ice concentration
         call restart_appendfield(filebase,'i_fice'//ccat,0,1,tmp,ierr)
c
         tmp=icestate%ice(hk)%hice        ! Ice thickness
         call restart_appendfield(filebase,'i_hice'//ccat,0,1,tmp,ierr)
c
         tmp=icestate%ice(hk)%hsnw        ! Snow thickness
         call restart_appendfield(filebase,'i_hsnw'//ccat,0,1,tmp,ierr)
c
         tmp=icestate%ice(hk)%rhosnw      ! Snow density
         call restart_appendfield(filebase,'i_rsnw'//ccat,0,1,tmp,ierr)
c
         tmp=icestate%ice(hk)%tsrf        ! Surface temperature
         call restart_appendfield(filebase,'i_tsrf'//ccat,0,1,tmp,ierr)
c
#if defined (SSNOWD)
         tmp=icestate%ice(hk)%hprcp        ! Accumulated depth
         call restart_appendfield(filebase,'i_hpcp'//ccat,0,1,tmp,ierr)

         tmp=icestate%ice(hk)%hmelt        ! Melt depth
         call restart_appendfield(filebase,'i_hmlt'//ccat,0,1,tmp,ierr)
#endif
#if defined (ICEAGE)
         tmp=icestate%ice(hk)%age      ! Ice age (years)
         call restart_appendfield(filebase,'i_iage'//ccat,0,1,tmp,ierr)
#endif          
         do hl=1,nlaymax   ! loop over vertical levels
            write(clev,'(i2.2)') hl
c ---       vertical temperature profile
            tmp=icestate%ice(hk)%vtp(hl)   
            call restart_appendfield(filebase,'i_t'//ccat//'_'//clev,
     &                            0,1,tmp,ierr)
         end do
      end do
      end function icestate_save_restart
#endif /*ICESTATE*/

c
c
c --- ================================================================
c --- ================================================================
c ---  ---------------- Auxillary routines --------------------------
c --- ================================================================
c --- ================================================================
c
c
      ! Base name of restart file
      function restartfile(rt)
         use mod_year_info, only: year_info
         use mod_hycom_nersc
         implicit none
         character(len=80) restartfile
         type(year_info) rt
         character(len=3) :: cmem
         integer indx
         restartfile=''
         restartfile(1:3)=rungen(1:3)
         restartfile(4:10)='restart'
         restartfile(11:14)=rt%cyy(1:4)
         restartfile(15:15)='_'
         restartfile(16:18)=rt%cdd(1:3)
         restartfile(19:19)='_'
         restartfile(20:21)=rt%chh(1:2)
         ! If imem > 1 - append member index ...
         if (ensflag) then
            write(cmem,'(i3.3)') imem
            restartfile=restartfile(1:len_trim(restartfile)) 
     &                  //'_mem'//cmem
         end if
      end function restartfile
c
      function fileice(rt)
      use mod_year_info, only: year_info
      use mod_hycom_nersc
      implicit none
         character(len=27) fileice
         type(year_info) rt
         fileice(1:3)=rungen(1:3)
         fileice(4:10)='restart'
         fileice(11:14)=rt%cyy(1:4)
         fileice(15:15)='_'
         fileice(16:18)=rt%cdd(1:3)
         fileice(19:19)='_'
         fileice(20:21)=rt%chh(1:2)
         fileice(22:27)='ICE.uf'
      end function fileice


c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c This retrieves the field corresponding to variable name, layer index and
c time level. Calls restart_getindex first
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine restart_getfield(filebase,varname,ilayer,itime,
     &   field,ierr)
      use mod_za
      implicit none
      character(len=*), intent(in)  :: filebase, varname
      integer,          intent(in)  :: ilayer, itime
      real,             intent(out) :: 
     &   field(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer,          intent(out) :: ierr

      integer :: indx,k
      real    :: bmin, bmax, amin, amax
      integer, parameter ::  nop=711

      ierr = 0
      call restart_getindex(filebase,varname,ilayer,itime,indx,
     &   bmin,bmax)
      if (indx<=0)  then
         if (mnproc==1)
     &      print '(a,2i5)','restart_getfield: Can not find field '// 
     &      trim(varname),ilayer, itime
         ierr=-1
         return
      end if
      if (mnproc==1) write(lp,'(a,i4,a,2i4,1x,a)') 
     &   'restart_getfield: index is ',indx,' for ',ilayer,itime,varname

c --- Safe to open
      call zaiopf(filebase(1:len_trim(filebase))//'.a','old', nop)
c
c --- Skip to indx "indx" and read
      do k=1,indx-1
         call zaiosk(nop)
      end do
      call zaiord(field, ip,.false., amin,amax, nop)
c
      if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. 
     &        abs(amax-bmax).gt.abs(bmax)*1.e-4 ) then
         if     (mnproc.eq.1) then
         write(lp,'(/ a / a,2i3 / a,2i3 /  a,1p3e14.6 / a,1p3e14.6 /)') 
     &         'error - .a and .b files not consistent:',               
     &         'iunit,indx = ',711,indx,                                
     &         'vname,itime,ilayer = '//trim(varname),itime,ilayer,     
     &         '.a,.b min = ',amin,bmin,amin-bmin,                      
     &         '.a,.b max = ',amax,bmax,amax-bmax 
         endif !1st tile
         ierr=-1
c ---    Deemed a critical error - something is way wrong
         call xcstop('(restart_getfield)')
         stop '(restart_getfield)'
      endif
      call zaiocl(nop)
      end subroutine
c  
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c This appends the field with varname, ilayer, itime to the file with base
c name "filebase".  Retrieves last index of .a-file first
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine restart_appendfield(filebase,varname,ilayer,itime,
     &   field,ierr)
      use mod_za
      implicit none
      character(len=*), intent(in)    :: filebase
      character(len=8), intent(in)    :: varname
      integer,          intent(in)    :: ilayer, itime
      real,             intent(inout) :: 
     &   field(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer,          intent(out)   :: ierr

      real    :: xmin, xmax
      integer :: indx,k
      integer, parameter :: nop=711

      ierr=0 ! Do something more meaningful later on
      indx=restart_lastindex(trim(filebase))
      call xcsync(flush_lp)
cdiag if (mnproc==1) print *,indx
      call zaiopf(filebase(1:len_trim(filebase))//'.a','unknown', nop)
      if (mnproc==1) then
         open (unit=nop,file=filebase(1:len_trim(filebase))//'.b', 
     &         status='unknown',action='readwrite',form='formatted')
      end if
      if (mnproc==1) write(lp,'(a,i4,a,2i4,1x,a)') 
     &   'restart_appendfield: index is ',indx, ' for ',
     &   ilayer,itime,varname
c
c --- skip first two lines + indx lines
      do k=1,2+indx
         if (mnproc==1) read(nop,*)
      end do
      do k=1,indx
         call zaiosk(nop)
      end do
c
      call zaiowr(field(1-nbdy,1-nbdy),ip,.false., 
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) write(nop,4100) varname,ilayer,itime,xmin,xmax
      if (mnproc==1) close(nop)
      call zaiocl(nop)
4100  format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
      end subroutine
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c This retrieves the number of the last index in .ab files. Input is file
c base (without .a or .b ending)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function restart_lastindex(filebase)
      use mod_za
      implicit none
      character(len=*), intent(in) :: filebase
      integer :: i1,i2,k,ios, indx
      character(len=8) :: char8
      real :: xmin, xmax, rli
      integer, parameter :: nop=711
      integer :: iunitb, iunita

      restart_lastindex=-1
      rli=-1.
      if (mnproc==1)  then
         open (unit=nop,file=filebase(1:len_trim(filebase))//'.b', 
     &                  status='unknown',action='read',form='formatted')
         ! Skip two first lines
         read(nop,*,iostat=ios) 
         read(nop,*,iostat=ios) 
         if (ios/=0) then
            close(nop)
            restart_lastindex=-1
            goto 100
         end if

         ! Read until EOF
         ios=0
         indx=0
         do while   (ios==0)
            read(nop,4100,iostat=ios) char8,i1,i2, xmin, xmax
            indx=indx+1
         end do
         close(nop)
         restart_lastindex=indx-1
         rli=restart_lastindex
      end if
100   continue
      call xcmaxr(rli) ; restart_lastindex=rli
4100  format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
      end function restart_lastindex
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c This retrieves the .ab-index corresponding to variable name, layer index and
c time level. One task reads.
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine restart_getindex(filebase,varname,ilayer,itime,fldindx,
     &                            bmin,bmax)
      use mod_za
      implicit none
      character(len=*), intent(in)  :: filebase, varname
      integer,          intent(in)  :: ilayer, itime
      integer,          intent(out) :: fldindx
      real,             intent(out) :: bmin, bmax

      integer :: i1,i2,k,ios, indx
      character(len=8) :: char8
      real :: xmin, xmax, rfi
      logical :: match
      integer, parameter :: nop=711

      fldindx=-1 ; rfi=-1.
      bmax=-1e80
      bmin=1e80
      if (mnproc==1) then
         open (unit=nop,file=filebase(1:len_trim(filebase))//'.b', 
     &         status='unknown',action='read',form='formatted')
c
c ---    Skip two first lines (header)
         read(nop,*,iostat=ios) 
         read(nop,*,iostat=ios) 
         if (ios/=0) then
            close(nop)
            goto 100
         end if
c
c ---    Read until match
         ios=0
         indx=0
         match=.false.
         do while   (ios==0)
            read(nop,4100,iostat=ios) char8,i1,i2, xmin, xmax
            indx=indx+1
            if (trim(char8)==trim(varname) .and. i1==ilayer .and. 
     &          i2==itime .and.  ios==0) then
               fldindx=indx; rfi=indx;
               bmin=xmin
               bmax=xmax
               match=.true.
            end if
         end do
         close(nop)
      end if ! mnproc==1
100   continue
      call xcmaxr(rfi) ; fldindx=rfi
      call xcmaxr(bmax) 
      call xcminr(bmin) 
c     if (mnproc==1) write(lp,'(a,i4,a,2i4,1x,a)'), &
c        'restart_getindex is ',fldindx, ' for ',ilayer,itime,varname
4100  format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
      end subroutine restart_getindex
c


      end module mod_restart

