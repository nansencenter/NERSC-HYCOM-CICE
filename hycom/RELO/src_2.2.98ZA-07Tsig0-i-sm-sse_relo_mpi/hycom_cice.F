      program hycom_cice
#if defined USE_ESMF 
#if !defined(USE_ESMF_5)
#error hycom_cice requires ESMF 5
#endif
c
c --- ESMF driver for HYCOM ocean model and CICE sea-ice model
c --- Cleanup. Now ESMF-5 only
c
c ---
c --- The routine handles cases with different ice flags in HYCOM. 
c --- This makes it possible to run standalone hycom using the same 
c --- executable.
c --- Note that some of the ice component calls are still called even
c --- though cice is not requested, but the main cice part is not run.
c ---
      use ESMF
      use mod_esmf_utils
      use mod_hycom, only :
     &      OCN_SetServices => HYCOM_SetServices

      use CICE_comp_esmf, only :
     &      ICE_SetServices => CICE_SetServices

      use mod_OICPL, only :
     &      CPL_i2o         => ice2ocn_phase,
     &      CPL_o2i         => ocn2ice_phase,
     &      CPL_SetServices => OICPL_SetServices
c
      implicit none
c
c --- Local variables
c
c --- Gridded Components
      type(ESMF_GridComp) :: ocnGridComp,     !HYCOM as an ESMF component
     &                       iceGridComp      ! CICE as an ESMF component
c
c --- Coupler Components
      type(ESMF_CplComp)  :: o2iCplComp
c
c --- States, Virtual Machines, and Layouts
      type(ESMF_State)    :: ocnImpState,     ! HYCOM import state
     &                       ocnExpState,     ! HYCOM export state
     &                       iceImpState,     ! CICE  import state
     &                       iceExpState,     ! CICE  export state
     &                       cplImpState,     ! OICPL import state
     &                       cplExpState      ! OICPL export state
c
      type(ESMF_VM) :: worldVM
      integer :: petCount, localPet, split
c
c --- Calendars and clocks
      type(ESMF_Clock) :: worldClock
      type(ESMF_Clock) :: ocnClock
      type(ESMF_Clock) :: iceClock
      type(ESMF_Time) :: startTime, stopTime
      type(ESMF_Time) :: ocnTime, iceTime
      type(ESMF_TimeInterval) :: timeStep
c
c --- Return codes for error checks
      integer :: rc,rc2
c
c --- Miscellaneous
      integer :: i,its,icpl,iday !,its_ocn,its_ice
c
c --- KAL - Moved from module acces to state access
      integer :: OCN_nts_cpl, OCN_iceflg
      logical :: OCN_put_export, OCN_get_import, OCN_end_of_run
      logical :: ICE_get_import, ICE_put_export
      integer :: ICE_nts_cpl
      logical :: ICE_restart, OCN_restart
      character(len=256) :: msg,tstr
      integer :: ICE_ktherm
      character(len=256) :: ICE_tfrz_option
c
c
c-------------------------------------------------------------------------------
c  Initialize the ESMF Framework
c-------------------------------------------------------------------------------
c
c --- Set default calendar and log type; get world VM
      rc = ESMF_Success
      call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN,
     &                      logkindflag=ESMF_LOGKIND_MULTI,
     &                                  vm=worldVM,
     &                                  rc=rc)
      if (rc .ne. ESMF_SUCCESS) stop 99
c
c --- Get VM info
      call ESMF_VMGet(worldVM, petCount=petCount, localPET=localPet,
     &                rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="ESMF_VMGet failed", rcToReturn=rc2))
     &   goto 10
c
c-------------------------------------------------------------------------------
c --- Create section
c-------------------------------------------------------------------------------
c
c --- Create the OCEAN gridded component
      ocnGridComp = ESMF_GridCompCreate(
     &  name="OCEAN Gridded Component",rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  GridCompCreate failed", rcToReturn=rc2))
     &   goto 10
c
c --- Create empty OCEAN  import/export states
      ocnImpState = ESMF_StateCreate(name="OCEAN Import",
     &    stateintent=ESMF_STATEINTENT_IMPORT,  rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  ImpState Create failed", rcToReturn=rc2))
     &   goto 10
c
      ocnExpState = ESMF_StateCreate(Name="OCEAN Export",
     &    stateIntent=ESMF_STATEINTENT_EXPORT,  rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  ExpState Create failed", rcToReturn=rc2))
     &   goto 10
c
c --- Create the SEAICE gridded component
      iceGridComp = ESMF_GridCompCreate(
     &  name='SEAICE Component', rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE GridCompCreate failed", rcToReturn=rc2))
     &   goto 10
c
c --- Create empty SEAICE import/export states
      iceImpState = ESMF_StateCreate(Name="SEAICE Import",
     &      stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE ImpState Create failed", rcToReturn=rc2))
     &   goto 10
c
      iceExpState = ESMF_StateCreate(Name="SEAICE Export",
     &      stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE ExpState Create failed", rcToReturn=rc2))
     &   goto 10
c
c --- Create the OICPL coupler component
      o2iCplComp = ESMF_CplCompCreate(
     &     name="OICPL Coupler Component", rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPLE CplCompCreate failed", rcToReturn=rc2))
     &   goto 10
c
c --- Create empty OICPL import/export states
      cplImpState = ESMF_StateCreate(Name="OICPL Import",
     &      stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL ImpState Create failed", rcToReturn=rc2))
     &   goto 10
c
      cplExpState = ESMF_StateCreate(Name="OICPL Export",
     &      stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL ExpState Create failed", rcToReturn=rc2))
     &   goto 10
c
c --- Add OCEAN and SEAICE states to OICPL states
      CALL ESMF_StateAdd(cplImpState, (/ocnImpState/), rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL: Add OCEAN  impState failed", rcToReturn=rc2))
     &   goto 10
c
      CALL ESMF_StateAdd(cplImpState, (/iceImpState/), rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL: Add SEAICE impState failed", rcToReturn=rc2))
     &   goto 10
c
      CALL ESMF_StateAdd(cplExpState, (/ocnExpState/), rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL: Add OCEAN  expState failed", rcToReturn=rc2))
     &   goto 10
c
      CALL ESMF_StateAdd(cplExpState, (/iceExpState/), rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL: Add SEAICE impState failed", rcToReturn=rc2))
     &   goto 10
c
c-------------------------------------------------------------------------------
c --- Register section
c-------------------------------------------------------------------------------
c
c --- Register the OCEAN  gridded component
      call ESMF_GridCompSetServices(ocnGridComp,
     &                              OCN_SetServices, rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  Registration failed", rcToReturn=rc2))
     &   goto 10
c
c --- Register the SEAICE gridded component
      call ESMF_GridCompSetServices(iceGridComp,
     &                              ICE_SetServices, rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE Registration failed", rcToReturn=rc2))
     &   goto 10
c
c --- Register the OICPL coupler component
      call ESMF_CplCompSetServices(o2iCplComp,
     &                             CPL_SetServices,rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL  Registration failed", rcToReturn=rc2))
     &   goto 10
c
c-------------------------------------------------------------------------------
c --- Initalize Section
c------------------------------------ -------------------------------------------
c --- Init clocks with dummy values
      call ESMF_TimeSet(startTime,yy=1901,mm=1,dd=1,h=0,s=0,rc=rc)
      call ESMF_TimeIntervalSet(timeStep,s_r8=3600.)
      ocnClock = ESMF_ClockCreate(name="HYCOM Clock", 
     &   startTime=startTime, timeStep=timeStep,  rc=rc)
      iceClock = ESMF_ClockCreate(name="CICE Clock",
     &   startTime=startTime, timeStep=timeStep,  rc=rc)
c
c --- Initialize OCEAN  gridded component
      call ESMF_GridCompInitialize(ocnGridComp,
     &                              importState=ocnImpState,
     &                              exportState=ocnExpState,
     &                                    phase=1,
     &                                 syncflag=ESMF_SYNC_BLOCKING,
     &                                    clock=ocnClock,
     &                                       rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN Initialize failed", rcToReturn=rc2))
     &   goto 10
c
c --- Get HYCOM ice flag. 
      call ESMF_AttributeGet(ocnGridComp, name="OCN_iceflg",
     &   value=OCN_iceflg, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &       msg="cice_setup_esmf: attributeget OCN_iceflg", 
     &       rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
c
c --- Get HYCOM clock properties 
      call ESMF_ClockGet(ocnClock,startTime=startTime,
     &  stopTime=stopTime,rc=rc)
      call ESMF_TimeGet(startTime,timestring=tstr)
      write(msg,'("hycom_cice startTime from ocnClock:",a)') trim(tstr)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
      if (localPet .eq. 0 ) print '(a)',trim(msg)
C
      call ESMF_TimeGet(stopTime,timestring=tstr)
      write(msg,'("hycom_cice stopTime from ocnClock :",a)') trim(tstr)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
      if (localPet .eq. 0 ) print '(a)',trim(msg)
c
c --- HYCOM determines the clock - here we copy HYCOMS ocnClock to CICEs iceClock
      iceClock = ESMF_ClockCreate(ocnClock)
c
c --- Initialize SEAICE gridded component
      call ESMF_GridCompInitialize(    gridComp=iceGridComp,
     &                              importState=iceImpState,
     &                              exportState=iceExpState,
     &                                    phase=1,
     &                                 syncflag=ESMF_SYNC_BLOCKING,
     &                                    clock=iceClock,
     &                                       rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE Initialize failed", rcToReturn=rc2))
     &   goto 10
c
c --- Consistency checks for time steps
      if (.not.check_gridcomp_timesteps(ocngridcomp,icegridcomp,
     &   ocnClock, iceClock,
     &   localPet,petCount,ocn_nts_cpl,ice_nts_cpl)) goto 10
c
c --- Consistency checks for grids
      if (.not.check_gridcomp_grids(ocngridcomp,icegridcomp,
     &   localPet,petCount)) goto 10
c
c --- Ice gridcomp carries information about thermodynamic model used.
c --- Copy to ocean grid comp to make consistent freezing point temp calculations
      call ESMF_AttributeGet(iceGridComp, name="CICE_ktherm",
     &   value=ICE_ktherm, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &   msg="cice_setup_esmf: attributeget ktherm", rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
      call ESMF_AttributeSet(ocnGridComp, name="CICE_ktherm",
     &   value=ICE_ktherm, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &   msg="cice_setup_esmf: attributeset ktherm", rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
      call ESMF_AttributeGet(iceGridComp, name="CICE_tfrz_option",
     &   value=ICE_tfrz_option, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &   msg="cice_setup_esmf: attributeget tfrz_option", 
     &   rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
      call ESMF_AttributeSet(ocnGridComp, name="CICE_tfrz_option",
     &   value=ICE_tfrz_option, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &   msg="cice_setup_esmf: attributeset tfrz_option", 
     &   rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
c
c --- Initialize OCEAN  gridded component. Phase 2 after getting info
C --- from CICE component
      call ESMF_GridCompInitialize(ocnGridComp,
     &                              importState=ocnImpState,
     &                              exportState=ocnExpState,
     &                                    phase=2,
     &                                 syncflag=ESMF_SYNC_BLOCKING,
     &                                    clock=ocnClock,
     &                                       rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN Initialize phase 2 failed", rcToReturn=rc2))
     &   goto 10
c
c --- Initialize OICPL coupler component
      call ESMF_CplCompInitialize(     cplComp=o2iCplComp,
     &                             importState=cplImpState,
     &                             exportState=cplExpState,
     &                                   phase=1,
     &                                syncflag=ESMF_SYNC_BLOCKING,
     &                                      rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL  Initialize failed", rcToReturn=rc2))
     &   goto 10
c
c --- Couple SEAICE to OCEAN
      call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                      importState=cplImpState,
     &                      exportState=cplExpState,
     &                            phase=CPL_i2o,
     &                         syncflag=ESMF_SYNC_BLOCKING,
     &                               rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL I2O Run failed", rcToReturn=rc2))
     &   goto 10
c
c --- Couple OCEAN to SEAICE
      call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                      importState=cplImpState,
     &                      exportState=cplExpState,
     &                            phase=CPL_o2i,
     &                         syncflag=ESMF_SYNC_BLOCKING,
     &                               rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL I2O Run failed", rcToReturn=rc2))
     &   goto 10
c
c-------------------------------------------------------------------------------
c --- Run Section
c-------------------------------------------------------------------------------
c
c     TODO: ESMF v 6 gives error in gridcomprun - figure out why...
c --- Run Ocean and SeaIce in lockstep, both looking backwards for imports
      do icpl=1,huge(iday)/2 !until end of run
c
c
c ----- OCEAN
        do its= 1,OCN_nts_cpl !couple period, OCEAN
          if     (localPet.eq.0) then !master
            write(6,'(a,3i8)') 
     &         'OCEAN  run - icpl,its,ocn_nts_cpl = ',
     &         icpl,its,ocn_nts_cpl
          endif
          OCN_get_import = its.eq.1           !import at start of period
          OCN_put_export = its.eq.OCN_nts_cpl !export at   end of period
c   
c ------- Set get_import on ocn import state
          call ESMF_AttributeSet(ocnImpState, 
     &       name="get_import",value=OCN_get_import,rc=rc)
          if (ESMF_LogFoundError(rc,msg=
     &      "hycom_cice: attributeset OCN_get_import", rcToReturn=rc2))
     &      call ESMF_Finalize(rc=rc)
c
c -------- Set put_export on ocn export state
           call ESMF_AttributeSet(ocnExpState, 
     &       name="put_export",value=OCN_put_export,rc=rc)
           if (ESMF_LogFoundError(rc, msg=
     &       "hycom_cice: attributeset OCN_put_export", rcToReturn=rc2))
     &       call ESMF_Finalize(rc=rc)
c
c -------- Run HYCOM model for one step
          call ESMF_GridCompRun(    gridComp=ocnGridComp,
     &                           importState=ocnImpState,
     &                           exportState=ocnExpState,
     &                                 phase=1,
     &                              syncflag=ESMF_SYNC_BLOCKING,
     &                                 clock=ocnClock,
     &                                    rc=rc)
          if (ESMF_LogFoundError(rc,
     &        msg="OCEAN Run failed", rcToReturn=rc2))
     &        goto 10
c
c ------- Advance ocean Clock TODO
          call ESMF_ClockAdvance(ocnClock,rc=rc)
        enddo !its; OCEAN
c
c ----- Get OCN_restart attribute. TODO: Will only work if restart at
c ----- end of coupling interval
        call ESMF_AttributeGet(ocnGridComp, 
     &      name="OCN_restart",value=OCN_restart,rc=rc)
        if (ESMF_LogFoundError(rc, msg=
     &     "hycom_cice: attributeget OCN_restart", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
c
c
c ----- SEAICE. Only run if esmf set in HYCOM blkdat.input
        if (OCN_iceflg .eq. 2) then
        do its= 1,ice_nts_cpl !couple period, SEAICE
          if     ( localPet.eq.0) then !master
            write(6,'(a,3i8)') 
     &         'SEAICE run - icpl,its,ice_nts_cpl = ',
     &         icpl,its,ice_nts_cpl
          endif
          ICE_get_import = its.eq.1           !import at start of period
          ICE_put_export = its.eq.ice_nts_cpl !export at   end of period
c ------- TODO: Will only work if restart at end of coupling interval
          ICE_restart = OCN_restart .and. ICE_put_export
c
c ------- Set get_import on ice import state
          call ESMF_AttributeSet(iceImpState, 
     &      name="get_import",value=ICE_get_import,rc=rc)
          if (ESMF_LogFoundError(rc, msg=
     &       "hycom_cice: attributeset ICE_get_import", rcToReturn=rc2))
     &      call ESMF_Finalize(rc=rc)
c
c ------- Set put_export on ice export state
          call ESMF_AttributeSet(iceExpState, 
     &      name="put_export",value=ICE_put_export,rc=rc)
          if (ESMF_LogFoundError(rc,msg=
     &       "hycom_cice: attributeset ICE_put_export", rcToReturn=rc2))
     &      call ESMF_Finalize(rc=rc)
c
c ------- Set ICE_restart on ice grid comp
          call ESMF_AttributeSet(iceGridComp, 
     &        name="ICE_restart",value=OCN_restart,rc=rc)
          if (ESMF_LogFoundError(rc, msg=
     &       "hycom_cice: attributeset ICE_restart", rcToReturn=rc2))
     &       call ESMF_Finalize(rc=rc)
c
c ------- Run CICE model one step
          call ESMF_GridCompRun(    gridComp=iceGridComp,
     &                           importState=iceImpState,
     &                           exportState=iceExpState,
     &                                 phase=1,
     &                              syncflag=ESMF_SYNC_BLOCKING,
     &                                 clock=iceClock,
     &                                    rc=rc)
          if (ESMF_LogFoundError(rc,
     &        msg="SEAICE Run failed (last half day)", rcToReturn=rc2))
     &         goto 10
c 
c ------- Advance CICE clock
          call ESMF_ClockAdvance(iceClock,rc=rc)
        enddo !its; SEAICE
        endif !if (OCN_iceflg .eq. 2) then
c
c ----- Get end_of_run attribute from ocean model
        call ESMF_AttributeGet(ocnGridComp, 
     &      name="end_of_run",value=OCN_end_of_run,rc=rc)
        if (ESMF_LogFoundError(rc, msg=
     &     "hycom_cice: attributeget OCN_end_of_run", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
c
c ----- use end_of_run, rather than a ESMF Clock (KAL: Why ?)
        if     (OCN_end_of_run) then
           exit !icpl loop
        endif !end_of_run
c
c ----- Couple SEAICE to OCEAN
        call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                        importState=cplImpState,
     &                        exportState=cplExpState,
     &                              phase=CPL_i2o,
     &                           syncflag=ESMF_SYNC_BLOCKING,
     &                                 rc=rc)
        if (ESMF_LogFoundError(rc,
     &      msg="OICPL I2O Run failed", rcToReturn=rc2))
     &     goto 10
c
c ---   Couple OCEAN to SEAICE
        call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                        importState=cplImpState,
     &                        exportState=cplExpState,
     &                              phase=CPL_o2i,
     &                           syncflag=ESMF_SYNC_BLOCKING,
     &                                 rc=rc)
        if (ESMF_LogFoundError(rc,
     &      msg="OICPL I2O Run failed", rcToReturn=rc2))
     &     goto 10
c
c ----- Check that the clock times match. Only if actual coupled model run
        call ESMF_ClockGet(ocnClock,currTime=ocnTime,rc=rc)
        call ESMF_ClockGet(iceClock,currTime=iceTime,rc=rc)
        if (ocnTime .ne. iceTime .and. OCN_iceflg.eq.2) then
c
           write(msg,'("hycom_cice clock mismatch")')
           call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
           if (localPet.eq.0) print '(a)',trim(msg)
c
c -------- Get HYCOM clock properties 
           call ESMF_TimeGet(ocnTime,timestring=tstr)
           write(msg,'("ocnClock:",a)') trim(tstr)
           call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
           if (localPet.eq.0) print '(a)',trim(msg)
c
c -------- Get CICE clock properties 
           call ESMF_TimeGet(iceTime,timestring=tstr)
           write(msg,'("iceClock:",a)') trim(tstr)
           call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
           if (localPet.eq.0) print '(a)',trim(msg)
c
           goto 10
        end if
c
c ----- Print ocnClock (same as ice clock)
        call ESMF_TimeGet(ocnTime,timestring=tstr)
        write(msg,'("ocnClock and iceClock:",a)') trim(tstr)
        call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
        if (localPet.eq.0) print '(a)',trim(msg)
      enddo !icpl
      call ESMF_VMBarrier(worldVM)
c
c-------------------------------------------------------------------------------
c  Finalize Section
c-------------------------------------------------------------------------------
c
c --- Finalize OCEAN gridded component
      call ESMF_GridCompFinalize(    gridComp=ocnGridComp,
     &                            importState=ocnImpState,
     &                            exportState=ocnExpState,
     &                                  phase=1,
     &                               syncflag=ESMF_SYNC_BLOCKING,
     &                                     rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  Finalize failed", rcToReturn=rc2))
     &   goto 10
c
c --- Finalize SEAICE gridded component
      call ESMF_GridCompFinalize(    gridComp=iceGridComp,
     &                            importState=iceImpState,
     &                            exportState=iceExpState,
     &                                  phase=1,
     &                               syncflag=ESMF_SYNC_BLOCKING,
     &                                     rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE Finalize failed", rcToReturn=rc2))
     &   goto 10
c
c --- Finalize OACPL coupler component
      call ESMF_CplCompFinalize(     cplComp=o2iCplComp,
     &                           importState=cplImpState,
     &                           exportState=cplExpState,
     &                                 phase=1,
     &                              syncflag=ESMF_SYNC_BLOCKING,
     &                                    rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL  Finalize failed", rcToReturn=rc2))
     &   goto 10
c
10    continue
      if (localPet .eq. 0) write(6,'(a)') 'Exiting hycom_cice'
      call ESMF_VMBarrier(worldVM)
      call ESMF_Finalize(rc=rc)
#endif
      end program hycom_cice
