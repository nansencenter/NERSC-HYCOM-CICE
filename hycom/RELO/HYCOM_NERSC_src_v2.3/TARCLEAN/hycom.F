#if defined(USE_ESMF)
      program hycom
c
c --- ESMF driver for stand-alone HYCOM ocean model
c
#if (USE_ESMF_5)
      use ESMF
#else
      use ESMF_Mod
#endif
      use mod_hycom, only : OCEAN_SetServices => HYCOM_SetServices
cKAL  use mod_hycom, only : end_of_run,
cKAL &                      OCEAN_SetServices => HYCOM_SetServices
c
      implicit none
c
c --- Local variables
c
c --- Gridded Components
      type(ESMF_GridComp) :: oceanGridComp
c
c --- States, Virtual Machines, and Layouts
      type(ESMF_VM) :: worldVM
      type(ESMF_State) :: oceanImpState, oceanExpState
      integer :: petCount, localPet, split
c
c --- Calendars and clocks
      type(ESMF_Clock) :: worldClock
      type(ESMF_Clock) :: oceanClock
c
c --- Return codes for error checks
      integer :: rc
c
c --- Miscellaneous
      integer :: i
c
c-------------------------------------------------------------------------------
c  Initialize the ESMF Framework
c-------------------------------------------------------------------------------
c
c --- Set default calendar and log type; get world VM
#if (USE_ESMF_5)
      call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN,
     &                     logKindflag=ESMF_LOGKIND_SINGLE,
     &                     vm=worldVM, rc=rc)
#else
      call ESMF_Initialize(defaultCalendar=ESMF_CAL_GREGORIAN,
     &                     defaultLogType=ESMF_LOG_SINGLE,
     &                     vm=worldVM, rc=rc)
#endif
      if (rc .ne. ESMF_SUCCESS) stop 99
c
c --- Get VM info
      call ESMF_VMGet(worldVM, petCount=petCount, localPET=localPet,
     &                rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="ESMF_VMGet failed", 
     &   rcToReturn=rc)) 
     &   goto 10
#else
      if (ESMF_LogMsgFoundError(rc, "ESMF_VMGet failed", rc)) 
     &   goto 10
#endif
c
c-------------------------------------------------------------------------------
c --- Create section
c-------------------------------------------------------------------------------
c
c --- Create the OCEAN gridded component
#if (USE_ESMF_5)
      oceanGridComp = ESMF_GridCompCreate(
     &                name="OCEAN Gridded Component",
     &                rc=rc)
#else
      oceanGridComp = ESMF_GridCompCreate(vm=worldVM,
     &                name="OCEAN Gridded Component",
     &                gridCompType=ESMF_OCEAN,
     &                rc=rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="OCEAN GridCompCreate failed", 
     &     rcToReturn=rc))
     &   goto 10
#else
      if (ESMF_LogMsgFoundError(rc, "OCEAN GridCompCreate failed", rc))
     &   goto 10
#endif
c
c --- Create empty OCEAN import/export states
#if (USE_ESMF_5)
      oceanImpState = ESMF_StateCreate(
     &   name="OCEAN Import State",
     &   stateIntent=ESMF_STATEINTENT_IMPORT,
     &   rc=rc)
#else
      oceanImpState = ESMF_StateCreate(stateName="OCEAN Import State",
     &                                 stateType=ESMF_STATE_IMPORT,
     &                                 rc=rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="OCEAN ImpState Create failed",
     &   rcToReturn=rc))
     &   goto 10
#else
      if (ESMF_LogMsgFoundError(rc, "OCEAN ImpState Create failed", rc))
     &   goto 10
#endif
#if (USE_ESMF_5)
      oceanExpState = ESMF_StateCreate(
     &   name="OCEAN Export State",
     &   stateIntent=ESMF_STATEINTENT_EXPORT,
     &   rc=rc)
#else
      oceanExpState = ESMF_StateCreate(stateName="OCEAN Export State",
     &                                 stateType=ESMF_STATE_EXPORT,
     &                                 rc=rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="OCEAN ExpState Create failed",
     &     rcToReturn=rc))
     &   goto 10
#else
      if (ESMF_LogMsgFoundError(rc, "OCEAN ExpState Create failed", rc))
     &   goto 10
#endif
c
c-------------------------------------------------------------------------------
c --- Register section
c-------------------------------------------------------------------------------
c
c --- Register the OCEAN gridded component
#if (USE_ESMF_5)
      call ESMF_GridCompSetServices(oceanGridComp,
     &                              OCEAN_SetServices, rc=rc)
#else
      call ESMF_GridCompSetServices(oceanGridComp,
     &                              OCEAN_SetServices, rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="OCEAN Registration failed", 
     &    rcToReturn=rc))
     &   goto 10
#else
      if (ESMF_LogMsgFoundError(rc, "OCEAN Registration failed", rc))
     &   goto 10
#endif
c
c-------------------------------------------------------------------------------
c --- Initalize Section
c-------------------------------------------------------------------------------
c
c --- Initialize OCEAN gridded component
#if (USE_ESMF_5)
      call ESMF_GridCompInitialize(oceanGridComp,
     &                             importState=oceanImpState,
     &                             exportState=oceanExpState,
     &                             clock=worldClock,
     &                             phase=1,
     &                             syncflag=ESMF_SYNC_NONBLOCKING,
     &                             rc=rc)
#else
      call ESMF_GridCompInitialize(gridComp=oceanGridComp,
     &                             importState=oceanImpState,
     &                             exportState=oceanExpState,
     &                             clock=worldClock,
     &                             phase=ESMF_SINGLEPHASE,
     &                             blockingflag=ESMF_NONBLOCKING,
     &                             rc=rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="OCEAN Initialize failed", 
     &      rcToReturn=rc))
     &   goto 10
#else
      if (ESMF_LogMsgFoundError(rc, "OCEAN Initialize failed", rc))
     &   goto 10
#endif
*c
*c --- Get copy of OCEAN clock
*      call ESMF_GridCompGet(oceanGridComp, clock=oceanClock, rc=rc)
*c
*c --- Initialize WORLD clock using OCEAN clock
*      worldClock = ESMF_ClockCreate(clock=oceanClock, rc=rc)
c
c-------------------------------------------------------------------------------
c --- Run Section
c-------------------------------------------------------------------------------
c
      do !until end of run
#if (USE_ESMF_5)
        call ESMF_GridCompRun(oceanGridComp,
     &                        importState=oceanImpState,
     &                        exportState=oceanExpState,
     &                        clock=worldClock,
     &                        phase=1,
     &                        syncflag=ESMF_SYNC_NONBLOCKING,
     &                        rc=rc)
#else
        call ESMF_GridCompRun(gridComp=oceanGridComp,
     &                        importState=oceanImpState,
     &                        exportState=oceanExpState,
     &                        clock=worldClock,
     &                        phase=ESMF_SINGLEPHASE,
     &                        blockingflag=ESMF_NONBLOCKING,
     &                        rc=rc)
#endif
c
c ---   use end_of_run, rather than a ESMF Clock
c ---   KAL - end_of_run is attached to model comp, so we should query it here
        if     (end_of_run) then
          exit
        endif
      enddo
c
c-------------------------------------------------------------------------------
c  Finalize Section
c-------------------------------------------------------------------------------
c
c --- Finalize OCEAN gridded component
#if (USE_ESMF_5)
      call ESMF_GridCompFinalize(oceanGridComp,
     &                           importState=oceanImpState,
     &                           exportState=oceanExpState,
     &                           clock=worldClock,
     &                           phase=1,
     &                           syncflag=ESMF_SYNC_NONBLOCKING,
     &                           rc=rc)
#else
      call ESMF_GridCompFinalize(gridComp=oceanGridComp,
     &                           importState=oceanImpState,
     &                           exportState=oceanExpState,
     &                           clock=worldClock,
     &                           phase=ESMF_SINGLEPHASE,
     &                           blockingflag=ESMF_NONBLOCKING,
     &                           rc=rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="OCEAN Finalize failed", 
     &      rcToReturn=rc)) 
     &   goto 10
#else
      if (ESMF_LogMsgFoundError(rc, "OCEAN Finalize failed", rc)) 
     &   goto 10
#endif
c
10    continue
      call ESMF_VMBarrier(worldVM)
      call ESMF_Finalize(rc=rc)
c
      stop
      end program hycom
#else
      program hycom
c
c --- Non-ESMF driver for stand-alone HYCOM ocean model
c
      use mod_hycom, only : end_of_run,
     &                      HYCOM_Init,
     &                      HYCOM_Run,
     &                      HYCOM_Final
c
      implicit none
c
c --- Initialize HYCOM.
      call HYCOM_Init

c --- Run HYCOM.
      do !until end of run
        call HYCOM_Run
        if     (end_of_run) then
          exit
        endif
      enddo
c
c --- Finalize HYCOM.
      call HYCOM_Final
c
      stop
      end program hycom
#endif
