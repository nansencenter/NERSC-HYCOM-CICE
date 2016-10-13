      program hycom_cice
c
c --- ESMF driver for HYCOM ocean model and CICE sea-ice model
c
#if (USE_ESMF_5)
      use ESMF
#else
      use ESMF_Mod
#endif
      use mod_esmf_utils
      use mod_hycom, only :
cKAL     &      OCN_put_export  => put_export,
cKAL     &      OCN_get_import  => get_import,
!KAL     &      OCN_end_of_run  => end_of_run,
!KAL     &      OCN_nts_day     => nts_day,
!KAL     &      OCN_nts_cpl     => nts_ice,
     &      OCN_SetServices => HYCOM_SetServices


!KAL - uncommented - pass as attributes if necessary
!KAL      use ice_kinds_mod

!KAL      use CICE_ComponentMod, only :
!KAL     &      ICE_SetServices => CICE_SetServices
!KAL  Replaced with the below
      use CICE_comp_esmf, only :
     &      ICE_SetServices => CICE_SetServices

!KAL      use CICE_InitMod, only :
!KAL     &      ICE_nts_day     => nts_day

cKAL  Passed on import / export states
cKAL  use CICE_RunMod, only :
cKAL &      ICE_put_export  => put_export,
cKAL &      ICE_get_import  => get_import
cKAL &      ICE_end_of_run  => end_of_run

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
      type(ESMF_Time) :: startTime, stopTime, currTime
      type(ESMF_TimeInterval) :: timeStep
!KAL  type(ESMF_Grid) :: ocnGrid
!KAL  type(ESMF_Grid) :: iceGrid
!KAL  type(ESMF_TimeInterval) :: iceTimeStep
!KAL  type(ESMF_TimeInterval) :: ocnTimeStep
c
c --- Return codes for error checks
      integer :: rc,rc2
c
c --- ICE coupling frequency
CKAL  integer :: ice_nts_cpl,ocn_cpl_day
c --- Miscellaneous
      integer :: i,its,icpl,iday !,its_ocn,its_ice
c
c --- KAL - Moved from module acces to state access
!KAL  integer :: OCN_nts_day, OCN_nts_cpl
      integer :: OCN_nts_cpl
!KAL  real    :: rOCN_nts_day, rOCN_nts_cpl, rICE_nts_day
      logical :: OCN_put_export, OCN_get_import, OCN_end_of_run
!KAL  integer :: ICE_nts_day, ICE_get_import, ICE_put_export
      logical :: ICE_get_import, ICE_put_export
      integer :: ICE_nts_cpl
      logical :: ICE_restart, OCN_restart
!KAL  integer :: ICE_nx, ICE_ny ! New
!KAL  integer :: OCN_nx, OCN_ny ! New
!KAL  real    :: iceTimeStep_sr8
!KAL  real    :: ocnTimeStep_sr8
      character(len=256) :: msg
      integer :: ICE_ktherm
      character(len=256) :: ICE_tfrz_option

!KAL  integer,dimension(2) :: ocnLBound,ocnUBound
!KAL  integer,dimension(2) :: iceLBound,iceUBound
!KAL  integer :: ocnTileCount, iceTileCount
!KAL  integer :: ocnLocalDeCount, iceLocalDeCount
!KAL  integer :: ocnCoordDimCount(2), iceCoordDimCount(2)
!KAL  integer, target :: ocnMaxIndex(2), iceMaxIndex(2)
!KAL  integer, target :: ocnMinIndex(2), iceMinIndex(2)
!KAL  TYPE(ESMF_Array) :: iceXCoord, iceYCoord, iceMask
!KAL  TYPE(ESMF_Array) :: ocnXCoord, ocnYCoord, ocnMask
!KAL  real(KIND=ESMF_KIND_R4),  pointer ::
!KAL &   xc_ocn(:,:),yc_ocn(:,:),msk_ocn(:,:)



c-------------------------------------------------------------------------------
c  Initialize the ESMF Framework
c-------------------------------------------------------------------------------
c
c --- Set default calendar and log type; get world VM
      rc = ESMF_Success
#if (USE_ESMF_5)
      call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN,
     &                      logkindflag=ESMF_LOGKIND_MULTI,
     &                                  vm=worldVM,
     &                                  rc=rc)
#else
      call ESMF_Initialize(defaultCalendar=ESMF_CAL_GREGORIAN,
     &                      defaultLogType=ESMF_LOG_MULTI,
     &                                  vm=worldVM,
     &                                  rc=rc)
#endif
      if (rc .ne. ESMF_SUCCESS) stop 99
c
c --- Get VM info
      call ESMF_VMGet(worldVM, petCount=petCount, localPET=localPet,
     &                rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="ESMF_VMGet failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "ESMF_VMGet failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c-------------------------------------------------------------------------------
c --- Create section
c-------------------------------------------------------------------------------
c
c --- Create the OCEAN gridded component
      ocnGridComp = ESMF_GridCompCreate(
     &                                name="OCEAN Gridded Component",
#if (USE_ESMF_5)
#else
     &                        gridCompType=ESMF_OCEAN,
#endif
     &                                  rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  GridCompCreate failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "OCEAN  GridCompCreate failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Create empty OCEAN  import/export states
#if (USE_ESMF_5)
      ocnImpState = ESMF_StateCreate(name="OCEAN Import",
     &                              stateintent=ESMF_STATEINTENT_IMPORT,
     &                                      rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  ImpState Create failed", rcToReturn=rc2))
#else
      ocnImpState = ESMF_StateCreate(stateName="OCEAN Import",
     &                               stateType=ESMF_STATE_IMPORT,
     &                                      rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OCEAN  ImpState Create failed", rcToReturn=rc2))
#endif
     &   goto 10

#if (USE_ESMF_5)
      ocnExpState = ESMF_StateCreate(Name="OCEAN Export",
     &                              stateIntent=ESMF_STATEINTENT_EXPORT,
     &                                      rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  ExpState Create failed", rcToReturn=rc2))
#else
      ocnExpState = ESMF_StateCreate(stateName="OCEAN Export",
     &                               stateType=ESMF_STATE_EXPORT,
     &                                      rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OCEAN  ExpState Create failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Create the SEAICE gridded component
      iceGridComp = ESMF_GridCompCreate(
     &                                name='SEAICE Component',
#if (USE_ESMF_5)
#else
     &                        gridcomptype=ESMF_SEAICE,
#endif
     &                                  rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE GridCompCreate failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "SEAICE GridCompCreate failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Create empty SEAICE import/export states
#if (USE_ESMF_5)
      iceImpState = ESMF_StateCreate(Name="SEAICE Import",
     &                              stateintent=ESMF_STATEINTENT_IMPORT,
     &                                      rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE ImpState Create failed", rcToReturn=rc2))
#else
      iceImpState = ESMF_StateCreate(stateName="SEAICE Import",
     &                               stateType=ESMF_STATE_IMPORT,
     &                                      rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "SEAICE ImpState Create failed", rcToReturn=rc2))
#endif
     &   goto 10

#if (USE_ESMF_5)
      iceExpState = ESMF_StateCreate(Name="SEAICE Export",
     &                              stateintent=ESMF_STATEINTENT_EXPORT,
     &                                      rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE ExpState Create failed", rcToReturn=rc2))
#else
      iceExpState = ESMF_StateCreate(stateName="SEAICE Export",
     &                               stateType=ESMF_STATE_EXPORT,
     &                                      rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "SEAICE ExpState Create failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Create the OICPL coupler component
      o2iCplComp = ESMF_CplCompCreate(
     &             name="OICPL Coupler Component",
     &             rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPLE CplCompCreate failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPLE CplCompCreate failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Create empty OICPL import/export states
#if (USE_ESMF_5)
      cplImpState = ESMF_StateCreate(Name="OICPL Import",
     &                              stateintent=ESMF_STATEINTENT_IMPORT,
     &                                      rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL ImpState Create failed", rcToReturn=rc2))
#else
      cplImpState = ESMF_StateCreate(stateName="OICPL Import",
     &                               stateType=ESMF_STATE_IMPORT,
     &                                      rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL ImpState Create failed", rcToReturn=rc2))
#endif
     &   goto 10

#if (USE_ESMF_5)
      cplExpState = ESMF_StateCreate(Name="OICPL Export",
     &                              stateintent=ESMF_STATEINTENT_EXPORT,
     &                                      rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL ExpState Create failed", rcToReturn=rc2))
#else
      cplExpState = ESMF_StateCreate(stateName="OICPL Export",
     &                               stateType=ESMF_STATE_EXPORT,
     &                                      rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL ExpState Create failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Add OCEAN and SEAICE states to OICPL states
#if (USE_ESMF_5)
      CALL ESMF_StateAdd(cplImpState, (/ocnImpState/), rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL: Add OCEAN  impState failed", rcToReturn=rc2))
#else
      CALL ESMF_StateAdd(cplImpState, ocnImpState, rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL: Add OCEAN  impState failed", rcToReturn=rc2))
#endif
     &   goto 10

#if (USE_ESMF_5)
      CALL ESMF_StateAdd(cplImpState, (/iceImpState/), rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL: Add SEAICE impState failed", rcToReturn=rc2))
#else
      CALL ESMF_StateAdd(cplImpState, iceImpState, rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL: Add SEAICE impState failed", rcToReturn=rc2))
#endif
     &   goto 10

#if (USE_ESMF_5)
      CALL ESMF_StateAdd(cplExpState, (/ocnExpState/), rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL: Add OCEAN  expState failed", rcToReturn=rc2))
#else
      CALL ESMF_StateAdd(cplExpState, ocnExpState, rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL: Add OCEAN  expState failed", rcToReturn=rc2))
#endif
     &   goto 10

#if (USE_ESMF_5)
      CALL ESMF_StateAdd(cplExpState, (/iceExpState/), rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL: Add SEAICE impState failed", rcToReturn=rc2))
#else
      CALL ESMF_StateAdd(cplExpState, iceExpState, rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL: Add SEAICE impState failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c-------------------------------------------------------------------------------
c --- Register section
c-------------------------------------------------------------------------------
c
c --- Register the OCEAN  gridded component
      call ESMF_GridCompSetServices(ocnGridComp,
     &                              OCN_SetServices, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  Registration failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "OCEAN  Registration failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Register the SEAICE gridded component
      call ESMF_GridCompSetServices(iceGridComp,
     &                              ICE_SetServices, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE Registration failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "SEAICE Registration failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Register the OICPL coupler component
      call ESMF_CplCompSetServices(o2iCplComp,
     &                             CPL_SetServices,rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL  Registration failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL  Registration failed", rcToReturn=rc2))
#endif
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
#if (USE_ESMF_5)
      call ESMF_GridCompInitialize(ocnGridComp,
     &                              importState=ocnImpState,
     &                              exportState=ocnExpState,
     &                                    phase=1,
     &                                 syncflag=ESMF_SYNC_NONBLOCKING,
     &                                    clock=ocnClock,
     &                                       rc=rc)
#else
      call ESMF_GridCompInitialize(    gridComp=ocnGridComp,
     &                              importState=ocnImpState,
     &                              exportState=ocnExpState,
     &                                    phase=ESMF_SINGLEPHASE,
     &                             blockingflag=ESMF_NONBLOCKING,
     &                                    clock=ocnClock,
     &                                       rc=rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN Initialize failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "OCEAN Initialize failed", rcToReturn=rc2))
#endif
     &   goto 10

      ! Get HYCOM clock properties 
      call ESMF_ClockGet(ocnClock,startTime=startTime,
     &  stopTime=stopTime,rc=rc)
      call ESMF_TimeGet(startTime,timestring=msg)
      write(msg,'("hycom_cice startTime from ocnClock:",a)') trim(msg)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
      if (localPet .eq. 0 ) print '(a)',msg
      call ESMF_TimeGet(stopTime,timestring=msg)
      write(msg,'("hycom_cice stopTime from ocnClock :",a)') trim(msg)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
      if (localPet .eq. 0 ) print '(a)',msg

C --- HYCOM determines the clock - here we copy HYCOMS ocnClock to CICEs iceClock
      iceClock = ESMF_ClockCreate(ocnClock)

c
c --- Initialize SEAICE gridded component
#if (USE_ESMF_5)
      call ESMF_GridCompInitialize(    gridComp=iceGridComp,
     &                              importState=iceImpState,
     &                              exportState=iceExpState,
     &                                    phase=1,
     &                                 syncflag=ESMF_SYNC_NONBLOCKING,
     &                                    clock=iceClock,
     &                                       rc=rc)
#else
      call ESMF_GridCompInitialize(    gridComp=iceGridComp,
     &                              importState=iceImpState,
     &                              exportState=iceExpState,
     &                                    phase=ESMF_SINGLEPHASE,
     &                             blockingflag=ESMF_NONBLOCKING,
     &                                    clock=iceClock,
     &                                       rc=rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE Initialize failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "SEAICE Initialize failed", rcToReturn=rc2))
#endif
     &   goto 10

      if (.not.check_gridcomp_timesteps(ocngridcomp,icegridcomp,
     &   ocnClock, iceClock,
     &   localPet,petCount,ocn_nts_cpl,ice_nts_cpl)) goto 10

      if (.not.check_gridcomp_grids(ocngridcomp,icegridcomp,
     &   localPet,petCount)) goto 10
c
c --- Ice gridcomp carries information about thermodynamic model used.
c --- Copy to ocean grid comp to make freezing point temp calculations
#if (USE_ESMF_5)
      call ESMF_AttributeGet(iceGridComp, name="CICE_ktherm",
     &   value=ICE_ktherm, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &   msg="cice_setup_esmf: attributeget ktherm", rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
#else
      stop '(hycom_cice: ESMF5 only)'
#endif
#if (USE_ESMF_5)
      call ESMF_AttributeSet(ocnGridComp, name="CICE_ktherm",
     &   value=ICE_ktherm, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &   msg="cice_setup_esmf: attributeset ktherm", rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
#else
      stop '(hycom_cice: ESMF5 only)'
#endif
#if (USE_ESMF_5)
      call ESMF_AttributeGet(iceGridComp, name="CICE_tfrz_option",
     &   value=ICE_tfrz_option, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &   msg="cice_setup_esmf: attributeget tfrz_option", 
     &   rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
#else
      stop '(hycom_cice: ESMF5 only)'
#endif
#if (USE_ESMF_5)
      call ESMF_AttributeSet(ocnGridComp, name="CICE_tfrz_option",
     &   value=ICE_tfrz_option, rc=rc)
      if (ESMF_LogFoundError(rc, 
     &   msg="cice_setup_esmf: attributeset tfrz_option", 
     &   rcToReturn=rc2)) 
     &   call ESMF_Finalize(rc=rc)
#else
      stop '(hycom_cice: ESMF5 only)'
#endif


c --- Initialize OCEAN  gridded component. Phase 2 after getting info
C --- from CICE component
#if (USE_ESMF_5)
      call ESMF_GridCompInitialize(ocnGridComp,
     &                              importState=ocnImpState,
     &                              exportState=ocnExpState,
     &                                    phase=2,
     &                                 syncflag=ESMF_SYNC_NONBLOCKING,
     &                                    clock=ocnClock,
     &                                       rc=rc)
#else
      call ESMF_GridCompInitialize(    gridComp=ocnGridComp,
     &                              importState=ocnImpState,
     &                              exportState=ocnExpState,
     &                                    phase=2,
     &                             blockingflag=ESMF_NONBLOCKING,
     &                                    clock=ocnClock,
     &                                       rc=rc)
#endif
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN Initialize phase 2 failed", rcToReturn=rc2))
#else
      if (ESMF_LogMsgFoundError(rc,
     &    "OCEAN Initialize phase 2 failed", rcToReturn=rc2))
#endif
     &   goto 10



c
c --- Initialize OICPL coupler component
#if (USE_ESMF_5)
      call ESMF_CplCompInitialize(     cplComp=o2iCplComp,
     &                             importState=cplImpState,
     &                             exportState=cplExpState,
     &                                   phase=1,
     &                                syncflag=ESMF_SYNC_BLOCKING,
     &                                      rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL  Initialize failed", rcToReturn=rc2))
#else
      call ESMF_CplCompInitialize(     cplComp=o2iCplComp,
     &                             importState=cplImpState,
     &                             exportState=cplExpState,
     &                                   phase=ESMF_SINGLEPHASE,
     &                            blockingflag=ESMF_BLOCKING,
     &                                      rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL  Initialize failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Couple SEAICE to OCEAN
#if (USE_ESMF_5)
      call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                      importState=cplImpState,
     &                      exportState=cplExpState,
     &                            phase=CPL_i2o,
     &                         syncflag=ESMF_SYNC_BLOCKING,
     &                               rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL I2O Run failed", rcToReturn=rc2))
#else
      call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                      importState=cplImpState,
     &                      exportState=cplExpState,
     &                            phase=CPL_i2o,
     &                     blockingflag=ESMF_BLOCKING,
     &                               rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL I2O Run failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Couple OCEAN to SEAICE
#if (USE_ESMF_5)
      call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                      importState=cplImpState,
     &                      exportState=cplExpState,
     &                            phase=CPL_o2i,
     &                         syncflag=ESMF_SYNC_BLOCKING,
     &                               rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL I2O Run failed", rcToReturn=rc2))
#else
      call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                      importState=cplImpState,
     &                      exportState=cplExpState,
     &                            phase=CPL_o2i,
     &                     blockingflag=ESMF_BLOCKING,
     &                               rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL I2O Run failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c-------------------------------------------------------------------------------
c --- Run Section
c-------------------------------------------------------------------------------
c
c --- Run Ocean and SeaIce in lockstep, both looking backwards for imports
      do icpl=1,huge(iday)/2 !until end of run

        !TODO: ESMF v 6 gives error in gridcomprun - figure out why...

c ---   OCEAN
        do its= 1,OCN_nts_cpl !couple period, OCEAN
          if     (localPet.eq.0) then !master
            write(6,'(a,3i8)') 
     &         'OCEAN run - icpl,its,ocn_nts_cpl = ',
     &         icpl,its,ocn_nts_cpl
          endif
          OCN_get_import = its.eq.1           !import at start of period
          OCN_put_export = its.eq.OCN_nts_cpl !export at   end of period
c   
c ---     Set get_import on ocn import state
          call ESMF_AttributeSet(ocnImpState, 
     &       name="get_import",value=OCN_get_import,rc=rc)
#if (USE_ESMF_5)
          if (ESMF_LogFoundError(rc,msg=
     &      "hycom_cice: attributeset OCN_get_import", rcToReturn=rc2))
     &      call ESMF_Finalize(rc=rc)
#else
          if (ESMF_LogMsgFoundError(rc,
     &     "hycom_cice: attributeset OCN_get_import", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
#endif
c
c ---      Set put_export on ocn export state
           call ESMF_AttributeSet(ocnExpState, 
     &       name="put_export",value=OCN_put_export,rc=rc)
#if (USE_ESMF_5)
           if (ESMF_LogFoundError(rc, msg=
     &       "hycom_cice: attributeset OCN_put_export", rcToReturn=rc2))
     &       call ESMF_Finalize(rc=rc)
#else
           if (ESMF_LogMsgFoundError(rc,
     &       "hycom_cice: attributeset OCN_put_export", rcToReturn=rc2))
     &       call ESMF_Finalize(rc=rc)
#endif
c
c ---      Run HYCOM model for one step
#if (USE_ESMF_5)
          call ESMF_GridCompRun(    gridComp=ocnGridComp,
     &                           importState=ocnImpState,
     &                           exportState=ocnExpState,
     &                                 phase=1,
     &                              syncflag=ESMF_SYNC_NONBLOCKING,
     &                                 clock=ocnClock,
     &                                    rc=rc)
          if (ESMF_LogFoundError(rc,
     &        msg="OCEAN Run failed", rcToReturn=rc2))
#else
          call ESMF_GridCompRun(    gridComp=ocnGridComp,
     &                           importState=ocnImpState,
     &                           exportState=ocnExpState,
     &                                 phase=ESMF_SINGLEPHASE,
     &                          blockingflag=ESMF_NONBLOCKING,
     &                                 clock=ocnClock,
     &                                    rc=rc)
          if (ESMF_LogMsgFoundError(rc,
     &        "OCEAN Run failed", rcToReturn=rc2))
#endif
     &        goto 10
c
c ---     ! Advance ocean Clock TODO
          call ESMF_ClockAdvance(ocnClock,rc=rc)
        enddo !its; OCEAN

c---    Get OCN_restart attribute. TODO: Will only work if restart at
c---    end of coupling interval
        call ESMF_AttributeGet(ocnGridComp, 
     &      name="OCN_restart",value=OCN_restart,rc=rc)
#if (USE_ESMF_5)
        if (ESMF_LogFoundError(rc, msg=
     &     "hycom_cice: attributeget OCN_restart", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
#else
        if (ESMF_LogMsgFoundError(rc,
     &     "hycom_cice: attributeget OCN_restart", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
#endif
c
c
c ---   SEAICE. TODO: use icflg 
        do its= 1,ice_nts_cpl !couple period, SEAICE
          if     ( localPet.eq.0) then !master
            write(6,'(a,3i4)') 
     &         'SEAICE run - icpl,its,ice_nts_cpl = ',
     &         icpl,its,ice_nts_cpl
          endif
          ICE_get_import = its.eq.1           !import at start of period
          ICE_put_export = its.eq.ice_nts_cpl !export at   end of period
c---      TODO: Will only work if restart at end of coupling interval
          ICE_restart = OCN_restart .and. ICE_put_export
c
c ---     Set get_import on ice import state
          call ESMF_AttributeSet(iceImpState, 
     &      name="get_import",value=ICE_get_import,rc=rc)
#if (USE_ESMF_5)
          if (ESMF_LogFoundError(rc, msg=
     &       "hycom_cice: attributeset ICE_get_import", rcToReturn=rc2))
     &      call ESMF_Finalize(rc=rc)
#else
          if (ESMF_LogMsgFoundError(rc,
     &     "hycom_cice: attributeset ICE_get_import", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
#endif
c
c ---     Set put_export on ice export state
          call ESMF_AttributeSet(iceExpState, 
     &      name="put_export",value=ICE_put_export,rc=rc)
#if (USE_ESMF_5)
          if (ESMF_LogFoundError(rc,msg=
     &       "hycom_cice: attributeset ICE_put_export", rcToReturn=rc2))
     &      call ESMF_Finalize(rc=rc)
#else
          if (ESMF_LogMsgFoundError(rc,
     &     "hycom_cice: attributeset ICE_put_export", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
#endif
c---      Set ICE_restart on ice grid comp
          call ESMF_AttributeSet(iceGridComp, 
     &        name="ICE_restart",value=OCN_restart,rc=rc)
#if (USE_ESMF_5)
          if (ESMF_LogFoundError(rc, msg=
     &       "hycom_cice: attributeset ICE_restart", rcToReturn=rc2))
     &       call ESMF_Finalize(rc=rc)
#else
          if (ESMF_LogMsgFoundError(rc,
     &       "hycom_cice: attributeset ICE_restart", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
#endif
c
c
c ---     Run CICE model one step
#if (USE_ESMF_5)
          call ESMF_GridCompRun(    gridComp=iceGridComp,
     &                           importState=iceImpState,
     &                           exportState=iceExpState,
     &                                 phase=1,
     &                              syncflag=ESMF_SYNC_NONBLOCKING,
     &                                 clock=iceClock,
     &                                    rc=rc)
          if (ESMF_LogFoundError(rc,
     &        msg="SEAICE Run failed (last half day)", rcToReturn=rc2))
#else
          call ESMF_GridCompRun(    gridComp=iceGridComp,
     &                           importState=iceImpState,
     &                           exportState=iceExpState,
     &                                 phase=ESMF_SINGLEPHASE,
     &                          blockingflag=ESMF_NONBLOCKING,
     &                                 clock=iceClock,
     &                                    rc=rc)
          if (ESMF_LogMsgFoundError(rc,
     &        "SEAICE Run failed (last half day)", rcToReturn=rc2))
#endif
     &         goto 10
c 
c ---     Advance CICE clock
          call ESMF_ClockAdvance(iceClock,rc=rc)
        enddo !its; SEAICE

c
c---    Get end_of_run attribute
        call ESMF_AttributeGet(ocnGridComp, 
     &      name="end_of_run",value=OCN_end_of_run,rc=rc)
#if (USE_ESMF_5)
        if (ESMF_LogFoundError(rc, msg=
     &     "hycom_cice: attributeget OCN_end_of_run", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
#else
        if (ESMF_LogMsgFoundError(rc,
     &     "hycom_cice: attributeget OCN_end_of_run", rcToReturn=rc2))
     &     call ESMF_Finalize(rc=rc)
#endif
cKAL
c ---   use end_of_run, rather than a ESMF Clock (KAL: Why ?)
        if     (OCN_end_of_run) then
           exit !icpl loop
        endif !end_of_run
c
c ---   Couple SEAICE to OCEAN
#if (USE_ESMF_5)
        call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                        importState=cplImpState,
     &                        exportState=cplExpState,
     &                              phase=CPL_i2o,
     &                           syncflag=ESMF_SYNC_BLOCKING,
     &                                 rc=rc)
        if (ESMF_LogFoundError(rc,
     &      msg="OICPL I2O Run failed", rcToReturn=rc2))
#else
        call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                        importState=cplImpState,
     &                        exportState=cplExpState,
     &                              phase=CPL_i2o,
     &                       blockingflag=ESMF_BLOCKING,
     &                                 rc=rc)
        if (ESMF_LogMsgFoundError(rc,
     &      "OICPL I2O Run failed", rcToReturn=rc2))
#endif
     &     goto 10
c
c ---   Couple OCEAN to SEAICE
#if (USE_ESMF_5)
        call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                        importState=cplImpState,
     &                        exportState=cplExpState,
     &                              phase=CPL_o2i,
     &                           syncflag=ESMF_SYNC_BLOCKING,
     &                                 rc=rc)
        if (ESMF_LogFoundError(rc,
     &      msg="OICPL I2O Run failed", rcToReturn=rc2))
#else
        call ESMF_CplCompRun(     cplComp=o2iCplComp,
     &                        importState=cplImpState,
     &                        exportState=cplExpState,
     &                              phase=CPL_o2i,
     &                       blockingflag=ESMF_BLOCKING,
     &                                 rc=rc)
        if (ESMF_LogMsgFoundError(rc,
     &      "OICPL I2O Run failed", rcToReturn=rc2))
#endif
     &     goto 10
c
c ---   Check that the clocks match
        if (ocnClock .ne. iceClock) then
c
           write(msg,'("hycom_cice clock mismatch")')
           call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
           if (localPet.eq.0) print '(a)',msg

c ---      Get HYCOM clock properties 
           call ESMF_ClockGet(ocnClock,currTime=currTime,rc=rc)
           call ESMF_TimeGet(currTime,timestring=msg)
           write(msg,'("ocnClock:",a)') trim(msg)
           call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
           if (localPet.eq.0) print '(a)',msg
c
c ---      Get CICE clock properties 
           call ESMF_ClockGet(iceClock,currTime=currTime,rc=rc)
           call ESMF_TimeGet(currTime,timestring=msg)
           write(msg,'("iceClock:",a)') trim(msg)
           call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
           if (localPet.eq.0) print '(a)',msg
c
           goto 10
        end if
c
c ---   Print ocnClock (same as ice clock)
        call ESMF_ClockGet(ocnClock,currTime=currTime,rc=rc)
        call ESMF_TimeGet(currTime,timestring=msg)
        write(msg,'("ocnClock and iceClock:",a)') trim(msg)
        call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
        if (localPet.eq.0) print '(a)',msg
c

      enddo !icpl
c
      call ESMF_VMBarrier(worldVM)
c
c-------------------------------------------------------------------------------
c  Finalize Section
c-------------------------------------------------------------------------------
c
c --- Finalize OCEAN gridded component
#if (USE_ESMF_5)
      call ESMF_GridCompFinalize(    gridComp=ocnGridComp,
     &                            importState=ocnImpState,
     &                            exportState=ocnExpState,
     &                                  phase=1,
     &                               syncflag=ESMF_SYNC_NONBLOCKING,
     &                                     rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OCEAN  Finalize failed", rcToReturn=rc2))
#else
      call ESMF_GridCompFinalize(    gridComp=ocnGridComp,
     &                            importState=ocnImpState,
     &                            exportState=ocnExpState,
     &                                  phase=ESMF_SINGLEPHASE,
     &                           blockingflag=ESMF_NONBLOCKING,
     &                                     rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OCEAN  Finalize failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Finalize SEAICE gridded component
#if (USE_ESMF_5)
      call ESMF_GridCompFinalize(    gridComp=iceGridComp,
     &                            importState=iceImpState,
     &                            exportState=iceExpState,
     &                                  phase=1,
     &                               syncflag=ESMF_SYNC_NONBLOCKING,
     &                                     rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="SEAICE Finalize failed", rcToReturn=rc2))
#else
      call ESMF_GridCompFinalize(    gridComp=iceGridComp,
     &                            importState=iceImpState,
     &                            exportState=iceExpState,
     &                                  phase=ESMF_SINGLEPHASE,
     &                           blockingflag=ESMF_NONBLOCKING,
     &                                     rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "SEAICE Finalize failed", rcToReturn=rc2))
#endif
     &   goto 10
c
c --- Finalize OACPL coupler component
#if (USE_ESMF_5)
      call ESMF_CplCompFinalize(     cplComp=o2iCplComp,
     &                           importState=cplImpState,
     &                           exportState=cplExpState,
     &                                 phase=1,
     &                              syncflag=ESMF_SYNC_BLOCKING,
     &                                    rc=rc)
      if (ESMF_LogFoundError(rc,
     &    msg="OICPL  Finalize failed", rcToReturn=rc2))
#else
      call ESMF_CplCompFinalize(     cplComp=o2iCplComp,
     &                           importState=cplImpState,
     &                           exportState=cplExpState,
     &                                 phase=ESMF_SINGLEPHASE,
     &                          blockingflag=ESMF_BLOCKING,
     &                                    rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &    "OICPL  Finalize failed", rcToReturn=rc2))
#endif
     &   goto 10

c
10    continue
      if (localPet .eq. 0) write(6,'(a)') 'Exiting hycom_cice'
      call ESMF_VMBarrier(worldVM)
      call ESMF_Finalize(rc=rc)
c
      !KAL stop
      end program hycom_cice
