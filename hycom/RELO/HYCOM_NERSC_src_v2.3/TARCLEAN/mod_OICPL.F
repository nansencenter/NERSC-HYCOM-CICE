      module mod_OICPL  !ocean-ice coupler
c
c --- ESMF Framework module
#if (USE_ESMF_5)
      use ESMF
#else
      use ESMF_Mod
#endif
c
      implicit none
      private
c
      public OICPL_SetServices
c
c --- phase info
      integer, parameter, public :: ice2ocn_phase = 1
      integer, parameter, public :: ocn2ice_phase = 2
c
c --- VM and PET info
      type(ESMF_VM), save :: vm
      integer,       save :: petCount,localPet
c
c --- Route handles for regridding and/or redistribution
      type(ESMF_RouteHandle), save ::
     &    i2oRouteHandle, o2iRouteHandle

      contains

      subroutine OICPL_SetServices(cplComp, rc)
c
      type(ESMF_CplComp)   :: cplComp
      integer, intent(out) :: rc
c
#if (USE_ESMF_5)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_METHOD_INITIALIZE,
     &     OICPL_Init,
     &     phase=1,
     &     rc=rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_METHOD_RUN,
     &     OICPL_Run_I2O,
     &     phase=ice2ocn_phase,
     &     rc=rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_METHOD_RUN,
     &     OICPL_Run_O2I,
     &     phase=ocn2ice_phase,
     &     rc=rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_METHOD_FINALIZE,
     &     OICPL_Final,
     &     phase=1,
     &     rc=rc)
#else
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_SETINIT,
     &     OICPL_Init,
     &     ESMF_SINGLEPHASE,
     &     rc=rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_SETRUN,
     &     OICPL_Run_I2O,
     &     ice2ocn_phase,
     &     rc=rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_SETRUN,
     &     OICPL_Run_O2I,
     &     ocn2ice_phase,
     &     rc=rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_SETFINAL,
     &     OICPL_Final,
     &     ESMF_SINGLEPHASE,
     &     rc=rc)
#endif
c
      end subroutine OICPL_SetServices


      subroutine OICPL_Init(cplComp, impState, expState, extClock, rc)
c
c --- Calling parameters
      type(ESMF_CplComp)   :: cplComp
      type(ESMF_State)     :: impState
      type(ESMF_State)     :: expState
      type(ESMF_Clock)     :: extClock
      integer, intent(out) :: rc
c
c --- Locals
      integer :: rc2
      type(ESMF_State)       :: oiState, oeState, iiState, ieState
      type(ESMF_FieldBundle) :: ocnBundle,        iceBundle
c
c --- Report
#if (USE_ESMF_5)
      call ESMF_LogWrite("OICPL initialize routine called",
     &                   ESMF_LOGMSG_INFO, rc=rc)
#else
      call ESMF_LogWrite("OICPL initialize routine called",
     &                   ESMF_LOG_INFO, rc=rc)
#endif
      call ESMF_LogFlush(rc=rc)
c
c --- Get VM
      call ESMF_CplCompGet(cplComp, vm=vm, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get VM failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get VM failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get PET info
      call ESMF_VMGet(vm, petCount=petCount, localPET=localPet, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get VM info failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get VM info failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get OCEAN and SEAICE import states
      call ESMF_StateGet(impState, "OCEAN Import", oiState, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get OCEAN impState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get OCEAN impState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
      call ESMF_StateGet(impState, "SEAICE Import", iiState, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get SEAICE impState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get SEAICE impState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get OCEAN and SEAICE export states
      call ESMF_StateGet(expState, "OCEAN Export", oeState, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get OCEAN expState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get OCEAN expState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
      call ESMF_StateGet(expState, "SEAICE Export", ieState, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get SEAICE expState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get SEAICE expState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Initialize I2O
c
c --- Get bundle for ocn
      call ESMF_StateGet(oiState, "HYCOM Import", ocnBundle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get HYCOM Import failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get HYCOM Import failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get bundle for ice
      call ESMF_StateGet(ieState, "CICE Export",  iceBundle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get CICE Export failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get CICE Export failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Transfer fields from ice state to ocn state
c --- TODO: MAke redist ignore unmapped areas (possible with fieldrediststore)
      call ESMF_FieldBundleRedistStore(iceBundle, ocnBundle,
     &                                 i2oRouteHandle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, 
     &   msg="FieldBundleRedistStore i2o failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, 
     &   "FieldBundleRedistStore i2o failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c  Initialize O2I
c
c --- Get bundle for ice
      call ESMF_StateGet(iiState, "CICE Import", iceBundle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get CICE Import failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get CICE Import failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get bundle for ocn
      call ESMF_StateGet(oeState, "HYCOM Export", ocnBundle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get HYCOM Export failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get HYCOM Export failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Setup O2I route handle
      call ESMF_FieldBundleRedistStore(ocnBundle, iceBundle,
     &                                 o2iRouteHandle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc,
     &   msg="FieldBundleRedistStore o2i failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc,
     &   "FieldBundleRedistStore o2i failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
      return
      end subroutine OICPL_Init

      subroutine OICPL_Run_I2O(cplComp, impState, expState, extClock,
     &                         rc)
c
c --- Calling parameters
      type(ESMF_CplComp)   :: cplComp
      type(ESMF_State)     :: impState
      type(ESMF_State)     :: expState
      type(ESMF_Clock)     :: extClock
      integer, intent(out) :: rc
c
c --- Locals
      integer :: rc2
      type(ESMF_State)       :: oiState,   ieState
      type(ESMF_FieldBundle) :: ocnBundle, iceBundle
c
c --- Report
#if (USE_ESMF_5)
      call ESMF_LogWrite("OICPL I2O run routine called",
     &                   ESMF_LOGMSG_INFO, rc=rc)
#else
      call ESMF_LogWrite("OICPL I2O run routine called",
     &                   ESMF_LOG_INFO, rc=rc)
#endif
      call ESMF_LogFlush(rc=rc)
c
c --- Get OCEAN import state
      call ESMF_StateGet(impState, "OCEAN Import", oiState, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get OCEAN impState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get OCEAN impState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get SEAICE export state
      call ESMF_StateGet(expState, "SEAICE Export", ieState, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get SEAICE expState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get SEAICE expState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get bundle for ocn
      call ESMF_StateGet(oiState, "HYCOM Import", ocnBundle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get HYCOM Import failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get HYCOM Import failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get bundle for ice
      call ESMF_StateGet(ieState, "CICE Export",  iceBundle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get CICE Export failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get CICE Export failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Transfer fields from ice state to ocn state
      call ESMF_FieldBundleRedist(iceBundle, ocnBundle,
     &                            i2oRouteHandle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="FieldBundleRedist i2o failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "FieldBundleRedist i2o failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
      return
      end subroutine OICPL_Run_I2O

      subroutine OICPL_Run_O2I(cplComp, impState, expState, extClock,
     &                         rc)
c
c --- Calling parameters
      type(ESMF_CplComp)   :: cplComp
      type(ESMF_State)     :: impState
      type(ESMF_State)     :: expState
      type(ESMF_Clock)     :: extClock
      integer, intent(out) :: rc
c
c --- Locals
      integer :: rc2
      type(ESMF_State)       :: oeState,   iiState
      type(ESMF_FieldBundle) :: ocnBundle, iceBundle
c
c --- Report
#if (USE_ESMF_5)
      call ESMF_LogWrite( "OICPL O2I run routine called",
     &                   ESMF_LOGMSG_INFO, rc=rc)
#else
      call ESMF_LogWrite( "OICPL O2I run routine called",
     &                   ESMF_LOG_INFO, rc=rc)
#endif
      call ESMF_LogFlush(rc=rc)
c
c --- Get SEAICE import state
      call ESMF_StateGet(impState, "SEAICE Import", iiState, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get SEAICE impState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get SEAICE impState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get OCEAN export state
      call ESMF_StateGet(expState, "OCEAN Export", oeState, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get OCEAN expState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get OCEAN expState failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get bundle for ice
      call ESMF_StateGet(iiState, "CICE Import", iceBundle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get CICE Import failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get CICE Import failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Get bundle for ocn
      call ESMF_StateGet(oeState, "HYCOM Export", ocnBundle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="Get HYCOM Export failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "Get HYCOM Export failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
c --- Transfer fields from ocn state to ice state
      call ESMF_FieldBundleRedist(ocnBundle, iceBundle,
     &                            o2iRouteHandle, rc=rc)
#if (USE_ESMF_5)
      if (ESMF_LogFoundError(rc, msg="FieldBundleRedist o2i failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#else
      if (ESMF_LogMsgFoundError(rc, "FieldBundleRedist o2i failed",
     &   rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
#endif
c
      return
      end subroutine OICPL_Run_O2I

      subroutine OICPL_Final(cplComp, impState, expState, extClock, rc)
c
c --- Calling parameters
      type(ESMF_CplComp)   :: cplComp
      type(ESMF_State)     :: impState
      type(ESMF_State)     :: expState
      type(ESMF_Clock)     :: extClock
      integer, intent(out) :: rc
c
c --- Locals
c
c --- Report
#if (USE_ESMF_5)
      call ESMF_LogWrite("OICPL finalize routine called", 
     &                   ESMF_LOGMSG_INFO, rc=rc)
#else
      call ESMF_LogWrite("OICPL finalize routine called", 
     &                   ESMF_LOG_INFO, rc=rc)
#endif
      call ESMF_LogFlush(rc=rc)
c
c --- Release i2o regrid/redist route handle
      call ESMF_FieldBundleRedistRelease(i2oRouteHandle, rc=rc)
c
c --- Release o2i regrid/redist route handle
      call ESMF_FieldBundleRedistRelease(o2iRouteHandle, rc=rc)
c
      return
      end subroutine OICPL_Final

      end module mod_OICPL
