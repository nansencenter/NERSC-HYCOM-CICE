      module mod_esmf_utils
#if defined(NERSC_USE_ESMF)
      use ESMF
#else
      use ESMF_Mod
#endif
      implicit none
      public
      integer :: my_dummy_var
      contains

!!!#if (USE_ESMF) There is no end to this
      !Function for reducing clutter 
      logical function ESMF_LogFoundErrorWrapper(msg,rc)
      implicit none
      integer,         intent(in) :: rc
      character(len=*),intent(in) :: msg
      integer rc2
#if defined(NERSC_USE_ESMF)
      ESMF_LogFoundErrorWrapper=ESMF_LogFoundError(rc,msg="tst", &
         rcToReturn=rc2)
#else
      ESMF_LogFoundErrorWrapper=ESMF_LogMsgFoundError(rc,msg,    &
         rcToReturn=rc2)
#endif 
      end function  ESMF_LogFoundErrorWrapper


      subroutine ESMF_LogWrite_Info_Wrapper(msg)
      implicit none
      character(len=*),intent(in) :: msg
      integer :: rc
#if defined(NERSC_USE_ESMF)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
#else
      call ESMF_LogWrite(msg, ESMF_LOG_INFO, rc=rc)
#endif
      call ESMF_LogFlush(rc=rc)
      end subroutine
!!!!#endif

      subroutine ESMF_LogWrite_Error_Wrapper(msg)
      implicit none
      character(len=*),intent(in) :: msg
      integer :: rc
#if defined(NERSC_USE_ESMF)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
#else
      call ESMF_LogWrite(msg, ESMF_LOG_ERROR, rc=rc)
#endif
      call ESMF_LogFlush(rc=rc)
      end subroutine



      ! Check if two gridcomps match in size, coordinate and mask
      logical function check_gridcomp_grids(ocngridcomp,icegridcomp, &
         localPet, petCount)
      type(ESMF_GridComp) :: ocnGridComp,      & !HYCOM as an ESMF component
                             iceGridComp      ! CICE as an ESMF component
      integer :: petCount, localPet, split
      type(ESMF_Grid) :: ocnGrid
      type(ESMF_Grid) :: iceGrid
      integer,dimension(2) :: ocnLBound,ocnUBound
      integer,dimension(2) :: iceLBound,iceUBound
      integer :: ocnTileCount, iceTileCount
      integer :: ocnLocalDeCount, iceLocalDeCount
      integer :: ocnCoordDimCount(2), iceCoordDimCount(2)
      integer, target :: ocnMaxIndex(2), iceMaxIndex(2)
      integer, target :: ocnMinIndex(2), iceMinIndex(2)
      TYPE(ESMF_Array) :: iceXCoord, iceYCoord, iceMask
      TYPE(ESMF_Array) :: ocnXCoord, ocnYCoord, ocnMask
      real(KIND=ESMF_KIND_R4),  pointer :: &
         xc_ocn(:,:),yc_ocn(:,:)
      real(KIND=ESMF_KIND_R4),  pointer :: &
         xc_ice(:,:),yc_ice(:,:)
      integer,  pointer ::msk_ocn(:,:), msk_ice(:,:)
      integer :: ICE_nx, ICE_ny ! New
      integer :: OCN_nx, OCN_ny ! New
      type(ESMF_TypeKind_Flag) :: mytypekind

      integer :: i,j,rc
      character(len=256) :: msg
      real :: tmpr
      integer :: tmpi

!KAL  Guilty until proven innocent
      check_gridcomp_grids = .false.


!KAL-------------------------------------------------------------------
!KAL  Now  get dims from gridcomp
!KAL-------------------------------------------------------------------
!KAL  Get grids from gridcomp
      call ESMF_GridCompGet(ocnGridComp,grid=ocnGrid,rc=rc)
      if (ESMF_LogFoundErrorWrapper("hycom_cice: get ocnGrid",rc)) &
          return 
      call ESMF_GridCompGet(iceGridComp,grid=iceGrid,rc=rc)
      if (ESMF_LogFoundErrorWrapper("hycom_cice: get iceGrid",rc))  &
          return 
!KAL
!KAL  Get tile counts (Ocean)
      call ESMF_GridGet(ocnGrid, tileCount=ocnTileCount, &
         localDeCount=ocnLocalDeCount,rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ocn grid tilecount",rc))  &
          return 
!KAL  Get tile counts (Ice)
      call ESMF_GridGet(iceGrid, tileCount=iceTileCount, &
         localDeCount=iceLocalDeCount,rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ice grid tilecount",rc))  &
          return 
!KAL
!KAL  Make sure grids are single-tile. Otherwise bail out
      if (ocnTileCount .ne.1 .or. iceTileCount.ne.1) then
         write(msg,'("tileCounts of ice and ocn grids must be 1")')
         call ESMF_Logwrite_Error_wrapper(msg)
         return 
      end if
!KAL
!KAL  get Index ranges (Ocean)
      call ESMF_GridGet(ocnGrid,tile=1, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         minIndex=ocnminIndex, &
         maxIndex=ocnMaxIndex,rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ocn grid indexrange",rc))  &
          return 
!KAL  get Index ranges (Ocean)
      call ESMF_GridGet(iceGrid,tile=1, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         minIndex=iceminIndex, &
         maxIndex=iceMaxIndex,rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ice grid tindexrangeilecount",rc))  &
          return 

! --- Diag
!     if (localPet == 0 ) then
!        print *,"check_gridcomp_grids:ocnTileCount   ",ocnTileCount
!        print *,"check_gridcomp_grids:ocnLocalDeCount",ocnLocalDeCount
!        print *,"check_gridcomp_grids:ocnMinIndex    ",ocnMinIndex
!        print *,"check_gridcomp_grids:ocnMaxIndex    ",ocnMaxIndex
!        print *,"check_gridcomp_grids:iceTileCount   ",iceTileCount
!        print *,"check_gridcomp_grids:iceLocalDeCount",iceLocalDeCount
!        print *,"check_gridcomp_grids:iceMinIndex    ",iceMinIndex
!        print *,"check_gridcomp_grids:iceMaxIndex    ",iceMaxIndex
!     end if
      OCN_nx = ocnMaxIndex(1)
      OCN_ny = ocnMaxIndex(2)
      ICE_nx = iceMaxIndex(1)
      ICE_ny = iceMaxIndex(2)

! -- Check for grid size. For now they must match
! -- Grid decomposition may differ, however
      if     (OCN_nx.ne.ICE_nx) then
         call ESMF_LogWrite("grid size mismatch ", &
              ESMF_LOGMSG_ERROR, rc=rc)
         if     (localPet.eq.0) then !master
           write(6,'(a,i5)') 'grid size mismatch: OCN_nx = ',OCN_nx
           write(6,'(a,i5)') 'grid size mismatch: ICE_nx = ',ICE_nx
        end if
        return 
      end if
      if     (OCN_ny.ne.ICE_ny) then
         call ESMF_LogWrite("grid size mismatch ", &
              ESMF_LOGMSG_ERROR, rc=rc)
         if     (localPet.eq.0) then !master
           write(6,'(a,i5)') 'grid size mismatch: OCN_ny = ',OCN_ny
           write(6,'(a,i5)') 'grid size mismatch: ICE_ny = ',ICE_ny
        end if
        return 
      end if
      if     (localPet.eq.0) then !master
        write(6,'(a)') 'check_gridcomp_grids: grid sizes ok:'
        write(6,'(a,i5)') 'grid size : OCN_nx = ',OCN_nx
        write(6,'(a,i5)') 'grid size : ICE_nx = ',ICE_nx
        write(6,'(a,i5)') 'grid size : OCN_ny = ',OCN_ny
        write(6,'(a,i5)') 'grid size : ICE_ny = ',ICE_ny
      end if

! --- Now gather array and compare xcoord, ycoord and mask (OCEAN)
      call ESMF_GridGetCoord(ocnGrid,1, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         array=ocnXCoord, &
         rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ocn grid XCoord",rc))  &
          return 
! --- OCN Y Coord
      call ESMF_GridGetCoord(ocnGrid,2, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         array=ocnYCoord, &
         rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ocn grid YCoord",rc))  &
          return 
! --- OCN MAsk
      call ESMF_GridGetItem(ocnGrid, &
         ESMF_GRIDITEM_MASK, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         array=ocnMask, &
         rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ocn grid Mask",rc))  &
          return 
! --- Now gather array and compare xcoord, ycoord and mask (ICE)
      call ESMF_GridGetCoord(iceGrid,1, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         array=iceXCoord, &
         rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ice grid XCoord",rc))  &
          return 
! --- ICE Y Coord
      call ESMF_GridGetCoord(iceGrid,2, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         array=iceYCoord, &
         rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ice grid YCoord",rc))  &
          return 
! --- ICE mask
      call ESMF_GridGetItem(iceGrid, &
         ESMF_GRIDITEM_MASK, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         array=iceMask, &
         rc=rc)
      if (ESMF_LogFoundErrorWrapper( &
         "hycom_cice: get ice grid Mask",rc))  &
          return 

! --- Gather coords an masks in fortran arrays on al pets
      allocate(xc_ocn(OCN_nx,OCN_ny))
      allocate(yc_ocn(OCN_nx,OCN_ny))
      allocate(msk_ocn(OCN_nx,OCN_ny))
      allocate(xc_ice(OCN_nx,OCN_ny))
      allocate(yc_ice(OCN_nx,OCN_ny))
      allocate(msk_ice(OCN_nx,OCN_ny))
      do i=1,petCount

! ---    Gather OCN x coord
         call ESMF_ArrayGather(ocnXCoord,xc_ocn,rootPet=i-1,rc=rc)
         if (ESMF_LogFoundErrorWrapper( &
            "hycom_cice: gather ocn grid xc",rc))  &
             return 

! ---    Gather OCN y coord
         call ESMF_ArrayGather(ocnYCoord,yc_ocn,rootPet=i-1,rc=rc)
         if (ESMF_LogFoundErrorWrapper( &
            "hycom_cice: gather ocn grid yc",rc))  &
             return 

! ---    Gather OCN mask
!        call ESMF_ArrayGet(ocnMask,typekind=mytypekind,rc=rc)
         call ESMF_ArrayGather(ocnMask,msk_ocn,rootPet=i-1,rc=rc)
         if (ESMF_LogFoundErrorWrapper( &
            "hycom_cice: gather ocn grid mask",rc))  &
             return 

! ---    Gather ICE x coord
         call ESMF_ArrayGather(iceXCoord,xc_ice,rootPet=i-1,rc=rc)
         if (ESMF_LogFoundErrorWrapper( &
            "hycom_cice: gather ice grid xc",rc))  &
             return 

! ---    Gather ICE y coord
         call ESMF_ArrayGather(iceYCoord,yc_ice,rootPet=i-1,rc=rc)
         if (ESMF_LogFoundErrorWrapper( &
            "hycom_cice: gather ice grid yc",rc))  &
             return 

! ---    Gather ICE mask
!        call ESMF_ArrayGet(iceMask,typekind=mytypekind,rc=rc)
         call ESMF_ArrayGather(iceMask,msk_ice,rootPet=i-1,rc=rc)
         if (ESMF_LogFoundErrorWrapper( &
            "hycom_cice: gather ice grid mask",rc))  &
             return 
      end do

      ! Now compare masks, etc
!     print *,"localPet,xc_ocn,xc_ice:",localPet,
!    &   xc_ocn(100,100),xc_ice(100,100),xc_ice(99,99),xc_ice(101,101)
!     print *,"localPet,yc_ocn,yc_ice:",localPet,
!    &   yc_ocn(100,100),yc_ice(100,100)
!     print *,"localPet,msk_ocn,msk_ice:",localPet,
!    &   msk_ocn(1,1),msk_ice(1,1)
!     print *,"localPet,msk_ocn100,msk_ice:",localPet,
!    &   msk_ocn(100,100),msk_ice(100,100)


! --- Compare masks and coordinates
      do j=1,OCN_NY
        do i=1,OCN_NX
           tmpi=msk_ice(i,j)-msk_ocn(i,j)
! ---      Check that masks agree
           if (abs(tmpi) > 0) then
              if (localPet==0) write(6,'(a,2i5,a,3i10)')  &
                 "check_gridcomp_grids:ice and ocn msjk differ in ", &
                 i,j," ice,ocn,diff=",msk_ice(i,j),msk_ocn(i,j),tmpi
              return ! Slightly dangerous ...
           end if

! ---      Only compare coords where masks gree
           if (msk_ocn(i,j) > 0) then

              tmpr=xc_ocn(i,j)-xc_ice(i,j)
              if (abs(tmpr) > 1e-4) then
                 if (localPet==0) write(6,'(a,2i5,a,3f12.4)')  &
                    "check_gridcomp_grids:ice and ocn xc differ in ", &
                    i,j," ice,ocn,diff=",xc_ice(i,j),xc_ocn(i,j),tmpr
                 return ! Slightly dangerous ...
              end if

              tmpr=yc_ocn(i,j)-yc_ice(i,j)
              if (abs(tmpr) > 1e-4) then
                 if (localPet==0) write(6,'(a,2i5,a,3f12.4)')  &
                    "check_gridcomp_grids:ice and ocn yc differ in ", &
                    i,j," ice,ocn,diff=",yc_ice(i,j),yc_ocn(i,j),tmpr
                 return ! Slightly dangerous ...
              end if
           end if

        end do
      end do

      if (localPet==0) then
         write(6,'(a)') &
          "check_gridcomp_grids:All grid checks passed"
      end if

      check_gridcomp_grids = .true.

      end function check_gridcomp_grids


      logical function check_gridcomp_timesteps(ocngridcomp,icegridcomp, &
         ocnClock, iceClock, &
         localPet, petCount,OCN_nts_cpl,ICE_nts_cpl)
      implicit none
      type(ESMF_GridComp) :: ocnGridComp,      & !HYCOM as an ESMF component
                             iceGridComp      ! CICE as an ESMF component
      type(ESMF_Clock), intent(in)    :: ocnClock,      & !HYCOM as an ESMF component
                                         iceClock      ! CICE as an ESMF component
      integer :: petCount, localPet
      integer, intent(out) :: OCN_nts_cpl, ICE_nts_cpl
      
      type(ESMF_Grid) :: ocnGrid
      type(ESMF_Grid) :: iceGrid
      integer,dimension(2) :: ocnLBound,ocnUBound
      integer,dimension(2) :: iceLBound,iceUBound
      integer :: ocnTileCount, iceTileCount
      integer :: ocnLocalDeCount, iceLocalDeCount
      integer :: ocnCoordDimCount(2), iceCoordDimCount(2)
      integer, target :: ocnMaxIndex(2), iceMaxIndex(2)
      integer, target :: ocnMinIndex(2), iceMinIndex(2)
      TYPE(ESMF_Array) :: iceXCoord, iceYCoord, iceMask
      TYPE(ESMF_Array) :: ocnXCoord, ocnYCoord, ocnMask
      real(KIND=ESMF_KIND_R4),  pointer :: &
         xc_ocn(:,:),yc_ocn(:,:)
      real(KIND=ESMF_KIND_R4),  pointer :: &
         xc_ice(:,:),yc_ice(:,:)
      integer,  pointer ::msk_ocn(:,:), msk_ice(:,:)
      integer :: ICE_nx, ICE_ny ! New
      integer :: OCN_nx, OCN_ny ! New
      type(ESMF_TypeKind_Flag) :: mytypekind
!
! --- KAL - Moved from module acces to state access
      integer :: OCN_nts_day
      real    :: rOCN_nts_day, rOCN_nts_cpl, rICE_nts_day
      logical :: OCN_put_export, OCN_get_import, OCN_end_of_run
      integer :: ICE_nts_day, ICE_get_import, ICE_put_export
!KAL  integer :: ICE_nx, ICE_ny ! New
!KAL  integer :: OCN_nx, OCN_ny ! New
      real    :: iceTimeStep_sr8
      real    :: ocnTimeStep_sr8
      type(ESMF_TimeInterval) :: iceTimeStep
      type(ESMF_TimeInterval) :: ocnTimeStep
      character(len=256) :: msg
      integer :: rc,rc2
! --- ICE coupling frequency
      integer :: ocn_cpl_day
! --- Miscellaneous
      integer :: i,its,its_ocn,its_ice,icpl,iday

! --- Guilty until proven innocent...
      check_gridcomp_timesteps=.false.

! --- Get OCN_nts_cpl from ocean gridComp attribute
       call ESMF_AttributeGet(ocnGridComp,  &
          name="nts_ice",value=OCN_nts_cpl,rc=rc)
#if defined(NERSC_USE_ESMF)
      if (ESMF_LogFoundError(rc, &
         msg="hycom_cice: attribute get OCN_nts_cpl", rcToReturn=rc2)) &
#else
      if (ESMF_LogMsgFoundError(rc, &
         "hycom_cice: attributeget OCN_nts_cpl", rcToReturn=rc2)) &
#endif 
       return

! --- Here we use clock data to get info
      call ESMF_ClockGet(ocnClock, timeStep=ocnTimeStep ,rc=rc)
      if (ESMF_LogFoundErrorWrapper("hycom_cice: get ocnTS",rc))  &
          return 
      call ESMF_TimeIntervalGet(ocnTimeStep , s_r8=ocnTimeStep_sr8, &
                                rc=rc)
      if (ESMF_LogFoundErrorWrapper("hycom_cice: get ocnTS_sr8",rc))  &
          return 
      call ESMF_ClockGet(iceClock, timeStep=iceTimeStep ,rc=rc)
      if (ESMF_LogFoundErrorWrapper("hycom_cice: get iceTS",rc))  &
          return 
      call ESMF_TimeIntervalGet(iceTimeStep ,  &
         s_r8=iceTimeStep_sr8,rc=rc)
      if (ESMF_LogFoundErrorWrapper("hycom_cice: get iceTS_sr8",rc)) &
          return 
!      
      rocn_NTS_day = 86400./ocnTimeStep_sr8
      rice_NTS_day = 86400./iceTimeStep_sr8
      ocn_NTS_day  = 86400./ocnTimeStep_sr8
      ice_NTS_day  = 86400./iceTimeStep_sr8
      if (abs(ice_nts_day-rice_nts_day)>1e-6) then 
         write(msg,'("secprday not div. with ice time step",g14.6)')  &
            rice_nts_day
         call ESMF_Logwrite_Error_wrapper(msg)
         write(6,'(a)') msg
         return 
      end if
      if (abs(ocn_nts_day-rocn_nts_day)>1e-6) then 
         write(msg,'("secprday not div. with ocn time step",g14.6)')  &
            rocn_nts_day
         call ESMF_Logwrite_Error_wrapper(msg)
         write(6,'(a)') msg
         return 
      end if
!
      ocn_cpl_day = OCN_nts_day/OCN_nts_cpl
      ice_nts_cpl = ICE_nts_day/ocn_cpl_day
      if     (localPet.eq.0) then !master
        write(6,'(a,f10.2)') 'OCN_ts      = ',ocnTimeStep_sr8
        write(6,'(a,f10.2)') 'ICE_ts      = ',iceTimeStep_sr8
        write(6,'(a,i5)')    'OCN_nts_day = ',OCN_nts_day
        write(6,'(a,i5)')    'ICE_nts_day = ',ICE_nts_day
        write(6,'(a,i5)')    'OCN_nts_cpl = ',OCN_nts_cpl
        write(6,'(a,i5)')    'ICE_nts_cpl = ',ice_nts_cpl
      endif
      if     (OCN_nts_day.ne.ocn_cpl_day*OCN_nts_cpl) then
        if     (localPet.eq.0) then !master
          write(6,*) 'ERROR OCN_nts_cpl not a divisor of OCN_nts_day'
        endif
        return 
      endif
      if     (ICE_nts_day.ne.ocn_cpl_day*ice_nts_cpl) then
        if     (localPet.eq.0) then !master
          write(6,*) 'ERROR ice_nts_cpl not a divisor of ICE_nts_day'
        endif
        return 
      endif

      if (localPet==0) then
         write(6,'(a)')  &
            "check_gridcomp_timesteps:All clock checks passed"
      end if

      check_gridcomp_timesteps=.true.
      return
      end function check_gridcomp_timesteps

      end module
