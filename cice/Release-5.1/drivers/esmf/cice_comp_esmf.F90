module cice_comp_esmf
USE ESMF
use ice_kinds_mod,   only : int_kind, dbl_kind, char_len_long, log_kind
use ice_constants,   only : pi, c180, rad_to_deg, c0
use ice_blocks,      only : nx_block, ny_block
use ice_domain_size, only : max_blocks
IMPLICIT NONE
   !--Import Fields
   integer, parameter :: numImpFields=23 !KAL
   character(ESMF_MAXSTR), save :: impFieldName(    numImpFields),    &
                                   impFieldLongName(numImpFields), &
                                   impFieldStdName( numImpFields), &
                                   impFieldUnits(   numImpFields)
   real(ESMF_KIND_R4),     save :: impFieldSclFac(  numImpFields), &
                                   impFieldAddOff(  numImpFields)
!   integer,                save :: impFieldHalo(    numImpFields)
!
!--Export Fields
   integer, parameter :: numExpFields=11
   character(ESMF_MAXSTR), save :: expFieldName(    numExpFields), &
                                   expFieldLongName(numExpFields), &
                                   expFieldStdName( numExpFields), &
                                   expFieldUnits(   numExpFields)
   real(ESMF_KIND_R4),     save :: expFieldSclFac(  numExpFields), &
                                   expFieldAddOff(  numExpFields)
!   integer,                save :: expFieldHalo(    numExpFields)

   !--Data types for Import/Export array pointers
   type ArrayPtrReal2D
     real(ESMF_KIND_R4), dimension(:,:), pointer :: p
   end type ArrayPtrReal2D
!--ESMF related variables
   type(ESMF_FieldBundle), save :: expBundle, &
                                   impBundle
   type(ESMF_Field),       save :: expField(numExpFields), &
                                   impField(numImpFields)

   !KAL  type(ArrayPtrReal2D),   save :: expData( numExpFields), &
   !KAL                                  impData( numImpFields)
   !KAL use blocks as well - allocated below
   type(ArrayPtrReal2D), allocatable,   save :: expData(:,:), &
                                                impData(:,:)

   integer,                save :: petCount, localPet,  mpiCommunicator
   type(ESMF_DELayout),    save :: deLayout
   type(ESMF_VM),          save :: vm
   type(ESMF_Grid),        save :: grid2D
   type(ESMF_DistGrid),    save :: distgrid2D
   type(ESMF_ArraySpec),   save :: arraySpec2Dr

   ! Used to create deLayout
   integer, allocatable ::  deBlockList(:,:,:), petMap(:)
   integer, allocatable ::  deBlockList_land(:,:,:), petMap_land(:)

   ! Used to keep accumulated fluxes 
   real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), private :: &
      ave_sic          , &
      ave_uvel         , &
      ave_vvel         , & 
      accum_fresh_ai   , & ! fresh water flux to ocean (kg/m^2/s)
      accum_fsalt_ai   , & ! salt flux to ocean (kg/m^2/s)
      accum_fhocn_ai   , & ! net heat flux to ocean (W/m^2)
      accum_fswthru_ai     ! shortwave penetrating to ocean (W/m^2)
   real (kind=dbl_kind), private :: &
      accum_time


!=======================================================================

contains



!=======================================================================
subroutine cice_setservices(comp, rc)
    use ice_communicate, only : my_task
    implicit none
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc
    rc = ESMF_SUCCESS
    !if (my_task==0) print *, "In cice_setservices"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      ice_init_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      ice_run_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      ice_final_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine

  function ice_distgrid_esmf(gsize)
    use ice_blocks,      only : block, get_block, nx_block, ny_block, nblocks_tot
    use ice_domain,      only : blocks_ice, nblocks
    use ice_domain_size, only : nx_global, ny_global
    implicit none
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(out)    :: gsize
    !
    ! Return
    type(esmf_distgrid)     :: ice_distgrid_esmf
    !
    ! Local variables
    !
    integer,allocatable :: gindex(:)
    integer     :: lat
    integer     :: lon
    integer     :: i, j, iblk, n, gi
    integer     :: lsize
    integer     :: rc
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer     :: rc2
    !-------------------------------------------------------------------

    ! number the local grid

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
          enddo !i
       enddo    !j
    enddo        !iblk
    lsize = n

    ! not valid for padded decomps
    !    lsize = block_size_x*block_size_y*nblocks
    gsize = nx_global*ny_global

    allocate(gindex(lsize))
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             lon = this_block%i_glob(i)
             lat = this_block%j_glob(j)
             gi = (lat-1)*nx_global + lon
             gindex(n) = gi
          enddo !i
       enddo    !j
    enddo        !iblk
   
    print *,"nblocks=",nblocks
    print *,"nblocks_tot=",nblocks_tot
    !print *,"gindex=",gindex
    ice_distgrid_esmf = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rc, msg="cice_Setup_ESMF: ESMF_DistGridCreate", rcToReturn=rc2)) call ESMF_Finalize(rc=rc)

    deallocate(gindex)

  end function ice_DistGrid_esmf


!KAL - my implementation of ice_distgrid_esmf ++, that avoids arbSeqIndexList logic
  subroutine ice_grid_esmf()
    use ice_blocks,      only : block, get_block_parameter, nblocks_tot, get_block
    use ice_domain,      only : distrb_info, blocks_ice, nblocks
    use ice_distribution,only : ice_DistributionGet
    use ice_domain_size, only : nx_global, ny_global
    use ice_grid,        only : tlon,tlat,tmask
    use ice_communicate, only : my_task
    implicit none
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    ! Return
    !
    ! Local variables
    !
    integer(int_kind)     :: i, j,ig,jg
    integer(int_kind)  :: ilo, ihi, jlo, jhi                   ! beginning and end of physical domain
    integer(int_kind)  :: iblk, iblk_lnd, iblk_ocn, count_landblocks
    integer(int_kind), dimension(:), pointer :: i_glob, j_glob ! global domain location for each point
    integer(int_kind), dimension(:), pointer :: blockLocation
    type(block) :: this_block                                  ! block information
    integer(int_kind)     :: rc,rc2
    integer                    :: lbnd(2),ubnd(2)
    real(ESMF_KIND_R4),pointer :: Xcoord(:,:),Ycoord(:,:)
    integer,           pointer :: mask_ptr(:,:)
    character(ESMF_MAXSTR)     :: msg
    integer                    :: localDECount
    character(10)              :: dimNames(2),dimUnits(2)
    type(ESMF_Logical)         :: periodic(2)
    integer                    :: allocstat
    !-------------------------------------------------------------------

!   ! Create an ESMF grid that matches the HYCOM 2D grid
!   dimNames(1)="longitude";    dimNames(2)="latitude";
!   dimUnits(1)="degrees_east"; dimUnits(2)="degrees_north";
!   periodic(1)=ESMF_TRUE
!   periodic(2)=ESMF_FALSE

    !Get number of active blocks from distribution
    call ESMF_LogWrite("Setting up ESMF grid for CICE",ESMF_LOGMSG_INFO, rc=rc)
    call ice_DistributionGet(distrb_info, blockLocation=blockLocation)
    !print *,"task ",my_task," blocLocation: ",blockLocation
    !print *,"task ",my_task," blocLocalId:  ",distrb_info%blockLocalId
    !print *,"task ",my_task," blocGlobalId: ",distrb_info%blockGlobalId
    !print *,"task ",my_task," size(blocLocalId) : ",size(distrb_info%blockLocalId)
    !print *,"task ",my_task," size(blocGlobalId): ",size(distrb_info%blockGlobalId)
    !print *,"task ",my_task," size(blocLocation):",size(blockLocation)

    !Calculate deBlocklist, which gives block domain mapping to domain grid
    !Calculate petMap,      which gives block assignment to PETs
    allocate(deBlockList(2, 2, nblocks_tot),stat=allocstat)
    if (ESMF_LogFoundAllocError(allocstat, msg="deBlockList", rcToReturn=rc2))  call ESMF_Finalize(rc=rc)
    allocate(deBlockList_land(2, 2, nblocks_tot),stat=allocstat)
    if (ESMF_LogFoundAllocError(allocstat, msg="deBlockList_land", rcToReturn=rc2))  call ESMF_Finalize(rc=rc)
    allocate(petMap(nblocks_tot),stat=allocstat)
    if (ESMF_LogFoundAllocError(allocstat, msg="petMap", rcToReturn=rc2))  call ESMF_Finalize(rc=rc)
    allocate(petMap_land(nblocks_tot),stat=allocstat)
    if (ESMF_LogFoundAllocError(allocstat, msg="petMap_land", rcToReturn=rc2))  call ESMF_Finalize(rc=rc)
    iblk_ocn=0
    iblk_lnd =0
    do iblk = 1, nblocks_tot ! loob global block id
       call get_block_parameter(iblk,ilo=ilo,jlo=jlo,ihi=ihi,jhi=jhi, &
                                i_glob=i_glob,j_glob=j_glob)

       if (blockLocation(iblk)/=0) then

          iblk_ocn=iblk_ocn + 1

          ! Assign block location on distributed grid
          deBlockList(1,:,iblk_ocn) = (/i_glob(ilo),i_glob(ihi)/)
          deBlockList(2,:,iblk_ocn) = (/j_glob(jlo),j_glob(jhi)/)

          ! Assign block to persistent execution thread. 
          petMap(iblk_ocn) = blockLocation(iblk) - 1 !Starts from zero

          ! Diag
          !print '("my_task ",i4, ":glob ID =",i4, " ocn ID=",i4, " location=",i4)', &
          !   my_task, iblk, iblk_ocn,blockLocation(iblk)
          !print *,"iblk_ocn,deblocks i:",my_task,iblk_ocn,deBlockList(1,:,iblk_ocn)
          !print *,"iblk_ocn,deblocks j:",my_task,iblk_ocn,deBlockList(2,:,iblk_ocn)

       else 
          iblk_lnd=iblk_lnd + 1

          ! Assign block location on distributed grid (LAND POINTS)
          deBlockList_land(1,:,iblk_lnd) = (/i_glob(ilo),i_glob(ihi)/)
          deBlockList_land(2,:,iblk_lnd) = (/j_glob(jlo),j_glob(jhi)/)

          ! Assign block to persistent execution thread (LAND POINTS are placed on PET 0)
          ! NB: These points are not touched by CICE code, but looks like it must be present
          ! for redistrubiton
          petMap_land(iblk_lnd) = 0

          ! Diag
          !print '("my_task ",i4, ":glob ID =",i4, " lnd ID=",i4, " location=",i4)', &
          !   my_task, iblk, iblk_lnd,blockLocation(iblk)
          !print *,"iblk_lnd,deblocks i:",my_task,iblk_lnd,deBlockList_land(1,:,iblk_lnd)
          !print *,"iblk_lnd,deblocks j:",my_task,iblk_lnd,deBlockList_land(2,:,iblk_lnd)
       end if

    enddo        !iblk

    ! Append land blocks to end of deBlocklist. This way nothing touches them
    ! TODO: See if masking/redist ignore can be used to skip this
    count_landblocks=iblk_lnd
    deBlockList(:,:,iblk_ocn+1:nblocks_tot) = deBlockList_land(:,:,1:iblk_lnd)
    petMap(iblk_ocn+1:nblocks_tot) = petMap_land(1:iblk_lnd)


    call ESMF_VMBarrier(vm, rc=rc)
    print '("my_task ",i4, " deComp info: size petmap, my petmap count, nblocks =",i4,i4,i4)', &
       my_task,size(petmap), count(petmap.eq.my_task), nblocks




    ! The following for safety. Make sure nblocks matches count of tasks on this pet
    !TODO: Create clean exit
    !TODO: Assumption is that block ordering is increasing from 1
    if (my_task .eq. 0) then
       if (count(petmap.eq.my_task) .ne. nblocks+count_landblocks) then
          write(msg,'("ERROR my_task ",i4, ":my petmap and nblocks+count_landblocks dont match. count=",i4, " nblocks=",i4)') &
             my_task,count(petmap.eq.my_task),nblocks+count_landblocks
          print '(a)',msg
          call ESMF_LogWrite(msg,ESMF_LOGMSG_ERROR, rc=rc)
          !TODO: Create clean exit
          call ESMF_Finalize(rc=rc,endflag=ESMF_END_ABORT)
       end if
    elseif (count(petmap.eq.my_task) .ne. nblocks) then
       write(msg,'("ERROR my_task ",i4, ":my petmap and nblocks dont match. count=",i4, " nblocks=",i4)') &
          my_task,count(petmap.eq.my_task),nblocks
       print '(a)',msg
       call ESMF_LogWrite(msg,ESMF_LOGMSG_ERROR, rc=rc)
       !TODO: Create clean exit
       call ESMF_Finalize(rc=rc,endflag=ESMF_END_ABORT)
    end if


    ! Create ESMF mapping from de/block region to PETs
    call ESMF_LogWrite("Creating ESMF deLayout for CICE",ESMF_LOGMSG_INFO, rc=rc)
    deLayout = ESMF_DELayoutCreate(petMap=petMap, rc=rc)
    if (ESMF_LogFoundError(rc, msg="ice_grid_esmf: DELayoutCreate failed", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

    ! Create ESMF mapping from de/block to grid locations
    call ESMF_LogWrite("Creating ESMF distgrid for CICE",ESMF_LOGMSG_INFO, rc=rc)
    distgrid2D = ESMF_DistGridCreate(   &
        minIndex=(/1,1/),                          &
        maxIndex=(/nx_global,ny_global/),          &
        indexflag=ESMF_INDEX_GLOBAL,               &
        deBlockList=deBlockList,&
        deLayout=deLayout,&
        rc = rc)
    if (ESMF_LogFoundError(rc, msg="ice_grid_esmf: DistGridCreate failed", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

   !create the 2D grid
   call ESMF_LogWrite("Creating ESMF grid for CICE",ESMF_LOGMSG_INFO, rc=rc)
   grid2D=ESMF_GridCreate(distGrid=distGrid2D,            &
                          coordTypeKind=ESMF_TYPEKIND_R4, &
                          coordDimCount=(/2,2/),          &
                          indexflag=ESMF_INDEX_GLOBAL,    &
                          rc=rc)
   if (ESMF_LogFoundError(rc, msg="cice_comp_esmf:GridCreate failed",rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

   ! Add coordinates for cell centers
   call ESMF_GridAddCoord(grid=grid2D,                       &
                          staggerloc=ESMF_STAGGERLOC_CENTER, &
                          rc=rc)
   if (ESMF_LogFoundError(rc, msg="ice_grid_esmf:GridAddCoord failed",rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

   ! Assing coordinate values ( x-coordinate )
   call ESMF_LogWrite("Setting up ESMF grid x-coords for CICE",ESMF_LOGMSG_INFO, rc=rc)
   !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,ig,jg,this_block, &
   !$OMP                     lbnd,ubnd,Xcoord)
   do iblk=1,nblocks

      ! Set grid coordinates from CICE coordinates. 
      ! TODO: Not sure about the thread-safety of all this ...
      ! Alternative is pointers to pointer 
      call ESMF_GridGetCoord(grid=grid2D,                       &
                             CoordDim=1,                        &
                             localDe=iblk-1,                    &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             computationalLbound=lbnd,          &
                             computationalUbound=ubnd,          &
                             fArrayptr=Xcoord,                  &
                             rc=rc)
      if (ESMF_LogFoundError(rc,msg="ice_grid_esmf:GridGetCoord-1 failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      call check_block_bounds (iblk,"x-coord",lbnd,ubnd)
!      print '("my_task:",i4," xcoord lbnd=",2i4,"CICE i_glob(ilo),j_glob(jlo):",i4,i4)',my_task,lbnd, &
!         this_block%i_glob(ilo), this_block%j_glob(jlo)
!      print '("my_task:",i4," xcoord lbnd=",2i4,"CICE i_glob(ihi),j_glob(jhi):",i4,i4)',my_task,ubnd, &
!         this_block%i_glob(ihi), this_block%j_glob(jhi)
#if defined(ARCTIC)
!---- Arctic (tripole) domain, top row is replicated (ignore it)
      write(6,*) "cice_grid_esmf: ARCTIC tripole not supported yet"
      call ESMF_Finalize(rc=rc)
      !KAL - TODO - replace with cice equivalent
      !jja = min( jj, jtdm-1-j0 )
      ! jhi = min ( ...) ! (CICE equivalent, may be supported in block ! functions)
#endif
      do j= jlo,jhi
        jg  = this_block%j_glob(j)
        do i= ilo,ihi
          ig  = this_block%i_glob(i)
          Xcoord(ig,jg) = TLON(i,j,iblk) * rad_to_deg
        enddo
      enddo
   enddo
   !$OMP PARALLEL END DO
   call ESMF_VMBarrier(vm, rc=rc)
   call ESMF_LogFlush(rc=rc)

   ! Assing coordinate values ( y-coordinate )
   call ESMF_LogWrite("Setting up ESMF grid y-coords for CICE",ESMF_LOGMSG_INFO, rc=rc)
   !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,ig,jg,this_block, &
   !$OMP                     lbnd,ubnd,Ycoord)
   do iblk=1,nblocks
      call ESMF_GridGetCoord(grid=grid2D,                       &
                             CoordDim=2,                        &
                             localDe=iblk-1,                    &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             computationalLbound=lbnd,          &
                             computationalUbound=ubnd,          &
                             fArrayptr=Ycoord,                  &
                             rc=rc)
      if (ESMF_LogFoundError(rc,  msg="GridGetCoord-2 failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      call check_block_bounds (iblk,"y-coord",lbnd,ubnd)
!      print '("my_task:",i4," ycoord lbnd=",2i4,"CICE i_glob(ilo),j_glob(jlo):",i4,i4)',my_task,lbnd, &
!         this_block%i_glob(ilo), this_block%j_glob(jlo)
!      print '("my_task:",i4," ycoord lbnd=",2i4,"CICE i_glob(ihi),j_glob(jhi):",i4,i4)',my_task,ubnd, &
!         this_block%i_glob(ihi), this_block%j_glob(jhi)
      do j= jlo,jhi
        jg  = this_block%j_glob(j)
        do i= ilo,ihi
          ig  = this_block%i_glob(i)
          Ycoord(ig,jg) = TLAT(i,j,iblk) * rad_to_deg
        enddo
      enddo
   enddo
   !$OMP PARALLEL END DO
   call ESMF_VMBarrier(vm, rc=rc)

   ! Set up grid2D mask
   call ESMF_LogWrite("Setting up ESMF grid mask for CICE",ESMF_LOGMSG_INFO, rc=rc)
   CALL ESMF_GridAddItem(grid2D,                            &
                         ESMF_GRIDITEM_MASK,                &
                         staggerloc=ESMF_STAGGERLOC_CENTER, &
                         rc=rc)
   if (ESMF_LogFoundError(rc, msg="GridAddItem failed",rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

   ! Assing mask values
   mask_ptr(:,:) = 0  !all land, outside active tile
   !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,ig,jg,this_block, &
   !$OMP                     lbnd,ubnd,mask_ptr)
   do iblk=1,nblocks
      CALL ESMF_GridGetItem(grid2D,                             &
                            ESMF_GRIDITEM_MASK,                 &
                            localDE=iblk-1,                     &
                            staggerloc=ESMF_STAGGERLOC_CENTER,  &
                            computationalLbound=lbnd,          &
                            computationalUbound=ubnd,          &
                            fArrayptr=mask_ptr,                 &
                            rc=rc)
      if (ESMF_LogFoundError(rc, msg="GridGetItem failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
       
      ! Set grid mask from CICE mask
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      call check_block_bounds (iblk,"mask",lbnd,ubnd)
!      print '("my_task:",i4," mask lbnd=",2i4,"CICE i_glob(ilo),j_glob(jlo):",i4,i4)',my_task,lbnd, &
!         this_block%i_glob(ilo), this_block%j_glob(jlo)
!      print '("my_task:",i4," mask lbnd=",2i4,"CICE i_glob(ihi),j_glob(jhi):",i4,i4)',my_task,ubnd, &
!         this_block%i_glob(ihi), this_block%j_glob(jhi)


      do j= jlo,jhi
        jg  = this_block%j_glob(j)
        do i= ilo,ihi
          ig  = this_block%i_glob(i)
          if (TMASK(i,j,iblk)) then 
             mask_ptr(ig,jg) = 1
          else
             mask_ptr(ig,jg) = 0
          end if
        enddo
      enddo
   enddo
   !$OMP PARALLEL END DO
   call ESMF_VMBarrier(vm, rc=rc)
   call ESMF_LogFlush(rc=rc)


   ! Set up grid2D mask on land-only blocks. These are only defined in ESMF
   ! Array, and only set on first Pet
   if (localPet==0) then
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,ig,jg,this_block, &
      !$OMP                     lbnd,ubnd,mask_ptr)
      do iblk=nblocks+1,nblocks+count_landblocks
         CALL ESMF_GridGetItem(grid2D,                             &
                               ESMF_GRIDITEM_MASK,                 &
                               localDE=iblk-1,                     &
                               staggerloc=ESMF_STAGGERLOC_CENTER,  &
                               computationalLbound=lbnd,          &
                               computationalUbound=ubnd,          &
                               fArrayptr=mask_ptr,                 &
                               rc=rc)
         if (ESMF_LogFoundError(rc, msg="GridGetItem failed",rcToReturn=rc2)) &
            call ESMF_Finalize(rc=rc)

          
        ! Set grid mask from CICE mask
         mask_ptr(:,:)=0

         ! The rest (xc,yc) is not set
      enddo

      !$OMP PARALLEL END DO
   end if
   call ESMF_VMBarrier(vm, rc=rc)
   call ESMF_LogFlush(rc=rc)

  end subroutine ice_Grid_esmf



subroutine cice_setup_esmf(gridComp, impState, expState, extClock, rc)
use ice_calendar   , only : dt, days_per_year, use_leap_years, daymo, &
   nyr, month, mday, hour, sec, npt, calendar, time, year_init, time2sec, idate, &
   basis_seconds
use ice_domain_size, only : nx_global, ny_global
use ice_blocks,      only : block, get_block, nx_block, ny_block
use ice_domain,      only : blocks_ice, nblocks
use ice_communicate, only : my_task, master_task
use ice_restart_shared, only : runtype, use_restart_time
use ice_therm_shared,   only : ktherm
use ice_ocean, only: tfrz_option
use CICE_InitMod
implicit none
!KAL : based on hycom setup_esmf
!
!--Calling parameters
   type(ESMF_GridComp)  :: gridComp
   type(ESMF_State)     :: impState
   type(ESMF_State)     :: expState
   type(ESMF_Clock)     :: extClock
   integer, intent(out) :: rc
!
!--set up ESMF data structures for HYCOM.
!
   integer                    :: i,j
!  integer(ESMF_KIND_I4)      :: year,month,day,hour,minute
!  integer(ESMF_KIND_I4)      :: sec,msec,usec,nsec
!  real(8)                    :: dsec,dmsec,dusec,dnsec
!  type(ESMF_TimeInterval)    :: timeStep, runDuration
!  type(ESMF_Time)            :: startTime
!  character(ESMF_MAXSTR)     :: msg
   integer :: nts_day,rc2, gsize, my_sec
   integer :: iblk,ilo,jlo,ihi,jhi,ig,jg
   integer :: localDeCount, deCount
   real(kind=dbl_kind) :: fnts_day
   character(ESMF_MAXSTR)     :: msg
   type(ESMF_Calendar)     :: e_calendar
   type(ESMF_TimeInterval)     :: timeStep, runDuration
   type(ESMF_Time)     :: startTime, stopTime, refTime, tmpTime
   integer :: start_year, start_month, start_mday, start_sec
   integer :: stop_year, stop_month, stop_mday, stop_sec
   integer :: tmp_year, tmp_month, tmp_mday, tmp_sec
   real*8  :: tmp1,tmp2
   integer :: allocstat 

   !--Report
   call ESMF_LogWrite("CICE ESMF Setup routine called",ESMF_LOGMSG_INFO, rc=rc)
   call ESMF_LogFlush(rc=rc)


   !---------------------------------------------------------------------
   !--Set up attributes for import and export fields
   !---------------------------------------------------------------------
   expFieldAddOff(:) = 0.0     !default is no offset
   expFieldSclFac(:) = 1.0     !default is no scale factor
!   expFieldHalo(  :) = halo_ps !default is scalar p-grid
   !
   expFieldName(     1) = "sic"
   expFieldLongName( 1) = "Sea Ice Concentration"
   expFieldStdName(  1) = "sea_ice_area_fraction"
   expFieldUnits(    1) = "1"
   expFieldName(     2) = "sitx"
   expFieldLongName( 2) = "Sea Ice X-Stress"
   expFieldStdName(  2) = "downward_x_stress_at_sea_ice_base"
   expFieldSclFac(   2) = -1.0  !field is upward
   expFieldUnits(    2) = "Pa"
!   expFieldHalo(     2) = halo_pv !vector p-grid
   expFieldName(     3) = "sity"
   expFieldLongName( 3) = "Sea Ice Y-Stress"
   expFieldStdName(  3) = "downward_y_stress_at_sea_ice_base"
   expFieldSclFac(   3) = -1.0  !field is upward
   expFieldUnits(    3) = "Pa"
!   expFieldHalo(     3) = halo_pv !vector p-grid
   expFieldName(     4) = "siqs"
   expFieldLongName( 4) = "Solar Heat Flux thru Ice to Ocean"
   expFieldStdName(  4) = "downward_sea_ice_basal_solar_heat_flux"
   expFieldUnits(    4) = "W m-2"
   expFieldName(     5) = "sifh"
   expFieldLongName( 5) = "Ice Freezing/Melting Heat Flux"
   expFieldStdName(  5) = "upward_sea_ice_basal_heat_flux"
   expFieldSclFac(   5) = -1.0  !field is downward
   expFieldUnits(    5) = "W m-2"
   expFieldName(     6) = "sifs"
   expFieldLongName( 6) = "Ice Freezing/Melting Salt Flux"
   expFieldStdName(  6) = "downward_sea_ice_basal_salt_flux"
   expFieldUnits(    6) = "kg m-2 s-1"
   expFieldName(     7) = "sifw"
   expFieldLongName( 7) = "Ice Net Water Flux"
   expFieldStdName(  7) = "downward_sea_ice_basal_water_flux"
   expFieldUnits(    7) = "kg m-2 s-1"
   expFieldName(     8) = "sit" !diagnostic
   expFieldLongName( 8) = "Sea Ice Temperature"
   expFieldStdName(  8) = "sea_ice_temperature"
   expFieldAddOff(   8) = +273.15 !field is in degC
   expFieldUnits(    8) = "K"
   expFieldName(     9) = "sih" !diagnostic
   expFieldLongName( 9) = "Sea Ice Thickness"
   expFieldStdName(  9) = "sea_ice_thickness"
   expFieldUnits(    9) = "m"
   expFieldName(    10) = "siu" !diagnostic
   expFieldLongName(10) = "Sea Ice X-Velocity"
   expFieldStdName( 10) = "sea_ice_x_velocity"
   expFieldUnits(   10) = "m s-1"
!   expFieldHalo(    10) = halo_pv !vector p-grid
   expFieldName(    11) = "siv" !diagnostic
   expFieldLongName(11) = "Sea Ice Y-Velocity"
   expFieldStdName( 11) = "sea_ice_y_velocity"
   expFieldUnits(   11) = "m s-1"
!   expFieldHalo(    11) = halo_pv !vector p-grid
!
!  expFieldName(    12) = "patm"
!  expFieldLongName(12) = "Surface Air Pressure"
!  expFieldStdName( 12) = "surface_air_pressure"
!  expFieldUnits(   12) = "Pa"
!  expFieldName(    13) = "xwnd"
!  expFieldLongName(13) = "X-Wind"
!  expFieldStdName( 13) = "x_wind"
!  expFieldUnits(   13) = "m s-1"
!  expFieldHalo(    13) = halo_pv !vector p-grid
!  expFieldName(    14) = "ywnd"
!  expFieldLongName(14) = "Y-Wind"
!  expFieldStdName( 14) = "y_wind"
!  expFieldUnits(   14) = "m s-1"
!  expFieldHalo(    14) = halo_pv !vector p-grid
!
   !  Attributes for export fields, identical to CICE import fields
   impFieldAddOff(:) = 0.0 !default is no offset
   impFieldSclFac(:) = 1.0 !default is no scale factor
!   impFieldHalo(  :) = halo_ps !default is scalar p-grid
   !
   impFieldName(     1) = "sst"
   impFieldLongName( 1) = "Sea Surface Temperature"
   impFieldStdName(  1) = "sea_surface_temperature"
   impFieldAddOff(   1) = +273.15 !field is in degC
   impFieldUnits(    1) = "K"
   impFieldName(     2) = "sss"
   impFieldLongName( 2) = "Sea Surface Salinity"
   impFieldStdName(  2) = "sea_surface_salinity"
   impFieldUnits(    2) = "1e-3"
   impFieldName(     3) = "ssu"
   impFieldLongName( 3) = "Sea Surface X-Current"
   impFieldStdName(  3) = "sea_water_x_velocity"
   impFieldUnits(    3) = "m s-1"
!   impFieldHalo(     3) = halo_pv !vector p-grid
   impFieldName(     4) = "ssv"
   impFieldLongName( 4) = "Sea Surface Y-Current"
   impFieldStdName(  4) = "sea_water_y_velocity"
   impFieldUnits(    4) = "m s-1"
!   impFieldHalo(     4) = halo_pv !vector p-grid
   impFieldName(     5) = "ssh"
   impFieldLongName( 5) = "Sea Surface Height"
   impFieldStdName(  5) = "sea_surface_height_above_sea_level"
   impFieldUnits(    5) = "m"
   impFieldName(     6) = "ssfi"
   impFieldLongName( 6) = "Oceanic Heat Flux Available to Sea Ice"
   impFieldStdName(  6) = "upward_sea_ice_basal_available_heat_flux"
   impFieldSclFac(   6) = -1.0  !field is downward
   impFieldUnits(    6) = "W m-2"
   impFieldName(     7) = "mlt"  !diagnostic
   impFieldLongName( 7) = "Ocean Mixed Layer Thickness"
   impFieldStdName(  7) = "ocean_mixed_layer_thickness"
   impFieldUnits(    7) = "m"
!KAL - new. Imported from hycom -> cice
   impFieldName(     8) = "tair"  !diagnostic
   impFieldLongName( 8) = ""
   impFieldStdName(  8) = ""
   impFieldUnits(    8) = "K"
   impFieldName(     9) = "uair"  !diagnostic
   impFieldLongName( 9) = ""
   impFieldStdName(  9) = ""
   impFieldUnits(    9) = "m s-1"
   impFieldName(    10) = "vair"  !diagnostic
   impFieldLongName(10) = ""
   impFieldStdName( 10) = ""
   impFieldUnits(   10) = "m s-1"
   impFieldName(    11) = "zair"  !diagnostic
   impFieldLongName(11) = ""
   impFieldStdName( 11) = ""
   impFieldUnits(   11) = "m"
   impFieldName(    12) = "pott"  !diagnostic
   impFieldLongName(12) = ""
   impFieldStdName( 12) = ""
   impFieldUnits(   12) = "K"
   impFieldName(    13) = "ss_tltx"
   impFieldLongName(13) = "Sea Surface slope x"
   impFieldStdName( 13) = ""
   impFieldUnits(   13) = "[]"
   impFieldName(    14) = "ss_tlty"
   impFieldLongName(14) = "Sea Surface Slope y"
   impFieldStdName( 14) = ""
   impFieldUnits(   14) = "[]"
   impFieldName(    15) = "Qair"
   impFieldLongName(15) = "Specific humidity of air"
   impFieldStdName( 15) = ""
   impFieldUnits(   15) = "kg kg**-1"
   impFieldName(    16) = "rhoair"
   impFieldLongName(16) = "Specific humidity of air"
   impFieldStdName( 16) = ""
   impFieldUnits(   16) = "kg m**-3"
   impFieldName(    17) = "swvdr"
   impFieldLongName(17) = "sw down, visible, direct"
   impFieldStdName( 17) = ""
   impFieldUnits(   17) = "W m**-3"
   impFieldName(    18) = "swvdf"
   impFieldLongName(18) = "sw down, visible, diffuse"
   impFieldStdName( 18) = ""
   impFieldUnits(   18) = "W m**-3"
   impFieldName(    19) = "swidr"
   impFieldLongName(19) = "sw down, infrared, direct"
   impFieldStdName( 19) = ""
   impFieldUnits(   19) = "W m**-3"
   impFieldName(    20) = "swidf"
   impFieldLongName(20) = "sw down, infrared, diffuse"
   impFieldStdName( 20) = ""
   impFieldUnits(   20) = "W m**-2"
   impFieldName(    21) = "flw"
   impFieldLongName(21) = "Longwave radiation, downward"
   impFieldStdName( 21) = ""
   impFieldUnits(   21) = "W m**-3"
   impFieldName(    22) = "rrate"
   impFieldLongName(22) = "Rainfall Rate"
   impFieldStdName( 22) = ""
   impFieldUnits(   22) = "kg m**-2 s*-1"
   impFieldName(    23) = "srate"
   impFieldLongName(23) = "Snowfall Rate"
   impFieldStdName( 23) = ""
   impFieldUnits(   23) = "kg m**-2 s*-1"
!KAL - new


   !---------------------------------------------------------------------
   !--Initial operations on VM
   !---------------------------------------------------------------------

   ! Get VM from gridComp
   call ESMF_GridCompGet(gridComp, vm=vm, rc=rc)
   if (ESMF_LogFoundError(rc, msg="CICE Init: GridCompGet failed", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

   ! Get VM info (local pe and pe count etc
   call ESMF_VMGet(vm, petCount=petCount, localPET=localPet,  &
        mpiCommunicator=mpiCommunicator, rc=rc)
   if (ESMF_LogFoundError(rc,msg="CICE Init: VMGet failed", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)
   write(msg,'(a,i4)') "CICE Init: petCount = ",petCount
   call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
   call ESMF_LogFlush(rc=rc)

   !---------------------------------------------------------------------
   !-- Set up ESMF clock            
   !---------------------------------------------------------------------
   ! Get input clock properties 
   call ESMF_ClockGet(extClock,startTime=startTime,stopTime=stopTime,rc=rc)
   if (ESMF_LogFoundError(rc, msg="cice_setup_Esmf: unable to get extClock properties", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

   call ESMF_TimeGet(startTime,timestring=msg)
   write(msg,'("CICE startTime from extClock:",a)') trim(msg)
   call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
   if (my_task==master_task) print '(a)',msg

   call ESMF_TimeGet(stopTime,timestring=msg)
   write(msg,'("CICE stopTime from extClock :",a)') trim(msg)
   call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
   if (my_task==master_task) print '(a)',msg


   ! Set calendar
   !if (use_leap_years) then
   !   calendar = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,name="Gregorian", rc=rc)
   !else 
   !   if (days_per_year==360) then
   !      calendar = ESMF_CalendarCreate(ESMF_CALKIND_360DAY,name="360DAY", rc=rc)
   !   elseif (days_per_year==365) then
   !      calendar = ESMF_CalendarCreate(name="365DAYS", daysPerMonth=daymo, rc=rc)
   !   elseif (days_per_year==366) then
   !      calendar = ESMF_CalendarCreate(name="366DAYS", daysPerMonth=daymo, rc=rc)
   !   else 
   !      write(msg,'("Unknown days_per_year:",i4)') days_per_year
   !      call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
   !      call ESMF_Finalize(rc=rc)
   !   end if
   !end if
   ! KAL - our restriction
   ! TODO: This test does not work since days_per_year is set in set_calendar,
   ! TODO: and overrides namelist value. For now we only test for leap year
   ! if (use_leap_years .and. days_per_year==365) then
   if (use_leap_years) then
      e_calendar = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,name="Gregorian", rc=rc)
   else 
      write(msg,'("Must use leap years and days_per_year=365")')
      call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
      call ESMF_Finalize(rc=rc)
   end if

   call ESMF_TimeGet(startTime,yy=start_year,mm=start_month, dd=start_mday,s=start_sec,rc=rc)
   if (ESMF_LogFoundError(rc,                               &
      msg="cice_setup_Esmf: Error getting time from extClock startTime", &
      rcToReturn=rc2))                                      & 
      call ESMF_Finalize(rc=rc)

   call ESMF_TimeGet(stopTime,yy=stop_year,mm=stop_month, dd=stop_mday,s=stop_sec,rc=rc)
   if (ESMF_LogFoundError(rc,                               &
      msg="cice_setup_Esmf: Error getting time from extClock startTime", &
      rcToReturn=rc2))                                      & 
      call ESMF_Finalize(rc=rc)

   ! Set timeStep (thermo)
   call ESMF_TimeIntervalSet(timeStep,s_r8=dt)
   if (ESMF_LogFoundError(rc, msg="cice_setup_Esmf: unable to set ESMF timeStep", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

   ! Set CICE internal time state to match that of extClock for initialization
   if (localPet==0) then 
      print '(a,5i6)',"cice_setup_esmf: CICE time[ice_in]:",nyr,year_init+nyr-1,month,mday,sec
      print '(a,5i6)',"cice_setup_esmf: Ext time         :",start_year,start_month,start_mday,start_sec
      !print '(a,5i6)',"cice_setup_esmf: Time             :",start_year,start_month,start_mday,start_sec
      !print *,"time ",time
   end if
   if (runtype=='initial') then 
      call time2sec(start_year, start_month, start_mday,time)
      time = time + start_sec
      time = time - basis_seconds
      call calendar(time)
   !KALelseif (runtype=='restart') then 
   elseif (runtype=='continue') then 
      if (use_restart_time) then
         !Need to check that clock and restart time info actually match.
         !Note that at this time, the model has jumped one time step in the init routine ..
         tmpTime = startTime + timeStep
         call ESMF_TimeGet(tmpTime,yy=tmp_year,mm=tmp_month, dd=tmp_mday,s=tmp_sec,rc=rc)
         if (localPet == 0 ) &
            print '(a,5i6)',"cice_setup_esmf: Ext  time + dt   :",tmp_year,tmp_month,tmp_mday,tmp_sec
         if (tmp_year <> year_init+nyr-1 .or.  tmp_month <> month .or. tmp_mday <> mday .or.  sec <> tmp_sec) then
            tmp1=start_year*10000+start_mday*100+start_mday + start_sec/86400.
            tmp2=(nyr+year_init-1)*10000+month*100+mday + sec/86400.
            write(msg,"('cice_setup_esmf: Time from restart do not match external clock',2f14.4)") tmp1,tmp2
            if (my_task==master_task) print *,msg
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
            call ESMF_Finalize(rc=rc)
         else 
            !KAL - yes, reset back to startTime
            call time2sec(start_year, start_month, start_mday,time)
            time = time + start_sec
            time = time - basis_seconds
            call calendar(time)
         end if
      else 
         if (localPet == 0 ) &
            print '(a,5i6)',"cice_setup_esmf: Using cice time  :",nyr,year_init+nyr-1,month,mday,sec
         !call time2sec(start_year, start_month, start_mday,time)
         !time = time + start_sec
         !time = time - basis_seconds
         !call calendar(time)
      end if
      write(msg,"('Not checking init clocks on restartfor runtype =',a)") runtype
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
   else 
      write(msg,"('Dont know how to handle runtype =',a)") runtype
      if (my_task==master_task) print *,msg
      call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
      call ESMF_Finalize(rc=rc)
   end if

   ! Apply settings to clock
   call ESMF_ClockSet(extClock,name="CICE Clock", timeStep=timeStep, rc=rc)
   if (ESMF_LogFoundError(rc, msg="hycom_init: unable to set extClock", rcToReturn=rc2)) &
     call ESMF_Finalize(rc=rc)

   ! Diag
   call ESMF_TimeGet(startTime,timestring=msg)
   write(msg,'("CICE start  Time:",a)') trim(msg)
   call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
   if (my_task==master_task) print '(a)',msg
   call ESMF_TimeGet(stopTime,timestring=msg)
   write(msg,'("CICE stop   Time:",a)') trim(msg)
   call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
   if (my_task==master_task) print '(a)',msg
   call ESMF_TimeIntervalGet(timeStep,timestring=msg)
   write(msg,'("CICE Time Step  :",a)') trim(msg)
   call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
   if (my_task==master_task) print '(a)',msg


   !---------------------------------------------------------------------
   !-- Set up grid using CICE domain info
   !---------------------------------------------------------------------

   ! Create array specifications
   call ESMF_ArraySpecSet(arraySpec2Dr, rank=2, typekind=ESMF_TYPEKIND_R4, rc=rc)
   if (ESMF_LogFoundError(rc, msg="Setup_ESMF: ArraySpecSet failed", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)

   ! Create distributed grid
   call ice_grid_esmf()

   ! Associate grid with ESMF gridded component
   call ESMF_GridCompSet(gridComp, grid=grid2D, rc=rc)
   if (ESMF_LogFoundError(rc, msg="cice_setup_Esmf: GridCompSet", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)
   call ESMF_LogFlush(rc=rc)


   !---------------------------------------------------------------------
   !-- Get local DE count from Grid
   !---------------------------------------------------------------------
   call ESMF_GridGet(grid2D,localDeCount=localDeCount,rc=rc) 
   if (ESMF_LogFoundError(rc, msg="cice_setup_Esmf: Gridget localDeCount export", rcToReturn=rc2)) &
      call ESMF_Finalize(rc=rc)


   !---------------------------------------------------------------------
   !-- Set up export fields, initialize variables
   !---------------------------------------------------------------------
!  !Get local DE count
!  call ESMF_FieldGet(expField(i),localDeCount=localDeCount,rc=rc) 
!  if (ESMF_LogFoundError(rc, msg="cice_setup_Esmf: Fieldget localDeCount export", rcToReturn=rc2)) &
!     call ESMF_Finalize(rc=rc)

   ! Setup export fields, bundles & state
   allocate(expData(numExpFields,localDeCount),stat=allocstat) ! Allocate using localDeCount, not nblocks
   if (ESMF_LogFoundAllocError(allocstat, msg="expData", rcToReturn=rc2)) &
        call ESMF_Finalize(rc=rc)
   do i=1,numExpFields
     expField(i)=ESMF_FieldCreate(grid=grid2D,                       &
                                  arrayspec=arraySpec2Dr,            &
                                  indexflag=ESMF_INDEX_GLOBAL,       &
                                  staggerLoc=ESMF_STAGGERLOC_CENTER, &
                                  name=trim(expFieldName(i)),        &
                                  rc=rc)

     !TODO: OpenMP?
     ! Initialize here, using ESMF de info. In actual transfers we will use CICE
     ! block info
     do deCount = 0, localDeCount-1
        call ESMF_FieldGet(expField(i),localDe=deCount,farrayptr=expData(i,deCount+1)%p,rc=rc)
        if (ESMF_LogFoundError(rc, msg="cice_setup_Esmf: Fieldget export", rcToReturn=rc2)) &
           call ESMF_Finalize(rc=rc)
        expData(i,deCount+1)%p(:,:) = 0.0
        !print '("my_task ",i4, "deCount=",i4, "localDeCount=",i4, ": init exp fieldname=",a)', &
        !  my_task,deCount, localDeCount,trim(expFieldName(i))
     end do

   enddo

   ! Create bundle from list of export fields
   expBundle=ESMF_FieldBundleCreate(                    &
                    fieldList=expField(1:numExpFields), &
                    name='CICE Export',                 &
                    rc=rc)
   if (ESMF_LogFoundError(rc,                               &
      msg="cice_setup_Esmf: FieldBundleCreate CICE Export", &
      rcToReturn=rc2))                                      & 
      call ESMF_Finalize(rc=rc)

   ! Add bundle to the export state
   call ESMF_StateAdd(expState,fieldBundleList=(/expBundle/),rc=rc)
   if (ESMF_LogFoundError(rc,                                 &
      msg="cice_setup_Esmf: StateAdd bundle to export state", &
      rcToReturn=rc2))                                        & 
      call ESMF_Finalize(rc=rc)
   call ESMF_LogFlush(rc=rc)


   !---------------------------------------------------------------------
   !-- Set up import fields, initialize variables
   !---------------------------------------------------------------------

   !Setup import fields, bundles & state
   allocate(impData(numImpFields,localDeCount),stat=allocstat) ! Allocate using localDeCount, not nblocks
   if (ESMF_LogFoundAllocError(allocstat, msg="impData", rcToReturn=rc2)) &
        call ESMF_Finalize(rc=rc)
   do i = 1,numImpFields
     impField(i)=ESMF_FieldCreate(grid2D,arraySpec2Dr,               &
                                  StaggerLoc=ESMF_STAGGERLOC_CENTER, &
                                  name=trim(impFieldName(i)),        &
                                  rc=rc) 

!    !Get local DE count
!    call ESMF_FieldGet(impField(i),localDeCount=localDeCount,rc=rc) 
!    if (ESMF_LogFoundError(rc, msg="cice_setup_Esmf: Fieldget localDeCount import", rcToReturn=rc2)) &
!       call ESMF_Finalize(rc=rc)

     !TODO: OpenMP?
     ! Initialize here, using ESMF de info. In actual transfers we will use CICE
     ! block info
     do deCount = 0, localDeCount-1
        call ESMF_FieldGet(impField(i),localDe=deCount,farrayptr=impData(i,deCount+1)%p,rc=rc)
        if (ESMF_LogFoundError(rc, msg="cice_setup_Esmf: Fieldget import", rcToReturn=rc2)) &
           call ESMF_Finalize(rc=rc)
        impData(i,deCount+1)%p(:,:) = 0.0
!        print '("my_task ",i4, "deCount=",i4, "localDeCount=",i4, ": init imp fieldname=",a)', &
!             my_task,deCount, localDeCount,trim(impFieldName(i))
     end do
   enddo
   call ESMF_LogFlush(rc=rc)

   ! Create bundle from list of fields
   impBundle=ESMF_FieldBundleCreate(                      &
                      fieldList=impField(1:numImpFields), &
                      name='CICE Import',                 &
                      rc=rc)
   if (ESMF_LogFoundError(rc,                               &
      msg="cice_setup_Esmf: FieldBundleCreate CICE Import", &
      rcToReturn=rc2))                                      & 
      call ESMF_Finalize(rc=rc)

   ! Add bundle to the import state
   call ESMF_StateAdd(impState,fieldBundleList=(/impBundle/),rc=rc)
   if (ESMF_LogFoundError(rc,                                 &
      msg="cice_setup_Esmf: StateAdd bundle to import state", &
      rcToReturn=rc2))                                        & 
      call ESMF_Finalize(rc=rc)
   call ESMF_LogFlush(rc=rc)


   !---------------------------------------------------------------------
   !-- Add some relevant metadata to gridComp
   !---------------------------------------------------------------------
   call ESMF_AttributeSet(gridComp, name="CICE_ktherm",value=ktherm, rc=rc)
   if (ESMF_LogFoundError(rc, msg="cice_setup_esmf: attributeset CICE_ktherm", rcToReturn=rc2)) call ESMF_Finalize(rc=rc)
   
   call ESMF_AttributeSet(gridComp, name="CICE_tfrz_option",value=tfrz_option, rc=rc)
   if (ESMF_LogFoundError(rc, msg="cice_setup_esmf: attributeset tfrz_option", rcToReturn=rc2)) call ESMF_Finalize(rc=rc)


end subroutine cice_setup_esmf

subroutine cice_put_export(export_state)
   use ice_communicate, only: my_task, master_task
   use ice_blocks,      only : block, get_block
   use ice_domain,      only : nblocks, blocks_ice
   use ice_flux,        only : strocnxT, strocnyT, fhocn, fsalt, fresh, fswthru, frzmlt
   use ice_flux,        only : fhocn_ai, fsalt_ai, fresh_ai, fswthru_ai
   use ice_state,       only : aice, vice, trcr,  uvel, vvel
   use ice_constants,   only : Tffresh
   use ice_grid,        only : u2tgrid_vector
   implicit none
   type(ESMF_State)       :: export_state
   integer rc
   type(block)            :: this_block
   integer i, j, ifld, iblk, ilo, ihi, jlo, jhi, ig, jg

   real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: uwork,vwork

   ! KAL: Note that all fluxes are grid cell averages. This is whats assumed by
   ! HYCOM. For interation with CESM, for instance - it is assumed that fluxes
   ! are ice area averages (!)

   ! Transfer ice velocities to grid cell centers
   ! TODO: Use averaged ice conc and velocity to hycom?
   uwork=uvel
   vwork=vvel
   call u2tgrid_vector(uwork)
   call u2tgrid_vector(vwork)

   call ESMF_LogWrite("CICE Put Export routine called", ESMF_LOGMSG_INFO, rc=rc)
   call ESMF_LogFlush(rc=rc)

   ! Get import fields - Import fields must match export fields of "the other" model
   ! TODO: OMP
   do ifld=1,numExpFields
      !print *,ifld,trim(impFieldName(ifld))
      do iblk=1,nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         ! Total ice concentration
         if (trim(expFieldName(ifld)) == "sic") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               expData(ifld,iblk)%p(ig,jg) = aice(i,j,iblk)
            end do
            end do
         ! Ice-ocean stress
         elseif (trim(expFieldName(ifld)) == "sitx") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               expData(ifld,iblk)%p(ig,jg) = strocnxT(i,j,iblk)*aice(i,j) ! NB: Check T to U and sign
            end do
            end do
         ! Ice-ocean stress
         elseif (trim(expFieldName(ifld)) == "sity") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               expData(ifld,iblk)%p(ig,jg) = strocnyT(i,j,iblk)*aice(i,j) ! NB: Check T to U and sign
            end do
            end do
         ! Penetrating shortwave radiation
         elseif (trim(expFieldName(ifld)) == "siqs") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !expData(ifld,iblk)%p(ig,jg) = fswthru_ai(i,j,iblk) 
               expData(ifld,iblk)%p(ig,jg) = accum_fswthru_ai(i,j,iblk) 
            end do
            end do
         ! Ice Freezing/Melting heat flux
         elseif (trim(expFieldName(ifld)) == "sifh") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               ! KAL Frazil ice creation. Heat ocean with energy used to create ice.
               if (frzmlt(i,j,iblk) > 0.) then
                  expData(ifld,iblk)%p(ig,jg) = frzmlt(i,j,iblk)
               ! KAL Melting from below. fhocn returns part of frzmlt used to
               ! heat ice (cools ocean)
               else 
                  !expData(ifld,iblk)%p(ig,jg) = fhocn_ai(i,j,iblk)
                  expData(ifld,iblk)%p(ig,jg) = accum_fhocn_ai(i,j,iblk)
               end if
            end do
            end do
         ! Ice Freezing/Melting salt flux
         elseif (trim(expFieldName(ifld)) == "sifs") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !expData(ifld,iblk)%p(ig,jg) = fsalt_ai(i,j,iblk) 
               expData(ifld,iblk)%p(ig,jg) = accum_fsalt_ai(i,j,iblk) 
            end do
            end do
         ! Ice Freezing/Melting mass (water) flux
         elseif (trim(expFieldName(ifld)) == "sifw") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !expData(ifld,iblk)%p(ig,jg) = fresh_ai(i,j,iblk) 
               expData(ifld,iblk)%p(ig,jg) = accum_fresh_ai(i,j,iblk) 
            end do
            end do
         ! Sea ice Temperature (here: srf)
         elseif (trim(expFieldName(ifld)) == "sit") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !create Tsrf from the first tracer (trcr) in ice_state.F
               !expData(ifld,iblk)%p(ig,jg) = Tffresh + trcr(i,j,1,iblk)     !Kelvin (original ???)
               expData(ifld,iblk)%p(ig,jg) = trcr(i,j,1,iblk)     !Kelvin (original ???)
            end do
            end do
         ! Sea ice Thickness
         elseif (trim(expFieldName(ifld)) == "sih") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !from vice and aice
               expData(ifld,iblk)%p(ig,jg) = vice(i,j,iblk) / max(aice(i,j,iblk),1e-6)
            end do
            end do
         ! Sea ice U-velocity
         elseif (trim(expFieldName(ifld)) == "siu") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !expData(ifld,iblk)%p(ig,jg) = uvel(i,j,iblk) ! NB: CICE U-point
               ! CICE P-point
               expData(ifld,iblk)%p(ig,jg) = uwork(i,j,iblk)
                 
            end do
            end do
         ! Sea ice V-velocity
         elseif (trim(expFieldName(ifld)) == "siv") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !expData(ifld,iblk)%p(ig,jg) = vvel(i,j,iblk) ! NB: CICE U-point
               ! CICE P-point
               expData(ifld,iblk)%p(ig,jg) = vwork(i,j,iblk)
                 
            end do
            end do
         else 
            if (my_task==master_task .and. iblk==1) then
               print '(a)',"WARN: Unhandled CICE export field "//trim(impFieldName(ifld))
            end if
         end if
      end do
   end do
  !stop '(cice_put_export)'
end subroutine cice_put_export

! --- Extract import state.
subroutine cice_get_import(import_state)
   use ice_communicate, only: my_task, master_task
   use ice_blocks,      only : block, get_block
   use ice_domain,      only : nblocks, blocks_ice
   use ice_flux,        only : frzmlt, uocn, vocn, sss, sst, hmix, &
                               uatm, vatm, Tair, zlvl, potT, &
                               ss_tltx, ss_tlty, Qa, flw, &
                               fsnow, frain, swvdr, swvdf, swidr, swidf
   use ice_grid,        only : t2ugrid_vector
   implicit none
   type(ESMF_State)       :: import_state
   type(block)            :: this_block
   character(len=100) :: msg
   integer rc, rc2, itemCount, ifld, iblk
   integer i, j, ilo, ihi, jlo, jhi, ig, jg

   call ESMF_LogWrite("CICE Get Import routine called", ESMF_LOGMSG_INFO, rc=rc)
   call ESMF_LogFlush(rc=rc)

   ! Get import fields - Import fields must match export fields of "the other" model
   ! TODO: OMP
   do ifld=1,numImpFields
      !KAL !print *,ifld,trim(impFieldName(ifld))
      !KAL if (my_task==master_task .and. iblk==1) print '(a)',"CICE:importing "//trim(impFieldName(ifld))
      do iblk=1,nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         ! Assign array values to relevant target arrays. For now 
         !print *,my_task,iblk,trim(impFieldName(ifld)),impData(ifld,iblk)%p(10,10)
         if (trim(impFieldName(ifld)) == "ssfi") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               frzmlt(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "ssu") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !NB: Must be in T/P-point
               uocn(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "ssv") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !NB: Must be in T/P-point
               vocn(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "sst") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               sst(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "sss") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               sss(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "mlt") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               hmix(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "uair") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !NB: Must be in T/P-point
               uatm(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "vair") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !NB: Must be in T/P-point
               vatm(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "tair") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               Tair(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "zair") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               zlvl(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "pott") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               potT(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "ss_tltx") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !NB: Must be in T/P-point
               ss_tltx(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "ss_tlty") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               !NB: Must be in T/P-point
               ss_tlty(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "Qair") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               Qa(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "swvdr") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               swvdr(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
               swvdf(i,j,iblk) = 0._dbl_kind    ! shortwave radiation (W/m^2)
               swidr(i,j,iblk) = 0._dbl_kind    ! shortwave radiation (W/m^2)
               swidf(i,j,iblk) = 0._dbl_kind    ! shortwave radiation (W/m^2)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "flw") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               flw(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
               !print *,"Setting flw",flw(i,j,iblk)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "rrate") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               frain(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
               !print *,"Setting flw",flw(i,j,iblk)
            end do
            end do
         elseif (trim(impFieldName(ifld)) == "srate") then
            do j=jlo,jhi
            do i=ilo,ihi
               ig  = this_block%i_glob(i)
               jg  = this_block%j_glob(j)
               fsnow(i,j,iblk) = impData(ifld,iblk)%p(ig,jg)
               !print *,"Setting flw",flw(i,j,iblk)
            end do
            end do
         ! TODO: Add stop for unknown fields(?)
         else 
            if (my_task==master_task .and. iblk==1) then
               print '(a)',"WARN: Unhandled CICE import field "//trim(impFieldName(ifld))
            end if
         end if
      end do
   end do
   !stop '(cice_get_import)'

   ! Go from p-points to CICE u-point for tilt and ocean current components
   call t2ugrid_vector(ss_tltx)
   call t2ugrid_vector(ss_tlty)
   call t2ugrid_vector(uocn)
   call t2ugrid_vector(vocn)

end subroutine cice_get_import


subroutine ice_init_esmf(comp, import_state, export_state, EClock, rc)
use CICE_InitMod
use ice_restart_shared, only : runtype, use_restart_time
use ice_calendar   , only : dt, days_per_year, use_leap_years, daymo, &
   nyr, month, mday, hour, sec, npt, calendar, time, year_init, time2sec, idate, &
   basis_seconds
use ice_communicate, only : my_task, master_task
implicit none
! !ARGUMENTS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    !type(ESMF_Time)     :: startTime, stopTime, refTime
    !integer :: start_year, start_month, start_mday, start_sec
    !integer :: stop_year, stop_month, stop_mday, stop_sec
    !real*8  :: tmp1,tmp2
    !character(ESMF_MAXSTR)     :: msg
    integer, intent(out)         :: rc

    call ESMF_LogWrite("CICE init routine started",ESMF_LOGMSG_INFO, rc=rc)
    call CICE_Init() 
    call cice_setup_esmf(comp,import_state,export_state,EClock,rc)
    !call ESMF_LogWrite("CICE init routine ended",ESMF_LOGMSG_INFO, rc=rc)

end subroutine


subroutine ice_run_esmf(comp, import_state, export_state, EClock, rc)
use CICE_RunMod
use ice_calendar, only: istep,npt,write_restart
implicit none
! !ARGUMENTS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

    integer :: rc2
    logical :: get_import, put_export, ice_restart
    character(len=100) :: msg

    call ESMF_LogWrite("CICE run routine started",ESMF_LOGMSG_INFO, rc=rc)
    call ESMF_LogFlush(rc=rc)

    ! Use import state
    call ESMF_AttributeGet(import_state, name="get_import",value=get_import,rc=rc)
    if (ESMF_LogFoundError(rc,  msg="ice_run_esmf: attributeget get_import", rcToReturn=rc2)) & 
       call ESMF_Finalize(rc=rc)
    if (get_import) then
       call cice_get_import(import_state)
    end if

    ! Flag from gridcomp that we should save restart
    call ESMF_AttributeGet(comp, name="ICE_restart",value=ICE_restart,rc=rc)
    if (ESMF_LogFoundError(rc,  msg="ice_run_esmf: attributeget ICE_restart", rcToReturn=rc2)) & 
       call ESMF_Finalize(rc=rc)
    write(msg,'("ice_run_esmf: ICE_restart:",l7)') ice_restart
    call ESMF_LogWrite(msg,ESMF_LOGMSG_INFO, rc=rc)

    ! Run single time step
    call CICE_Run(ice_restart)

    ! Update accumulated fields
    call update_accumulated_fields()

    ! Assign to export state
    call ESMF_AttributeGet(export_state, name="put_export",value=put_export,rc=rc)
    if (ESMF_LogFoundError(rc,  msg="ice_run_esmf: attributeget put_export", rcToReturn=rc2)) & 
       call ESMF_Finalize(rc=rc)
    if (put_export) then
       call average_accumulated_fields()
       call cice_put_export(export_state)
       call reset_accumulated_fields()
    end if

    write(msg,'("istep, npt, write_restart:",i5,i5,i5)') istep,npt,write_restart
    call ESMF_LogWrite(msg,ESMF_LOGMSG_INFO, rc=rc)


    call ESMF_LogWrite("CICE run routine ended",ESMF_LOGMSG_INFO, rc=rc)
    call ESMF_LogFlush(rc=rc)
end subroutine

subroutine ice_final_esmf(comp, import_state, export_state, EClock, rc)
use CICE_FinalMod
implicit none
! !ARGUMENTS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
    ! Use import state
    call ESMF_LogWrite("CICE finalize routine started",ESMF_LOGMSG_INFO, rc=rc)
    call ESMF_LogFlush(rc=rc)
    call CICE_Finalize()
    !call ESMF_LogWrite("CICE finalize routine ended",ESMF_LOGMSG_INFO, rc=rc)
    !call ESMF_LogFlush(rc=rc)
end subroutine

subroutine check_block_bounds (iblk,fldname,lbnd,ubnd)
    use ESMF
    use ice_communicate, only : my_task
    use ice_blocks,      only : block,get_block
    use ice_domain,      only : blocks_ice
    implicit none

    integer, intent(in) :: iblk
    character(len=*), intent(in) :: fldname
    character(ESMF_MAXSTR)     :: msg
    integer,intent(in)         :: lbnd(2),ubnd(2)

    type(block) :: this_block  
    integer     :: ilo,ihi,jlo,jhi
    integer     :: rc,rc2
          
    ! Set grid mask from CICE mask
    this_block = get_block(blocks_ice(iblk),iblk)         
    ilo = this_block%ilo
    ihi = this_block%ihi
    jlo = this_block%jlo
    jhi = this_block%jhi

    ! Another test....
    if ( this_block%i_glob(ilo) .ne. lbnd(1)) then
        write(msg,'("ERROR my_task ",i4, a10, ":lower bound mismatch. lbnd(1)=",i4, "i_glob(ilo)=",i4)') &
           my_task,fldname,lbnd(1),this_block%i_glob(ilo)
        call ESMF_LogWrite(msg,ESMF_LOGMSG_ERROR, rc=rc)
        print '(a)',msg
        call ESMF_Finalize(rc=rc,endflag=ESMF_END_ABORT)
    end if
    if ( this_block%j_glob(jlo) .ne. lbnd(2)) then
        write(msg,'("ERROR my_task ",i4,a10, ":lower bound mismatch. lbnd(2)=",i4, "j_glob(jlo)=",i4)') &
           my_task,fldname,lbnd(2),this_block%j_glob(jlo)
        call ESMF_LogWrite(msg,ESMF_LOGMSG_ERROR, rc=rc)
        print '(a)',msg
        call ESMF_Finalize(rc=rc,endflag=ESMF_END_ABORT)
    end if
    if ( this_block%i_glob(ihi) .ne. ubnd(1)) then
        write(msg,'("ERROR my_task ",i4, a10 ": higher bound mismatch. ubnd(1)=",i4, "i_glob(ihi)=",i4)') &
           my_task,fldname,ubnd(1),this_block%i_glob(ihi)
        call ESMF_LogWrite(msg,ESMF_LOGMSG_ERROR, rc=rc)
        print '(a)',msg
        call ESMF_Finalize(rc=rc,endflag=ESMF_END_ABORT)
    end if
    if ( this_block%j_glob(jhi) .ne. ubnd(2)) then
        write(msg,'("ERROR my_task ",i4, a10, ": higher bound mismatch. ubnd(2)=",i4, "j_glob(jhi)=",i4)') &
           my_task,fldname,ubnd(2),this_block%j_glob(jhi)
        call ESMF_LogWrite(msg,ESMF_LOGMSG_ERROR, rc=rc)
        print '(a)',msg
        call ESMF_Finalize(rc=rc,endflag=ESMF_END_ABORT)
    end if

end subroutine check_block_bounds


subroutine init_accumulated_fields()
   implicit none
   accum_fresh_ai   = c0 ! fresh water flux to ocean (kg/m^2/s)
   accum_fsalt_ai   = c0 ! salt flux to ocean (kg/m^2/s)
   accum_fhocn_ai   = c0 ! net heat flux to ocean (W/m^2)
   accum_fswthru_ai = c0 ! shortwave penetrating to ocean (W/m^2)
   ave_sic          = c0
   ave_uvel         = c0
   ave_vvel         = c0
   accum_time=0.
end subroutine

subroutine reset_accumulated_fields()
   use ice_blocks,      only : block, get_block
   use ice_domain,      only : nblocks, blocks_ice
   use ice_flux,        only : fresh_ai, fsalt_ai, fswthru_ai, fhocn_ai
   implicit none
   type(block)            :: this_block
   character(len=100) :: msg
   integer iblk, i, j, ilo, ihi, jlo, jhi
   do iblk=1,nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      do j=jlo,jhi
      do i=ilo,ihi
         accum_fresh_ai  (i,j,iblk) = c0
         accum_fsalt_ai  (i,j,iblk) = c0
         accum_fhocn_ai  (i,j,iblk) = c0
         accum_fswthru_ai(i,j,iblk) = c0
         ave_sic         (i,j,iblk) = c0
         ave_uvel        (i,j,iblk) = c0
         ave_vvel        (i,j,iblk) = c0
      end do
      end do
   end do
   accum_time = c0
end subroutine

subroutine update_accumulated_fields()
   use ice_blocks,      only : block, get_block
   use ice_domain,      only : nblocks, blocks_ice
   use ice_flux,        only : fresh_ai, fsalt_ai, fswthru_ai, fhocn_ai
   use ice_calendar,    only : dt
   use ice_state,       only : aice, uvel, vvel
   implicit none
   type(block)            :: this_block
   character(len=100) :: msg
   integer iblk, i, j, ilo, ihi, jlo, jhi
   do iblk=1,nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      do j=jlo,jhi
      do i=ilo,ihi
         ! Flux
         accum_fresh_ai  (i,j,iblk) = accum_fresh_ai  (i,j,iblk) + fresh_ai(i,j,iblk)   * dt
         accum_fsalt_ai  (i,j,iblk) = accum_fsalt_ai  (i,j,iblk) + fsalt_ai(i,j,iblk)   * dt
         accum_fhocn_ai  (i,j,iblk) = accum_fhocn_ai  (i,j,iblk) + fhocn_ai(i,j,iblk)   * dt
         accum_fswthru_ai(i,j,iblk) = accum_fswthru_ai(i,j,iblk) + fswthru_ai(i,j,iblk) * dt
         ! Average
         ave_sic         (i,j,iblk) = ave_sic         (i,j,iblk) + aice      (i,j,iblk) * dt
         ave_uvel        (i,j,iblk) = ave_uvel        (i,j,iblk) + uvel      (i,j,iblk) * dt
         ave_vvel        (i,j,iblk) = ave_vvel        (i,j,iblk) + vvel      (i,j,iblk) * dt
      end do
      end do
   end do
   accum_time = accum_time + dt
end subroutine

subroutine average_accumulated_fields()
   use ice_blocks,      only : block, get_block
   use ice_domain,      only : nblocks, blocks_ice
   use ice_flux,        only : fresh_ai, fsalt_ai, fswthru_ai, fhocn_ai
   implicit none
   type(block)            :: this_block
   character(len=100) :: msg
   integer iblk, i, j, ilo, ihi, jlo, jhi
   real (kind=dbl_kind) ::  i_accum_time
   i_accum_time = 1./max(accum_time,1e-4)
   do iblk=1,nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      do j=jlo,jhi
      do i=ilo,ihi
         ! Flux
         accum_fresh_ai  (i,j,iblk) = accum_fresh_ai  (i,j,iblk) * i_accum_time
         accum_fsalt_ai  (i,j,iblk) = accum_fsalt_ai  (i,j,iblk) * i_accum_time
         accum_fhocn_ai  (i,j,iblk) = accum_fhocn_ai  (i,j,iblk) * i_accum_time
         accum_fswthru_ai(i,j,iblk) = accum_fswthru_ai(i,j,iblk) * i_accum_time
         ! Average
         ave_sic         (i,j,iblk) = ave_sic         (i,j,iblk) * i_accum_time
         ave_uvel        (i,j,iblk) = ave_uvel        (i,j,iblk) * i_accum_time
         ave_vvel        (i,j,iblk) = ave_vvel        (i,j,iblk) * i_accum_time
      end do
      end do
   end do
end subroutine



end module
