      MODULE mod_cpl_oasis
      !!----------------------------------------------------------------------
#ifdef CPL_OASIS_HYCOM
      use mod_oasis                    ! OASIS3-MCT module
#endif
      use mod_cpl_oasis_init
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_dimensions
      use mod_xc         ! HYCOM communication interface

      IMPLICIT NONE

      PUBLIC  cpl_write_grid
      PUBLIC  cpl_enddef
      PUBLIC  cpl_define, cpl_cleanup
      PUBLIC  cpl_snd, cpl_rcv
      
      CONTAINS

!AS subroutine for writing grid inforamtion needed by oasis.
      SUBROUTINE cpl_write_grid(flag, part_id)
      use mod_oasis

      IMPLICIT NONE
      !
      INTEGER, INTENT(out) :: flag
      INTEGER, INTENT(in) :: part_id 
      INTEGER, parameter :: nc=4
      REAL, dimension(itdm,jtdm) :: gfld1, gfld2,gip 
      !
      ! Mask inversion to follow (historical) OASIS convention (0=not masked;1=masked)
      ALLOCATE(indice_mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), STAT=ierror )
      IF ( ierror /= 0 ) print*, 'Error allocating indice_mask'

! Grid in pressure coordinate
      WHERE(ip == 1) 
          indice_mask=0
      ELSEWHERE
          indice_mask=1
      END WHERE
      ! 
!AS OBS! routine does noe overwrite so make sure you move files from SCRATCH-Area
      CALL oasis_start_grids_writing(flag)
      CALL xcaget(gfld1,plon,1)
      CALL xcaget(gfld2,plat,1)
      if (mnproc==1) then
         CALL oasis_write_grid('oc_p', itdm, jtdm,gfld1,gfld2)
      end if
!      CALL oasis_write_corner('torc', idm+2*nbdy, jdm+2*nbdy, nc, clon, clat) !AS only needed if we use SCRIPR/CONSERV sec 4.3
      CALL xcaget(gip,real(indice_mask),1)
      CALL xcaget(gfld1,scp2,1)
      if (mnproc==1) then
         CALL oasis_write_mask('oc_p', itdm, jtdm, int(gip))
         CALL oasis_write_area('oc_p', itdm, jtdm, gfld1) 
      end if
      CALL oasis_terminate_grids_writing()

      END SUBROUTINE cpl_write_grid

      SUBROUTINE cpl_define( krcv, ksnd, kcplmodel, part_id )
      use mod_oasis
      IMPLICIT NONE
      !!-------------------------------------------------------------------
      !! ** Purpose :   Define grid and field information for ocean
      !!    exchange between HYCOM, NEXTSIM and coupler.
      !!
      !!              * Read namsbc_cpl namelist 
      !!              * define the receive interface
      !!              * define the send    interface
      !!              * initialise the OASIS coupler

      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   krcv, ksnd     ! Number of received and sent coupling fields
      INTEGER, INTENT(in) ::   kcplmodel      ! Maximum number of models to/from which NEMO is potentialy sending/receiving data
      !
      INTEGER, INTENT(out) :: part_id
      INTEGER :: il_paral(5)       ! OASIS3 box partition
      INTEGER :: ishape(2,2)    ! shape of arrays passed to PSMILe
      INTEGER :: ji,jc,jm,ifld       ! local loop indicees
      CHARACTER(LEN=64) :: zclname
      CHARACTER(LEN=2) :: cli2
      !!--------------------------------------------------------------------
      ! HERE I NEED TO SPECIFY PROPERLY FOR FOLLOWING STEPS:
      ! (1) some checks 
      ! (2) partioning used in oasis for exchange
      ! (3) preparing & allocation of send & reciewve fileds
      !     here, I define the variables in a pointer-wise manner. 
      ! -----------------------------------------------------------------
      !!
      INTEGER ::   ios  ! Local integer output status for namelist read
      INTEGER ::   lwp  ! Local integer output status for namelist read
      INTEGER ::   maxrcv 


      !  Define the partition 
      ! -----------------------------------------------------------------
      ! TODO: please set here the information related to the box partitioning 
      il_paral(1)=2             !2 (indicates a Box partition)
      il_paral(2)=j0*itdm+i0    !the upper left corner global offset
      il_paral(3)=ii           !the local extent in x
      il_paral(4)=jj          !the local extent in y
      il_paral(5)=itdm          !the global extent in x.
 
      print*, 'Partit I', i0,ii,idm
      print*, 'Partit J', j0,jj,jdm
      call oasis_def_partition(part_id,il_paral,ierror,itdm*jtdm)

      print*, 'IL_PARAL', il_paral(2:5), mnproc
      print*, 'PAR_ID', part_id, ierror

      ncplmodel = kcplmodel
      IF( kcplmodel > nmaxcpl ) THEN
        CALL oasis_abort ( ncomp_id,'cpl_define','ncplmodel is larger & 
          than nmaxcpl, increase nmaxcpl')   ;   RETURN
      ENDIF

      IF( max(nrcv,nsnd) > nmaxfld ) THEN
         CALL oasis_abort ( ncomp_id,'cpl_define','nrcv or nsend &
         is larger than nmaxfld, increase nmaxfld')   ;   RETURN
      ENDIF

      do ifld=1,nsnd 
         flds_send(ifld)%var_actual_shape=(/1,il_paral(4),1,il_paral(4)/)
         CALL oasis_def_var(flds_send(ifld)%var_id,flds_send(ifld)%var_name,il_part_id,flds_send(ifld)%var_nodims, &
                            OASIS_Out,flds_send(ifld)%var_actual_shape,OASIS_Real,ierror)     
         allocate(flds_send(ifld)%var(idm,jdm))
      end do

      maxrcv = nrcv

      ! One coupling less (MSLP) if no BGC model
      IF ( ntracr == 0 ) maxrcv = nrcv - 1

      do ifld=1,maxrcv 
         flds_recv(ifld)%var_actual_shape=(/1,il_paral(4),1,il_paral(4)/)
         CALL oasis_def_var(flds_recv(ifld)%var_id,flds_recv(ifld)%var_name,il_part_id,flds_recv(ifld)%var_nodims, &
                            OASIS_In,flds_recv(ifld)%var_actual_shape,OASIS_Real,ierror)     
         allocate(flds_recv(ifld)%var(idm,jdm))
      end do

      !EM Copy of incoming field from previous coupling time step
      allocate(cplts_recv(idm,jdm,nrcv))

      END SUBROUTINE cpl_define

     
      SUBROUTINE cpl_cleanup
      use mod_oasis
      IMPLICIT NONE
      INTEGER :: ifld       ! local loop indicees 
      do ifld=1,nsnd
         deallocate(flds_send(ifld)%var)
      end do
  
      do ifld=1,nrcv
         deallocate(flds_recv(ifld)%var)
      end do

      !EM Copy of incoming field from previous coupling time step
      deallocate(cplts_recv)

      CALL oasis_terminate (ierror)

      END SUBROUTINE cpl_cleanup


      SUBROUTINE cpl_snd(time, n, write_res, info )
      use mod_oasis
      IMPLICIT NONE

      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_snd  ***
      !!
      !! ** Purpose : - Sending fields at each coupling time-step
      !
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(out) ::   info     ! OASIS3 info argument
      INTEGER                   , INTENT(in ) :: time,n     ! ocean time-step in seconds
      LOGICAL                   , INTENT(in ) :: write_res  ! write restart or not 
      !!
      INTEGER    :: ifld
      !!--------------------------------------------------------------------

      flds_send(o2i_sst)%var=temp(1:idm,1:jdm,1,n)
      flds_send(o2i_sss)%var=saln(1:idm,1:jdm,1,n)

      util1 = u(:,:,1,n)+ubavg(:,:,n)  ! total velocity
      call cpl_u2p(util1)
      flds_send(o2i_uve)%var=util1(1:idm,1:jdm)


      util1 = v(:,:,1,n)+vbavg(:,:,n)  ! total velocity
      call cpl_v2p(util1)
      flds_send(o2i_vve)%var=util1(1:idm,1:jdm)

      flds_send(o2i_ssh)%var=srfhgt(1:idm,1:jdm)/g
      flds_send(o2i_1st)%var=dp(1:idm,1:jdm,1,n)/onem

!AS      flds_send(o2i_frr)%var is set to the correct value inside mxkprf.F

      do ifld=1,nsnd 
         CALL oasis_put(flds_send(ifld)%var_id,time,flds_send(ifld)%var(1:ii,1:jj),&
         flds_send(ifld)%info, write_restart=write_res)     
!EM bug ?         flds_recv(ifld)%info)     
      end do

      END SUBROUTINE cpl_snd


      SUBROUTINE cpl_rcv(time, info )
      use mod_oasis
      IMPLICIT NONE

      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_snd  ***
      !!
      !! ** Purpose : - Sending fields at each coupling time-step
      !
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(out) ::   info     ! OASIS3 info argument
      INTEGER                   , INTENT(in ) :: time     ! ocean time-step in seconds
      !!
      INTEGER    :: ifld, maxrcv
      !!--------------------------------------------------------------------

      maxrcv = nrcv

      ! One coupling less (MSLP) if no BGC model
      IF ( ntracr == 0 ) maxrcv = nrcv - 1

      do ifld=1, maxrcv

         CALL oasis_get(flds_recv(ifld)%var_id,time,flds_recv(ifld)%var(1:ii,1:jj),&
         flds_recv(ifld)%info)     

         ! EM: if coupling field actually received at current time step,
         !     save into cplts_recv array to be use in the model until
         !     the next update by oasis_get call
         if (flds_recv(ifld)%info>0) &
            cplts_recv(1:ii,1:jj,ifld)=flds_recv(ifld)%var(1:ii,1:jj)

         ! halo filling

      end do
!AS If receives put variables to the hycom variables if they were received from nextsim
!      if (flds_recv(i2o_taux)%info>0) then
!          surtx(1:ii,1:jj)=flds_recv(i2o_taux)%var(1:ii,1:jj)
!      end if
!      if (flds_recv(i2o_tauy)%info==3.or.flds_recv(i2o_tauy)%info==12) then
!          surty(1:ii,1:jj)=flds_recv(i2o_tauy)%var(1:ii,1:jj)
!      end if

! for filling the halos: - look at the cice coupling:
!     call xctilr(surtx,1,1,nbdy,nbdy,halo_pv) 
!     call xctilr(surty,1,1,nbdy,nbdy,halo_pv) 

      END SUBROUTINE cpl_rcv


      SUBROUTINE cpl_enddef
!AS  May move this routine into the define routine later.
      use mod_oasis
      IMPLICIT NONE

      call oasis_enddef(ierror)
      print*, 'Called routine to finalize definition in oasis, ierror', ierror
      END SUBROUTINE cpl_enddef

      SUBROUTINE cpl_u2p(fld)
      ! Subroutine to regrid fileds from the u-grid to the p-grid
          implicit none
          real, intent(inout) :: fld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
          real :: tmpfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

          tmpfld(idm+nbdy,:)=fld(idm+nbdy,:)
          tmpfld(1-nbdy:idm+nbdy-1,:)  = 0.5*(fld(1-nbdy:idm+nbdy-1,:) + fld(2-nbdy:idm+nbdy,:))
          
          fld=tmpfld

      END SUBROUTINE cpl_u2p

      SUBROUTINE cpl_p2u(fld)
      ! Subroutine to regrid fileds from the p-grid to the u-grid                                                                                    
          implicit none
          real, intent(inout) :: fld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
          real :: tmpfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

          tmpfld(1-nbdy,:)=fld(1-nbdy,:)
          tmpfld(2-nbdy:idm+nbdy,:)  = 0.5*(fld(1-nbdy:idm+nbdy-1,:) + fld(2-nbdy:idm+nbdy,:))

          fld=tmpfld

      END SUBROUTINE cpl_p2u

      SUBROUTINE cpl_v2p(fld)
      ! Subroutine to regrid fileds from the v-grid to the p-grid                                                                                  
          implicit none
          real, intent(inout) :: fld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
          real :: tmpfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

          tmpfld(:,jdm+nbdy)=fld(:,jdm+nbdy)
          tmpfld(:,1-nbdy:jdm+nbdy-1)  = 0.5*(fld(:,1-nbdy:jdm+nbdy-1) + fld(:,2-nbdy:jdm+nbdy))
          
          fld=tmpfld

      END SUBROUTINE cpl_v2p

      SUBROUTINE cpl_p2v(fld)
      ! Subroutine to regrid fileds from the v-grid to the p-grid                                                                                  
          implicit none
          real, intent(inout) :: fld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
          real :: tmpfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

          tmpfld(:,1-nbdy)=fld(:,1-nbdy)
          tmpfld(:,2-nbdy:jdm+nbdy)  = 0.5*(fld(:,1-nbdy:jdm+nbdy-1) + fld(:,2-nbdy:jdm+nbdy))
          
          fld=tmpfld

      END SUBROUTINE cpl_p2v

      END MODULE mod_cpl_oasis
