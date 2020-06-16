      MODULE mod_cpl_oasis_init
      !!----------------------------------------------------------------------
#ifdef CPL_OASIS_HYCOM
      use mod_oasis                    ! OASIS3-MCT module
#endif
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_dimensions

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      PRIVATE

      INTEGER, PUBLIC :: localComm  ! local MPI communicator for a particular grid
      !AS - from SHOM:INTEGER :: localComm_save ! local MPI communicator for all processes
      CHARACTER(len=6), PUBLIC :: comp_name='hycomm'
      INTEGER, PUBLIC          :: ncomp_id          ! component identification returned by oasis_init_comp
      INTEGER, PUBLIC          :: ierror            ! return error code
! Variable necessary for writing grid
      INTEGER, PUBLIC          :: il_flag  ! Flag for grid writing by proc 0
      INTEGER, PUBLIC, DIMENSION(:,:), POINTER     :: indice_mask ! mask, 0 == valid point, 1 == masked point 
      INTEGER, PUBLIC          :: il_part_id

      PUBLIC  cpl_init
      INTEGER, PUBLIC::   oasis_info

      INTEGER, PARAMETER, PUBLIC  ::   o2i_sst  =  1            ! ocean SST on p-grid
      INTEGER, PARAMETER, PUBLIC  ::   o2i_sss  =  2            ! ocean SSS on p-grid
      INTEGER, PARAMETER, PUBLIC  ::   o2i_uve  =  3            ! ocean U-velocity on p-grid
      INTEGER, PARAMETER, PUBLIC  ::   o2i_vve  =  4            ! ocean V-velocity on p-grid
      INTEGER, PARAMETER, PUBLIC  ::   o2i_ssh  =  5            ! ocean Sea Surface Height
      INTEGER, PARAMETER, PUBLIC  ::   o2i_1st  =  7            ! Depth of the first layer
      INTEGER, PARAMETER, PUBLIC ::    o2i_frr  =  6    ! Ocean Fraction of solar net radiation absorbed in the first layer, 
                                                       !  needs field from mixing routine

      INTEGER, PARAMETER, PUBLIC  ::   i2o_taux   =  1            ! Eastward wind stress 
      INTEGER, PARAMETER, PUBLIC  ::   i2o_tauy   =  2            ! Northward wind stress
      INTEGER, PARAMETER, PUBLIC  ::   i2o_fwfl   =  3            ! Freshwater flux [kg/m^2]
      INTEGER, PARAMETER, PUBLIC  ::   i2o_lwlh   =  4            ! Total non soler heat flux
      INTEGER, PARAMETER, PUBLIC  ::   i2o_swra   =  5            ! Totol short wave radiation
      INTEGER, PARAMETER, PUBLIC  ::   i2o_saln   =  6            ! Salt flux [kg/m^2] 
      INTEGER, PARAMETER, PUBLIC  ::   i2o_taut   =  7            ! Wind stress module
      INTEGER, PARAMETER, PUBLIC  ::   i2o_sico   =  8            ! Sea ice cover
 
      INTEGER, PUBLIC, PARAMETER ::         nrcv=8         ! total number of fields received 
      INTEGER, PUBLIC, parameter ::         nsnd=6         ! total number of fields sent 
      INTEGER, PUBLIC                    ::   ncplmodel    ! Maximum number of models to/from which HYCOM is potentialy sending/receiving data
      INTEGER, PUBLIC, PARAMETER ::   nmaxfld=10   ! Maximum number of coupling fields
      INTEGER, PUBLIC, PARAMETER ::   nmaxcat=5    ! Maximum number of coupling categories
      INTEGER, PUBLIC, PARAMETER ::   nmaxcpl=5    ! Maximum number of coupling components
      INTEGER, PUBLIC, PARAMETER ::   numnam_ref=100  
   
!AS define new type of coupling field:
      TYPE, PUBLIC ::   FLD_CPL               !: Type for coupling field information
      CHARACTER(len = 8)    ::   var_name    ! Name of the coupling field   
      INTEGER               ::   var_id
      CHARACTER(len = 1)    ::   var_grid='p'    ! Grid type  - for couple ing to Nextsim - p_grid by default
      INTEGER, dimension(2) ::   var_nodims = (/1,1/)
      INTEGER, dimension(4) ::   var_actual_shape      
      REAL*8, allocatable, dimension(:,:) ::   var
      INTEGER               :: info
      END TYPE FLD_CPL

      TYPE(FLD_CPL), DIMENSION(nmaxfld), PUBLIC ::   flds_recv, flds_send   !: Coupling fields

      !!----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE cpl_init( cd_modname, kl_comm )
      use mod_oasis
      IMPLICIT NONE
      !!--------------------------------------------------------------------
      CHARACTER(len = *), INTENT(in) ::   cd_modname   ! model name as set in namcouple file
      INTEGER          , INTENT(out) ::   kl_comm      ! local communicator of the model
      !!--------------------------------------------------------------------

      !------------------------------------------------------------------
      ! First Initialize the OASIS system for the application
      !------------------------------------------------------------------
      !AS-CHECKING      
      print*, 'Before call initialization of oasis, ierror', ierror, 'ncomp_id', ncomp_id
      !AS-CHECKING      

      CALL oasis_init_comp ( ncomp_id, TRIM(cd_modname), ierror )

      !AS-CHECKING      
      print*, 'Called initialization of oasis, ierror', ierror, 'ncomp_id', ncomp_id

      !AS-CHECKING      
      IF ( ierror /= OASIS_Ok ) &
        CALL oasis_abort (ncomp_id,'cpl_init','Failure oasis_init_comp')

      !------------------------------------------------------------------
      ! 3rd Get an MPI communicator for OPA local communication
      !------------------------------------------------------------------
      
      CALL oasis_get_localcomm ( kl_comm, ierror )
      !AS-CHECKING      

      print*, 'Called get localcomm, ierror', ierror, 'localComm', kl_comm
      !AS-CHECKING      
      IF ( ierror /= OASIS_Ok ) &
         CALL oasis_abort (ncomp_id,'cpl_init',&
              'Failure oasis_get_localcomm' )


      !! DEFINE ALL FIELDS RELATED TO SENDING FIELDS - names here must correspond to namcouple
      flds_send(o2i_sst)%var_name = 'O_SSTSST'  
      flds_send(o2i_sss)%var_name = 'O_SSSal'  
      flds_send(o2i_uve)%var_name = 'O_OCurx1'  
      flds_send(o2i_vve)%var_name = 'O_OCury1'  
      flds_send(o2i_ssh)%var_name = 'O_SSHght'  
      flds_send(o2i_1st)%var_name = 'O_E3T1st'  
      flds_send(o2i_frr)%var_name = 'O_FraQsr'  

      !! DEFINE ALL FIELDS RELATED TO RECIEVING FIELDS  - names here must correspond to namcouple

      flds_recv(i2o_taux)%var_name = 'O_OTaux1'  
      flds_recv(i2o_tauy)%var_name = 'O_OTauy1'  
      flds_recv(i2o_fwfl)%var_name = 'OOEvaMPr'  
      flds_recv(i2o_lwlh)%var_name = 'O_QnsOce'  
      flds_recv(i2o_swra)%var_name = 'O_QsrOce'  
      flds_recv(i2o_saln)%var_name = 'O_SFLX'  
      flds_recv(i2o_taut)%var_name = 'O_TauMod'  
      flds_recv(i2o_sico)%var_name = 'RIceFrc'  
       
      print*, 'Finnished setting the field names: ' 
      print*,  flds_send(1:nsnd)%var_name
      print*,  flds_recv(1:nrcv)%var_name

      END SUBROUTINE cpl_init

      END MODULE mod_cpl_oasis_init
