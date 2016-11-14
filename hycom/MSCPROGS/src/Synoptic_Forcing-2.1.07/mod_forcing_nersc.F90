module mod_forcing_nersc
   use mod_xc
   implicit none

   ! Temporal arrays for reading forcing fields - synoptic
   real, allocatable, dimension(:,:)  ::                synuwind, &
     synvwind, synwndspd, synairtmp, synrelhum, synprecip, &
     synclouds, syntaux, syntauy, synvapmix,&
     synradflx, synshwflx, synslp, synssr

   ! Temporal variables - used when reading climatology
   real, allocatable, dimension(:,:,:)  :: &
      clmuwind, clmvwind, clmwndspd, clmairtmp, clmrelhum, clmprecip, &
      clmclouds, clmsst,  clmsss, clmtaux, clmtauy, clmvapmix, &
      clmradflx, clmshwflx, clmslp

   


   ! Logical vars, denote forcing fields read
   logical,save ::         &
     lsynwnd    = .false., &
     lsynairtmp = .false., &
     lsynrelhum = .false., &
     lsynprecip = .false., &
     lsynclouds = .false., &
     lsynvapmix = .false., &
     lsynradflx = .false., &
     lsynshwflx = .false., &
     lsynslp    = .false., &
     lsynssr    = .false.


! ! Forcing parameters -- particular to nersc (in addition to 
! ! those in blkdat)
! real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::   &
!    uwind  ,    &! Winds            [m/s]
!    vwind  ,    &! Winds            [m/s]
!    clouds ,    &! Cloud Cover      [0-1]
!    slp    ,    &! Cloud Cover      [0-1]
!    relhum       ! Relative humidity[0-1]
!
!  ! For testing original hycom thermf
!  real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::   &
!    tstsurflx,tstsalflx
!
!  ! For diagnosing relaxation fields
!  real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::   &
!    relsurflx,relsalflx
!
  ! Units to connect to for reading fields above
  integer, parameter :: &
     unit_uwind  =111, &
     unit_vwind  =112, &
     unit_clouds =113, &
     unit_relhum =114, &
     unit_slp    =115   

  ! Time info -- Also used to catch last time
  ! Data was read in randforc
  real*8,save ::  &
   dtime_uwind   ,&
   dtime_vwind   ,&
   dtime_relhum  ,&
   dtime_clouds  ,&
   dtime_slp

  ! 
  real, parameter :: cd=0.0012
  real,parameter :: airdns  =    1.2
  real,parameter :: slp0=1012.5


  real, save :: time_srelax
  real, save :: time_trelax
  real, save :: dtime_read(4)



     

end module mod_forcing_nersc
