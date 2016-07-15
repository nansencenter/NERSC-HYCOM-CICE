module mod_necessary_ecovars

  use mod_xc
  implicit none 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  * * * I M P O R T A N T * * *
!
!  The variables in this module must be defined for ANY biochemical model coupled to HYCOM, 
!  i.e., same names, and with nbio set to the correct number of compartments in your model
!  and eco_tag set to a 5 character tag for your model.
!
!  Then, in mxkpp, tsadvc, regrid_tracer, read_write_restart_bio, and hycom 
!  (+ your specific ecosystem routines), 
!  make a use statement to the correct version of this module...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef EVA85
  integer, PARAMETER :: nbio=3         ! number of biochemical compartments
  character(len=5) :: eco_tag='EVA85'      ! tag for the current ecomodel
#endif  
#ifdef FDM02
  integer, PARAMETER :: nbio=11         ! number of biochemical compartments
  character(len=5) :: eco_tag='FDM02'      ! tag for the current ecomodel
#endif  
#ifdef SCH02
  integer, PARAMETER :: nbio=15         ! number of biochemical compartments
  character(len=5) :: eco_tag='SCH02'      ! tag for the current ecomodel
#endif  
#if defined (ECOSM)
  integer, PARAMETER :: nsed=3         ! number of sediment compartments
  integer, PARAMETER :: nbiox=14            ! number of biochemical compartments
  integer, PARAMETER :: ndbiox=6        ! diagnostics from biology
  character(len=5) :: eco_tag='ECOSM'      ! tag for the current ecomodel
  integer, PARAMETER :: ifla=1,idia=2,imicro=3,imeso=4,idet = 5, inh4= 6,  &
                        idom =7,init=8,ipho=9,isil =10,ioxy =11,  iopa=12, &
                        ibg=13,ifish=14   !field: bio
  integer, PARAMETER ::ised1=1,ised2=2,ised3=3   !field: bot_layer
  integer, PARAMETER ::idiap=1,iflap=2, inpp=3,ispp=4, &
                       izsp=5,izlp=6   !field: bio_diagn
                      
#if defined (ECO2)
  integer, PARAMETER :: nbio=nbiox+2        ! number of biochemical compartments+DIC+ALK
  integer, PARAMETER :: ndbio=ndbiox+3        ! diagnostics from biology+ carbon diagnostics
  integer, PARAMETER :: idic=nbiox+1,ialk=nbiox+2  !field: bio
  integer, PARAMETER :: ipco=ndbiox+1,iph=ndbiox+2,iatc=ndbiox+3   !field: bio_diagn
#else
  integer, PARAMETER :: nbio=nbiox        ! number of biochemical compartments+DIC+ALK
  integer, PARAMETER :: ndbio=ndbiox        ! diagnostics from biology+ carbon diagnostics

#endif 
#endif 

#if defined (NOR05)
  integer, PARAMETER :: ndbio=9        ! diagnostics from biology
#if defined (ZOOPL)
  integer, PARAMETER :: nbio=11         ! number of biochemical compartments
#elif defined (DETPHO)
  integer, PARAMETER :: nbio=9         ! number of biochemical compartments
#else
  integer, PARAMETER :: nbio=8         ! number of biochemical compartments
#endif /* DETPHO-ZOOPL */
  character(len=5) :: eco_tag='NOR05'      ! tag for the current ecomodel
  integer, PARAMETER :: init= 1, ipho = 2, isil = 3,  idet  = 4,  isis= 5, &
                        ifla= 6, idia = 7, ioxy = 8,  ised  = 12, iyel=13, &
                        idetp=9, imeso=10, imicro=11, icha=14;
!AS15062010: 
!The mortality cc(3) maqybe a bit high when ZOOPL is not defined
!   CC(1)    : FRACTION OF PHOSPHATE AND NITRATE IN A CELL
!   CC(2)    : FRACTION OF SILICATE AND NITRATE IN A CELL
!   CC(3)    : DEATH RATE
!   CC(4)    : RATE OF DECOMPOSITION OF DETRITUS
#if defined (ZOOPL)
  real, PARAMETER :: CC(1:4)=(/0.138,1.15,4.0E-7,1.52E-7/)
#else
  real, PARAMETER :: cc(1:4)=(/0.138,1.15,1.6E-6,1.52E-7/)
#endif
#endif

!previously defined in blkdat_eco
  integer, PARAMETER :: sinkfl=1
  integer, PARAMETER :: rivflg=1  !if 1 add river utrients
  integer, PARAMETER :: sedflg=1


!!! just for debugging
  integer, PARAMETER :: susp_i=104
  integer, PARAMETER :: susp_j=6
  integer, PARAMETER :: susp_k=14
  integer, PARAMETER :: susp_var=idia

! time step-used for the biochem. sources and sinks: = baclin*nphys_sou:
  real :: biodt_sou                           


  real pp_m2                              ! volume of box-area for primp comp.
  
#if defined (NOR05)
  real :: bot_layer(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2,nbio  )     ! Sediment layer
  real :: grosspp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)            ! gross primary production      
  real :: netpp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)              ! net primary production      
#endif

#if defined (ECOSM)
  real :: bot_layer(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2,nsed  )     ! Sediment layer
#endif
  real :: bio_diagn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,ndbio)    ! diagnostic biochemical compartments
  real :: bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2,nbio)         ! prognostic biochemical compartments
  real :: rad(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)         ! biochemical compartments
  real :: rivnit(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                 ! river nitrogen
  real :: rivpho(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                 ! river phosphate
  real :: rivsil(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                 ! river silicate
  real :: rivdon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                 ! river DON
  real :: wvel(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)               ! veritcal velocity      
!  logical pp_area(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)               ! Area for primp   
           
end module mod_necessary_ecovars
