module mod_forcing_nersc
! -- Module contains time info, forcing info, climate/synaptic data holders,
! -- parameters flags and various routines for processing data which is common
! -- to all forcing options
!
! -- init_forcing_nersc - parces input data (file, command line and environment)
!                         sets synoptic flags for different forcing options
!                         (hardcoded). 
! -- init_fvars_nersc   - Allocates and initializes climate and synoptic vars
!                         (private routine called by init_forcing_nersc)
! -  opnfrcwr           - Routine for opening output forcing files - reduces 
!                         clutter in forfun_nersc. Opened for WRITING
!
! --------------------------------------------------------------------------
!
! NB - call xcspmd() before using this module - otherwise array dimensions are
!      undefined



   use mod_xc
   implicit none

   ! Time and forcing flags
   character(len=5), save :: rforce   ! Forcing (ecmwf/ncep ..)
   character(len=5), save :: clmflag  ! Climate forcing to use (era40 or "old")
   real*8, save :: time1,time2        ! Start and end times of processing
   integer, save :: yrflag            ! yrflag from blkdat.input
   integer, save :: wndflag           ! wndflag from blkdat.input(1=u-pt,2=p-pt)
   real   , save :: delt1             ! Time step of model
   integer, save :: nday1,nhour1      ! First day and hour of integration
   integer, save :: nday2,nhour2      ! Last  day and hour of integration
   real, save :: rday1, rday2         ! First and last julian day of integration
   integer, save  :: refyear          ! Reference year of integration
   real,save ::  w_switch=0.          ! fanf: "weight of the linear combinaison after
   real,save ::  g_switch=1e10        ! fanf: "gap before the switch"
   real,save ::  l_switch             ! fanf: "length of the smooth transition"
   real*8, save :: rdtime=0.            ! Time step of forcing update

   integer, parameter :: rotateflag=1! 1=lonlat rotate, 2=sphere rotate
                                     ! For now this is hardcoded to 1

! Reconfigurable  paths


   character(len=80), save ::   &     ! headers for forcing files,
     forc_header_syn, forc_header_clm ! synoptic and climate headers
   character(len=1),parameter :: flnmfor=' '


   ! Temporal arrays for reading forcing fields - synoptic
   real, allocatable, dimension(:,:)  ::                synuwind, &
     synvwind, synwndspd, synairtmp, synrelhum, synprecip, &
     synclouds, syntaux, syntauy, synvapmix,&
     synradflx, synshwflx, synslp, synssr, cawdir

   ! Logical vars, denote forcing fields read
   logical,save ::         &
     lsynuwnd    = .false., &
     lsynvwnd    = .false., &
     lsynairtmp = .false., &
     lsynprecip = .false., &
     lsynvapmix = .false., &
     lsynradflx = .false., &
     lsynshwflx = .false., &
     lsynslp    = .false.

  ! Parameters used in forcing calc
  real, parameter :: cd=0.0012
  real,parameter :: airdns  =    1.2
  real,parameter :: slp0=1012.5
  real, parameter :: stefanb=5.67e-8
  real, parameter :: t0deg  =273.15


  ! Units to connect to for reading fields above
  integer, parameter :: &
     unit_uwind  =921, &
     unit_vwind  =922, &
     unit_slp    =925 

!   integer, parameter :: nrunits=14
!   integer, dimension(nrunits) :: units=                 &
!       (/901,901,902,903,904,905,906,907,908,unit_relhum,&
!         unit_uwind,unit_vwind,unit_clouds,unit_slp/)

   private init_fvars_nersc

contains


   subroutine init_fvars_nersc()
   implicit none
   ! Allocate fields - some are not used...
   allocate(synuwind (idm,jdm))
   allocate(synvwind (idm,jdm))
   allocate(synwndspd(idm,jdm))
   allocate(syntaux  (idm,jdm))
   allocate(syntauy  (idm,jdm))
   allocate(synvapmix(idm,jdm))
   allocate(synairtmp(idm,jdm))
   allocate(synrelhum(idm,jdm))
   allocate(synprecip(idm,jdm))
   allocate(synclouds(idm,jdm))
   allocate(synradflx(idm,jdm))
   allocate(synshwflx(idm,jdm))
   allocate(synslp   (idm,jdm))
   allocate(cawdir   (idm,jdm))
   synuwind (:,:)=0.
   synvwind (:,:)=0.
   synwndspd(:,:)=0.
   syntaux  (:,:)=0.
   syntauy  (:,:)=0.
   synvapmix(:,:)=0.
   synairtmp(:,:)=0.
   synrelhum(:,:)=0.
   synprecip(:,:)=0.
   synclouds(:,:)=0.
   synradflx(:,:)=0.
   synshwflx(:,:)=0.
   synslp   (:,:)=0.
   cawdir   (:,:)=0.
   end subroutine


   subroutine init_forcing_nersc()
   use mod_year_info22
   use m_parse_blkdat
   implicit none
   integer :: ivar, iday, ihr
   real ::rvar, fver
   real*8 dtime,dtime2
   integer*4, external ::iargc
   logical  :: lflux, ex, lclim=.true.
   character(len=80) :: tmparg,cenv
   type(year_info) :: rtst


   !Transition flag - for now an argument to this routine
   if (iargc()==2) then
      call getarg(1,rforce)
      call getarg(2,clmflag)
      l_switch=1.
      g_switch=1e10 ! Ie; switch will never happen...
   elseif (iargc()==4) then
      call getarg(1,rforce)
      call getarg(2,clmflag)
      call getarg(3,tmparg)
      read(tmparg,*) g_switch
      call getarg(4,tmparg)
      read(tmparg,*) l_switch
      print *
      print '(a)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print '(a)','!!!!!!!!!  Using Climatology transition !!!!!!!!!'
      print '(a,f10.3,x,f10.3,a)', &
             '!!!!!!!!! params:',g_switch,l_switch, &
              '  !!!!!!!!!'  
      print '(a)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *
   else
      print *
      print '(a)','This routine will only calculate perturbed fields used by NERSC-HYCOM'
      print '(a)','You must have the synoptic fields done by Python scripts'
      print '(a)','datasets to use.'
      print *
      print '(a)','Usage :'
      print '(a)', '  force_perturb-2.2 rforce clmflag '
!      print '(a)', '  forfun_nersc rforce clmflag g_switch l_switch'
!      print '(a)','rforce is forcing option, clmflag is climate flag'
      print *
      call exit(1)
   end if

   ! Header flag used in forcing files, set flags on which files to
   ! read from synoptic fields, and which to base on climatology.
   ! This only applies if synoptic forcing is used
   forc_header_syn=''
   lclim=.false.
   if (trim(rforce)=='era-i') then
      forc_header_syn='ERA40 reanalyzed fields after Python2.7'
      lsynslp   =.true.
      lsynuwnd   =.true.
      lsynvwnd   =.true.
      lsynvapmix=.true.
      lsynairtmp=.true.
      lsynprecip=.true.
      lsynshwflx=.false.
      lsynradflx=.false.
      write(lp,'(a)')'Perturbing ERA-I Forcing'
      rdtime=6.d0/24.d0

   elseif (trim(rforce)=='month') then
      forc_header_syn=''
      write(lp,'(a)')'Reading Climatology Forcing (forcing option month)'
      lclim=.true.
   else
      write(lp,'(a)') 'Unrecognized forcing flag '//trim(rforce)
      print *, '(forfun_nersc)'
      call exit(1)
   end if


   ! This environment variable will only calculate synoptic wind forcing
   lflux=.true.
   call getenv('FORFUN_WIND',cenv)
   if (trim(cenv)=='yes') then
      print *,'wind flag detected - will only calc synoptic wind '
      lflux=.false.
   end if
   lsynslp   = lsynslp    .and. lflux
   lsynuwnd   = lsynuwnd    .and. .true.
   lsynvwnd   = lsynvwnd    .and. .true.
   lsynairtmp= lsynairtmp .and. lflux
   lsynprecip= lsynprecip .and. lflux
   lsynradflx= lsynradflx .and. lflux
   lsynshwflx= lsynshwflx .and. lflux
   lsynvapmix= lsynvapmix .and. lflux


   if (trim(clmflag)=='old') then
      forc_header_clm='"Old" NERSC climatology'
   elseif (trim(clmflag)=='era40') then
      forc_header_clm='ERA40-Based climatology'
   elseif (trim(clmflag)=='ncepr') then
      forc_header_clm='NCEP-Based climatology'
   else
      write(lp,'(a)') 'Unrecognized climate forcing flag '//trim(clmflag)
      print *, '(forfun_nersc)'
      call exit(1)
   end if
   write(lp,'(a)') 'Climatology forcing flag '//trim(clmflag)


   ! Allocates variables and initializes them
   call init_fvars_nersc()

   end subroutine


   ! Helper routine, opens forcing file for writing and returns unit
   subroutine opnfrcwr(iunit,cvar,synflg,cstr1,cstr2,prefix)
   use mod_za
   implicit none
   character(len=*), intent(in) :: cvar,cstr1,cstr2
   logical,          intent(in) :: synflg
   integer,          intent(in) :: iunit
   character(len=*), intent(in),optional :: prefix

   integer lgth
   character(len=80) :: forc_header,prfx
   prfx=''
   if(present(prefix)) prfx=trim(prefix)
   if (synflg) then
      forc_header=forc_header_syn
   else
      forc_header=forc_header_clm
   end if
   lgth = len_trim(flnmfor)
   !print *,iunit,prfx
   open (unit=iunit,file=flnmfor(1:lgth)//trim(prfx)//'forcing.'//cvar//'.b', &
         status='replace', action='write')
   write(iunit,'(a)') cstr1//trim(forc_header)
   write(iunit,'(a)') cstr2
   write(iunit,'(a)') ''
   write(iunit,'(a)') ''
   write(iunit,'(a,2i5)') 'i/jdm = ',idm,jdm
   call zaiopf(flnmfor(1:lgth)//trim(prfx)//'forcing.'//cvar//'.a', 'replace', iunit)
   end subroutine opnfrcwr


   ! Helper routine, Opens forcing file for reading and returns unit
   subroutine opnfrcrd(iunit,cvar)
   use mod_za
   implicit none
   character(len=*), intent(in) :: cvar
   integer,          intent(in) :: iunit
   integer lgth
   !print *,iunit
   lgth = len_trim(flnmfor)
   open (unit=iunit,file=flnmfor(1:lgth)//'forcing.'//cvar//'.b', &
         status='old', action='read')
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   read(iunit,*)
   call zaiopf(flnmfor(1:lgth)//'forcing.'//cvar//'.a', 'old', iunit)
   end subroutine opnfrcrd




! ---    Write-out of forcing field to .[ab] files
   subroutine writeforc(fld,iunit,cfld,lsynfld,dtime,forcesyn)
   use mod_grid
   use mod_za
   implicit none
   !real, dimension(idm,jdm), intent(in) :: fld
   real, dimension(idm,jdm), intent(inout) :: fld ! zaiowr has inout
   integer,          intent(in) :: iunit
   character(len=7), intent(in) :: cfld
   logical,          intent(in) :: lsynfld
   real*8 ,          intent(in) :: dtime
   logical,optional, intent(in) :: forcesyn

   logical :: forcesyn2
   character(len=5) :: charid
   real :: hmin,hmax

   forcesyn2=.false.
   if (present(forcesyn)) forcesyn2=.true.
   if (lsynfld) then
      charid='(syn)'
   else
      charid='(clm)'
   end if
   call zaiowr(fld,ip,.false.,hmin,hmax,iunit,.true.)
   if (lsynfld .or. forcesyn2) then
      write(iunit,'(a7,a5,":dtime1,range = ",2f12.4,2e14.6)') &
            cfld,charid,dtime,rdtime,hmin,hmax
   else
      write(iunit,'(a7,":month,range = ",i2.2,2e16.8)') &
            cfld,int(dtime),hmin,hmax
   end if
   call flush(iunit)
   end subroutine writeforc



   ! Routine for time interpolation using either synoptic fields 
   ! (lsynfld is true) or climatology fields
   subroutine interptime(fldout,synfld,clmfld, &
      tw0,tw1,lwslot,upslot,lsynfld)
   implicit none
   real, dimension(idm,jdm  ), intent(out) :: fldout
   real, dimension(idm,jdm  ), intent(in ) :: synfld
   real, dimension(idm,jdm,2), intent(in ) :: clmfld
   real   , intent(in) :: tw0,tw1
   integer, intent(in) :: lwslot,upslot
   logical, intent(in) :: lsynfld
   integer :: i,j

   ! Time interpolation of forcing. Two situations:
   ! 1) Synoptic forcing was always set for a field 
   ! 2) Synoptic forcing was never  set for a field
   ! Transition period is caught by case 1
   if (lsynfld) then
     !$OMP PARALLEL DO PRIVATE(i,j)SCHEDULE(STATIC,jblk)
     do j=1,jdm
     do i=1,idm
        fldout(i,j) = w_switch*(tw0*clmfld(i,j,lwslot) +  &
                                tw1*clmfld(i,j,upslot))   &
                    +(1.-w_switch)*synfld(i,j)
     end do
     end do
     !$OMP END PARALLEL DO
  else 
     !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(STATIC,jblk)
     do j=1,jdm
     do i=1,idm
        fldout(i,j) = tw0*clmfld(i,j,lwslot) +  &
                      tw1*clmfld(i,j,upslot)
     end do
     end do
     !$OMP END PARALLEL DO
  end if
  end subroutine



end module mod_forcing_nersc
