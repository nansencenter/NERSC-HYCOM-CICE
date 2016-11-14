! --- generate monthly forcing fields
! --- Now also High-frequency forcing


! --- When this routine is finished, The following fields should be 
! --- prepared with these variable names and units. Note that the fields
! --- are either read all at once (climatology), or sequentially in HYCOM
! --- (synoptic forcing).
! --- -------------------------------------------------------------
! --- Field description:     Variable name:     Unit:
! --- ------------------     --------------     -----              
! --- Atmosphere temp        airtmp             Celsius
! --- Relative humidity      relhum             [] (Fraction)
! --- Precipitation          precip             m/s (water equivalent)
! --- Total Cloud Cover      cloud              [] (fraction)
! --- Wind speed             wndspd             m/s
! --- U Wind component       uwnd               m/s (component on grid)
! --- V Wind component       uwnd               m/s (component on grid
! --- U Drag component       taux               N/m^2
! --- V Drag component       tauy               N/m^2
! --- SLP                    slp                mBar
! --- -------------------------------------------------------------



      program forfun_nersc
      use mod_xc
      use mod_za
      use mod_grid
      use mod_year_info
      use mod_forcing_nersc
      use mod_storm
      use m_read_clim
      use m_read_ecmwf
      use m_read_era40
      use m_read_ecmwf_nc
      use m_read_ncep
      use m_qsw0
      implicit none

      integer*4, external :: iargc

      logical, parameter :: fdiag=.false.
      integer, parameter :: ifdiag=12
      integer :: refyear
      real*8 :: time1,time2
      character(len=5) :: rforce

      real fac,w4,wfact,span
      logical ex
      integer i,j,k,id,mo,l,ip1,im1,jp1,jm1,lgth

      integer :: lstyear,fstyear,lstmonth,fstmonth
      integer :: iyear,cmonth,ihour,juld
      integer :: upjuld,lwjuld,upmonth,lwmonth,upyear,lwyear
      real    :: tw0,tw1

      real, dimension(:,:), allocatable ::   tmptile, gfld, cawdir,
     &   util1, radfl0
      type (year_info) :: rttmp,rtinit
      character(len=5) :: charid

      real*8 :: time0,rdtime,dtime, time0keep, timetmp
      integer, parameter :: hyc_base_year = 1901
      integer, parameter :: hyc_base_month= 1
      integer, parameter :: hyc_base_dom  = 1

      real :: vapp, shum
      character(len=80) :: syntstf,tmparg

      logical, parameter :: lprrelax=.false.

      integer, parameter :: nrunits=14
      integer :: iunit
      integer, dimension(nrunits) :: units=
     & (/901,901,902,903,904,905,906,907,908,unit_relhum,
     &   unit_uwind,unit_vwind,unit_clouds,unit_slp/)
      logical :: lopen, lopenall

      logical :: uwindex, vwindex, wndabex, airtpex, 
     &   humidex, precpex, cloudex, sssex, sstex

      ! Turns on diagnostic for synoptic forcing
      logical,parameter  :: syndiagon=.false.
      logical :: lsyntst
      real ::  w_switch=1.   !fanf: "weight of the linear combinaison after
      real ::  g_switch   !fanf: "gap before the switch"
      real ::  l_switch   !fanf: "length of the smooth transition"

      integer :: imargin

      !character(len=*), parameter :: clmflag='era40'
      !character(len=*), parameter :: clmflag='old'
      character(len=5) :: clmflag
      character(len=80) :: forc_header_syn, forc_header_clm,
     &   forc_header,cline, cenv
      character(len=200) :: path_era40_clim

      integer :: ticker

      integer, dimension(2) :: clmslot,slot_read, slot_place 
      integer :: lwslot, upslot
      integer :: islot, nslots_read, iplace

      real :: fversion,rvar
      real,parameter :: version30=3.0
      real,parameter :: version31=3.1
      integer nday1,nhour1
      integer nday2,nhour2
      real rday1
      real rday2
      real delt1
      integer yrflag
      real, parameter :: stefanb=5.67e-8
      real, parameter :: t0deg  =273.15
      real, parameter :: radian =57.2957795
      real, parameter :: pi     =3.1415926536

      real :: hmin,hmax,hminb,hmaxb

      logical :: bailout, frcinput, lflux
      character(len=1) :: flnmfor=' '
      character(len=6) :: cvarin
      character(len=20) :: tmpchar

      bailout=.false.



         


C ---------------------------------------------------------------------------------------------
      inquire(file='blkdat.input',exist=ex)
      if (.not.ex) then
         write(lp,*) 'ERROR: blkdat.input does not exist'
         print *, '(forfun_nersc:blkdat.input parse)'
         call exit(1)
      endif
      open(10,file='blkdat.input')
      cvarin=''
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*)
      do while (cvarin/='idm   ')
         read(10,*) rvar,cvarin
      end do
      if (cvarin/='idm   ') then
         write(lp,*) 'ERROR: Could not get idm'
         print *, '(forfun_nersc:blkdat.input parse)'
         call exit(1)
      else 
         idm=nint(rvar)
      end if

      do while (cvarin/='jdm   ')
         read(10,*) rvar,cvarin
      end do
      if (cvarin/='jdm   ') then
         write(lp,*) 'ERROR: Could not get jdm'
         print *, '(forfun_nersc:blkdat.input parse)'
         call exit(1)
      else 
         jdm=nint(rvar)
      end if

      do while (cvarin/='yrflag')
         read(10,*) rvar,cvarin
      end do
      if (cvarin/='yrflag') then
         write(lp,*) 'ERROR: Could not get year flag'
         print *, '(forfun_nersc:blkdat.input parse)'
         call exit(1)
      else 
         yrflag=nint(rvar)
      end if
      do while (cvarin/='baclin')
         read(10,*) delt1,cvarin
      end do
      if (cvarin/='baclin') then
         write(lp,*) 'ERROR: Could not get baclin'
         print *, '(forfun_nersc:blkdat.input parse)'
         call exit(1)
      end if
      close(10)
C ---------------------------------------------------------------------------------------------
      ! This environment variable will reset the yrflag if encountered -
      ! useful for ignoring yearflag when we just want to produce a
      ! climatology.
      call getenv('FORFUN_FORCE_CLIMATOLOGY',cenv)
      print *,cenv
      if (trim(cenv)=='yes') then
         print *,'Overriding yearflag from blkdat.input -- '//
     &           'will produce climatology '//trim(clmflag)
         yrflag=0
      end if






      ! Open infile.in - get start/end times
C ---------------------------------------------------------------------------------------------
      inquire(file='infile.in',exist=ex)
      if (.not.ex) then
         write(lp,*) 'ERROR: infile.in does not exist'
         print *, '(forfun_nersc:infile.in parse)'
         call exit(1)
      endif
      write(*,'(a)')'### infile.in'
      open(10,file='infile.in')
      read(10,*)fversion
      if (version30 /= fversion .and. version31 /= fversion) then
         write(lp,*) 'Wrong version of input file infile.in'
         write(lp,*) 'Expected version is: ',version30 ,' or ',version31
         if (fversion == 2.8) then
            write(lp,*) 'have now introduced:'
            write(lp,*) 'lslowtide'
            write(lp,*)  'day hour rather than day.dec for diag days'
         endif
         print *, '(forfun_nersc:infile.in parse)'
         call exit(1)
      endif
      read(10,*)
 

      ! parse infile.in - these are not used when choosing 
      ! climatology forcing (here yrflag <3 )
      read(10,*)refyear
      write(*,'(a,i4)') 'infile.in:  refyear = ',refyear
      read(10,*)nday1,nhour1
      rday1=float(nday1)+float(nhour1)/24.0
      write(*,'(a,i4,i3,a,f8.2)')
     &   'infile.in: start at = ',nday1,nhour1,' or ',rday1
      read(10,*)nday2,nhour2
      rday2=float(nday2)+float(nhour2)/24.0
      write(*,'(a,i4,i3,a,f8.2)') 
     &   'infile.in:  end at  = ',nday2,nhour2,' or ',rday2
      close(10)


      ! Transition flag - for now an argument to this routine
      time1=max(0.,nday1+nhour1/24.-.25)
      time2=max(0.,nday2+nhour2/24.+.25)
      if (iargc()==2) then
         call getarg(1,rforce)
         call getarg(2,clmflag)
         l_switch=1.
         g_switch=time2+10. ! Ie; switch will never happen...
      elseif (iargc()==4) then
         call getarg(1,rforce)
         call getarg(2,clmflag)
         call getarg(3,tmparg)
         read(tmparg,*) g_switch
         call getarg(4,tmparg)
         read(tmparg,*) l_switch
         !l_switch=5.
         !g_switch=17. ! GOM  setup
         !
         !l_switch=6.
         !g_switch=2. ! Test  setup
         print *
         print '(a)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         print '(a)','!!!!!!!!!  Using Climatology transition !!!!!!!!!'
         print '(a,f10.3,x,f10.3,a)',
     &           '!!!!!!!!! params:',g_switch,l_switch,
     &           '  !!!!!!!!!' 
         print '(a)','!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         print *
      else
         print *
         print '(a)','Routine creates forcing files which can be  '
         print '(a)','read by hycom models. The forcing files are '
         print '(a)','tuned towards the NERSC version of HYCOM    '
         print '(a)','You will need the file infile.in, grid and  '
         print '(a)','bathy files (regional.*), blkdat.input  and '
         print '(a)','raw forcing fields available in current     '
         print '(a)','directory (typically linked in).            '
         print *
         print '(a)','Usage :'
         print '(a)', '  forfun_nersc rforce clmflag '
         print '(a)', '  forfun_nersc rforce clmflag g_switch l_switch'
         print '(a)','rforce is forcing option, clmflag is climate flag'
         print '(a)','Example:                                         '
         print '(a)', '  forfun_nersc era40 era40 '
         print '(a)','Creates forcing using ERA40 synoptic and '
         print '(a)','climatological forcing                   '
         print '(a)', '  forfun_nersc ecmwf era40 '
         print '(a)','Creates forcing using ECMWF synoptic and '
         print '(a)','ERA40 climatological forcing. ERA40 climatology'
         print '(a)','will substitute for any missing ECMWF fields   '
         call exit(1)
      end if
      if(rforce=='storm') time1=max(0.,nday1+nhour1/24.)


      if (trim(clmflag)/='old' .and. trim(clmflag)/='era40') then
         ! Revert to "old"
         !clmflag='old'
         print *,'Climate flag either old or era40'
         print *,'Climate flag is ',clmflag
         call exit(1)
      end if




      ! Header flag used in forcing files, set flags on which files to
      ! read from synoptic fields, and which to base on climatology.
      ! This only applies if synoptic forcing is used
      forc_header_syn=''
      if (trim(rforce)=='ecmwf' .or. trim(rforce)=='metno') then
         forc_header_syn='ECMWF operational fields'
         lsynslp   =.true.
         lsynwnd   =.true.
         lsynrelhum=.true.
         lsynairtmp=.true.
         lsynvapmix=lsynairtmp .and. lsynrelhum
         write(lp,*)'Reading ECMWF Forcing'
      elseif (trim(rforce)=='era40') then
         forc_header_syn='ERA40 reanalyzed fields'
         lsynslp   =.true.
         lsynwnd   =.true.
         lsynrelhum=.true.
         lsynairtmp=.true.
         lsynprecip=.true.
         lsynssr   =.false.
         lsynclouds=.true.
         write(lp,*)'Reading ERA40 Forcing'
      elseif (trim(rforce)=='ecnc') then
         forc_header_syn='Ecmwf reanalyzed fields'
         lsynslp   =.true.
         lsynwnd   =.true.
         lsynrelhum=.true.
         lsynairtmp=.true.
         lsynprecip=.false. ! Precipitation not yet available
         lsynssr   =.false. ! ssr not yet available
         lsynclouds=.true.
         write(lp,*)'Reading Ecmwf Forcing (NetCDF)'
      elseif (trim(rforce)=='storm') then
         forc_header_syn='STORM reanalyzed fields'
         lsynslp   =.true.
         lsynwnd   =.true.
         lsynrelhum=.true.
         lsynairtmp=.true.
         lsynprecip=.true.
         lsynssr   =.true.
         lsynclouds=.true.
         write(lp,*)'Reading STORM Forcing'
      elseif (trim(rforce)=='ncepr') then
         forc_header_syn='NCEP/NCAR reanalyzed fields'
         lsynslp   =.true.
         lsynwnd   =.true.
         lsynrelhum=.true.
         lsynprecip=.true.
         lsynclouds=.true.
         lsynairtmp=.true.
         lsynvapmix=lsynairtmp .and. lsynrelhum
         write(lp,*)'Reading NCEP Forcing'
      elseif (trim(rforce)=='month') then
         forc_header_syn=''
         write(lp,*)'Reading Climatology Forcing '
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
      lsynwnd   = lsynwnd    .and. .true.
      lsynrelhum= lsynrelhum .and. lflux
      lsynairtmp= lsynairtmp .and. lflux
      lsynprecip= lsynprecip .and. lflux
      lsynssr   = lsynssr    .and. lflux
      lsynclouds= lsynclouds .and. lflux
      lsynvapmix= lsynvapmix .and. lflux



      if (trim(clmflag)=='old') then
         forc_header_clm='"Old" NERSC climatology'
      elseif (trim(clmflag)=='era40') then
         forc_header_clm='ERA40-Based climatology'
      else
         write(lp,'(a)') 'Unrecognized climate forcing flag '//
     &         trim(clmflag)
         print *, '(forfun_nersc)'
         call exit(1)
      end if

      write(lp,*) 'Reading Climatology Forcing '//trim(clmflag)





C ---------------------------------------------------------------------------------------------



      call xcspmd()
      call zaiost()
      call get_grid()

      allocate(clmuwind (idm,jdm,2))
      allocate(clmvwind (idm,jdm,2))
      allocate(clmwndspd(idm,jdm,2))
      allocate(clmtaux  (idm,jdm,2))
      allocate(clmtauy  (idm,jdm,2))
      allocate(clmvapmix(idm,jdm,2))
      allocate(clmairtmp(idm,jdm,2))
      allocate(clmrelhum(idm,jdm,2))
      allocate(clmprecip(idm,jdm,2))
      allocate(clmclouds(idm,jdm,2))
      allocate(clmradflx(idm,jdm,2))
      allocate(clmshwflx(idm,jdm,2))
      allocate(clmslp   (idm,jdm,2))
      allocate(clmsss   (idm,jdm,2))
      allocate(clmsst   (idm,jdm,2))

      allocate(gfld(idm,jdm))
      allocate(tmptile(idm,jdm))
      allocate(cawdir(idm,jdm))
      allocate(util1(idm,jdm))
      allocate(radfl0(idm,jdm))
      allocate(ip(idm,jdm))

      clmuwind=0.
      clmvwind=0.
      clmwndspd=0.
      clmtaux=0.
      clmtauy=0.
      clmvapmix=0.
      clmairtmp=0.
      clmrelhum=0.
      clmprecip=0.
      clmclouds=0.
      clmradflx=0.
      clmshwflx=0.
      clmsss=0.
      clmsst=0.









! --- ------------------------------------------------------------------------
! --- 
! ---   Open NERSC Climatology files -- Actual data & Header
! --- 
! --- ------------------------------------------------------------------------


      write(lp,'(a)') '#############################################'
      write(lp,'(a)') 'Reading climatological forcing fields forfun'


! --- ------------------------------------------------------------------------
! --- 
! ---   Open HYCOM forcing files -- Actual data & Header
! --- 
! --- ------------------------------------------------------------------------
      if (lsynwnd) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      lgth = len_trim(flnmfor)
      open (unit=901,file=flnmfor(1:lgth)//'forcing.tauewd.b',
     &      status='replace', action='write')
      write(901,'(a)') 'Wind Forcing, '//trim(forc_header)
      write(901,'(a)') 'Stresses on u grid [N/m^2]'
      write(901,'(a)') ''
      write(901,'(a)') ''
      write(901,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.tauewd.a', 'replace', 901)

      open (unit=902,file=flnmfor(1:lgth)//'forcing.taunwd.b',
     &      status='replace', action='write')
      write(902,'(a)') 'Wind Forcing, '//trim(forc_header)
      write(902,'(a)') 'Stresses on v grid [N/m^2]'
      write(902,'(a)') ''
      write(902,'(a)') ''
      write(902,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.taunwd.a', 'replace', 902)


      open (unit=903,file=flnmfor(1:lgth)//'forcing.wndspd.b',
     &      status='replace', action='write')
      write(903,'(a)') 'Wind Forcing, '//trim(forc_header)
      write(903,'(a)') 'Windspeed [m/s]'
      write(903,'(a)') ''
      write(903,'(a)') ''
      write(903,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.wndspd.a', 'replace', 903)


      if (lsynairtmp) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      open (unit=904,file=flnmfor(1:lgth)//'forcing.airtmp.b',
     &      status='replace', action='write')
      write(904,'(a)') 'Atmospheric Forcing, '//trim(forc_header)
      write(904,'(a)') 'Airtemperature [C]'
      write(904,'(a)') ''
      write(904,'(a)') ''
      write(904,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.airtmp.a', 'replace', 904)


      if (lsynvapmix) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      open (unit=905,file=flnmfor(1:lgth)//'forcing.vapmix.b',
     &      status='replace', action='write')
      write(905,'(a)') 'Atmospheric Forcing, '//trim(forc_header)
      write(905,'(a)') 'Vapmix [] (0-1)'
      write(905,'(a)') ''
      write(905,'(a)') ''
      write(905,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.vapmix.a', 'replace', 905)

      if (lsynprecip) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      open (unit=906,file=flnmfor(1:lgth)//'forcing.precip.b',
     &      status='replace', action='write')
      write(906,'(a)') 'Atmospheric Forcing, '//trim(forc_header)
      write(906,'(a)') 'Precipitation [m^3/m^2s = m/s]'
      write(906,'(a)') ''
      write(906,'(a)') ''
      write(906,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.precip.a', 'replace', 906)

      forc_header=forc_header_clm
      open (unit=907,file=flnmfor(1:lgth)//'forcing.radflx.b',
     &      status='replace', action='write')
      write(907,'(a)') 'Atmospheric Forcing, '//trim(forc_header)
      write(907,'(a)') 'Radiation flux [W/m^2]'
      write(907,'(a)') ''
      write(907,'(a)') ''
      write(907,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.radflx.a', 'replace', 907)

      if (lsynssr) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      open (unit=908,file=flnmfor(1:lgth)//'forcing.shwflx.b',
     &      status='replace', action='write')
      write(908,'(a)') 'Atmospheric Forcing, '//trim(forc_header)
      write(908,'(a)') 'Short wave flux [W/m^2]'
      write(908,'(a)') ''
      write(908,'(a)') ''
      write(908,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.shwflx.a', 'replace', 908)



C --  NERSC specific fields

      if (lsynwnd) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      open (unit=unit_uwind,file=flnmfor(1:lgth)//'forcing.uwind.b',
     &      status='replace', action='write')
      write(unit_uwind,'(a)') 'Wind Forcing, '//trim(forc_header)
      write(unit_uwind,'(a)') 'Wind p grid [m/s]'
      write(unit_uwind,'(a)') ''
      write(unit_uwind,'(a)') ''
      write(unit_uwind,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.uwind.a', 'replace',
     &unit_uwind)

      open (unit=unit_vwind,file=flnmfor(1:lgth)//'forcing.vwind.b',
     &      status='replace', action='write')
      write(unit_vwind,'(a)') 'Wind Forcing, '//trim(forc_header)
      write(unit_vwind,'(a)') 'Wind p grid [m/s]'
      write(unit_vwind,'(a)') ''
      write(unit_vwind,'(a)') ''
      write(unit_vwind,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.vwind.a', 'replace',
     &unit_vwind)

      if (lsynclouds) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      open (unit=unit_clouds,file=flnmfor(1:lgth)//'forcing.clouds.b',
     &       status='replace', action='write')
      write(unit_clouds,'(a)')'Atmospheric Forcing, '//trim(forc_header)
      write(unit_clouds,'(a)') 'Clouds [] (0-1)'
      write(unit_clouds,'(a)') ''
      write(unit_clouds,'(a)') ''
      write(unit_clouds,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.clouds.a', 'replace',
     &unit_clouds)

      if (lsynrelhum) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      open (unit=unit_relhum,file=flnmfor(1:lgth)//'forcing.relhum.b',
     &      status='replace', action='write')
      write(unit_relhum,'(a)')'Atmospheric Forcing, '//trim(forc_header)
      write(unit_relhum,'(a)') 'Relative humidity [] (0-1)'
      write(unit_relhum,'(a)') ''
      write(unit_relhum,'(a)') ''
      write(unit_relhum,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.relhum.a', 'replace',
     &unit_relhum)

      if (lsynslp) then
         forc_header=forc_header_syn
      else
         forc_header=forc_header_clm
      end if
      open (unit=unit_slp,file=flnmfor(1:lgth)//'forcing.slp.b',
     &      status='replace', action='write')
      write(unit_slp,'(a)') 'Atmospheric Forcing, '//trim(forc_header)
      write(unit_slp,'(a)') 'Slp [mBar]'
      write(unit_slp,'(a)') ''
      write(unit_slp,'(a)') ''
      write(unit_slp,'(a,2i5)') 'i/jdm = ',idm,jdm
      call zaiopf(flnmfor(1:lgth)//'forcing.slp.a', 'replace',
     &unit_slp)

C --  END NERSC specific fields


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Climatology forcing First time only !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (yrflag<3) then  ! Climatology forcing



         do mo=1,12
         id =0                   !=1 dump matrix

         if (clmflag=='old') then

            CALL read_clim(mo,clmuwind (:,:,1),'./Data/uwind.forc',
     &         id,clmflag,plon,plat,depths)       !m/s

            CALL read_clim(mo,clmvwind (:,:,1),'./Data/vwind.forc',
     &         id,clmflag,plon,plat,depths)       !m/s

            CALL read_clim(mo,clmwndspd(:,:,1),'./Data/wndab.forc',
     &         id,clmflag,plon,plat,depths)       !m/s

            CALL read_clim(mo,clmairtmp(:,:,1),'./Data/airtp.forc',
     &         id,clmflag,plon,plat,depths)       !kelvin

            CALL read_clim(mo,clmrelhum(:,:,1),'./Data/humid.forc',
     &         id,clmflag,plon,plat,depths)       !rh fraction 0-1

            CALL read_clim(mo,clmprecip(:,:,1),'./Data/precp.forc',
     &         id,clmflag,plon,plat,depths)       !mm/month precepitation

            CALL read_clim(mo,clmclouds(:,:,1),'./Data/cloud.forc',
     &         id,clmflag,plon,plat,depths)       !cc fraction 0-1

            call sss_sst_nodc(clmsst(:,:,1),clmsss(:,:,1),mo)

            write (lp,'(a,i3)')     'Read old clim for month ',mo

            ! Produce wind stress etc for the old nersc climatology...
            call nersc_clm_calc(1,mo,refyear)

         else if (clmflag=='era40') then

            ! Default path to ERA40 climatology
            path_era40_clim='./ERA40-Clim/'
            call getenv('NERSC_ERA40_CLIM',cenv)
            if (trim(cenv)/='') then
               path_era40_clim=trim(cenv)
            end if

            CALL read_clim(mo,clmuwind (:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_U10M_sfc.nc',id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmvwind (:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_V10M_sfc.nc',id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmwndspd(:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_wndspd.nc'  ,id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmairtmp(:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_T2M_sfc.nc' ,id,clmflag,
     &         plon,plat,depths)       !kelvin

            CALL read_clim(mo,clmrelhum(:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_relhum.nc'  ,id,clmflag,
     &         plon,plat,depths)       !rh fraction 0-1

            CALL read_clim(mo,clmprecip(:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_TP.nc'      ,id,clmflag,
     &         plon,plat,depths)       !mm/month precepitation

            CALL read_clim(mo,clmclouds(:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_TCC_sfc.nc' ,id,clmflag,
     &         plon,plat,depths)       !cc fraction 0-1

            ! Climatologies were calculated for these fields as well 
            CALL read_clim(mo,clmtaux (:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_taux.nc'    ,id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmtauy (:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_tauy.nc'    ,id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmslp (:,:,1),trim(path_era40_clim)//
     &         'era40_climatology_MSL_sfc.nc' ,id,clmflag,
     &         plon,plat,depths)       !m/s

            call sss_sst_nodc(clmsst(:,:,1),clmsss(:,:,1),mo)

            call era40_clm_calc(1,mo,refyear)

            write (lp,'(a,i3)') 
     &         'Read ERA40 clim for month ',mo

         else 
            write(lp,*) 'Unknown Clim flag '//clmflag
            print *, '(m_forfun_nersc)'
            call exit(1)
         end if




C ---    ----------------------------------------------
C ---    File IO to HYCOM-type files
C ---    ----------------------------------------------
         call zaiowr(clmtaux(:,:,1),ip,.false.,hmin,hmax,901,.true.)
         write(901,'("tau_ewd:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(901)

         call zaiowr(clmtauy(:,:,1),ip,.false.,hmin,hmax,902,.true.)
         write(902,'("tau_nwd:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(902)

         call zaiowr(clmwndspd(:,:,1),ip,.false.,hmin,hmax,903,.true.)
         write(903,'(" wndspd:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(903)

         call zaiowr(clmairtmp(:,:,1),ip,.false.,hmin,hmax,904,.true.)
         write(904,'(" airtmp:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(904)

         call zaiowr(clmvapmix(:,:,1),ip,.false.,hmin,hmax,905,.true.)
         write(905,'(" vapmix:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(905)

         call zaiowr(clmprecip(:,:,1),ip,.false.,hmin,hmax,906,.true.)
         write(906,'(" precip:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(906)


C --  NERSC specific fields

         call zaiowr(clmclouds(:,:,1),ip,.false.,hmin,hmax,
     &               unit_clouds,.true.)
         write(unit_clouds,'(" clouds:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(unit_clouds)

         call zaiowr(clmrelhum(:,:,1),ip,.false.,hmin,hmax,
     &               unit_relhum,.true.)
         write(unit_relhum,'(" relhum:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(unit_relhum)

         call zaiowr(clmuwind(:,:,1),ip,.false.,hmin,hmax,
     &               unit_uwind,.true.)
         write(unit_uwind,'("  uwind:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(unit_uwind)

         call zaiowr(clmvwind(:,:,1),ip,.false.,hmin,hmax,
     &               unit_vwind,.true.)
         write(unit_vwind,'("  vwind:month,range = ",i2.2,2e16.8)')
     &         mo,hmin,hmax
         call flush(unit_vwind)

         ! No slp climatology
         write(unit_slp,'("  slp  :month,range = ",i2.2,2e16.8)')
     &         mo,slp0,slp0
         call flush(unit_slp)



C --  END NERSC specific fields

         call zaiowr(clmshwflx(:,:,1),ip,.false.,hmin,hmax,908,.true.)
         write(908,'(" shwflx:month,range = ",i2.2,2e16.8)')
     &         mo, hmin, hmax
         call flush(908)


         ! Very simple par. of radiative flux
         call zaiowr(clmradflx,ip,.false.,hmin,hmax,907,.true.)
         write(907,'(" radflx:month,range = ",i2.2,2e16.8)')
     &         mo, hmin, hmax
         call flush(907)
C ---    ----------------------------------------------
C ---    End File IO to HYCOM-type files
C ---    ----------------------------------------------


         call clim_diag(mo,fdiag,ifdiag)

         end do ! mo
      else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!A)  High-frequency forcing         !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
      ! ECMWF/Met.no/NCEP files -- 6hr forcing...
      if (rforce=='ecmwf' .or. rforce=='ncepr'.or. 
     &    rforce=='metno' .or. rforce=='era40'.or.
     &    rforce=='storm' .or. rforce=='ecnc' ) then
         rdtime=6./24.
      else
         print *,'Unknown forcing ID:'//rforce
         print *,'If you want to create  climatological forcing,'
         print *,'make sure that yrflag=0 in blkdat.input'
         print *

         print *, '(forfun_nersc)'
         call exit(1)
      end if
      span=rdtime
      
      ! Lowest record needed
      ihour=floor((time1-floor(time1))*24)
      ihour=(ihour/6)*6
      time0keep=floor(time1) + float(ihour)/24. 

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

      clmslot=-1


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



      ! Start cycling records
      time0=time0keep
      call year_day(time0,refyear,rtinit,rforce)
      ticker=0
      do while (time0<time2+rdtime .and. .not. bailout)

         call year_day(time0,refyear,rttmp,rforce)
         ticker=ticker+1


         ! Read forcing fields
         !print *,time0,time0keep+g_switch
         if (time0<time0keep+g_switch) then

            if (rforce=='ecmwf' .or. rforce=='metno') then
                call read_ecmwf(rttmp,plon,plat,depths,rforce)
                !print *,'read_ecmwf skipped ...'
            elseif (rforce=='ncepr') then
               call read_ncep(rttmp,plon,plat,depths)
            elseif (rforce=='era40') then
               call read_era40(rttmp,plon,plat,depths)
            elseif (rforce=='ecnc ') then
               call read_ecmwf_nc(rttmp,plon,plat,depths)
            elseif (rforce=='storm') then
            call read_storm(rttmp,rtinit,plon,plat,depths,bailout)
            if (bailout) cycle
            else 
               print *,'Unknown forcing ID:'//rforce
               print *, '(forfun_nersc)'
               call exit(1)
            end if

         else

            !after g_switch days we switch smoothly during  
            !l_switch days to the clim : fanf

            ! KAL -- switch everything synoptic off
            lsynslp   =.false.
            lsynwnd   =.false.
            lsynrelhum=.false.
            lsynairtmp=.false.
            lsynprecip=.false.
            lsynssr   =.false.
            lsynclouds=.false.

            ! Switch weight 1 after l_switch period is over
Cdiag       print *,'max/min synairtmp:', 
Cdiag&         maxval(synairtmp),minval(synairtmp)
            w_switch=min(1.,(time0-time0keep-g_switch)/l_switch)
            write(lp,'(a,3i6,a,f6.2)')
     &         'rt= ',rttmp%iyy,rttmp%idd,rttmp%ihh,
     &         ' '//clmflag//' Climatology transition weight is:',
     &         w_switch

CKAL --     TODO - Change this unpredictable behaviour.....
C           add t0deg that it substact each time fror getting Celsius
C           instead of kelvin (here always on the
C           same variable ......(its really tricky) fanf
!$OMP PARALLEL DO PRIVATE(i,j,vapp,shum)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm               
               synairtmp(i,j)=synairtmp(i,j)+t0deg
            end do
            end do
!$OMP END PARALLEL DO
         endif


         ! Convert rel hum to vapmix
!$OMP PARALLEL DO PRIVATE(i,j,vapp,shum)
!$OMP&SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm

            ! From Relative humidity to mixing ratio
            vapp = satvappw(synairtmp(i,j))
            shum = synrelhum(i,j)*humid(1e5,vapp)
            !print *,vapp,shum
            synvapmix(i,j) = shum/(1.-shum)

            ! From Kelvin to Celsius (NB: NEW )
            synairtmp(i,j)= synairtmp(i,j)-t0deg
         end do
         end do
!$OMP END PARALLEL DO
         
         ! To "HYCOM Date"
         juld= datetojulian(rttmp%iyy,rttmp%imm,rttmp%idm+1, 
     &                      hyc_base_year,hyc_base_month,hyc_base_dom)
         dtime=juld+(rttmp%ihh+rttmp%iss/3600.)/24. 


         ! Linear interpolation of climatology fields
         lwmonth=rttmp%imm -1  + min(floor(rttmp%idm/14.),1)
         lwyear =rttmp%iyy+min(lwmonth-1,0)
         lwmonth=mod(lwmonth+12-1,12)+1
         lwjuld=datetojulian(lwyear,lwmonth,15,
     &                      hyc_base_year,hyc_base_month,hyc_base_dom)

         upmonth=lwmonth+1
         upyear = rttmp%iyy+(upmonth/12)
         upmonth=mod(upmonth-1,12)+1
         upjuld=datetojulian(upyear,upmonth,15,
     &                      hyc_base_year,hyc_base_month,hyc_base_dom)


         ! Check if climate forcing slots are filled with correct data
         ! Find missing month(s)
         lwslot=-1 ; upslot=-1
         do islot=1,2
            if(clmslot(islot)==lwmonth) lwslot=islot
            if(clmslot(islot)==upmonth) upslot=islot
         end do

         
         nslots_read=0
         if (lwslot==-1 .and. upslot == -1) then

            ! Both slots must be read
            nslots_read=2
            slot_read(1)=lwmonth
            slot_read(2)=upmonth
            lwslot=1
            upslot=2
            slot_place(1)=lwslot
            slot_place(2)=upslot
         else if (lwslot==-1) then

            ! forcing for lwmonth must be read
            nslots_read=1
            lwslot=mod(upslot,2)+1
            slot_read(1)=lwmonth
            slot_place(1)=lwslot
         else if (upslot==-1) then

            ! forcing for upmonth must be read 
            nslots_read=1
            upslot=mod(lwslot,2)+1
            slot_read(1)=upmonth
            slot_place(1)=upslot
         end if

         !if (mnproc==1) then
         !   print *,'lwslot,upslot',lwslot,upslot
         !   print *,'weights      ',tw0,tw1
         !   print *,'clmslot      ',clmslot
         !   print *
         !   print *,'nslot_read   ',nslots_read
         !   print *,'slot_read    ',slot_read
         !   print *,'slot_place   ',slot_place
         !end if

         do islot=1,nslots_read
            
            print *
            mo=slot_read(islot)
            iplace=slot_place(islot)
            print '(a,i4,a,i4)','Climatology slot:',iplace,
     &                         '   month: ',mo
            if (clmflag=='old') then

            CALL read_clim(mo,clmuwind (:,:,iplace),
     &         './Data/uwind.forc',id,clmflag,plon,plat,depths)       !m/s

            CALL read_clim(mo,clmvwind (:,:,iplace),
     &         './Data/vwind.forc',id,clmflag,plon,plat,depths)       !m/s

            CALL read_clim(mo,clmwndspd(:,:,iplace),
     &         './Data/wndab.forc',id,clmflag,plon,plat,depths)       !m/s

            CALL read_clim(mo,clmairtmp(:,:,iplace),
     &         './Data/airtp.forc',id,clmflag,plon,plat,depths)       !kelvin

            CALL read_clim(mo,clmrelhum(:,:,iplace),
     &         './Data/humid.forc',id,clmflag,plon,plat,depths)       !rh fraction 0-1

            CALL read_clim(mo,clmprecip(:,:,iplace),
     &         './Data/precp.forc',id,clmflag,plon,plat,depths)       !mm/month precepitation

            CALL read_clim(mo,clmclouds(:,:,iplace),
     &         './Data/cloud.forc',id,clmflag,plon,plat,depths)       !mm/month precepitation

            call sss_sst_nodc(clmsst(:,:,iplace),
     &         clmsss(:,:,iplace),mo)


            call nersc_clm_calc(iplace,mo,refyear)
            print *


            else if (clmflag=='era40') then

            CALL read_clim(mo,clmuwind (:,:,iplace),
     &         './ERA40-Clim/era40_climatology_U10M_sfc.nc',id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmvwind (:,:,iplace),
     &         './ERA40-Clim/era40_climatology_V10M_sfc.nc',id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmwndspd(:,:,iplace),
     &         './ERA40-Clim/era40_climatology_wndspd.nc'  ,id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmairtmp(:,:,iplace),
     &         './ERA40-Clim/era40_climatology_T2M_sfc.nc' ,id,clmflag,
     &         plon,plat,depths)       !kelvin

            CALL read_clim(mo,clmrelhum(:,:,iplace),
     &         './ERA40-Clim/era40_climatology_relhum.nc'  ,id,clmflag,
     &         plon,plat,depths)       !rh fraction 0-1

            CALL read_clim(mo,clmprecip(:,:,iplace),
     &         './ERA40-Clim/era40_climatology_TP.nc'      ,id,clmflag,
     &         plon,plat,depths)       !mm/month precepitation

            CALL read_clim(mo,clmclouds(:,:,iplace),
     &         './ERA40-Clim/era40_climatology_TCC_sfc.nc' ,id,clmflag,
     &         plon,plat,depths)       !cc fraction 0-1

            ! Climatologies were calculated for these fields as well 
            CALL read_clim(mo,clmtaux (:,:,iplace),
     &         './ERA40-Clim/era40_climatology_taux.nc'    ,id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmtauy (:,:,iplace),
     &         './ERA40-Clim/era40_climatology_tauy.nc'    ,id,clmflag,
     &         plon,plat,depths)       !m/s

            CALL read_clim(mo,clmslp (:,:,iplace),
     &         './ERA40-Clim/era40_climatology_MSL_sfc.nc' ,id,clmflag,
     &         plon,plat,depths)       !m/s

            call sss_sst_nodc(clmsst(:,:,iplace),clmsss(:,:,iplace),mo)

            call era40_clm_calc(iplace,mo,refyear)

            write (lp,'(a,i3)') 'Read ERA40 clim for month ',mo

            else 
               write(lp,*) 'Unknown Clim flag '//clmflag
               print *, '(m_forfun_nersc)'
               call exit(1)
            end if
            ! Mark slot filled
            clmslot(iplace)=mo

C ---    ----------------------------------------------
         end do
         if (nslots_read>0)  then
            print *
            do islot=1,2
               print '(a,i4,a,i4)','Climatology slot:',islot,
     &                             '   month: ',clmslot(islot)
            end do
            print *
         end if
Cdiag       print *,'max/min clmairtmp:', 
Cdiag&         maxval(clmairtmp),minval(clmairtmp)
Cdiag    call zaiopf ('tst.a','replace',11)
Cdiag    call zaiowr(clmairtmp(:,:,1),ip,.false.,hmin,hmax,11,.true.)
Cdiag    call zaiowr(clmairtmp(:,:,2),ip,.false.,hmin,hmax,11,.true.)
Cdiag    call zaiocl(11)







         tw1=(dtime-float(lwjuld))/float(upjuld-lwjuld)
         tw0=(1.-tw1)

         !lsyntst =mod(rttmp%idd,2)==0.and.rttmp%ihh==12.and.syndiagon
         lsyntst =mod(rttmp%idd,2)==0.and.syndiagon

         ! Zzanity check
         if (tw1<0. .or. tw1>1. ) then
            print *,'ups...tw0,tw1',tw0,tw1
            print *,'lwjuld,lwmonth,lwyear : ',lwjuld,lwmonth,lwyear
            print *,'dtime,juld  : ',dtime,juld
            print *,'upjuld,upmonth,upyear : ',upjuld,upmonth,upyear
            print *, '(forfun_nersc)'
            call exit(1)
         end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- tau_ewd (taux)
         if (lsynwnd) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               syntaux(i,j)=w_switch*(tw0*clmtaux(i,j,lwslot) + 
     &                                tw1*clmtaux(i,j,upslot))
     &                     +(1.-w_switch)*syntaux(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(syntaux,ip,.false.,hmin,hmax,901,.true.)
         write(901,'("tau_ewd",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &         charid,dtime,span,hmin,hmax
         call flush(901)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- tau_nwd (tauy)
         if (lsynwnd) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               syntauy(i,j)=w_switch*(tw0*clmtauy(i,j,lwslot) + 
     &                                tw1*clmtauy(i,j,upslot) )
     &                     +(1.-w_switch)*syntauy(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(syntauy,ip,.false.,hmin,hmax,902,.true.)
         write(902,'("tau_nwd",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &         charid,dtime,span,hmin,hmax
         call flush(902)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- wind speed
         if (lsynwnd) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synwndspd(i,j)=w_switch*(tw0*clmwndspd(i,j,lwslot) + 
     &                                  tw1*clmwndspd(i,j,upslot) )
     &                       +(1.-w_switch)*synwndspd(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(synwndspd,ip,.false.,hmin,hmax,903,.true.)
         if (mnproc==1) then
         write(903,'(" wndspd",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &         charid,dtime,span,hmin,hmax
         call flush(903)
         end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- air temp
         if (lsynairtmp) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synairtmp(i,j)=w_switch*(tw0*clmairtmp(i,j,lwslot) + 
     &                                  tw1*clmairtmp(i,j,upslot) )
     &                       +(1.-w_switch)*synairtmp(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(synairtmp,ip,.false.,hmin,hmax,904,.true.)
         write(904,'(" airtmp",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &         charid,dtime,span,hmin,hmax
         call flush(904)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- mixratio
         if (lsynvapmix) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synvapmix(i,j)=w_switch*(tw0*clmvapmix(i,j,lwslot) + 
     &                                  tw1*clmvapmix(i,j,upslot) )
     &                       +(1.-w_switch)*synvapmix(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(synvapmix,ip,.false.,hmin,hmax,905,.true.)
         write(905,'(" vapmix",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &         charid,dtime,span,hmin,hmax
         call flush(905)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- precipitation
         if (lsynprecip) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synprecip(i,j)=w_switch*(tw0*clmprecip(i,j,lwslot) + 
     &                                  tw1*clmprecip(i,j,upslot) )
     &                       +(1.-w_switch)*synprecip(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(synprecip,ip,.false.,hmin,hmax,906,.true.)
         write(906,'(" precip",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &         charid,dtime,span,hmin,hmax
         call flush(906)

!##########################
! --  NERSC specific fields
!##########################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- clouds
         if (lsynclouds) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synclouds(i,j)=w_switch*(tw0*clmclouds(i,j,lwslot) + 
     &                                  tw1*clmclouds(i,j,upslot) )
     &                       +(1.-w_switch)*synclouds(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            !if (mnproc==1) then
            !   print *,'max/min clouds'
            !   print *,maxval(clmclouds)
            !   print *,minval(clmclouds)
            !   print *,maxval(synclouds)
            !   print *,minval(synclouds)
            !end if
            charid="(clm)"
         end if
         call zaiowr(synclouds,ip,.false.,hmin,hmax,unit_clouds,.true.)
         write(unit_clouds,
     &      '(" clouds",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &      charid,dtime,span, hmin, hmax
         call flush(unit_clouds)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- relative humidity
         if (lsynrelhum) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synrelhum(i,j)=w_switch*(tw0*clmrelhum(i,j,lwslot) + 
     &                                  tw1*clmrelhum(i,j,upslot) )
     &                       +(1.-w_switch)*synrelhum(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(synrelhum,ip,.false.,hmin,hmax,unit_relhum,.true.)
         write(unit_relhum,
     &      '(" relhum",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &      charid,dtime,span, hmin, hmax
         call flush(unit_relhum)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- u wind component
         if (lsynwnd) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synuwind(i,j)=w_switch*(tw0*clmuwind(i,j,lwslot) + 
     &                                 tw1*clmuwind(i,j,upslot) )
     &                      +(1.-w_switch)*synuwind(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(synuwind,ip,.false.,hmin,hmax,unit_uwind,.true.)
         write(unit_uwind,
     &      '(" uwind ",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &      charid,dtime,span, hmin, hmax
         call flush(unit_uwind)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- v wind component
         if (lsynwnd) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synvwind(i,j)=w_switch*(tw0*clmvwind(i,j,lwslot)  +
     &                                 tw1*clmvwind(i,j,upslot) )
     &                      +(1.-w_switch)*synvwind(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(synvwind,ip,.false.,hmin,hmax,unit_vwind,.true.)
         write(unit_vwind,
     &      '(" vwind ",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &      charid,dtime,span, hmin, hmax
         call flush(unit_vwind)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- sea level pressure
         if (lsynslp) then
            charid='(syn)'
         else 
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               synslp(i,j)=w_switch*slp0+(1.-w_switch)*synslp(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid="(clm)"
         end if
         call zaiowr(synslp,ip,.false.,hmin,hmax,unit_slp,.true.)
         write(unit_slp,
     &      '(" slp   ",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &      charid,dtime,span, hmin, hmax
         call flush(unit_slp)

!##############################
! --  END NERSC specific fields
!##############################



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- shortwave flux
          if (.not.lsynssr) then
            ! Calc sw rad forcing
            if (lsynclouds) then
               tmptile=synclouds
            else
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
               do j=1,jdm
               do i=1,idm
                  tmptile(i,j)=tw0*clmclouds(i,j,lwslot) + 
     &                         tw1*clmclouds(i,j,upslot) 
               end do
               end do
!$OMP END PARALLEL DO
            end if
            timetmp=rttmp%idd+rttmp%ihh/24.+delt1/86400.
            call qsw0(radfl0,cawdir,tmptile,plat/radian,
     &                timetmp,rttmp%daysinyear)

!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               util1(i,j)=radfl0(i,j)*tmptile(i,j) 
               synshwflx(i,j)=util1(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid='(clm)'
         else ! Use shortwave flux from atm model (NCEPR / ERA40) 
            charid='(syn)'
         end if
         call zaiowr(synshwflx,ip,.false.,hmin,hmax,908,.true.)
         write(908,'(" shwflx",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &         charid,dtime,span, hmin, hmax
         call flush(908)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Synoptic forcing -- Radiative flux
         ! Very simple par.
         if (lsynairtmp) then
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               tmptile(i,j)= (synairtmp(i,j)+t0deg)**4 -
     &                    (tw0*clmsst(i,j,lwslot) + 
     &                     tw1*clmsst(i,j,upslot) +t0deg)**4
               tmptile(i,j)=stefanb*0.95*tmptile(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid='(mix)'
         else
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
            do j=1,jdm
            do i=1,idm
               tmptile(i,j)=(tw0*clmairtmp(i,j,lwslot) + 
     &                   tw1*clmairtmp(i,j,upslot) +t0deg)**4 -
     &                  (tw0*clmsst(i,j,lwslot) + 
     &                   tw1*clmsst(i,j,upslot) +t0deg)**4
               tmptile(i,j)=stefanb*0.95*tmptile(i,j)
            end do
            end do
!$OMP END PARALLEL DO
            charid='(clm)'
         end if
         !print *,'rad flux unset ...',mnproc
         tmptile=0.

!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm
            util1(i,j)=tmptile(i,j)+synshwflx(i,j)
         end do
         end do
!$OMP END PARALLEL DO

         !Quick fix
         where(depths<.1) util1=0.
         call zaiowr(util1,ip,.false.,hmin,hmax,907,.true.)
         write(907,'(" radflx",a5,":dtime1,range = ",2f12.4,2e14.6)')
     &         charid,dtime,span, hmin, hmax
         call flush(907)


         time0=time0+rdtime
         !print *,time0,dtime


      end do
      deallocate(synslp   )
      deallocate(synuwind )
      deallocate(synvwind )
      deallocate(synwndspd)
      deallocate(syntaux  )
      deallocate(syntauy  )
      deallocate(synvapmix)
      deallocate(synairtmp)
      deallocate(synrelhum)
      deallocate(synprecip)
      deallocate(synclouds)
      deallocate(synradflx)
      deallocate(synshwflx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!B) Climatology and  High-frequency forcing         !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#endif
         
      end if

      deallocate(clmuwind )
      deallocate(clmvwind )
      deallocate(clmwndspd)
      deallocate(clmtaux  )
      deallocate(clmtauy  )
      deallocate(clmvapmix)
      deallocate(clmairtmp)
      deallocate(clmrelhum)
      deallocate(clmprecip)
      deallocate(clmclouds)
      deallocate(clmradflx)
      deallocate(clmshwflx)
      deallocate(clmslp   )
      deallocate(clmsss   )
      deallocate(clmsst   )



      


! --- ------------------------------------------------------------------------
! --- 
! ---   Close  HYCOM forcing files -- Actual data & Header
! --- 
! --- ------------------------------------------------------------------------
      call zaiocl(901)
      call zaiocl(902)
      call zaiocl(903)
      call zaiocl(904)
      call zaiocl(905)
      call zaiocl(906)
      call zaiocl(907)
      call zaiocl(908)
      call zaiocl(unit_relhum)
      call zaiocl(unit_clouds)
      call zaiocl(unit_uwind)
      call zaiocl(unit_vwind)
      call zaiocl(unit_slp)

      if (mnproc==1) then
         close(901)
         close(902)
         close(903)
         close(904)
         close(905)
         close(906)
         close(907)
         close(908)
         close(unit_relhum)
         close(unit_clouds)
         close(unit_uwind)
         close(unit_vwind)
         close(unit_slp)
      end if


  99  FORMAT(30I4) 
 100  FORMAT(10(1x,e10.4)) 
      write(lp, '(a)') '###############################################'
      write(lp, '(a)') 'Done producing '//rforce//' Forcing fields'
      write(lp, '(a)') '###############################################'
      call exit(0)

      contains
      



C #######################################################################
C #######################################################################
C #######################################################################
C #######################################################################
C #######################################################################
C  END MAIN ROUTINE -- AUXILLARY ROUTINES FOLLOW
C #######################################################################
C #######################################################################
C #######################################################################
C #######################################################################
C #######################################################################






      subroutine nersc_clm_calc(islot,month,refyear)
      use mod_xc
      use mod_forcing_nersc
      use mod_year_info
      use m_qsw0
      implicit none
      integer, intent(in) :: islot,refyear,month

      integer :: imargin,jm1,jp1,im1,ip1,i,j
      real :: wfact,vapp,shum,w4, fac
      real*8 :: dtime
      type(year_info) :: rttmp



C ---    ----------------------------------------------
C ---    Calculate properties dependent on what we read
C ---    ----------------------------------------------
         ! Calc sw rad forcing
         call year_day(real(0.,kind=8),0,rttmp,'month')
         dtime=float(rttmp%daysinyear)*(month/12.-1./24.)
         call year_day(dtime,refyear,rttmp,'month')
         dtime=rttmp%idd+rttmp%ihh/24.+delt1/86400.
         call qsw0(clmshwflx(:,:,islot),cawdir,clmclouds(:,:,islot),
     &      plat/radian,dtime,rttmp%daysinyear)

! --- ------------------------------------------------------------------------
! --- compute air-water wind stresses  tau_x  and  tau_y  from the surface wind 
! --- components. surface wind components have dim [m/s], and are defined on a 
! --- c-grid;  airdns  has dim [kg/m^3], and  cd  is dimensionless. 
! ---  tau_x  and  tau_y  have dim [kg/(m s^2) = N/m^2]
! --- ------------------------------------------------------------------------

         fac=airdns*cd
         imargin=0
         ! KAL: This calculates monthly average stress from monthly average
         !      winds, which is .. kinda wrong..
!$OMP PARALLEL DO PRIVATE(j,l,w4,im1,ip1,jm1,jp1,wfact,vapp,shum)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm


            im1=max(1,i-1)
            jm1=max(1,j-1)
            ip1=min(i+1,idm)
            jp1=min(j+1,jdm)

            !print *,i,j,idm,jdm


            !windstress in u-dir
            w4=.25*(clmvwind(im1,jp1,islot)+clmvwind(i,jp1,islot)+ 
     &              clmvwind(im1,j  ,islot)+clmvwind(i,j  ,islot))     
C           print '(10i6,2f10.2)',mproc,nproc,i,j,im1,jm1,ip1,jp1,ii,jj,
C    &         clmuwind(i,j,islot),w4
            wfact=sqrt(clmuwind(i,j,islot)*clmuwind(i,j,islot)+w4*w4)
     &             *fac 
            clmtaux(i,j,islot)=clmuwind(i,j,islot)*wfact               

            !windstress in v-dir  
            w4=.25*(clmuwind(i,jm1,islot)+clmuwind(ip1,jm1,islot) 
     &             +clmuwind(ip1,j,islot)+clmuwind(i,j,islot    ))
            !call xcsync(flush_lp)
C           print '(11i5,4f15.2)',
C    &         mnproc,i0,j0,ii,jj,i,j,im1,jm1,ip1,jp1,
C    &         clmuwind(i,jm1,islot),clmuwind(ip1,jm1,islot),
C    &         clmvwind(ip1,j,islot),clmvwind(i,j,islot)
C           print '(6i6,2f10.2)',mproc,nproc,i,j,im1,jm1,
C    &         clmvwind(i,j,islot),w4
            wfact=sqrt(clmvwind(i,j,islot)*clmvwind(i,j,islot)+w4*w4)
     &            *fac
            clmtauy(i,j,islot)=clmvwind(i,j,islot)*wfact               

            !from mm/month to m/month 
            clmprecip(i,j,islot)=clmprecip(i,j,islot)*1.e-3            

            !from m/month to m/s
            clmprecip(i,j,islot)=clmprecip(i,j,islot)/(30*86400)

            ! From Relative humidity to mixing ratio
            vapp = satvappw(clmairtmp(i,j,islot))
            !print *,mnproc,i,j,islot,vapp,clmairtmp(i,j,islot)
            shum = clmrelhum(i,j,islot)*humid(1e5,vapp)
            clmvapmix(i,j,islot) = shum/(1.-shum)

            ! From Kelvin to Celsius (NB: NEW )
            clmairtmp(i,j,islot)= clmairtmp(i,j,islot)-t0deg

            ! Calculate Radiative flux
            clmradflx(i,j,islot)= (clmairtmp(i,j,islot)+t0deg)**4    -
     &                   (clmsst   (i,j,islot)+t0deg)**4
            clmradflx(i,j,islot)=clmradflx(i,j,islot)*stefanb*0.95
            clmradflx(i,j,islot)=clmradflx(i,j,islot)
     &                          +clmshwflx(i,j,islot)
         enddo
         enddo
!$OMP END PARALLEL DO


      end subroutine nersc_clm_calc



      subroutine era40_clm_calc(islot,month,refyear)
      use mod_xc
      use mod_forcing_nersc
      use mod_year_info
      use m_qsw0
      use m_rotate
      implicit none
      integer, intent(in) :: islot,refyear,month

      integer :: imargin,jm1,jp1,im1,ip1,i,j
      real :: wfact,vapp,shum,w4, fac
      real*8 :: dtime
      type(year_info) :: rttmp



C ---    ----------------------------------------------
C ---    Calculate properties dependent on what we read
C ---    ----------------------------------------------
         ! Calc sw rad forcing
         call year_day(real(0.,kind=8),0,rttmp,'month')
         dtime=float(rttmp%daysinyear)*(month/12.-1./24.)
         call year_day(dtime,refyear,rttmp,'month')
         dtime=rttmp%idd+rttmp%ihh/24.+delt1/86400.
         call qsw0(clmshwflx(:,:,islot),cawdir,clmclouds(:,:,islot),
     &      plat/radian,dtime,rttmp%daysinyear)


         fac=airdns*cd
         ! KAL: This calculates monthly average stress from monthly average
         !      winds, which is .. kinda wrong..
!$OMP PARALLEL DO PRIVATE(j,l,w4,im1,ip1,jm1,jp1,wfact,vapp,shum)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm


            ! From Relative humidity to mixing ratio
            vapp = satvappw(clmairtmp(i,j,islot))
            shum = clmrelhum(i,j,islot)*humid(1e5,vapp)
            clmvapmix(i,j,islot) = shum/(1.-shum)

            ! From Kelvin to Celsius (NB: NEW )
            clmairtmp(i,j,islot)= clmairtmp(i,j,islot)-t0deg

            ! ERA40 Precipitation is accumulated over 6hrs
            clmprecip(i,j,islot)= clmprecip(i,j,islot)/(6*3600.)

            ! Avoid neg precip (very small, scaling errors?)
            clmprecip(i,j,islot)=max(0.,clmprecip(i,j,islot))


            ! ERA40 SLP is in Pa
            clmslp(i,j,islot)= clmslp(i,j,islot)*0.01

            ! Calculate Radiative flux
            clmradflx(i,j,islot)= (clmairtmp(i,j,islot)+t0deg)**4    -
     &                   (clmsst   (i,j,islot)+t0deg)**4
            clmradflx(i,j,islot)=clmradflx(i,j,islot)*stefanb*0.95
            clmradflx(i,j,islot)=clmradflx(i,j,islot)
     &                          +clmshwflx(i,j,islot)
         enddo
         enddo
!$OMP END PARALLEL DO

         ! Rotate velocity model grid
         call rotate(clmuwind,
     &               clmvwind,
     &               plat    ,
     &               plon    ,
     &               idm,jdm,'l2m')


         ! Rotate stresses to model grid
         call rotate(clmtaux, 
     &               clmtauy, 
     &               plat   , 
     &               plon   , 
     &               idm,jdm,'l2m')

      end subroutine era40_clm_calc




C #######################################################################
C Saturation vapor påressure over ice                                    
C #######################################################################
      real function satvappi(tx)
      implicit none
      real, intent(in) :: tx
      satvappi=611.*10.**(9.5*(tx-273.16)/(tx-7.66)) ! tx in K
      end function satvappi




C #######################################################################
C Saturation vapor påressure over water
C #######################################################################
      real function satvappw(tx)
      implicit none
      real, intent(in) :: tx
      satvappw=611.*10.**(7.5*(tx-273.16)/(tx-7.50)) ! tx in K
      end function satvappw




C #######################################################################
C Specific humidity
C #######################################################################
      real function humid(px,es)
      real,intent(in) :: px,es
      humid=.622*es/(px-.378*es)
      end function humid






C #######################################################################
C Diagnostic for climatology forcing (tecplot file)
C #######################################################################
      subroutine clim_diag(mo,fdiag,ifdiag)
      use mod_xc
      use mod_forcing_nersc
      implicit none
      integer, intent(in) :: mo,ifdiag
      logical, intent(in) ::  fdiag

      integer :: i,j


C ---    --------------------------------------
C ---    Diagnostic ouput to tecforce.dat
C ---    --------------------------------------
         if (fdiag)THEN
            if (mo==1) then
               print *,'Diagnostics to tecforce.dat'
               OPEN(10,FILE='tecforce.dat',STATUS='replace')
               WRITE(10,'(''TITLE= "Forcing fields"'')')
               write(10,'(a)')'VARIABLES="i" "j" "lon" "lat" "depths[m]"
     &             "wndspd[m/s]" "airtmp[C]" "clouds[0-1]" "relhum[0-1]"
     &             "precip[mm/d]" "taux[N/m^2]" "tauy[N/m^2]" "sss[psu]"
     &             "sst[C]"'
               WRITE(10,'(''ZONE I='',I3,'', J='',I3,'', F=BLOCK'')')
     &                  idm,jdm
               WRITE(10,99)((i,i=1,idm),j=1,jdm)
               WRITE(10,99)((j,i=1,idm),j=1,jdm)
               WRITE(10,100)((plon(i,j),i=1,idm),j=1,jdm)
               WRITE(10,100)((plat(i,j),i=1,idm),j=1,jdm)
               WRITE(10,100)((depths(i,j),i=1,idm),j=1,jdm)
            else if (mo>1) then
               OPEN(10,FILE='tecforce.dat',STATUS='old',
     &               position='append')
            end if

            if (mo <= ifdiag) then

               if (mo > 1) then
                  WRITE(10,'(''ZONE I='',I3,'', J='',I3,'', F=BLOCK'')')
     &                  idm,jdm
                  write(10,'(a)')'D=(1,2,3,4,5)'
               endif

               WRITE(10,100)((clmwndspd(i,j,1),i=1,idm),j=1,jdm)
               WRITE(10,100)((clmairtmp(i,j,1),i=1,idm),j=1,jdm)
               WRITE(10,100)((clmclouds(i,j,1),i=1,idm),j=1,jdm)
               WRITE(10,100)((clmrelhum(i,j,1),i=1,idm),j=1,jdm)
               WRITE(10,100)((clmprecip(i,j,1)*86400000.,
     &            i=1,idm),j=1,jdm)
               WRITE(10,100)((clmtaux(i,j,1),i=1,idm),j=1,jdm)
               WRITE(10,100)((clmtauy(i,j,1),i=1,idm),j=1,jdm)
               WRITE(10,100)((clmsss(i,j,1),i=1,idm),j=1,jdm)
               WRITE(10,100)((clmsst(i,j,1),i=1,idm),j=1,jdm)
            ENDif
            CLOSE(10)
         end if
C ---    --------------------------------------
C ---    End Diagnostic ouput to tecforce.dat
C ---    --------------------------------------
  99  FORMAT(30I4) 
 100  FORMAT(10(1x,e10.4)) 
      end subroutine clim_diag








C #######################################################################
C Diagnostic for synoptic forcing (tecplot file)
C #######################################################################
      subroutine syn_diag(lsyntst,rttmp)
      use mod_xc
      use mod_forcing_nersc
      use mod_year_info
      implicit none
      logical, intent(in) :: lsyntst
      type(year_info) :: rttmp

      character(len=80) :: syntstf
      integer :: i,j


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnose synoptic forcing
      IF(lsyntst)THEN
         syntstf='synforce_y'//rttmp%cyy//'d'//rttmp%cdd//'h'
     &   //rttmp%chh//'.tec'
         print *,'Diagnostics to '//trim(syntstf)
         OPEN(10,FILE=trim(syntstf),STATUS='UNKNOWN')
         WRITE(10,'(''TITLE= "Forcing fields"'')')
         write(10,'(a)')
     &   'VARIABLES="i" "j" "lon" "lat" "depths[m]" "wndspd[m/s]" 
     &   "airtmp[C]" "clouds[0-1]" "relhum[0-1]" "vapmix[0-1]" 
     &   "precip[mm/d]" "taux[N/m^2]" "tauy[N/m^2]"
     &   "shwflx[W/m^2]" "radfl[W/m^2]x" "slp[mBar]"'
         WRITE(10,'(''ZONE I='',I3,'', J='',I3,'', F=BLOCK'')')
     &   idm,jdm
         WRITE(10,99)((i,i=1,idm),j=1,jdm)
         WRITE(10,99)((j,i=1,idm),j=1,jdm)
         WRITE(10,100)((plon     (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((plat     (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((depths   (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synwndspd(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synairtmp(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synclouds(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synrelhum(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synvapmix(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synprecip(i,j)*86400.*1000.,i=1,idm),j=1,jdm)
         WRITE(10,100)((syntaux  (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((syntauy  (i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synshwflx(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synradflx(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((synslp   (i,j),i=1,idm),j=1,jdm)

          CLOSE(10)
      end if
  99  FORMAT(30I4) 
 100  FORMAT(10(1x,e10.4)) 
      end subroutine syn_diag









C #######################################################################
C Reading - manipulating sss_nodc and sst_nodc files 
C #######################################################################

      subroutine sss_sst_nodc(sstfld,sssfld,month)
      use mod_xc
      use mod_za
      implicit none
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: sstfld,sssfld
      integer, intent(in) :: month

      logical :: sssex,sstex
      integer :: i,j,k
      real    :: hmin, hmax, hmin2, hmax2
      real, dimension(idm,jdm) :: gfld
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: lfld


      Print '(a,i2)','Reading sst/sss from old climatology - '//
     &            'used for radflx estimate. Month=',month

      ! Look for existing style files
      inquire(file='./Data/sst_nodc.dat',exist=sstex)
      inquire(file='./Data/sss_nodc.dat',exist=sssex)

      if (.not.sstex) then
         if (mnproc==1) then
            write(lp,*) 'Can not find file sst_nodc.dat'
            call flush(lp)
         end if
         print *, 'sss_sst_nodc)'
         call exit(1)
      end if

      if (.not.sssex) then
         if (mnproc==1) then
            write(lp,*) 'Can not find file sss_nodc.dat'
            call flush(lp)
         end if
         print *, '(sss_sst_nodc)'
         call exit(1)
      end if

!KAL - Dropped generation of .[ab] files

      open(14,file='./Data/sst_nodc.dat',form='formatted')
      do k=1,12
         read(14,'(10f9.4)')gfld
         if (k==month) then
            sstfld=lfld
            exit
         end if
      end do
      close(14)

      open(14,file='./Data/sss_nodc.dat',form='formatted')
      do k=1,12
         read(14,'(10f9.4)')gfld
         if (k==month) then
            sssfld=lfld
            exit
         end if
      end do
      close(14)

      ! Strange arctic feature....
C$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&SCHEDULE(STATIC,jblk)
      do j=1,jdm
      do i=1,idm
         if (plat(i,j) > 71.0)
     &      sssfld(i,j)=min(35.0,sssfld(i,j))
      enddo
      enddo
!$OMP END PARALLEL DO

      end subroutine sss_sst_nodc
C ------- ------------ ----------- -------------- ------------ ---
C ------  Change notes:
C ------- ------------ ----------- -------------- ------------ ---
c --- Cleaned up main routine - created many subroutines            KAL -- 06022005
c --- Got rid of global arrays for climatology (reduces mem usage)  KAL -- 06022005
C --- Introduced preprocessing of clim. - placed in ".a" files      KAL -- 06022005
C --- Clim files now placed in .ab-files -- more robust             KAL -- 22022005
C ------- ------------ ----------- -------------- ------------ ---
      end program forfun_nersc
