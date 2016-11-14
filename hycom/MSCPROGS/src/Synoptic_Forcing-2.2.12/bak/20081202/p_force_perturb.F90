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



      program force_perturb
      use mod_xc
      use mod_za
      use mod_grid
      use mod_forcing_nersc
      use mod_random_forcing
      implicit none

      integer,parameter :: ioff=30 ! offset for input  files
      integer,parameter :: ooff=0  ! offset for output files
      integer ::  ios,off,k
      real    :: dtime,span
      character(len=*), parameter ::  &
         catm='Atm. forcing', &
         cwnd='Wind forcing', &
         prfx='tst.' 

      call xcspmd()
      call zaiost()
      call get_grid()           ! inits various stuff (grid, etc) 
      call init_forcing_nersc() ! Inits flags, constants and allocatable fields
      call init_rand_update()   ! initializes random forcing and allocates vars

      if (.not.randf) then
         print *,'randf option switched off in infile2.in '
         call exit(0)
      end if

      ! Open forcing files for reading
      call opnfrcrd(901+ioff,'tauewd')
      call opnfrcrd(902+ioff,'taunwd')
      call opnfrcrd(903+ioff,'wndspd')
      call opnfrcrd(904+ioff,'airtmp')
      call opnfrcrd(906+ioff,'precip')
      call opnfrcrd(unit_uwind +ioff,'uwind')
      call opnfrcrd(unit_vwind +ioff,'vwind')
      call opnfrcrd(unit_relhum+ioff,'relhum')
      call opnfrcrd(unit_slp   +ioff,'slp')
      call opnfrcrd(unit_clouds+ioff,'clouds')

      ! Open forcing files for writing - treat all perturbations as synoptic
      call opnfrcwr(901+ooff,'tauewd',.true.,cwnd,'Stresses on u grid [N/m^2]'   ,prefix=prfx)
      call opnfrcwr(902+ooff,'taunwd',.true.,cwnd,'Stresses on u grid [N/m^2]'   ,prefix=prfx)
      call opnfrcwr(903+ooff,'wndspd',.true.,cwnd,'Wind speed [m/s]'             ,prefix=prfx)
      call opnfrcwr(904+ooff,'airtmp',.true.,catm,'Air temperature [C]'          ,prefix=prfx)
      call opnfrcwr(906+ooff,'precip',.true.,catm,'precipitation [m/s]'          ,prefix=prfx)
      call opnfrcwr(unit_uwind +ooff,'uwind', .true.,cwnd,'Wind u-dir [m/s]'     ,prefix=prfx)
      call opnfrcwr(unit_vwind +ooff,'vwind', .true.,cwnd,'Wind v-dir [m/s]'     ,prefix=prfx)
      call opnfrcwr(unit_clouds+ooff,'clouds',.true.,catm,'Clouds []'            ,prefix=prfx)
      call opnfrcwr(unit_relhum+ooff,'relhum',.true.,catm,'Relative humidity  []',prefix=prfx)
      call opnfrcwr(unit_slp   +ooff,'slp',   .true.,catm,'Slp [mBar]'           ,prefix=prfx)

      ios=0
      do while(ios==0)
         ! Read input fields into syn fields
         call readfrcitem(901        +ioff,syntaux  ,dtime,span,ios)
         call readfrcitem(902        +ioff,syntauy  ,dtime,span,ios)
         call readfrcitem(903        +ioff,synwndspd,dtime,span,ios)
         call readfrcitem(904        +ioff,synairtmp,dtime,span,ios)
         call readfrcitem(906        +ioff,synprecip,dtime,span,ios)
         call readfrcitem(unit_uwind +ioff,synuwind ,dtime,span,ios)
         call readfrcitem(unit_vwind +ioff,synvwind ,dtime,span,ios)
         call readfrcitem(unit_relhum+ioff,synrelhum,dtime,span,ios)
         call readfrcitem(unit_slp   +ioff,synslp   ,dtime,span,ios)
         call readfrcitem(unit_clouds+ioff,synclouds,dtime,span,ios)

         rdtime= span

         if (ios==0) then
            if ( mod(int(dtime/rdtime),20)==0) then
               write(lp,*)
               print *,'Random forcing update - dtime=',dtime
            else
               write(lp,'(a1)',advance='no')  '.'
               call flush(lp)
            end if

            call rand_update()

            ! Dump to new forcing fields
            call writeforc(syntaux  ,901+ooff,        ' tauewd',.true.,dtime)
            call writeforc(syntauy  ,902+ooff,        ' taunwd',.true.,dtime)
            call writeforc(synwndspd,903+ooff,        ' wndspd',.true.,dtime)
            call writeforc(synairtmp,904+ooff,        ' airtmp',.true.,dtime)
            call writeforc(synprecip,906+ooff,        ' precip',.true.,dtime)
            call writeforc(synuwind ,unit_uwind +ooff,'  uwind',.true.,dtime)
            call writeforc(synvwind ,unit_vwind +ooff,'  vwind',.true.,dtime)
            call writeforc(synrelhum,unit_relhum+ooff,' relhum',.true.,dtime)
            call writeforc(synslp   ,unit_slp   +ooff,'    slp',.true.,dtime)
            call writeforc(synclouds,unit_clouds+ooff,'  clouds',.true.,dtime)
         end if
      end do


      do k=1,2
         if (k==1) then
            off=ioff
         else
            off=ooff
         end if
         close(901+off) ; call zaiocl(901+off)
         close(902+off) ; call zaiocl(902+off)
         close(903+off) ; call zaiocl(903+off)
         close(904+off) ; call zaiocl(904+off)
         close(906+off) ; call zaiocl(906+off)
         close(unit_uwind +off) ; call zaiocl(unit_uwind +off)
         close(unit_vwind +off) ; call zaiocl(unit_vwind +off)
         close(unit_slp   +off) ; call zaiocl(unit_slp   +off)
         close(unit_relhum+off) ; call zaiocl(unit_relhum+off)
         close(unit_clouds+off) ; call zaiocl(unit_clouds+off)
      end do


      contains 
      subroutine readfrcitem(iunit,fld,dtime,span,ios)
      implicit none
      real   , intent(out) :: fld(idm,jdm),dtime,span
      integer, intent(out) :: ios
      integer, intent(in)  :: iunit

      real :: hmin,hminb,hmax,hmaxb
      character(len=80) :: cline
      integer :: i

      read(iunit,'(a80)',iostat=ios) cline

      if (ios/=0) return
      call zaiord(fld,ip,.false.,hmin,hmax,iunit)
      i = index(cline,'=')
      read (cline(i+1:),*) dtime,span, hminb,hmaxb
      if (abs(hmin-hminb).gt.abs(hminb)*1.e-4 .or. &
          abs(hmax-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
            'error - .a and .b files not consistent:', &
            '.a,.b min = ',hmin,hminb,hmin-hminb, &
            '.a,.b max = ',hmax,hmaxb,hmax-hmaxb
         print *,'unit=',iunit
         print *,'(p_force_perturb:readfrcitem)'
         call exit(1)
      end if

      if (span/=0.25) then
         print *,'unsupported span ',span
         call exit(1)
      end if
      end subroutine


      end program
