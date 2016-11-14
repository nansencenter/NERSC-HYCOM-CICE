C --- generate monthly forcing fields
C --- Now also High-frequency forcing


C --- When this routine is finished, The following fields should be 
C --- prepared with these variable names and units. Note that the fields
C --- are either read all at once (climatology), or sequentially in HYCOM
C --- (synoptic forcing).
C --- -------------------------------------------------------------
C --- Field description:     Variable name:     Unit:
C --- ------------------     --------------     -----              
C --- Atmosphere temp        airtmp             Celsius
C --- Relative humidity      relhum             [] (Fraction)
C --- Vapor mixing ratio     vapmix             [] (fraction)
C --- Precipitation          precip             m/s (water equivalent)
C --- Total Cloud Cover      cloud              [] (fraction)
C --- Wind speed             wndspd             m/s
C --- U Wind component       uwnd               m/s (component on grid)
C --- V Wind component       uwnd               m/s (component on grid
C --- U Drag component       taux               N/m^2
C --- V Drag component       tauy               N/m^2
C --- SLP                    slp                mBar
C --- -------------------------------------------------------------

C --- TODO uwind/vwind/taux/tauy on p-grid


      program forfun_nersc
      use mod_xc
      use mod_za
      use mod_grid
      use mod_parameters
      use mod_year_info22
      use mod_atm_func
      use mod_forcing_nersc
      use mod_storm
      use m_read_ecmwf
      use mod_read_era40
      use mod_read_erai
      use m_read_ecmwf_nc
      use m_read_ncep
      use m_qsw0
      use mod_clim_atm
      implicit none
#if defined (IARGC)
      integer*4, external :: iargc
#endif
      logical,parameter  :: syndiagon=.false.
      logical, parameter :: fdiag=.true.
      integer, parameter :: ifdiag=12
      integer, parameter :: hyc_base_year = 1901
      integer, parameter :: hyc_base_month= 1
      integer, parameter :: hyc_base_dom  = 1
      logical, save      :: lfirstsyn=.true.
      logical, parameter :: timer=.false.

      integer i,j,k,id,mo,l,ip1,im1,jp1,jm1,lgth
      integer :: ihour,juld
      integer :: upmonth,lwmonth,upyear,lwyear
      real    :: tw0,tw1
      real, dimension(:,:), allocatable ::   tmpcloud,  
     &   radfl0,tmplw, tmpradfl, tmpairtmp, fldout, tmpshwflx
      type (year_info) :: rttmp,rtinit,rttmp2
      real*8 :: dtime, timetmp, dtimelw, dtimeup
      character(len=80) :: syntstf,tmparg, forc_header
      integer :: iunit
      logical :: ex logical, lsyntst, lupdate
      integer :: imargin
      character(len=80) :: cline, cenv
      integer, dimension(2) :: clmslot
      integer :: lwslot, upslot, islot
      real :: hmin,hmax,hminb,hmaxb
      real*8 :: rmo
      logical :: bailout, frcinput
      character(len=6) :: cvarin
      character(len=20) :: tmpchar
      real*8 :: tmr0,tmr1
      real*8, external :: wtime

      bailout=.false.


      call xcspmd()
      call zaiost()
      call get_grid()
      call init_forcing_nersc() ! Inits flags, constants and allocatable fields
      call init_clmvars_nersc() ! Inits Atmospheric climatology variables

      allocate(tmpcloud  (idm,jdm))
      allocate(tmplw     (idm,jdm))
      allocate(tmpshwflx (idm,jdm))
      allocate(fldout    (idm,jdm))
      allocate(tmpairtmp (idm,jdm))
      allocate(tmpradfl  (idm,jdm))
      allocate(radfl0    (idm,jdm))


! --- ------------------------------------------------------------------------
! --- Open HYCOM forcing files -- Actual data & Header
! --- This is done for both the climate and synoptic forcing case
! --- ------------------------------------------------------------------------
      write(lp,'(a)') ' '
      write(lp,'(a)') ' '
      write(lp,'(a)') 'Writing climatological forcing fields'

      call opnfrcwr(901,'tauewd',lsynwnd,'Wind forcing',
     &   'Stresses on u grid [N/m^2]')

      call opnfrcwr(902,'taunwd',lsynwnd,'Wind forcing',
     &   'Stresses on u grid [N/m^2]')

      call opnfrcwr(903,'wndspd',lsynwnd,'Wind forcing',
     &   'Wind speed [m/s]')

      call opnfrcwr(904,'airtmp',lsynairtmp,'Atm. forcing',
     &   'Air temperature [C]')

      call opnfrcwr(905,'vapmix',lsynvapmix,'Atm. forcing',
     &   'Vapmix []')

      call opnfrcwr(906,'precip',lsynprecip,'Atm. forcing',
     &   'precipitation [m/s]')

      call opnfrcwr(907,'radflx',.false.,'Atm. forcing',
     &   'Radiative flux [W/m^2]')

      call opnfrcwr(908,'shwflx',lsynshwflx,'Atm. forcing',
     &   'Shortwave flux [W/m^2]')

      call opnfrcwr(unit_uwind,'uwind',lsynwnd,'Wind forcing',
     &   'Wind u-dir [m/s]')

      call opnfrcwr(unit_vwind,'vwind',lsynwnd,'Wind forcing',
     &   'Wind v-dir [m/s]')

      call opnfrcwr(unit_clouds,'clouds',lsynclouds,'Atm. forcing',
     &   'Clouds  []')

      call opnfrcwr(unit_relhum,'relhum',lsynclouds,'Atm. forcing',
     &   'Relative humidity  []')

      call opnfrcwr(unit_slp,'slp',lsynslp,'Atm. forcing',
     &   'Slp [mBar]')

Cdiag call opnfrcwr(998,'sss',lsynslp,'Ocn. forcing',
Cdiag&   'sst [C]')
Cdiag call opnfrcwr(999,'sst',lsynslp,'Ocn. forcing',
Cdiag&   'sst [C]')
C
C --- ----------------------------------------------------------------
C --- ----------------------------------------------------------------
C --- -------------------- Climatology forcing -----------------------
C --- ----------------------------------------------------------------
C --- ----------------------------------------------------------------

      ! TODO: yrflag=2 can be hf forcing too
      if (yrflag<3) then 



         do mo=1,12

c ---    Read all climatology data for month "mo" and place them in slot 1
         write (lp,'(a,i3)') 'Read '//clmflag//' clim for month ',mo
         call read_clim_all(1,mo)

         !print *,minval(clmuwind(:,:,1)),maxval(clmuwind(:,:,1))
         !print *,minval(clmvwind(:,:,1)),maxval(clmvwind(:,:,1))
         !print *,minval(clmairtmp(:,:,1)),maxval(clmairtmp(:,:,1))

C ---    ----------------------------------------------
C ---    File IO to HYCOM-type files
C ---    ----------------------------------------------
         if (timer) tmr0=wtime()
         rmo=mo
         call writeforc(clmtaux  (:,:,1),901,"tau_ewd",lsynwnd,rmo)
         call writeforc(clmtauy  (:,:,1),902,"tau_nwd",lsynwnd,rmo)
         call writeforc(clmwndspd(:,:,1),903," wndspd",lsynwnd,rmo)
         fldout=clmairtmp(:,:,1)-t0deg
         call writeforc(fldout,904," airtmp",lsynairtmp,rmo)
C        call writeforc(clmairtmp(:,:,1)-t0deg,
C    &      904," airtmp",lsynairtmp,rmo)
         call writeforc(clmvapmix(:,:,1),905," vapmix",lsynvapmix,rmo)
         call writeforc(clmprecip(:,:,1),906," precip",lsynprecip,rmo)
         call writeforc(clmshwflx(:,:,1),907," shwflx",lsynshwflx,rmo)
         call writeforc(clmshwflx(:,:,1),908," radflx",lsynradflx,rmo)
         iunit=unit_clouds
         call writeforc(clmclouds(:,:,1),iunit," clouds",lsynclouds,rmo)
         iunit=unit_relhum
         call writeforc(clmrelhum(:,:,1),iunit," relhum",lsynrelhum,rmo)
         iunit=unit_uwind
         call writeforc(clmuwind (:,:,1),iunit,"  uwind",lsynwnd,rmo)
         iunit=unit_vwind
         call writeforc(clmvwind (:,:,1),iunit,"  vwind",lsynwnd,rmo)
Cdiag    call writeforc(clmsss   (:,:,1),998,"    sss",.false.,rmo)
Cdiag    call writeforc(clmsst   (:,:,1),999,"    sst",.false.,rmo)
         if (timer) tmr1=wtime()
         if (timer) print *,'One clm field write took ',tmr1-tmr0

         ! No slp climatology for old
         if (trim(clmflag)=="old") then
            write(unit_slp,'("  slp  :month,range = ",i2.2,2e16.8)')
     &            mo,slp0,slp0
            call flush(unit_slp)
         else
            iunit=unit_slp
            call writeforc(clmslp (:,:,1),iunit,"    slp",lsynslp,rmo)
         end if
C ---    ----------------------------------------------
C ---    End File IO to HYCOM-type files
C ---    ----------------------------------------------

         if (timer) tmr0=wtime()
         call clim_diag(mo,fdiag,ifdiag)
         if (timer) tmr1=wtime()
         if (timer) print *,'One clm diag (tecplot) ',tmr1-tmr0

         end do ! mo
      else

C --- ----------------------------------------------------------------
C --- ----------------------------------------------------------------
C --- ---------- High Frequency (synpotic) forcing -------------------
C --- ----------------------------------------------------------------
C --- ----------------------------------------------------------------
      ! ECMWF/Met.no/NCEP files -- 6hr forcing...
      if (rforce=='ecmwf' .or. rforce=='ncepr'.or. 
     &    rforce=='metno' .or. rforce=='era40'.or.
     &    rforce=='storm' .or. rforce=='ecnc' .or.
     &    rforce=='era-i') then
      else
         print *,'Unknown forcing ID:'//rforce
         print *,'If you want to create  climatological forcing,'
         print *,'make sure that yrflag=0 in blkdat.input'
         print *

         print *, '(forfun_nersc)'
         call exit(1)
      end if
      

      ! Start cycling records
      clmslot=-1
      dtime=time1
      call year_day(dtime,refyear,rtinit,yrflag)
      do while (dtime<time2+rdtime .and. .not. bailout)

         
         call year_day(dtime,refyear,rttmp,yrflag)


         ! Read forcing fields
         if (dtime<time1+g_switch) then
            if (rforce=='ecmwf' .or. rforce=='metno') then
                call read_ecmwf(rttmp,plon,plat,depths)
            elseif (rforce=='ncepr') then
               call read_ncep(rttmp,plon,plat,depths)
            elseif (rforce=='era40') then
               call read_era40(rttmp,plon,plat,depths)
            elseif (rforce=='era-i') then
               call read_erai(rttmp,plon,plat,depths)
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
            !Fanf we take off erroneous values.
            synclouds(:,:)=min(1.0,max(synclouds(:,:),0.))
            synrelhum(:,:)=max(0.0,min(synrelhum(:,:),1.0))
            synprecip(:,:)=max(0.0,synprecip(:,:))

         else
C ---       after g_switch days we switch smoothly during  
C ---       l_switch days to the clim : fanf
            ! Switch weight 1 after l_switch period is over
            w_switch=min(1.,(dtime-time1-g_switch)/l_switch)
            write(lp,'(a,3i6,a,f6.2)')
     &         'rt= ',rttmp%iyy,rttmp%idd,rttmp%ihh,
     &         ' '//clmflag//' Climatology transition weight is:',
     &         w_switch
         endif


c
c --- -- ----------------------------------------------------------
c --- -- Time interpolation weights                                 
c --- -- ----------------------------------------------------------
c
c ---    Linear interpolation of climatology fields - set month
c ---    midpoints straddling current time (dtimeup, dtimelw)
         !print *,'setting months in main loop'
         lwmonth=rttmp%imm ; lwyear=rttmp%iyy
         if (rttmp%idm/14.<1) call previousmonth(rttmp,lwyear,lwmonth)
         call hycomday(dtimelw,lwyear,lwmonth,15,0,0,rttmp2,yrflag)
c    
         upmonth=rttmp%imm; upyear=rttmp%iyy
         if (rttmp%idm/14.>=1.) call nextmonth(rttmp,upyear,upmonth)
         call hycomday(dtimeup,upyear,upmonth,15,0,0,rttmp2,yrflag)
c
c ---    Check if climate forcing slots are filled with correct data
c ---    Find missing month(s)
         lwslot=-1 ; upslot=-1
         do islot=1,2
            if(clmslot(islot)==lwmonth) lwslot=islot
            if(clmslot(islot)==upmonth) upslot=islot
         end do
c
         lupdate=any((/lwslot,upslot/)==-1)
         if (lwslot==-1 .and. upslot == -1) then
c ---       Both slots must be read
            clmslot(1)=lwmonth ; lwslot=1
            clmslot(2)=upmonth ; upslot=2
            call read_clim_all(1,lwmonth)
            call read_clim_all(2,upmonth)
         else if (lwslot==-1) then
c ---       forcing for lwmonth must be read
            lwslot=mod(upslot,2)+1
            clmslot(lwslot)=lwmonth
            call read_clim_all(lwslot,lwmonth)
         else if (upslot==-1) then
c ---       forcing for upmonth must be read 
            upslot=mod(lwslot,2)+1
            clmslot(upslot)=upmonth
            call read_clim_all(upslot,upmonth) 
         end if
c ---    Final calculation of weights
         tw1=(dtime-dtimelw)/(dtimeup-dtimelw)
         tw0=(1.-tw1)

         if (lupdate) then
            do islot=1,2
               print '(a,i4,a,i4,a,2i4,a,2f6.2)',
     &         '-- Climatology slot:'
     &          ,islot,'  has data for  month: ', clmslot(islot) 
                                  
            end do
         end if
         if ( mod(int(dtime/rdtime),10) ==0 ) then
            print '(a,2i4,a,2f6.2,a,2i4,a,2i4)',
     &          '-- curent month/dom is ',rttmp%imm, rttmp%idm, 
     &          ' -- low/hi weights are ',tw0,tw1,
     &          ' -- low/hi slots are  ',lwslot,upslot, 
     &          ' -- low/hi months are  ',lwmonth,upmonth
         end if
         lsyntst =mod(rttmp%idd,2)==0.and.syndiagon

c ---    Sanity check
         if (tw1<0. .or. tw1>1. ) then
            print *,'ups...tw0,tw1',tw0,tw1
            print *,'dtimelw,lwmonth,lwyear : ',dtimelw,lwmonth,lwyear
            print *,'dtime,juld  : ',dtime,juld
            print *,'dtimeup,upmonth,upyear : ',dtimeup,upmonth,upyear
            print *, '(forfun_nersc)'
            call exit(1)
         end if
c
c ---    weight diagnostic
         if (lfirstsyn) then
            open(543,file='clmweight.asc',status='replace')
            lfirstsyn=.false.
         end if
         write(543,'(3f12.3,2i5)') dtime,tw0,tw1,lwmonth,upmonth
         call flush(543)
c
c
c --- -- ----------------------------------------------------------
c --- -- Post-processing of synoptic fields
c --- -- ----------------------------------------------------------
c
         call interptime(tmpcloud,synclouds,clmclouds, 
     &      tw0,tw1,lwslot,upslot,lsynclouds)
         call interptime(tmpairtmp,synairtmp,clmairtmp, 
     &      tw0,tw1,lwslot,upslot,lsynairtmp)
c ---
C$OMP PARALLEL DO PRIVATE(i,j)
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm
c ---       Convert rel hum to vapmix
            synvapmix(i,j) = rhtovpmix(synrelhum(i,j),
     &                       synairtmp(i,j),slp0*100)
c ---       Calculate longwave radiation for radiative fluxes
c --- TODO  Note that this doesnt take into account clouds (!)
c ---       For the moment this field isnt used in our thermo-
c ---       dynamical formulation. 
            tmplw(i,j)= tmpairtmp(i,j)**4 -
     &                  (tw0*clmsst(i,j,lwslot) + 
     &                   tw1*clmsst(i,j,upslot) +t0deg)**4
            tmplw(i,j)=stefanb*0.95*tmplw(i,j)
         end do
         end do
C$OMP END PARALLEL DO
c
c ---    Calculate shortwave flux. Clouds may be synoptic or
c ---    climatology
         timetmp=rttmp%idd+rttmp%ihh/24.
         call qsw0(radfl0,cawdir,tmpcloud,plat/radian,
     &             timetmp,rttmp%daysinyear)
C$OMP PARALLEL DO PRIVATE(i,j)
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1,jdm
         do i=1,idm
            tmpshwflx(i,j)=radfl0(i,j)*tmpcloud(i,j) 
            tmpradfl (i,j)=tmplw(i,j)+tmpshwflx(i,j)
         end do
         end do
C$OMP END PARALLEL DO
         if (lsynshwflx) then
c ---       Use standard time interpolation if data set provides a 
c ---       shwflx data set
            call interptime(tmpshwflx,synshwflx,clmshwflx, 
     &         tw0,tw1,lwslot,upslot,lsynshwflx)
         else
c ---       No time interpolation between climate shwflx and
c ---       synoptic shwflx if not. This is because clouds can
c ---       still be synoptic
         end if
c
c --- -- ----------------------------------------------------------
c --- -- Final write-out
c --- -- ----------------------------------------------------------
c
c ---    These fields were time - interpolated above - NB: airtmp output in celsius
         fldout=tmpairtmp-t0deg
         call writeforc(fldout,904," airtmp",lsynairtmp,dtime,.true.) 
C        call writeforc(tmpairtmp-t0deg,904," airtmp",lsynairtmp,
C    &                  dtime,.true.) 
         call writeforc(tmpcloud,unit_clouds," clouds",lsynclouds,
     &                  dtime,.true.)
         call writeforc(tmpshwflx,908," shwflx",lsynshwflx,dtime,.true.)
         call writeforc(tmpradfl,907," radflx",lsynradflx,dtime,.true.)

c ---    Do time interpolation and final write-out of remaining forcing fields
         call interptime(fldout,syntaux,clmtaux, 
     &      tw0,tw1,lwslot,upslot,lsynwnd)
         call writeforc(fldout,901," tauewd",lsynwnd,dtime,.true.)
c
         call interptime(fldout,syntauy,clmtauy, 
     &      tw0,tw1,lwslot,upslot,lsynwnd)
         call writeforc(fldout,902," taunwd",lsynwnd,dtime,.true.)
c
         call interptime(fldout,synwndspd,clmwndspd, 
     &      tw0,tw1,lwslot,upslot,lsynwnd)
         call writeforc(fldout,903," wndspd",lsynwnd,dtime,.true.)
c
         call interptime(fldout,synvapmix,clmvapmix, 
     &      tw0,tw1,lwslot,upslot,lsynvapmix)
         call writeforc(fldout,905," vapmix",lsynvapmix,dtime,.true.)
c
         call interptime(fldout,synprecip,clmprecip, 
     &      tw0,tw1,lwslot,upslot,lsynprecip)
         call writeforc(fldout,906," precip",lsynprecip,dtime,.true.)
c
         call interptime(fldout,synrelhum,clmrelhum, 
     &      tw0,tw1,lwslot,upslot,lsynrelhum)
         call writeforc(fldout,unit_relhum," relhum",lsynrelhum,dtime,
     &                  .true.)
c
         call interptime(fldout,synuwind,clmuwind, 
     &      tw0,tw1,lwslot,upslot,lsynwnd)
         call writeforc(fldout,unit_uwind,"  uwind",lsynwnd,dtime,
     &                  .true.)
c
         call interptime(fldout,synvwind,clmvwind, 
     &      tw0,tw1,lwslot,upslot,lsynwnd)
         call writeforc(fldout,unit_vwind,"  vwind",lsynwnd,dtime,
     &                  .true.)
c
         call interptime(fldout,synslp,clmslp, 
     &      tw0,tw1,lwslot,upslot,lsynslp)
         call writeforc(fldout,unit_slp,"    slp",lsynslp,dtime,.true.)
c

c ---    Prepare for next iteration
         dtime=dtime+rdtime
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
!!!!!!!!!   end climatology and high-frequency forcing         !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
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
! ---   Close  HYCOM forcing files -- Actual data & Header
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

      if (.not.lfirstsyn) close(543)


  99  FORMAT(30I4) 
 100  FORMAT(10(1x,e10.4)) 
      write(lp,*)
      write(lp,*)
      write(lp, '(a)') 'Done producing '//rforce//' Forcing fields'
      if (rotateflag==1) then
         write(lp, '(a)') 'Velocity rotation to model grid '//
     &                    'done in lon/lat space   '
      else if (rotateflag==2) then
         write(lp, '(a)') 'Velocity rotation to model grid '//
     &                    'done on sphere '
      end if
      !write(lp,'(a)') *,'These fields use synoptic forcing :'
      !if (lsynwnd    ) write(lp,'(a)') 'winds'
      !if (lsynairtmp ) write(lp,'(a)') 'airtmp'
      !if (lsynrelhum ) write(lp,'(a)') 'relative humidity'
      !if (lsynprecip ) write(lp,'(a)') 'precipitation'
      !if (lsynclouds ) write(lp,'(a)') 'clouds'
      !if (lsynvapmix ) write(lp,'(a)') 'vapor mixing ratio'
      !if (lsynradflx ) write(lp,'(a)') 'radiative flux'
      !if (lsynshwflx ) write(lp,'(a)') 'shwflx'
      !if (lsynslp    ) write(lp,'(a)') 'slp'
      !if (lsynshwflx)    ) write(lp,'(a)') 'show'
      call exit(0)
      end program forfun_nersc
