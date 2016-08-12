module m_limits

private :: limits_randf

contains
subroutine limits(yrflag, iniflg, relax, srelax, trelax,tidflg,lbflag,nestfq,bnstfq,flxflg)
   use mod_xc
   use mod_diagnostics
   use mod_forcing_nersc 
   use mod_average      , only : l_weekly_average_in_code, laverage, average_dt
   use mod_daily_average, only : l_daily_average_in_code, l_daily_average, l_daily_accum
   use mod_year_info, only: year_info,dtime_start,refyear
   use mod_nesting
   use mod_gridp
   use mod_tides_nersc  , only : ltides, ltideuv, lslowstart, ctide
   use mod_hycom_nersc
   implicit none
   real, parameter :: version=3.1     ! version of limits routine

   ! Input from blkdat - mainly sed for consistency checks
   integer, intent(in)  :: yrflag   ! Year flag from blkdat
   integer, intent(in)  :: iniflg   ! Init flag from blkdat (restart=3)
   logical, intent(in)  :: relax    ! Lateral T/S/rho nudging
   logical, intent(in)  :: srelax   ! Surface salinity nudging
   logical, intent(in)  :: trelax   ! Surface temperature nudging
   integer, intent(in)  :: tidflg   ! hycom tidal flag
   integer, intent(in)  :: lbflag   ! hycom nesting flag 
   real   , intent(in)  :: nestfq   ! hycom nest nudging (baroclinic)
   real   , intent(in)  :: bnstfq   ! hycom barotropic nesting 
   integer, intent(in)  :: flxflg   ! thermo forcing flag

   integer i,itmp
   logical ex

   integer nday1,nhour1    ! Current run starts at nday1,nhour1
   integer nday2,nhour2    ! Current run ends at nday2,nhour2
   real rday1    ! Current run starts at rday1
   real rday2    ! Current run ends at rday2
   real rday
   real fversion    ! version of input file

   ! KAL - version bump
   real :: baclin, batrop  !Local version - not set from input file any more
   integer :: seed            !Local version - not set from input file any more
   integer :: rstsgn
   integer :: ios, ios2

   character(len=80) :: inputtmp
   character(len=7 ) :: char7

   if (mnproc==1) print '(a)','####################################################'
   if (mnproc==1) print '(a)','Reading input parameters from infiles "limits"'
   average_dt=0


   inquire(file='infile.in',exist=ex)
   if (.not.ex) then
      if (mnproc==1) write(lp,*) 'ERROR: infile.in does not exist'
      call xcstop('(limits)')
      stop '(limits)'
   endif
   if (mnproc==1) write(*,'(a)')'### infile.in'
   open(10,file='infile.in',action='read')
   read(10,*)fversion
   if (version /= fversion) then
      if (mnproc==1) write(lp,*) 'Wrong version of input file infile.in'
      if (mnproc==1) write(lp,*) 'Expected version is: ',version
      if (fversion == 2.8) then
         if (mnproc==1) write(lp,*) 'have now introduced:'
         if (mnproc==1) write(lp,*) 'lslowtide'
         if (mnproc==1) write(lp,*) 'day hour rather than day.dec for diag days'
      else if (fversion == 2.9) then
         if (mnproc==1) then
            write(lp,*) 'Changes relative to 2.9:'
            write(lp,*) 'restart flag removed (to blkdat.input)'
            write(lp,*) 'Ensemble member moved (input argument)'
            write(lp,*) 'Added climatology option (after force option)'
            write(lp,*) 'baroclinic time step moved (to blkdat.input)'
            write(lp,*) 'barotropic time step moved (to blkdat.input)'
            write(lp,*) 'Added daily average option (after weekly average option)'
            write(lp,*) 'Accumulate sections removed (post processing routine)'
            write(lp,*) 'Movie option removed '
         end if
      else if (fversion == 3.0) then
         if (mnproc==1) then
            write(lp,*) 'Changes relative to 3.0:'
            write(lp,*) 'River flag removed (set in blkdat.input)'
            write(lp,*) 'Read latlon from file flag removed (always done)'
            write(lp,*) 'Print forcing functions flag removed'
         end if
      endif
      call xcstop('(limits)')
      stop '(limits)'
   endif
   read(10,'(a3)') rungen
   if (mnproc==1) write(*,'(a,a3)')   'infile.in:  rungen  = ',rungen

   read(10,*)refyear
   if (mnproc==1) write(*,'(a,i4)')   'infile.in:  refyear = ',refyear

   read(10,*)nday1,nhour1
   !rday1=float(nday1)+float(nhour1)/24.0
   call dayfor(rday1,yrflag,refyear,nday1,nhour1)
   if (mnproc==1) write(*,'(a,i4,i3,a,f8.2)') 'infile.in: start at = ',nday1,nhour1,' or ',rday1

   read(10,*)nday2,nhour2
   !rday2=float(nday2)+float(nhour2)/24.0
   call dayfor(rday2,yrflag,refyear,nday2,nhour2)
   if (mnproc==1) write(*,'(a,i4,i3,a,f8.2)') 'infile.in:  end at  = ',nday2,nhour2,' or ',rday2

   ! KAL - now accepts climate option
   read(10,'(a5,x,a5)')rforce,clmflag
   if (mnproc==1) write(*,'(a,a5,a,a5)')   'infile.in:  rforce  = ',rforce, &
      ', clmflag = ',clmflag
   if (rforce == 'ecmwf' .or. rforce == 'month' .or. rforce == 'waner' .or. &
       rforce == 'nwagr' .or. rforce == 'ncepr' .or. rforce == 'arc01' .or. &
       rforce == 'metno' .or. rforce == 'era40'.or. rforce == 'ecmo'   .or. &
       rforce=='storm' .or. rforce == 'era-i'                          .or. &
       rforce == 'ecnc'.or. rforce == 'ecnc2' ) then
   else
      if (mnproc==1)  &
      write(lp,*) 'limits: invalid option for rforce:'//rforce
      call xcstop('(limits)')
      stop '(limits)'
   endif

   ! Check climate option - can also be used in synoptic forcing
   if (trim(clmflag)/='old'.and.trim(clmflag)/='era40'.and.trim(clmflag)/='prep') then
      if (mnproc==1) then
         write(*,'(a)') 'infile.in:  WARNING: Unknown climate flag '//clmflag
         write(*,'(a)') 'infile.in:  Use era40 or old '
      end if
      call xcstop('(limits)')
      stop '(limits)'
   end if

#if defined (FORCING_INLINE) 
   if (trim(clmflag)=='prep') then
      if (mnproc==1) then
      write(lp,'(a)') 'To use the "prep" climate flag you must disable the'
      write(lp,'(a)') 'FORCING_INLINE flag when compiling HYCOM'
      end if
      call xcstop('(limits)')
   end if
#endif

   ! Consistency check vs yrflag
   if (yrflag==3) then 
      if (rforce == 'month') then 
         if (mnproc==1) then
            write(lp,*) ' For yrflag==3 use synoptic forcing'
            write(lp,*) ' Current forcing is '//rforce
         end if
         call xcstop('(limits)')
         stop '(limits)'
      end if
      if (refyear<1901) then
         if (mnproc==1) then
            write(lp,*) ' For yrflag==3 refyear>=1901'
         end if
         call xcstop('(hycom)')
         stop '(hycom)'
      end if
   else if (yrflag==0 .or. yrflag==1 .or. yrflag==2) then 
      if (rforce /= 'month') then 
         if (mnproc==1) then
            write(lp,*) ' For yrflag==0 use climatology forcing'
            write(lp,*) ' Current forcing is '//rforce
         end if
         call xcstop('(limits)')
      end if
   end if


   read(10,*)time_trelax
   if (mnproc==1) write(*,'(a,f9.3)') 'infile.in:  time_trelax    = ',time_trelax

   read(10,*)time_srelax
   if (mnproc==1) write(*,'(a,f9.3)') 'infile.in:  time_srelax    = ',time_srelax

   ! relaxation should be switched off in blkdat.input
   if ( time_trelax < 1.0   .or. time_srelax < 1.0 ) then
      if (mnproc==1) write(*,'(a)') 'Error: Surface relaxation times are too low!'
      if (mnproc==1) write(*,'(a)') '(If you want to turn off surface relaxation'
      if (mnproc==1) write(*,'(a)') ' do it in blkdat.input)'
      call xcstop('(limits)')
      stop '(limits)'
   end if

   read(10,'(l1,tr1,i2)')laverage,average_dt
   if (mnproc==1) write(*,'(a,l1,i3)')'infile.in:  ave_dt    = ',laverage,average_dt

   ! Test if weekly averages is actually compiled into model
   if (.not.l_weekly_average_in_code.and.laverage) then
      if (mnproc==1) then
         write(lp,'(a)') 'FATAL: Weekly averages is not compiled into '
         write(lp,'(a)') '       code, you need to recompile code with'
         write(lp,'(a)') '       CPP Flag "WEEKLY_AVERAGE" defined,  '
         write(lp,'(a)') '       or switch off averages altogether... '
      end if
      call xcstop('(limits)')
      stop '(limits)'
   end if

   !read(10,'(l1,x,l1)') l_daily_average, l_daily_accum
   read(10,*) l_daily_average, l_daily_accum
   if (mnproc==1) write(*,'(a,l1,x,l1)')'infile.in:  daily_ave = ',l_daily_average,l_daily_accum

   ! Test if weekly averages is actually compiled into model
   if (.not.l_daily_average_in_code.and.l_daily_average) then
      if (mnproc==1) then
         write(lp,'(a)') 'FATAL: Daily  averages is not compiled into '
         write(lp,'(a)') '       code, you need to recompile code with'
         write(lp,'(a)') '       CPP Flag "DAILY_AVERAGE" defined,  '
         write(lp,'(a)') '       or switch off daily averages altogether... '
      end if
      call xcstop('(limits)')
      stop '(limits)'
   end if
  
  

   read(10,'(l1,tr1,i2)')lnesto,nestdto
   if (mnproc==1) write(*,'(a,l1,tr1,i3)')'infile.in:  nesto     = ',lnesto,nestdto


   read(10,'(l1,tr1,i2)')lnesti,nestdti
   if (mnproc==1) write(*,'(a,l1,tr1,i3)')'infile.in:  nesti     = ',lnesti,nestdti
#ifndef NEST_INNER
   if (lnesti) then
      if (mnproc==1) then
         write(lp,'(a)') 'Nesting switched on at runtime but not at compile-time'
         write(lp,'(a)') 'Re-compile model with CPP-flag NEST_INNER defined'
      end if
      call xcstop('(limits)')
      stop '(limits)'
   end if
#endif
   ! Check against blkdat.input setup. For now barotropic and 3d nesting times must match
   ! TODO - this might be a good time to remove inner nesting flags entirely here..
   if (lnesti) then
      if (lbflag.ne.2) then
         if (mnproc==1) then
            write(lp,'(a)') 'Mismatch between blkdat.input and infile.in inner nesting flags'
         end if
         call xcstop('(limits)')
         stop '(limits)'
      end if
      if (nint(nestfq*24.).ne.nestdti .or. nint(bnstfq*24.).ne.nestdti) then
         if (mnproc==1) &
         write(lp,'(a)') 'Mismatch between blkdat.input and infile.in inner nesting times'
         call xcstop('(limits)')
         stop '(limits)'
      end if
   end if


   read(10,'(l1,tr1,a3,tr1,l1,tr1,l1)')ltides,ctide,ltideuv,lslowstart
   if (mnproc==1) write(*,'(a,l1,tr1,a3,tr1,l1,tr1,l1)')'infile.in:  tides     = ',&
                                       ltides,ctide,ltideuv,lslowstart

   read(10,'(l1)')lgridp
   if (mnproc==1) write(*,'(a,l1)')   'infile.in:  lgridp  = ',lgridp

   read(10,*)
   do i=0,nrdday
      read(10,'(t3,i4,t8,i2,t12,l1,t18,l1,t27,l1)',end=100,err=100)&
                  dday(i)%nday, dday(i)%nhour, dday(i)%ass, dday(i)%dia, dday(i)%res
      if (dday(i)%nhour > 23) then
         if (mnproc==1) write(lp,*) 'nhour must be integer in (0,..,23).  See line ',i+1
         call xcstop('(limits)')
         stop '(limits)'
      endif
      !dday(i)%day=float(dday(i)%nday)+float(dday(i)%nhour)/24.0
      call dayfor(dday(i)%day,yrflag,refyear,dday(i)%nday,dday(i)%nhour)
   enddo
   100 iday2=i-1
   close(10)

   if (laverage) then
      if (average_dt == 0) then
         if (mnproc==1) write(lp,*) 'infile.in:  laverage is true and average_dt =0'
         if (mnproc==1) write(lp,*) 'infile.in:  Correct parameter in infile.in'
         if (mnproc==1) write(lp,*) 'infile.in:  Check location in limits.F90'
         call xcstop('(limits)')
         stop '(limits)'
      endif
      if (mod(24,average_dt) /= 0) then
         if (mnproc==1) write(lp,*) 'infile.in:  laverage is true and average_dt not an integer in 24'
         if (mnproc==1) write(lp,*) 'infile.in:  valid numbers are 1,2,3,4,6,8,12,24'
         if (mnproc==1) write(lp,*) 'infile.in:  Correct parameter in infile.in'
         if (mnproc==1) write(lp,*) 'infile.in:  Check location in limits.F90'
      endif
   endif
   iday1=0
   !if (mnproc==1) print *,'iday1 ',iday1
   

   ! KAL - NEW: set a dtime to read upon first restart. infile.in gives time rel
   ! to 0, hycom time rel to 1.
   call dayfor(dtime_start,yrflag, refyear,nday1+1,nhour1)

   ! A diagnostic time equal to rday1 must also be given 
   ! - KAL but only if this is the first step !
   !if (iniflg==3) then
      iday1=-1
      do i=0,iday2
         if (dday(i)%day == rday1) then
            iday1=i
            exit
         endif
      enddo
      if (iday1 == -1 ) then
         if (mnproc==1) write(lp,*) 'infile.in: Problem with restart time vs diagno time'
         if (mnproc==1) write(lp,*) ' A diagnostic time equal to rday1 must be given'
         call xcstop('(limits)')
         stop '(limits)'
      endif
   !endif

!KAL - loosen this restriction
!KAL   if ((iniflg/=3).and.(rday1 /= dday(0)%day)) then
!KAL      if (mnproc==1) write(lp,*) 'infile.in: rday1=',rday1,' should be equal to dday(0)%day:',dday(0)%day 
!KAL      call xcstop('(limits)')
!KAL      stop '(limits)'
!KAL   endif

   ! Final diagnostic day
   itmp=iday2
   do i=iday1+1,itmp
      if (dday(i)%day > rday2) then
         iday2=i-1
         exit
      endif
   enddo
   if (iday1 == iday2) then
      if (mnproc==1) write(lp,*) 'infile.in: Warning: iday1=iday2=',iday1
      if (mnproc==1) write(lp,*) 'infile.in: No integration will be performed'
      call xcstop('(limits)')
      stop '(limits)'
   endif

   ! Check for increasing time coordinate
   do i=iday1+1,iday2
      if (dday(i)%day <= dday(i-1)%day) then
         if (mnproc==1) write(lp,*) 'infile.in: Inconsistency in sequence of diagnostic times'
         call xcstop('(limits)')
         stop '(limits)'
      endif
   enddo

   ! Check that a final diagnostic time is supplied
   if (dday(iday2)%day /= rday2) then
      if (mnproc==1) print *,'rday2      ',rday2
      if (mnproc==1) print *,'day(iday2) ',dday(iday2)
      if (mnproc==1) write(lp,*) 'infile.in: day(iday2) /= rday2'
      call xcstop('(limits)')
      stop '(limits)'
   endif

!KAL - Dont need these if we follow hycom time
!KAL   juld_offset = float(datetojulian(refyear,1,1,hrefyear,1,1))
!KAL   juld1=float(datetojulian(refyear,1,1,hrefyear,1,1))
!KAL   juld1=juld1+dday(iday1)%day
!KAL   juld2=float(datetojulian(refyear,1,1,hrefyear,1,1))
!KAL   juld2=juld2+dday(iday2)%day

!KAL - dont touch limits file
   !rstsgn=1
   !if (yrflag==3 .or. iniflg<3) rstsgn=-1
   !if (mnproc==1.and.imem==1) then
   !   open(10,file='limits',status='replace')
   !   write(10,*) rstsgn*juld1,juld2
   !   close(10)
   !end if

   ! New test - make sure times > 0  Advice what to do if not
   if (rday1<0.) then
      if (mnproc==1) then
         write(lp,'(a)') 'start day < 0 (in "hycom" time)'
         write(lp,'(a)') 'Make sure diagnostic dates satisfy this : '
         write(lp,'(a)') 'yrflag=0: older than day 16 year 1 '
         write(lp,'(a)') 'yrflag=1: older than day 16 year 1 '
         write(lp,'(a)') 'yrflag=2: older than day  1 year 1 '
         write(lp,'(a)') 'yrflag=3: older than day  1 year 1901 '
      end if
      call xcstop('(limits)')
      stop '(limits)'
   end if


   if (mnproc==1) write(*,'(a,i5)') 'limits:  iday1 = ',iday1
   if (mnproc==1) write(*,'(a,i5)') 'limits:  iday2=  ',iday2
   do i=iday1,iday2
      if (mnproc==1) print '(a,i3,i5,i3,f8.2,3l8)','limits:  ddays are:',i,dday(i)
   enddo
   if (mnproc==1) print '(a)','####################################################'
   if (mnproc==1) print '(a)','Randfom forcing parameters:'
   call limits_randf()
   if (mnproc==1) print '(a)','####################################################'
   if (mnproc==1) write(*,*)

   
   ! Final tests
   if (relax .and. lnesti) then
      if (mnproc==1) &
      write(lp,'(a)') 'Both inner nesting and relaxation is active - choose one ..'
      call xcstop('(limits)')
   end if

   if (ltides .and. tidflg>0) then
      if (mnproc==1) &
      write(lp,'(a)') 'Both NERSC and standard tidal forcing is active - choose one ..'
      call xcstop('(limits)')
   else if (ltides .and. lbflag.ne.2) then
      if (mnproc==1) &
      write(lp,'(a)') 'NERSC tides must have lbflag == 2 (for now)'
      call xcstop('(limits)')
   end if

   if (flxflg.ne.99 .and. flxflg.ne.0) then
      if (mnproc==1) then
         write(lp,'(a,i2)')   'flxflg in blkdat.input=',flxflg
         write(lp,'(a,i2,a)') 'Use flxflg ',99,' in nersc version'
      end if
      call xcstop('(limits)')
   end if



end subroutine limits


subroutine limits_randf()
   use mod_forcing_nersc
   use mod_random_forcing
   implicit none
   real, parameter :: version2=1.2     ! version of limits routine
   logical :: ex
   integer :: seed
   real    :: fversion

   inquire(file='infile2.in',exist=ex)
   if (.not.ex) then
      if (mnproc==1) write(lp,*) 'ERROR: infile2.in does not exist'
      call xcstop('(limits)')
      stop '(limits)'
   endif

   open(uoff+99,file='infile2.in',status='old',action='read')
   call blkinr(fversion, 'fversn','(a6," =",f10.4," ")')
   if (abs(version2-fversion)>1e-4) then
      if (mnproc==1) write(lp,*) 'Wrong version of input file infile2.in'
      if (mnproc==1) write(lp,*) 'Expected version is',version2
      if (mnproc==1) write(lp,*) 'file     version is',fversion
      call xcstop('(limits)')
      stop '(limits)'
   endif
   call blkinl(randf        ,'randf ')
   call blkini(seed         ,'seed  ')
   call blkinr(vars%slp     ,'vslp  ','(a6," =",g10.4," mBar**2")')
   call blkinr(vars%taux    ,'vtaux ','(a6," =",g10.4," N**2m**-4")')
   call blkinr(vars%tauy    ,'vtauy ','(a6," =",g10.4," N**2m**-4")')
   call blkinr(vars%wndspd  ,'vwspd ','(a6," =",g10.4," m**2s**-2")')
   call blkinr(vars%clouds  ,'vcloud','(a6," =",g10.4," []")')
   call blkinr(vars%airtmp  ,'vtair ','(a6," =",g10.4," C**2")')
   call blkinr(vars%precip  ,'vprcp ','(a6," =",g10.4," ")')
   call blkinr(vars%relhum  ,'vrlhum','(a6," =",g10.4," []")')
   call blkinr(rf_hradius   ,'scorr ','(a6," =",g10.4," m")')
   call blkinr(rf_tradius   ,'tcorr ','(a6," =",g10.4," days")')
   call blkini(rf_prsflg    ,'prsflg')

   if (rf_prsflg>2 .or. rf_prsflg<0) then
      if (mnproc==1) then
         write(lp,*) 'Pressure flag must be 0 or 2'
      end if
      call xcstop('(limits:limits_randf)')
      stop '(limits:limits_randf)'
   else if (rf_prsflg==1) then
      if (mnproc==1) then
         write(lp,*) 'Pressure flag = 1 should not be used'
      end if
      call xcstop('(limits:limits_randf)')
      stop '(limits:limits_randf)'
   end if
   close(99)

end subroutine limits_randf
end module
