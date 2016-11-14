      program p_copymem
      !------------------------------------------------
      ! Copies an ensemble member from member mem1 
      ! to member mem2. Modifies restart files for 
      ! EVP, Random forcing and ICE. 
      ! NB: Does not modify hycom restart  files 
      !------------------------------------------------
      use mod_xc
      use mod_za
      use mod_year_info
      use mod_rand_io
      use mod_ice_io
      use mod_evp_io
      use m_parse_blkdat
      implicit none


      character(len=40) :: arg(6)
      character(len=3)  :: rungen
      character(len=1)  :: answer
      integer :: i,mem1,mem2,year,day,hour,nmem
      type(year_info) :: rt

      integer*4,external :: iargc
      logical,  external :: LSAMEN
      real*8 :: juld_offset, dtime, juld1, time
      integer :: hrefyear, yrflag,kdm, nstep, ivar
      real :: rvar, baclin
      logical :: lice, levp, lrand

   
      ! Number of arguments supplied to routine
      if (iargc()==6) then
         do i=1,6
            call getarg(i,arg(i))
         end do
         
         read(arg(1),'(a3)') rungen   ; print *,'rungen=',rungen
         read(arg(2),'(i)') year     ; print *,'year  =',year
         read(arg(3),'(i)') day      ; print *,'day   =',day 
         read(arg(4),'(i)') hour     ; print *,'hour  =',hour
         read(arg(5),'(i)') mem1     ; print *,'mem1  =',mem1
         read(arg(6),'(i)') mem2     ; print *,'mem2  =',mem2
      else
         print *,'Copies an ensemble member from member mem1 '
         print *,'to member mem2. Modifies restart files for '
         print *,'EVP, Random forcing and ICE. '
         print *,'NB: Does not modify hycom restart  files '
         print *,'Usage:'
         print *,'copymem rungen year day hour mem1 mem2'
         print *,'     filename = name of file ...'
         print *,'     rungen   = 3 letter run ID    '
         print *,'     year     = year               '
         print *,'     day      = julian day-of-year '
         print *,'     hour     = hour of day        '
         print *,'     mem1     = member to copy from'
         print *,'     mem2     = member to copy to'
         print *
         stop '(copymem)'
      end if

      call parse_blkdat('yrflag','integer',rvar,yrflag)
      call parse_blkdat('kdm   ','integer',rvar,kdm)
      call parse_blkdat('baclin','real',baclin,ivar)



      time = day + hour/24.
      call year_day(time,year,rt,'ecmwf')


      ! initialize - sets idm and jdm from regional files
      call xcspmd()
      call zaiost()

      ! initalize reading modules
      call initrand()
      call initice()
      call initevp()

      ! Calculate nstep dtime
      if (yrflag==3) then
         hrefyear=1901
         juld_offset = float(datetojulian(year,1,1,hrefyear,1,1))
         juld1=float(datetojulian(year,1,1,hrefyear,1,1))
         juld1=juld1+time
         if (year<1901) then
            if (mnproc==1) then
               write(lp,*) ' For yrflag==3 year>=1901'
               call xcstop('(oldtonewrtestart)')
               stop '(oldtonewrtestart)'
            end if
         end if
      else if (yrflag==0) then
         hrefyear=0
         juld1=float(datetojulian(year,1,1,hrefyear,1,1))
         juld1=juld1+time
      else
         if (mnproc==1)
     &   write (lp,*) 'Only yrflag=3 or 0 supported at the moment..'
         call xcstop('(oldtonewrestart)')
         stop '(oldtonewrestart)'
      end if

      dtime=juld1
      nstep=nint(dtime*(86400.0d0/baclin))
      dtime=nstep/(86400.0d0/baclin)
      print *,'nstep and dtime= ',nstep,dtime
      print *

      lice=.false.
      levp=.false.
      lrand=.false.

      print '(a,i3,a,i3,a)','This wil copy RAND, ICE and EVP member ',
     &      mem1,' into member ',mem2,'.'
!      write (*,'(a)',advance='no') 
!     &      'Are you sure this is what you want ? '
!      read  (*,'(a)') answer




      ! Read mem1 from restart files -- HYCOM
      print '(a,i3,a,i3,a)','This wil copy RAND member ',mem1,
     &      ' into member ',mem2,'.'
      if (read_restart_rand(rungen,rt,mem1,nstep,dtime,kdm)) then
         if (save_restart_rand(rungen,rt,mem2,nstep,dtime,kdm)) then
            print *,'Copy OK'
            lrand=.true.
         else
            print *,'An error occured while writing to RAND member',mem2
               stop '(copymem)'
         end if
      else
         print *,'An error occured while reading the '//
     &           'RAND file -- Aborting'
      end if
      print *
      print *
      print *


      print '(a,i3,a,i3,a)','This wil copy ICE member ',mem1,
     &   ' into member ',mem2,'.'
      if (read_restart_ice(rungen,rt,mem1)) then
         if (save_restart_ice(rungen,rt,mem2)) then
            print *,'test min/max ficem',minval(ficem),maxval(ficem)
            print *,'test min/max hicem',minval(hicem),maxval(hicem)
            print *,'Copy OK'
            lice=.true.
         else
            print *,'An error occured while writing to ICE member ',mem2
            stop '(copymem)'
         end if
      else
         print *,'Could not read ICE -- Aborting'
      end if
      print *
      print *
      print *




      print '(a,i3,a,i3,a)','This wil copy EVP member ',mem1,
     &      ' into member ',mem2,'.'
      if(read_restart_evp(rungen,rt,mem1,nstep,dtime)==0) then
         if (save_restart_evp(rungen,rt,mem2,nstep,dtime)==0) then
            print *,'Copy OK'
            levp=.true.
         else
            print *,'An error occured while writing to EVP member ',mem2
            stop '(copymem)'
         end if
      else
         print *,'An error occured while reading EVP -- Aborting'
      end if
      print *
      print *
      print *

      print '(a)','Summary:'
      if (lice) then
         print '(a)', 'ICE  copy succeeded'
      else
         print '(a)', 'ICE  copy failed'
      end if
      if (levp) then
         print '(a)', 'EVP  copy succeeded'
      else
         print '(a)', 'EVP  copy failed'
      end if
      if (lrand) then
         print '(a)', 'RAND copy succeeded'
      else
         print '(a)', 'RAND copy failed'
      end if

 
 

      end program
