      subroutine blkdat(lhycom,layadd)
      use mod_xc   ! HYCOM communication interface
      use mod_rest ! HYCOM restart arrays
c
      logical lhycom
      integer layadd
c
      integer   k
      character sigfmt*26
c
c --- initialize common variables.
c
      open(unit=99,file='blkdat.input')
c
c --- 'lp' = logical unit number for printer output
      lp = 6
c
c --- 'iversn' = hycom version number x10
c --- 'iexpt'  = experiment number x10
      write(lp,*)
      call blkini(iversn,'iversn')
      call blkini(iexpt, 'iexpt ')
c
      if (iversn.lt.20 .or. iversn.gt.21) then
        write(lp,'(/ a,i3,a,i3 /)') 
     &    'error - iversn must be between',20,' and',21
        call flush(lp)
        stop '(blkdat)'
      endif
c
c --- 'kdm   ' = number of layers (including layadd)
      write(lp,*)
      call blkini(kdm,   'kdm   ')
c
      if     (layadd.ne.0) then
c
c --- 'layadd' = number of layers to add between layers 1 and 2
        write(lp,*)
        call blkini(layadd,'layadd')
c
        if (layadd.le.0) then
          write(lp,'(/ a /)') 
     &      'error - layadd much be positive'
          call flush(lp)
          stop '(blkdat)'
        endif
      endif
c
      if     (.not.lhycom .or. layadd.ne.0) then
        sigfmt = '(a6," =",f10.4," sigma-0")'
c
c ---   MICOM and layadd only.
c ---   'thbase' = reference density (sigma units)
        call blkinr(thbase,'thbase',sigfmt)
        write(lp,*)
c
c ---   layer densities (sigma units)
        do k=1,kdm
          call blkinr(sigma(k),'sigma ',sigfmt)
c
          if     (k.gt.1) then
            if      (sigma(k).le.sigma(k-1)) then
              write(lp,'(/ a,i3 /)')
     &          'error - sigma is not stabally stratified'
              call flush(lp)
              stop '(blkdat)'
            endif
          endif
        enddo
      endif !MICOM or layadd
c
c --- 'yrflag' = days in year flag (0=360,1=366,2=366Jan1,3=actual)
c --- 'baclin' = baroclinic time step (seconds), int. divisor of 86400
      write(lp,*)
      call blkini(yrflag,'yrflag')
      call blkinr(baclin,'baclin','(a6," =",f10.4," sec")')
c
      if (yrflag.lt.0 .or. yrflag.gt.3) then
        write(lp,'(/ a /)') 
     &    'error - yrflag must be between 0 and 3'
        call flush(lp)
        stop '(blkdat)'
      endif
      if     (yrflag.le.1) then  ! monthly forcing
        if     (abs(nint(86400.0/baclin)-86400.0/baclin).gt.0.01) then
          write(lp,'(/ a /)')
     &      'error - baclin not an integer divisor of 86400 secs'
          call flush(lp)
          stop '(blkdat)'
        endif
      else  ! high frequency forcing
        if     (abs(nint(21600.0/baclin)-21600.0/baclin).gt.0.01) then
          write(lp,'(/ a /)')
     &      'error - baclin not an integer divisor of 21600 secs'
          call flush(lp)
          stop '(blkdat)'
        endif
      endif
c
c --- 'iceflg' = ice model flag (0=none,1=energy loan model)
      write(lp,*)
      call blkini(iceflg,'iceflg')
      icegln = iceflg.eq.1
c
      if (iceflg.lt.0 .or. iceflg.gt.1) then
        write(lp,'(/ a /)') 
     &    'error - iceflg must be between 0 and 1'
        call flush(lp)
        stop '(blkdat)'
      elseif (iceflg.eq.1 .and. .not.lhycom) then
        write(lp,'(/ a /)') 
     &    'error - iceflg must be 0 for MICOM input'
        call flush(lp)
        stop '(blkdat)'
      endif
c
c --- use 'huge' to initialize array portions that the code should never access
      huge = 1.e33
c
c --- i/o file names
c
      flnmdep  = 'regional.depth'
      flnmrsa  = 'archive_in'
      flnmrsi  = 'restart_in'
      flnmrso  = 'restart_out'
c
      close(unit=99)
      return
      end
      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read(99,*) rvar,cvarin
      write(6,cfmt) cvarin,rvar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      write(6,6000) cvarin,ivar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
 6000 format(a6,' =',i6)
      end
      subroutine blkinl(lvar,cvar)
      implicit none
c
      logical     lvar
      character*6 cvar
c
c     read in one logical value
c     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
c
      character*6 cvarin
      integer     ivar
c
      read(99,*) ivar,cvarin
      lvar = ivar .ne. 0
      write(6,6000) cvarin,lvar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
 6000 format(a6,' =',l6)
      end
