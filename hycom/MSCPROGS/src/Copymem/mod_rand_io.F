      module mod_rand_io


      real, dimension(:,:), allocatable, save  ::  ! Hard-coded real*8 for file IO
     &   slp, taux, tauy, wndspd, airtmp, relhum, clouds,
     &   precip,sss,sst


      contains


      subroutine initrand()
      use mod_xc
      implicit none
      allocate(slp   (idm,jdm))
      allocate(taux  (idm,jdm))
      allocate(tauy  (idm,jdm))
      allocate(wndspd(idm,jdm))
      allocate(airtmp(idm,jdm))
      allocate(relhum(idm,jdm))
      allocate(clouds(idm,jdm))
      allocate(precip(idm,jdm))
      allocate(sss   (idm,jdm))
      allocate(sst   (idm,jdm))
      end subroutine



! Routine to read random forcing fields 
!22.11.2004 - KAL: Now uses hycom-type IO
!22.11.2004 - KAL: Random forcing vars are tiled by default
      function   read_restart_rand(rungen,rt,imem,nstep0,dtime0,kdm)
      use mod_xc
      use mod_za
      use mod_year_info
      implicit none
      character(len=3)    , intent(in)  :: rungen
      type(year_info)     , intent(in)  :: rt
      integer,              intent(in)  :: imem
      integer,              intent(in)  :: nstep0
      real*8,               intent(in)  :: dtime0
      integer,              intent(in) :: kdm
      logical :: read_restart_rand
c
      character(len=50) :: filename
      integer :: j,iostat
      logical :: ex
      character*80    l_ctitle(4)
      integer :: l_iversn, l_iexpt, l_yrflag,l_idm, 
     &           l_jdm, l_kdm, iyy, idd, ihh, iss
      character(len=50) :: filebase
      character(len=8) :: char8
      real :: xmin, xmax
      integer :: nop,ifld
      real :: xmin2, xmax2, l_time
      integer :: l_nstep,l_imem,nfields,offset,k,offsetb
      integer, dimension(idm,jdm) :: ip

      ! New output format -- uses hycom-style IO (zaio)
      nop=710
      filebase=rungen//'restRAND'//rt%cyy//'_'//rt%cdd//'_'//rt%chh
      inquire (exist=ex, file=trim(filebase)//'.a')
      if (.not.ex) then
         if (mnproc==1) write(lp,*)'No such file '//trim(filebase)//'.a'
         read_restart_rand=.false.
         return
      else
         if (mnproc==1) then
            write(lp,'(a,i4)')'read_restart_rand: reading '
     &         //trim(filebase)//'.[ab] record=',imem
         end if
      end if

      call zaiopf(trim(filebase)//'.a','old',nop)
      open(nop,file=trim(filebase)//'.b',status='old',
     &         form='formatted',access='direct',recl=80)

      ! Number of fields for one member 
      nfields=10

      ! Offset for IO
      offset=(imem-1)*nfields
      offsetb=offset+1 ! First record in .b  is model info

      ! Read "header"
      read(nop,118,rec=1) l_idm,l_jdm,l_kdm
      if (idm/=l_idm .or. jdm /= l_jdm .or. kdm/=l_kdm) then
         if (mnproc==1) write(lp,*) 'Dimension mismatch in RAND file'
         call xcstop('(read_restart_rand)')
         stop '(read_restart_rand)'
      end if


      ! Skip offset number of records
      do k=1,offset
         call zaiosk(nop)
      end do

      write(lp,'(a,i4,a,i4)')   'read_restart_rand: restart for ',
     &   imem,' starts at .a offset ',offset+1


      ! Start reading
      do ifld=1,nfields

         ! Same order as out...
         read(nop,117,rec=offsetb+ifld)
     &      char8,l_nstep,l_time,l_imem,xmin2,xmax2

         if (l_imem/=imem) then
            print *, 'Ensemble member mismatch..'
            read_restart_rand=.false.
            return
         end if

         if (mnproc==1) write(lp,*) char8,imem,xmin2,xmax2

         if     (trim(char8)=='slp') then
            call zaiord(slp(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='taux') then
            call zaiord(taux(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='tauy') then
            call zaiord(tauy(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='wndspd') then
            call zaiord(wndspd(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='airtmp') then
            call zaiord(airtmp(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='relhum') then
            call zaiord(relhum(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='clouds') then
            call zaiord(clouds(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='precip') then
            call zaiord(precip(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='sss') then
            call zaiord(sss(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='sst') then
            call zaiord(sst(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         else
            if (mnproc==1) print *, 'Unknown field '//char8
            call xcstop('(read_restart_rand')
            stop '(read_restart_rand'
         end if

         ! Check for coherency between min and maxes
         if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.
     &       abs(xmax-xmax2)>abs(xmax)*1e-4) then
            if     (mnproc.eq.1) then
               write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           'error - .a and .b files not consistent:',
     &           '.a,.b min = ',xmin2,xmin ,xmin2-xmin ,
     &           '.a,.b max = ',xmax2,xmax ,xmax2-xmax
               print *,'Variable :',char8
            endif
            call xcstop('(read_restart_rand)')
            stop '(read_restart_rand)'
         end if


      end do

      call zaiocl(nop)
      close(nop)
      read_restart_rand=.true.


 117  format (a8,' =',i11,f11.2,i3,1p2e16.7)
 118  format ('itdm= ',i5, ' jtdm= ',i5,' kdm= ',i5)
      end function read_restart_rand




! Routine to save random forcing fields
! 22.11.2004 - KAL: Changed to use HYCOM type io

      logical function save_restart_rand(rungen,rt,imem,
     &   nstep0,dtime0,kdm)
      use mod_xc
      use mod_za
      use mod_year_info
      implicit none
      character(len=3)    , intent(in) :: rungen
      type(year_info)     , intent(in) :: rt
      integer,              intent(in) :: imem
      integer,              intent(in) :: nstep0
      real*8,               intent(in) :: dtime0
      integer,              intent(in) :: kdm
      character(len=50) :: filename
      integer :: j,iostat,ios
      integer ioproc, ierr
      integer :: nop,nfields,offset,k,l_nstep,l_imem,offsetb
      real :: xmin,xmax, l_time
      character(len=50) :: filebase
      character(len= 1) :: char1,carret
      character(len= 8) :: char8
      integer, dimension(idm,jdm) :: ip

      ! New output format -- uses hycom-style IO (zaio)
      nop=710
      filebase=rungen//'restRAND'//rt%cyy//'_'//rt%cdd//'_'//rt%chh
      if (mnproc==1) then 
         write(lp,'(a,i4)',advance='yes') 'save_restart_rand: saving '//
     &      trim(filebase)//'.[ab] record=',imem
      end if
      call zaiopf(trim(filebase)//'.a','unknown',nop)
      if (mnproc==1) then
         open(nop,file=trim(filebase)//'.b',status='unknown',
     &        form='formatted',access='direct',recl=80)
      end if


      ! Number of fields for one member 
      nfields=10

      ! Offset for IO
      offset=(imem-1)*nfields
      offsetb=offset+1 ! First record is model info 
      if (mnproc==1) then
         write(nop,118,rec=1) idm,jdm,kdm
         call flush(nop)
      end if

      ! Skip offset number of records
      do k=1,offset
         call zaiosk(nop)
      end do

      write(lp,'(a,i4,a,i4)')   'save_restart_rand: restart for ',
     &   imem,' starts at .a offset ',offset+1


      ! The rest ...
      call zaiowr(slp(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+1)
     &   'slp     ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'slp     ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(taux(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+2)
     &   'taux    ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'taux    ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(tauy(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+3)
     &   'tauy    ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'tauy    ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(wndspd(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+4)
     &   'wndspd  ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'wndspd  ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(airtmp(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+5)
     &   'airtmp  ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'airtmp  ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(relhum(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+6)
     &   'relhum  ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'relhum  ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(clouds(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+7)
     &   'clouds  ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'clouds  ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(precip(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+8)
     &   'precip  ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'precip  ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(sss   (1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+9)
     &   'sss     ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'sss     ',nstep0,dtime0,imem,xmin,xmax

      call zaiowr(sst   (1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      write(nop,117,rec=offsetb+10)
     &   'sst     ',nstep0,dtime0,imem,xmin,xmax
      write(6,117) 'sst     ',nstep0,dtime0,imem,xmin,xmax

      call zaiocl(nop)
      close(nop)
      save_restart_rand=.true.

      return

 117  format (a8,' =',i11,f11.2,i3,1p2e16.7)
 118  format ('itdm= ',i5, ' jtdm= ',i5,' kdm= ',i5)
      end function save_restart_rand

      end module mod_rand_io

