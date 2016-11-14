      module mod_evp_io

#if defined (EVP_IO_TAUICE)
      integer, private,parameter :: nfields=4
#else
      integer, private,parameter :: nfields=2
#endif

      real, dimension(:,:), save, allocatable ::
     &   uvel, vvel, tauxice, tauyice

      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine initevp()
      use mod_xc
      implicit none
      allocate(uvel    (idm,jdm))
      allocate(vvel    (idm,jdm))
      allocate(tauxice (idm,jdm))
      allocate(tauyice (idm,jdm))
      end subroutine

! Routine to read EVP fields

      integer function read_restart_evp(rungen,rt,imem,nstep0,dtime0)
      use mod_year_info
      use mod_za
      implicit none
      character(len=3)    , intent(in)  :: rungen
      type(year_info)     , intent(in)  :: rt
      integer,              intent(in)  :: imem
      integer,              intent(in)  :: nstep0
      real*8,               intent(in)  :: dtime0
c
      character(len=50) :: filename
      integer :: j,iostat,nop,ifld
      logical :: ex, exb, first
      character*80    l_ctitle(4)
      integer :: l_iversn, l_iexpt, l_yrflag,l_itdm, 
     &           l_jtdm, iyy, idd, ihh, iss
      character(len=50) :: filebase
      character(len=8) :: char8
      real :: xmin2, xmax2, l_time, xmin, xmax
      integer :: l_nstep,l_imem,offset,k,offsetb,ios,ifldb,nfields_input
      integer, dimension(idm,jdm) :: ip


      ! New output format -- uses hycom-style IO (zaio)
      nop=710
      filebase=rungen//'restart'//rt%cyy//'_'//rt%cdd//'_'//rt%chh//
     &         'EVP'
      inquire (exist=ex , file=trim(filebase)//'.a')
      inquire (exist=exb, file=trim(filebase)//'.b')
      if (.not.ex .or. .not.exb) then
         write(lp,*)'No such file '//trim(filebase)//'.a'
         read_restart_evp=1
         return
      else
         write(lp,'(a,i4)')'read_restart_evp: reading '
     &      //trim(filebase)//'.[ab] record=',imem
      end if

      call zaiopf(trim(filebase)//'.a','old',nop)
      open(nop,file=trim(filebase)//'.b',status='old',
     &         form='formatted',access='direct',recl=80)

      ! Number of fields for one member (only u/v)
      !KAL - TODO - must fix the offset/nfields logic when introducing new fields

      ! Offset for IO
      ! KAL -- offset not used anymore on input - search for correct
      ! KAL -- imem on input in stead
      offset=(imem-1)*nfields
      offsetb=offset+1 ! First record in .b  is model info

      ! Read "header"
      read(nop,118,rec=1) l_itdm,l_jtdm
      if (idm/=l_itdm .or. jdm /= l_jtdm) then
         if (mnproc==1) write(lp,*) 'Dimension mismatch in EVP file'
         call xcstop('(read_restart_evp)')
         stop '(read_restart_evp)'
      end if


      ! Read .b file until we find correct imem
      first=.true.
      ios=0
      ifldb=1 ! First record was read above
      nfields_input=0
      do while (ios==0) 

         read(nop,117,rec=ifldb+1,iostat=ios)
     &      char8,l_nstep,l_time,l_imem,xmin2,xmax2
         ifldb=ifldb+1

         ! Locate first field belonging to imem
         if (first) then
            if (l_imem==imem .and. ios==0) then
               offsetb=ifldb
               first=.false.
            end if
         end if

         ! Update num fields belonging to imem
         if (l_imem==imem .and. ios==0) nfields_input=nfields_input+1
         !print *,ifldb,nfields_input,ios
      end do

      ! Here we never found member imem
      if (first) then
         if (mnproc==1) write(lp,'(a,i3)') 'Could not find member ',imem
         read_restart_evp=1
         return
      end if

      ! Offset into ".a" - file
      offset=offsetb-1

      ! Adjust offset so that offset+1 is first record
      offset =offset -1
      offsetb=offsetb-1

      ! Info 
      write(lp,'(a,i4,a,i4)')   'read_restart_evp: restart for ',
     &   imem,' starts at .a offset ',offset+1
      write(lp,'(a,i4,a,i4,a)') 'read_restart_evp: restart for ',
     &   imem,' contains ',nfields_input,' model fields'
            


      ! Skip offset number of records in .a-file
      do k=1,offset
         call zaiosk(nop)
      end do


      ! Start reading at correct offsets
      do ifld=1,nfields_input

         ! Same order as out...
         read(nop,117,rec=offsetb+ifld,iostat=ios)
     &      char8,l_nstep,l_time,l_imem,xmin2,xmax2

         if (l_imem/=imem) then
            if (mnproc==1) write(lp,*) 'Ensemble member mismatch..'
            call xcstop('(read_restart_evp)')
            stop '(read_restart_evp)'

            ! TODO Change this with "cycle" ?
         end if

         if (mnproc==1)
     &   write(lp,*)char8,imem,xmin2,xmax2

         if     (trim(char8)=='uvel') then
            call zaiord(uvel(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='vvel') then
            call zaiord(vvel(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif     (trim(char8)=='tauxice') then
            call zaiord(tauxice(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         elseif (trim(char8)=='tauyice') then
            call zaiord(tauyice(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax, nop)
         else
            if (mnproc==1) print *, 'Unknown field '//char8
            call xcstop('(read_restart_evp')
            stop '(read_restart_evp'
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
            call xcstop('(read_restart_evp)')
            stop '(read_restart_evp)'
         end if


      end do

      call zaiocl(nop)
      close(nop)
      read_restart_evp=0


 117  format (a8,' =',i11,f11.2,i3,1p2e16.7)
 118  format ('itdm= ',i5, ' jtdm= ',i5)
      end function 



! Routine to save EVP fields

      integer function save_restart_evp(rungen,rt,imem,nstep0,dtime0)
      use mod_xc
      use mod_za
      use mod_year_info
      implicit none
      character(len=3)    , intent(in) :: rungen
      type(year_info)     , intent(in) :: rt
      integer,              intent(in) :: imem
      integer,              intent(in) :: nstep0
      real*8,               intent(in) :: dtime0
      character(len=50) :: filename
      integer :: j,iostat,ios
      integer ioproc, ierr
      
      integer :: nop,offset,k,l_nstep,l_imem,offsetb
      real :: xmin,xmax, l_time
      character(len=50) :: filebase
      character(len= 1) :: char1
      character(len= 8) :: char8
      integer, dimension(idm,jdm) :: ip


      save_restart_evp=0

      ! New output format -- uses hycom-style IO (zaio)
      nop=710
      filebase=rungen//'restart'//rt%cyy//'_'//
     &         rt%cdd//'_'//rt%chh//'EVP'
      if (mnproc==1) then
         write(lp,'(a,i4)',advance='yes') 'save_restart_evp: saving '//
     &         trim(filebase)//'.[ab] record=',imem
      end if
      call zaiopf(trim(filebase)//'.a','unknown',nop)
      if (mnproc==1) then
         open(nop,file=trim(filebase)//'.b',status='unknown',
     &        form='formatted',access='direct',recl=80)
      end if


      ! Offset for IO
      offset=(imem-1)*nfields
      offsetb=offset+1 ! First record is model info 
      if (mnproc==1) then
         write(nop,118,rec=1) idm,jdm
         call flush(nop)
      end if

      write(lp,'(a,i4,a,i4)')   'save_restart_evp: restart for ',
     &   imem,' starts at .a offset ',offset+1
      write(lp,'(a,i4,a,i4,a)') 'save_restart_evp: restart for ',
     &   imem,' contains ',nfields,' model fields'

      ! Skip offset number of records
      do k=1,offset
         call zaiosk(nop)
      end do

c --- Ready to write
      call zaiowr(uvel(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+1)
     &      'uvel    ',nstep0,dtime0,imem,xmin,xmax
         write(6,117)
     &      'uvel    ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(vvel(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+2)
     &      'vvel    ',nstep0,dtime0,imem,xmin,xmax
         write(6,117)
     &      'vvel    ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if


c --- Ready to write extra fields if defined
#if defined (EVP_IO_TAUICE)
      call zaiowr(tauxice(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+3)
     &      'tauxice ',nstep0,dtime0,imem,xmin,xmax
         write(6,117)
     &      'tauxice ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(tauyice(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+4)
     &      'tauyice ',nstep0,dtime0,imem,xmin,xmax
         write(6,117)
     &      'tauyice ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if
#endif

      call zaiocl(nop)
      if (mnproc==1)  close(nop)
      save_restart_evp=0

      return

 117  format (a8,' =',i11,f11.2,i3,1p2e16.7)
 118  format ('itdm= ',i5, ' jtdm= ',i5)


      end function 

      end module mod_evp_io


