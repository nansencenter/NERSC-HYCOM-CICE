      module mod_ice_io


      real*8, dimension(:,:), allocatable, save  ::  ! Hard-coded real*8 for file IO
     &   ficem,hicem,hsnwm,ticem,tsrfm,iceU,iceV


      contains

      subroutine initice()
      use mod_xc
      implicit none
      allocate(ficem (idm,jdm))
      allocate(hicem (idm,jdm))
      allocate(hsnwm (idm,jdm))
      allocate(ticem (idm,jdm))
      allocate(tsrfm (idm,jdm))
      allocate(iceU  (idm,jdm))
      allocate(iceV  (idm,jdm))
      end subroutine

      function read_restart_ice(rungen,rt,imem)
      use mod_xc
      use mod_year_info
      implicit none
      logical read_restart_ice
      type(year_info) , intent(in) :: rt
      character(len=3), intent(in) :: rungen
      integer, intent(in) :: imem
c
      character(len=27) icename
      logical ex
      integer i,j,iostat
      real x(10)

      read_restart_ice=.true.
      icename=fileice(rungen,rt)
      if(mnproc==1) write(lp,'(a,i4)',advance='yes') 
     &   'reading: '//trim(icename)//
     & '                     record=',imem
      inquire(file=icename,exist=ex)
      if (ex) then
         inquire(iolength=j)
     &      ficem,hicem,hsnwm,ticem,tsrfm  !,iceU,iceV
         open(10,file=icename,status='old',form='unformatted',
     &           access='direct',recl=j)
            iceU=0.0; iceV=0.0
            read(10,rec=imem,err=200,iostat=iostat)ficem,hicem,
     &           hsnwm,ticem,tsrfm !,liceU,liceV
  200       continue
         close(10)
         if (iostat/=0) then
            if (mnproc==1) then
            write(lp,*) 'An error occured while reading '//trim(icename)
            write(lp,*) 'IOSTAT is :',iostat
            end if
            read_restart_ice=.false.
         end if
      else 
         if (mnproc==1) then
            write(lp,*) 'Restart file for ICE does not exist.'
         end if
         read_restart_ice=.false.
      endif
      end function read_restart_ice



      function save_restart_ice(rungen,rt,imem)
      use mod_xc
      use mod_year_info
      implicit none

      type(year_info),       intent(in)  :: rt
      character(len=3),      intent(in) :: rungen
      integer,               intent(in) :: imem
      logical :: save_restart_ice
      character(len=27) icename
      integer i,j
      integer iostat
      integer ioproc, ierr
      

      ioproc=mnproc 


      ! Only one node is alloved to access this file (only important for MPI)
      save_restart_ice=.true.
      if (mnproc==ioproc) then
         icename=fileice(rungen,rt)
         write(lp,'(a,i4)',advance='yes')'Saving : '//
     &      trim(icename)//'                     record=',imem
         inquire(iolength=j)
     &      ficem,hicem,hsnwm,ticem,tsrfm  !,ticeU,ticeV
         open(10,file=icename,status='unknown',form='unformatted',
     &      access='direct',recl=j)
         write(10,rec=imem,err=200,iostat=iostat)ficem,hicem,
     &      hsnwm,ticem,tsrfm   !,io_iceU,io_iceV
200      continue
         close(10)
         if (iostat==0) then
            save_restart_ice=.true.
         else
            write(lp,'(a,i4)')
     &         'An Error encountered while writing record ',imem
            save_restart_ice=.false.
            write(lp,'(a,i4)') 'IOSTAT is :',iostat
         end if
      end if

      ! KAL - new: bail out on save errors.. Prevent disaster in the
      ! making
      if (iostat/=0) then
         if (mnproc==1) then
            write(lp,*) 'Exiting on iostat save  errors...'
            write(lp,*) 'IOSTAT is ',iostat
         end if
         call xcstop('(save_restart_ice)')
         stop '(save_restart_ice)'
      end if

      end function save_restart_ice

      function fileice(rungen,rt)
      use mod_year_info
      implicit none
         character(len=27) fileice
         character(len=3) rungen
         type(year_info) rt
         fileice(1:3)=rungen(1:3)
         fileice(4:10)='restart'
         fileice(11:14)=rt%cyy(1:4)
         fileice(15:15)='_'
         fileice(16:18)=rt%cdd(1:3)
         fileice(19:19)='_'
         fileice(20:21)=rt%chh(1:2)
         fileice(22:27)='ICE.uf'
      end function fileice


      end module mod_ice_io
