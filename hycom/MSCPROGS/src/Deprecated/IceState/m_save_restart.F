      module m_save_restart
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function save_restart_mem(rungen,rt,imem,nstep0,dtime0)
      use mod_xc
      use mod_za
      use mod_hycomfuncs
      use m_year_day
      implicit none

      logical :: save_restart_mem
      type(year_info),  intent(in) :: rt
      character(len=3), intent(in) :: rungen
      integer,          intent(in) :: imem
      integer,          intent(in) :: nstep0
      real*8,           intent(in) :: dtime0

      character(len=24) filename 
      character(len=48) filenew
      character(len=40) botfile
      integer i,j,iostat,indx
      character(len=9), parameter :: cident='HYCOM_1.2'
      character(len=9) rident
      character(len=3) cmem
      real, dimension(nx,ny) :: gpbot,gthkk,gpsikk
      integer ioproc, ierr

      ! NB get stop - signals valid for MPI runs as well !

      include 'common_blocks.h'
      
#if defined(MPI)
      include 'mpif.h'
#endif


      filename=filemem(rungen,rt)

      ! use hycom restart routine. Modify restart filename
      indx=index(filename,'.uf')
      filenew=filename(1:indx-1)
      ! If imem > 1 - append member index ...
      if (imem>1) then
         write(cmem,'(i3.3)') imem
         filenew=filenew(1:len_trim(filenew))//'_mem'//cmem
      end if
      print *,filenew
      flnmrso=trim(filenew)
      call restart_out(nstep0,dtime0)

      ! Hmm .. not closed automatically ...
      call zaiocl(12)
      if (mnproc==1)  then
         write(lp,*) 'Closing up new restart files...'
         close(11)
      end if


      save_restart_mem=.true.
      return
      end function save_restart_mem





      function save_restart_ice(rungen,rt,imem)
      use mod_xc
      use mod_hycomfuncs
      use mod_year_info
      use mod_common_ice
      implicit none

      type(year_info),       intent(in)  :: rt
      character(len=3),      intent(in) :: rungen
      integer,               intent(in) :: imem


      ! Local arrays for file-saving
      real, dimension(nx,ny)  :: 
     &   tficem,thicem,thsnwm,tticem,ttsrfm,ticeU,ticeV
      real*8, dimension(nx,ny)  ::  ! Hard-coded real*8 for file IO
     &   io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm,io_iceU,io_iceV

      logical :: save_restart_ice

      character(len=27) icename
      integer i,j
#if ! defined (SHMEM) 
      integer iostat
#endif
      integer ioproc, ierr
      
#if defined(MPI)
      include 'mpif.h'
#elif defined (SHMEM)
      include 'mpp/shmem.fh'
      INTEGER, save :: PSYNC(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER, save :: PWRK(MAX(1/2+1,SHMEM_REDUCE_MIN_WRKDATA_SIZE))
      integer, save :: isaveiceall, isaveice, alliostat, iostat
#endif


      ioproc=mnproc 
#if defined(MPI)
      ! Broadcast "io" mnproc (ioproc) from node with rank 0
      call mpi_bcast(ioproc,1,MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
#endif

      ! NB get stop - signals valid for MPI runs as well !

      ! Put ice data into non-tiled arrays
      call xcaget(tficem,ficem,0) ; io_ficem=tficem
      call xcaget(thicem,hicem,0) ; io_hicem=thicem
      call xcaget(thsnwm,hsnwm,0) ; io_hsnwm=thsnwm
      call xcaget(tticem,ticem,0) ; io_ticem=tticem
      call xcaget(ttsrfm,tsrfm,0) ; io_tsrfm=ttsrfm
      call xcaget(ticeU ,iceU ,0) ; io_iceu =ticeu 
      call xcaget(ticeV ,iceV ,0) ; io_icev =ticev

      ! Only one node is alloved to access this file (only important for MPI)
      save_restart_ice=.true.
      if (mnproc==ioproc) then
         icename=fileice(rungen,rt)
         write(lp,'(a,i4)',advance='yes')'Saving : '//
     &      trim(icename)//'                     record=',imem
         !inquire(iolength=j)tficem,thicem,thsnwm,tticem,ttsrfm  !,ticeU,ticeV
         inquire(iolength=j)
     &      io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm  !,ticeU,ticeV
         !print *,'iolength=',j ! Not important
         open(10,file=icename,status='unknown',form='unformatted',
     &      access='direct',recl=j)
         write(10,rec=imem,err=200,iostat=iostat)io_ficem,io_hicem,
     &      io_hsnwm,io_ticem,io_tsrfm   !,io_iceU,io_iceV
200      continue
         close(10)
         if (iostat==0) then
            !write(lp,*) '  ICE done'
            save_restart_ice=.true.
         else
            write(lp,'(a,i4)')
     &         'An Error encountered while writing record ',imem
            save_restart_ice=.false.
            write(lp,'(a,i4)') 'IOSTAT is :',iostat
         end if
      end if

      ! KAL : Should communicate result to all processes in MPI...
#if defined(MPI)
      call mpi_bcast(save_restart_ice,1,MPI_LOGICAL, 0,
     &               MPI_COMM_WORLD,ierr)
      call mpi_bcast(iostat,1,MPI_INTEGER, 0,  MPI_COMM_WORLD,ierr)
#elif defined (SHMEM)
      isaveice=0
      if (save_restart_ice) isaveice=1
      call SHMEM_INT8_AND_TO_ALL(isaveiceall,isaveice,1,0,0,N$PES, 
     &                           PWRK, PSYNC)
      save_restart_ice=isaveiceall==1
      call shmem_barrier_all()
      call shmem_int8_sum_to_all(alliostat,iostat,0,0,N$PES,PWRK,PSYNC)
      iostat=alliostat
#endif

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





#if defined (WEEKLY_AVERAGE)
      function save_restart_ave(rungen,rt,imem)
      use mod_xc
      use mod_hycomfuncs
      use mod_year_info
      use mod_average
      implicit none

      logical :: save_restart_ave
      type(year_info),  intent(in) :: rt
      character(len=3), intent(in) :: rungen
      integer,          intent(in) :: imem

      !type(average_state_global) :: gave_week

      character(len=27) avename
      character(len=60) tmpchar
      integer i,j,find
      integer ierr
      
#if defined(MPI)
      include 'mpif.h'
#endif

      ! Dump with new version
      tmpchar=fileaveres(rungen,rt)
      if (mnproc==1) 
     &   write(lp,'(a)',advance='no') 'save_restart_ave: saving:'//
     &      trim(tmpchar)//'.[ab]'

      find=index(tmpchar,'.',back=.true.)
      call save_ave2(tmpchar(1:find-1),ave_week)
      if (mnproc==1)
     &    write(lp,'(a,i5,a)') '   counter= ',ave_week%counter, 
     &      '   AVE done'
         
      ! Always true .. for now ..
      save_restart_ave=.true.
      end function save_restart_ave
#endif


! Routine to save random forcing fields
! 22.11.2004 - KAL: Changed to use HYCOM type io

      logical function save_restart_rand(rungen,rt,ran,imem,
     &   nstep0,dtime0)
      use mod_forcing_nersc
      use mod_year_info
      use mod_za
      implicit none
      character(len=3)    , intent(in) :: rungen
      type(year_info)     , intent(in) :: rt
      type(forcing_fields), intent(inout) :: ran
      integer,              intent(in) :: imem
      integer,              intent(in) :: nstep0
      real*8,               intent(in) :: dtime0
      character(len=50) :: filename
      integer :: j,iostat,ios
      integer ioproc, ierr
      
#if defined(MPI)
      include 'mpif.h'
#endif

      integer :: nop,nfields,offset,k,l_nstep,l_imem,offsetb
      real :: xmin,xmax, l_time
      character(len=50) :: filebase
      character(len= 1) :: char1,carret
      character(len= 8) :: char8

      include 'common_blocks.h'


      ! New output format -- uses hycom-style IO (zaio)
      nop=710
      filebase=rungen//'restRAND'//rt%cyy//'_'//rt%cdd//'_'//rt%chh
      write(lp,'(a,i4)',advance='yes') 'save_restart_rand: saving '//
     &   trim(filebase)//'.[ab] record=',imem
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
         write(nop,118,rec=1) itdm,jtdm,kdm
         call flush(nop)
      end if

      ! Skip offset number of records
      do k=1,offset
         call zaiosk(nop)
      end do


      ! The rest ...
      call zaiowr(ran%slp(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+1)
     &      'slp     ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%taux(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+2)
     &      'taux    ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%tauy(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+3)
     &      'tauy    ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%wndspd(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+4)
     &      'wndspd  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%airtmp(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+5)
     &      'airtmp  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%relhum(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+6)
     &      'relhum  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%clouds(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+7)
     &      'clouds  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%precip(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+8)
     &      'precip  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%sss   (1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+9)
     &      'sss     ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%sst   (1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+10)
     &      'sst     ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiocl(nop)
      if (mnproc==1)  close(nop)
      save_restart_rand=.true.

      return

 117  format (a8,' =',i11,f11.2,i3,1p2e16.7)
 118  format ('itdm= ',i5, ' jtdm= ',i5,' kdm= ',i5)


      end function save_restart_rand


      end module m_save_restart
