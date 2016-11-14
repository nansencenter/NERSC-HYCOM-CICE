module mod_read_rstab
! This routine is slightly different to the one for daily fields.
! Modifications done May 2005 by Hanne Sagen
contains

   subroutine read_rstfield2d(fbase,char8,field,idm,jdm,coord,undef,maxcoord)
! Routine which read instantanous fields. Both a and b files are used in the reading process.
   use mod_za , only: zaiopf, zaiord, zaiocl, zaiosk
   implicit none
   integer,         intent(in)  :: idm,jdm,coord !i,j,k coordinates
   character(len=*),intent(in)  :: fbase
   real,            intent(out) :: field(idm,jdm)
   character(len=8), intent(in) :: char8
   integer, optional,intent(in) :: maxcoord
   real,             intent(in) :: undef

   character(len=5) :: char5
   character(len=7) :: char7
   character(len=8) :: cfld
   character(len=80) :: fullstring
   integer, dimension(idm,jdm) :: mask
   real :: xmin,xmax,rday,axmax,axmin,bxmin,bxmax
   integer :: lcoord,irec,nrec,nstep,dens
   integer :: nop,ios,maxco
   integer tlevel, layer
   integer count,numheadlines
   if (present(maxcoord)) then
      maxco=maxcoord
   else
      maxco=-1
   end if

   nop=777

   ! Open input file (b file)
       char5='' ; ios=0
       char7='RESTART'
       open(nop,file=trim(fbase)//'.b',status='old')
       count=0
!  Skip two first lines info is decoded in the subroutine which read the header
!Find number of headerlines and spool forward.
       do while (char7=='RESTART' .and. ios==0)
         read(nop,'(a7)',iostat=ios) char7
!         write(*,*) char7
         if (char7=='RESTART') then
         count=count+1
         endif
       end do
       numheadlines=count
!       print*,'Num header lines:', numheadlines
       close(nop)
       open(nop,file=trim(fbase)//'.b',status='old')
!       print *,'End b loop 1 ',ios,char7
       do count=1,numheadlines
       read(nop,'(a)',iostat=ios) fullstring 
       end do 
! Read until we get the field we want
       nrec=0
       cfld=''
       lcoord=-1

! --- CH170605:Problems in this one....       
   do while((cfld/=char8 .or. lcoord/=coord) .and. ios==0)
      read(nop,'(a)',iostat=ios) fullstring 
      !write(*,*) fullstring
      !read(fullstring,'(a8, 22x, i4,i3, 1pe16.7,1pe18.7)') cfld,lcoord,dens,bxmin,bxmax
      read(fullstring,'(a8, 22x, i4,i3)') cfld,lcoord,dens
      read(fullstring(38:80),*) bxmin,bxmax
!      print *,cfld,lcoord,dens,bxmin,bxmax
      !print *,'|'//fullstring(38:80)//'|'
      nrec=nrec+1 ! Counting number of lines not providing the field requested
 !      print *,nrec,ios
   end do

!Checkes for all 
!       select case(ios)
!       case (-1)
!            print*,'no ',char8,' data in this file !!'
!       case (0)
!            print*,char8,' data to be found after skipping  ',nrec,'data records'   
!!            print*, 'Level looked for  ',coord,'Level in file  ', lcoord  
!!  print *,'Number of lines not providing data ',nrec,ios,cfld,char8,lcoord,coord
!       end select

 if (cfld==char8 .and. lcoord==coord .and. ios==0) then
!      read(nop,'(a)',iostat=ios) fullstring  the first line with temp has already been read above...
!      read(fullstring,'(a8, 22x, i4,i3, 1pe16.7,1pe18.7)') cfld,lcoord,dens,bxmin,bxmax
      read(fullstring,'(a8, 22x, i4,i3)') cfld,lcoord,dens
      read(fullstring(38:80),*) bxmin,bxmax
      if (dens==1) then
          if (maxco>0) then
                if (coord==1) then
                   write(*,'(a,i3.3)',advance='no')  &
                      'Got field           "'//char8//'" at coordinate',coord
                else
                    write(*,'(a1)',advance='no') '.'
                end if

                if (coord==maxco) write(*,'(i3.3)') coord
           else
               write(*,'(a,i3.3)',advance='yes') &
               'Got field           "'//char8//'" at coordinate',coord
           end if !maxco

         call zaiopf(trim(fbase)//'.a','old',nop)
             do irec=1,nrec-1
              call zaiosk(nop)
             end do
         call zaiord(field,mask,.false.,axmin,axmax,nop)
!KAL         if (cfld == 'temp    ') then
!KAL      !     write(*,*)'Info from bfile - data available'
!KAL      !     write(*,*) fullstring
!KAL      !     write(*,'(a8, 22x, i4,i3, 1p2e16.7)')  cfld,lcoord,dens,bxmin,bxmax
!KAL      !     write(*,*) ' *.a file info:   cfld,lcoord,dens,axmin,axmax'
!KAL           !write(*,'(a8, 22x, i4,i3, 1p2e16.7)') cfld,lcoord,dens,axmin,axmax
!KAL           print*, cfld,lcoord,dens,axmin,axmax
!KAL      !      write(*,*)cfld,bxmin,bxmax,axmin,axmax,dens,lcoord
!KAL         endif
      call zaiocl(nop)
!      print *,minval(field),maxval(field)
  end if ! dens
 else
      print '(a,i3.3)', 'Could not get field "'//char8//'" at coordinate',coord
      field=undef
 end if ! cfld

   close(nop)
 !117  format (a8,': layer,tlevel,range = ',i3,f11.2,i3,f7.3,1p2e16.7)
! 117  format (a8,': layer,tlevel,range = ',i3,i1,1p2e16.7,1p2e16.7)
! 118  format ('member ',i5.5,' = ',i1,' Ensemble member flag')
   end subroutine read_rstfield2d



   subroutine read_rstfield3d(fbase,char8,field,idm,jdm,kdm,undef)
   implicit none
   integer, intent(in) :: idm,jdm,kdm
   character(len=*), intent(in) :: fbase
   character(len=8), intent(in) :: char8
   real, intent(out)    :: field(idm,jdm,kdm)
   real, intent(in) :: undef

   integer :: k

   do k=1,kdm
      call read_rstfield2d(fbase,char8,field(:,:,k),idm,jdm,k,undef,maxcoord=kdm)
   end do
   end subroutine read_rstfield3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rst_read_header(fbase,rtdump,nrmem,idm,jdm,kdm)
   use mod_year_info
   use m_year_day
   use m_datetojulian
   use m_read_blkdat
   implicit none
   character(len=*), intent(in)  :: fbase
   type(year_info) , intent(out) :: rtdump
   integer,          intent(out) :: nrmem,idm,jdm,kdm

   character(len=80) :: fullstring
   character(len=24) :: first
   character (len=1) :: third
   integer  nstep,nop,yrflag,iexpt,iversn,idtime
   integer year,day,month
   real  :: dtime,hhd24
   nrmem=1
   nop=777

! Open input file
   open(nop,file=trim(fbase)//'.b',status='old')

! Read  first line and decode.
     read(nop,'(31x,3i6)') iexpt,iversn,yrflag    
     print *,iexpt,iversn,yrflag

     read(nop,'(a)') fullstring !Whole line into one string    
     print *,fullstring
     read(fullstring(25:len_trim(fullstring)),*) nstep,dtime ! x skipper en posisjon
     !print *,'second fourth',nstep,dtime

!    write(*,'(a,f9.1)') 'Julian day from header file :',fourth
! Calculate year-date infor from Julian day.
    idtime=floor(dtime)
    hhd24=dtime-floor(dtime)
    !rtdump%ihh= floor((dtime-floor(dtime))*24)
    call juliantodate(idtime,year,month,day,1901,1,1)
! Set rtdump
   !print *,dtime,year,month,day
    rtdump%iyy = year
    rtdump%imm = month
    rtdump%idm = day
    rtdump%idd= datetojulian(rtdump%iyy,rtdump%imm,rtdump%idm,year,1,1)
! Retrive info about grid from blkdat.dat
    call read_gridindex(idm,jdm,kdm)

   if (yrflag==3) then
      call year_day(float(rtdump%idd)+hhd24,rtdump%iyy,rtdump,'ecmwf')
!      call year_day(float(rtinit%idd),rtinit%iyy,rtinit,'ecmwf')
   elseif (yrflag==0) then
      call year_day(float(rtdump%idd)+hhd24,rtdump%iyy,rtdump,'month')
!      call year_day(float(rtinit%idd),rtinit%iyy,rtinit,'month')
   else
      print *,'Unknow yrflag --- ',yrflag
      stop 
   end if

   print *,'Day and year info for the restart field:            ', &
      rtdump%cyy//'-'//rtdump%cmm//'-'//rtdump%cdm
   close (nop)


    end subroutine rst_read_header

      

end module
