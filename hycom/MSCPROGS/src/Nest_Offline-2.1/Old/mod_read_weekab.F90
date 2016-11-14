module mod_read_weekab
! This routine is slightly different to the one for daily fields.
! Modifications done May 2005 by Hanne Sagen
contains

   subroutine read_weekfield2d(fbase,char8,field,idm,jdm,coord,undef,maxcoord,lsilent)
! Routine which read instantanous fields. Both a and b files are used in the reading process.
   use mod_za , only: zaiopf, zaiord, zaiocl, zaiosk
   implicit none
   integer,         intent(in)  :: idm,jdm,coord !i,j,k coordinates
   character(len=*),intent(in)  :: fbase
   real,            intent(out) :: field(idm,jdm)
   character(len=8), intent(in) :: char8
   integer, optional,intent(in) :: maxcoord
   real,             intent(in) :: undef
   logical, optional,intent(in) :: lsilent

   character(len=5) :: char5
   character(len=7) :: char7
   character(len=8) :: cfld
   character(len=80) :: fullstring
   integer, dimension(idm,jdm) :: mask
   real :: xmin,xmax,rday,axmax,axmin,bmin,bmax
   integer :: irec,nrec,nstep
   integer :: nop,ios,maxco
   integer tlevel, layer
   integer numheadlines
   integer lkdm,lidm,ljdm,yrflag,i,avecount
   real :: ldens, dtime
   logical :: lsilent2

   if (present(lsilent)) then
      lsilent2=lsilent
   else
      lsilent2=.false.
   end if


   if (present(maxcoord)) then
      maxco=maxcoord
   else
      maxco=-1
   end if

   nop=777

   ! Open input file (b file)
   ios=0
   cfld=''          
   open(nop,file=trim(fbase)//'.b',status='old')

   ! Skip 6 lines 
   do i=1,6
      read(nop,*)
   end do

   ! Read year flag, idm, jdm, kdm
   read(nop,*) yrflag
   read(nop,*) lidm
   read(nop,*) ljdm
   read(nop,*) lkdm
   
   ! Check dim match
   if (idm/=lidm .or. jdm/=ljdm) then
      print *,'Dimension mismatch '
      print *,idm,lidm
      print *,jdm,ljdm
      stop '(read_weekfield2d)'
   end if

   ! Skip 2 lines, read ave count, skip again
   read(nop,*) 
   read(nop,*) 
   read(nop,*) avecount !; print *,'avecount',avecount
   read(nop,'(a80)') fullstring
   !print *,trim(fullstring)

   ! We are now at the beginning of the records. each line is one record in .a
   ! file
   irec=0
   cfld=''
   ios=0
   do while((cfld/=char8 .or. layer/=coord) .and. ios==0)
      irec=irec+1
      read(nop,117,iostat=ios) cfld,nstep,dtime,layer,ldens,bmin,bmax
      !if (trim(char8)=='pres') &
      !print *,irec,cfld,char8,layer,coord
   end do
   nrec=irec
   close(nop)

   if (trim(cfld)==trim(char8) .and. layer==coord .and. ios==0) then

      call zaiopf(trim(fbase)//'.a','old',nop)
      do irec=1,nrec-1
         call zaiosk(nop)
      end do
      call zaiord(field,mask,.false.,axmin,axmax,nop)
      field=field/avecount
      !print *,cfld,layer,maxval(field),minval(field)
      call zaiocl(nop)
      if     (abs(axmin-bmin).gt.abs(bmin)*1.e-4 .or. &
              abs(axmax-bmax).gt.abs(bmax)*1.e-4     ) then
          write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
            'error - .a and .b files not consistent:', &
            '.a,.b min = ',axmin,bmin,axmin-bmin, &
            '.a,.b max = ',axmax,bmax,axmax-bmax
          stop '(read_weekfield2d)'
      endif


      if (maxco>0) then
         if (coord==1) then
            write(*,'(a,i3.3)',advance='no')  &
               'Got field           "'//char8//'" at coordinate',coord
         else
            write(*,'(a1)',advance='no') '.'
         end if
         if (coord==maxco) write(*,'(i3.3)') coord
      else if (.not. lsilent2) then
         write(*,'(a,i3.3)',advance='yes') &
         'Got field           "'//char8//'" at coordinate',coord
      end if

   else
      print '(a,i3.3)', 'Could not get field "'//char8//'" at coordinate',coord
      field=undef
   end if ! cfld

 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
   end subroutine read_weekfield2d



   subroutine read_weekfield3d(fbase,char8,field,idm,jdm,kdm,undef)
   implicit none
   integer, intent(in) :: idm,jdm,kdm
   character(len=*), intent(in) :: fbase
   character(len=8), intent(in) :: char8
   real, intent(out)    :: field(idm,jdm,kdm)
   real, intent(in) :: undef
   real, dimension(idm,jdm) :: dpfield

   integer :: k

   do k=1,kdm
      call read_weekfield2d(fbase,char8,field(:,:,k),idm,jdm,k,undef,maxcoord=kdm)
   end do
   end subroutine read_weekfield3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine week_read_header(fbase,rtdump,idm,jdm,kdm)
   use mod_year_info
   use m_year_day
   use m_datetojulian
   implicit none
   character(len=*), intent(in)  :: fbase
   type(year_info) , intent(out) :: rtdump
   integer,          intent(out) :: idm,jdm,kdm

   character(len=80) :: fullstring
   character(len=24) :: fiweek
   character (len=1) :: third
   character (len=8) ::char8
   integer  nstep,nop,yrflag,iexpt,iversn
   integer year,day,month
   integer :: i,ios,avecount
   real  :: dtime,coord,dens,bmin,bmax

   nop=777

   ! Open input file (b file)
   ios=0
   open(nop,file=trim(fbase)//'.b',status='old')

   ! Skip 6 lines 
   do i=1,6
      read(nop,*)
   end do

   ! Read year flag, idm, jdm, kdm
   read(nop,*) yrflag
   read(nop,*) idm
   read(nop,*) jdm
   read(nop,*) kdm
   
   ! Skip 2 lines, read ave count, skip again
   read(nop,*) 
   read(nop,*) 
   read(nop,*) avecount
   read(nop,*) 

   ! We are now at the beginning of the records. each line is one record in .a
   ! file
   read(nop,117,iostat=ios) char8,nstep,dtime,coord,dens,bmin,bmax 
   ! Only need one read
   close(nop)


!    write(*,'(a,f9.1)') 'Julian day from header file :',fourth
! Calculate year-date infor from Julian day.
    rtdump%ihh= floor((dtime-floor(dtime))*24)
    call juliantodate(floor(dtime),year,month,day,1901,1,1)
! Set rtdump
    rtdump%iyy = year
    rtdump%imm = month
    rtdump%idm = day
    rtdump%idd= datetojulian(rtdump%iyy,rtdump%imm,rtdump%idm,year,1,1)

   if (yrflag==3) then
      call year_day(real(rtdump%idd,kind=8),rtdump%iyy,rtdump,'ecmwf')
   elseif (yrflag==0) then
      call year_day(real(rtdump%idd,kind=8),rtdump%iyy,rtdump,'month')
   else
      print *,'Unknow yrflag --- ',yrflag
      stop 
   end if

   print *,'Day and year info for the weekly field:            ', &
      rtdump%cyy//'-'//rtdump%cmm//'-'//rtdump%cdm
   close (nop)

 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)

    end subroutine week_read_header

   !Reads different info from blkdat.input 
   subroutine read_gridindex(idm,jdm,kdm)
      implicit none

      ! Dummy variables
      integer, intent(out) :: idm,jdm,kdm
      integer    dump1,dump2 
      integer iversn,iexpt
      character(len=80) :: file1
      logical :: ex
      character(len=80) :: ctitle(4)


      ! Check for existence of blkdat.input file
      file1='blkdat.input'
      inquire(file=trim(file1),exist=ex)
   if (.not.ex) then
      print *,'File blkdat.input does not exist'
      stop
   end if
 
    

open(10,file=file1,form='formatted',status='old')
 read(10,116) ctitle,iversn,iexpt,idm,jdm,dump1,dump2,kdm

close (10)
!116 format (a80/a80/a80/a80/ &
!      i4,4x,'''iversn'' = hycom version number x10'/  &
!      i4,4x,'''iexpt '' = experiment number x10'/  &
!      i4,4x,'''idm   '' = longitudinal array size'/  &
!      i4,4x,'''jdm   '' = latitudinal  array size'/  &
!      i4,4x,'''itest '' = Year of integration start '/  &
!      i4,4x,'''jtest '' = Day of integration start'/  &
!      i4,4x,'''kdm   '' = Vertical     array size')  
116  format (a80/a80/a80/a80/ &
       i4,4x,35x/ &
       i4,4x,32x/  &
       i4,4x,34x/  &
       i4,4x,34x/  &
       i4,4x,37x/  &
       i4,4x,35x/  &
       i4,4x,34x)  

   end subroutine read_gridindex

   subroutine read_thbase(thbase)
      implicit none

      ! Dummy variables
      real, intent(out) :: thbase
      integer    dump1,dump2 
      integer iversn,iexpt
      character(len=80) :: file1
      logical :: ex
      character(len=80) :: ctitle(22)


      ! Check for existence of blkdat.input file
      file1='blkdat.input'
      inquire(file=trim(file1),exist=ex)
   if (.not.ex) then
      print *,'File blkdat.input does not exist'
      stop
   end if
 
    

open(10,file=file1,form='formatted',status='old')
 read(10,116) ctitle,thbase
print*,'thbase inside read_thbase',thbase
close (10)
116  format (a80/a80/a80/a80/a80/a80/a80/a80/a80/a80/a80/ &
             a80/a80/a80/a80/a80/a80/a80/a80/a80/a80/a80/ &
       2x,f3.1,3x,'''thbase'' = reference density (sigma units)') 

   end subroutine read_thbase



end module
