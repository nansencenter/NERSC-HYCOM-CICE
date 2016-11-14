module mod_read_gp

   ! TODO - the below approach will allow only file to be processed at a time
   ! Fix by using types if necessary...

   character(len=12),allocatable, dimension(:), save :: varname
   integer          ,allocatable, dimension(:), save :: vardim
   integer, allocatable, save :: varstop(:),varstart(:)
   integer,save :: nvars
   integer,save :: numr4
   integer,save :: kdm

   ! Holders of data read
   real,save      :: fyear
   integer,save   :: iyear, imonth, iday, ihour
   integer*4,save :: ix,jx,itime
   real*4           ,allocatable, dimension(:), save :: r4s
   logical, save  :: read_success

   interface readvar
      module procedure readvar0d,readvar1d
   end interface
   
contains



   ! Read gpheader and read variable names and dimensions. Return
   ! These along with number of variables.
   !
   ! NB:The variables in the gp header constitute one record (time) in the gp files
   subroutine read_gpheader(rungen)
   implicit none
   character(len=3),intent(in) :: rungen

   integer :: ivar,ios,i
   logical :: ex
   character(len=12) :: c12dummy
   character(len=80) :: gphdr

   ! Get number of variables in one gp record
   gphdr=rungen//'gprecord.asc'
   inquire(exist=ex,file=rungen//'gprecord.asc')
   if (.not.ex) then
      print *,'No gp header file '//trim(gphdr)
      stop '(read_gpheader)'
   end if
   open(10,file=rungen//'gprecord.asc', status='old')
   read(10,*) ! Skip first line
   ios=0
   nvars=0
   do while (ios==0)
      read(10,'(a12,a3)',iostat=ios) c12dummy,i
      if (ios==0) nvars=nvars+1
   end do

   ! Read variable names and dimensions into module arrays
   allocate(varname(nvars))
   allocate(vardim (nvars))
   rewind(10)
   read(10,*) ! Skip first line
   kdm=-100
   do ivar=1,nvars
      read(10,'(a12,i3)') varname(ivar),vardim(ivar)
      kdm=max(kdm,vardim(ivar))
   end do
   close(10)


   ! Allocate index olders for where data from a variable stops and ends.
   ! (Points into values in real4 array read from the record)
   allocate(varstart(nvars))
   allocate(varstop (nvars))

   ! Count number of variables in one record. Three first is integer, the rest
   ! are real*4s
   numr4=0
   print *,'One gp record contains the following data:'
   print *,' ivar ',' var name    ',' start','  stop '
   do i=1,nvars
      if (trim(varname(i))/='itime' .and. &
          trim(varname(i))/='i'     .and. &
          trim(varname(i))/='j' ) then
         varstart(i)=numr4+1
         numr4=numr4+vardim(i)
         varstop (i)=numr4
         print '(i6,"  ",a12,2i6)',i,varname(i),varstart(i),varstop(i)
      end if
   end do

   ! Allocate variable holding real 4 vars
   allocate(r4s(numr4))
   end subroutine



   ! Getvarind returns index in varname which matches "vartofind"
   integer function getvarind(vartofind)
   implicit none
   character(len=*),intent(in) :: vartofind
   integer :: ivar
   getvarind=-1
   do ivar=1,nvars
      if (trim(vartofind)==trim(varname(ivar))) then
         getvarind=ivar
      end if
   end do
   end function

   ! Routine reads one gprecord into variable r4s, itime, ix,jx and time info
   ! fields. NB: No data is returned here
   subroutine read_gprecord(filename,irec,ios)
   use mod_year_info
   implicit none
   character(len=*), intent(in) :: filename
   integer, intent(in ) :: irec
   integer, intent(out) :: ios

   integer, parameter :: nop=710
   integer :: diy,diy_now, iol

   inquire(iolength=iol) itime, ix, jx, r4s
!  open(nop,file=trim(filename),status='old',          &
!       form='unformatted',access='direct',iostat=ios, &
!       recl=iol)
   open(nop,file=trim(filename),status='old',          &
        form='unformatted',access='direct',iostat=ios, &
        recl=iol)
   read(nop,rec=irec,iostat=ios) itime,ix,jx,r4s

   ! Succesful read, and record contains meaningful data
   if (ios==0.and.itime>0) then

      iyear=itime/1e6
      imonth=(itime-iyear*1000000)/10000
      iday=(itime-iyear*1000000-imonth*10000)/100
      ihour=(itime-iyear*1000000-imonth*10000-iday*100)

      ! "Floating point year"  from time info
      if (iyear<1000) then ! climatology I presume?
         fyear=(imonth-1)*30 + iday + ihour/24.
         fyear=iyear+fyear/360.
      else ! not climatology I presume?
         diy=datetojulian(iyear+1,1,1,iyear,1,1)
         diy_now=datetojulian(iyear,imonth,iday,iyear,1,1)
         fyear=float(diy_now)/float(diy)
         fyear=fyear+iyear+float(ihour)/(24.*diy)
      end if

      read_success=.true.

   ! Succesful read, but record contains no meaningful data
   elseif (ios==0) then
      read_success=.false.

   ! Unsuccesful read
   else
      itime=-1
      read_success=.false.
   end if

   ! TODO: always keep open..
   close(nop) 
   end subroutine


   ! Routine returns data kept in r4s - scalar version
   function readvar0d(varname)
   implicit none
   character(len=*), intent(in) :: varname
   real readvar0d
   integer :: istart, istop, varind

   ! Return index in variable list
   varind=getvarind(varname)
   if (varind==-1) then
      print '(a)','Can not find variable '//varname
      stop '(readvar0d)'
   end if
   if (vardim(varind)/=1) then
      print '(a)','Sought 0d variable '//varname//' but it is listed as 1d'
      stop '(readvar0d)'
   end if
   if (.not.read_success) then
      print '(a)','Last read was unsuccesful, yet you try to assign....'
      stop '(readvar0d)'
   end if

   ! Get start and stop index of variable data
   istart=varstart(varind)
   istop =varstop (varind)
   readvar0d=r4s(istart)
   end function


   ! Routine returns data kept in r4s - vector version
   function readvar1d(varname,nz)
   implicit none
   character(len=*), intent(in) :: varname
   integer         , intent(in) :: nz
   real readvar1d(nz)
   integer :: istart, istop, varind

   ! Return index in variable list
   varind=getvarind(varname)
   if (varind==-1) then
      print '(a)','Can not find variable '//varname
      stop '(readvar1d)'
   end if
   if (vardim(varind)/=nz) then
      print '(a)','Sought 1d variable '//varname//' has ',vardim(varind), &
                  ' entries, but call requested ',nz, ' entries'
      stop '(readvar1d)'
   end if
   if (nz/=kdm) then
      print '(a)','Sought 1d variable '//varname//' with ',nz , &
                  ' entries, but module indicates ',kdm, ' entries'
      stop '(bug? readvar1d)'
   end if
   if (.not.read_success) then
      print '(a)','Last read was unsuccesful, yet you try to assign....'
      stop '(readvar1d)'
   end if

   ! Get start and stop index of variable data
   istart=varstart(varind)
   istop =varstop (varind)
   readvar1d=r4s(istart:istop)
   end function

   end module
      

