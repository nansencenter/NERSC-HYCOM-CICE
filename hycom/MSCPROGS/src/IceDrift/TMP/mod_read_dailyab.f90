module mod_read_dailyab

contains

   subroutine read_field2d(fbase,char8,field,idm,jdm,coord,undef,maxcoord)
   use mod_za , only: zaiopf, zaiord, zaiocl, zaiosk
   implicit none
   integer,         intent(in)  :: idm,jdm,coord
   character(len=*),intent(in)  :: fbase
   real,            intent(out) :: field(idm,jdm)
   character(len=8), intent(in) :: char8
   integer, optional,intent(in) :: maxcoord
   real,             intent(in) :: undef

   character(len=5) :: char5
   character(len=8) :: cfld
   integer, dimension(idm,jdm) :: mask
   real :: xmin,dens,xmax,rday
   integer :: lcoord,irec,nrec,nstep
   integer :: nop,ios,maxco

   if (present(maxcoord)) then
      maxco=maxcoord
   else
      maxco=-1
   end if

   nop=777

   ! Open input file
   char5='' ; ios=0
   open(nop,file=trim(fbase)//'.b',status='old')
   do while (char5/='field' .and. ios==0)
      read(nop,'(a5)',iostat=ios) char5
      !write(*,*) char5
   end do
   !print *,'End b loop 1 ',ios,char5

   ! Read until we get the field we want
   nrec=0
   cfld=''
   lcoord=-1
   do while((cfld/=char8 .or. lcoord/=coord) .and. ios==0)
      read(nop,117,iostat=ios) cfld,nstep,rday,lcoord,dens,xmin,xmax
      !write(*,*) cfld,char8,lcoord,coord,ios
      nrec=nrec+1
      !print *,nrec,ios
   end do
   !print *,'End b loop 2 ',nrec,ios,cfld,char8,lcoord,coord

   if (cfld==char8 .and. lcoord==coord .and. ios==0) then
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
      end if

      call zaiopf(trim(fbase)//'.a','old',nop)
      do irec=1,nrec-1
         call zaiosk(nop)
      end do
      call zaiord(field,mask,.false.,xmin,xmax,nop)
      call zaiocl(nop)
      !print *,minval(field),maxval(field)
   else
      print '(a,i3.3)', 'Could not get field "'//char8//'" at coordinate',coord
      field=undef
   end if

   close(nop)
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
 118  format ('member ',i5.5,' = ',i1,' Ensemble member flag')
   end subroutine read_field2d



   subroutine read_field3d(fbase,char8,field,idm,jdm,kdm,undef)
   implicit none
   integer, intent(in) :: idm,jdm,kdm
   character(len=*), intent(in) :: fbase
   character(len=8), intent(in) :: char8
   real, intent(out)    :: field(idm,jdm,kdm)
   real, intent(in) :: undef

   integer :: k

   do k=1,kdm
      call read_field2d(fbase,char8,field(:,:,k),idm,jdm,k,undef,maxcoord=kdm)
   end do
   end subroutine read_field3d


   subroutine daily_average_read_header(fbase,rtinit,rtdump,nrmem,idm,jdm,kdm)
   use mod_year_info
   use m_year_day
   implicit none
   character(len=*), intent(in)  :: fbase
   type(year_info) , intent(out) :: rtinit,rtdump
   integer,          intent(out) :: nrmem,idm,jdm,kdm

   character(len=80) :: ctitle(4)
   integer :: nop,yrflag,iexpt,iversn
   

   nop=777

   ! Open input file
   open(nop,file=trim(fbase)//'.b',status='old')
   ! Read .b Header
   read(nop,116) ctitle,iversn,iexpt,yrflag, &
      idm,jdm,kdm,rtinit%iyy,rtinit%idd,  &
      rtdump%iyy,rtdump%idd,nrmem

   !print 116, ctitle,iversn,iexpt,yrflag, &
   !   idm,jdm,kdm,rtinit%iyy,rtinit%idd,  &
   !   rtdump%iyy,rtdump%idd,nrmem

      


   if (yrflag==3) then
      call year_day(float(rtdump%idd),rtdump%iyy,rtdump,'ecmwf')
      call year_day(float(rtinit%idd),rtinit%iyy,rtinit,'ecmwf')
   elseif (yrflag==0) then
      call year_day(float(rtdump%idd),rtdump%iyy,rtdump,'month')
      call year_day(float(rtinit%idd),rtinit%iyy,rtinit,'month')
   else
      print *,'Unknow yrflag --- ',yrflag
      stop 
   end if
   !print *,rtinit
   !print *,rtdump

   close (nop)



116  format (a80/a80/a80/a80/ &
       i5,4x,'''iversn'' = hycom version number x10'/  &
       i5,4x,'''iexpt '' = experiment number x10'/  &
       i5,4x,'''yrflag'' = days in year flag'/  &
       i5,4x,'''idm   '' = longitudinal array size'/  &
       i5,4x,'''jdm   '' = latitudinal  array size'/  &
       i5,4x,'''kdm   '' = Vertical     array size'/  &
       i5,4x,'''syear '' = Year of integration start '/  &
       i5,4x,'''sday  '' = Day of integration start'/  &
       i5,4x,'''dyear '' = Year of this dump      '/  &
       i5,4x,'''dday  '' = Day of this dump     '/  &
       i5,4x,'''count '' = Ensemble counter       ')
 

    end subroutine daily_average_read_header

      

end module
