program section_transport2
   use mod_xc
   use mod_za
   use mod_grid
   use mod_sections
   !use mod_fieldtypes
   use mod_year_info
   use mod_transport
   use mod_read_rstab
   use mod_read_weekab
   use mod_read_dailyab
   use mod_readpak
   use m_read_section_intersect
   use m_transport2
   use m_datetojulian
   use m_get_header_info
   use m_getfiletype
   implicit none

#if defined (IARGC)
   integer*4, external :: iargc
#endif

   character(len=50) :: char50
   character(len=80) :: fnamein,fbase
   character(len=24) :: char24,plot_time
   CHARACTER(len=3)  :: rungen

   integer :: n2d,n3d,ii1,jj1
   integer :: iday,iyear,ihour,imonth,iweek
   integer :: thisday,idayinmonth,endyear

   integer :: nrmem,hycfile,kdm

   integer, parameter :: maxch=70
   character(len=8), dimension(maxch) :: ch2d,ch3d
   real :: rtime

   type (year_info) :: rti, rtd


   if (iargc()==1) then
      call getarg(1,fnamein)
   else
      print *,'One input argument only'
      print *, '(section_transport)'
      call exit(1)
   end if


   ! What file is this? (daily? weekly? restart? pak?)
   call getfiletype(trim(fnamein),fbase,hycfile)

   ! Read different fields with header info - returns undef if not defined
   if (hycfile==2) then
      call daily_average_read_header(trim(fbase),rti,rtd,nrmem,idm,jdm,kdm,ch2d,n2d,ch3d,n3d,maxch)
   else if (hycfile==1) then
      call rst_read_header(trim(fbase),rti,nrmem,idm,jdm,kdm,ch2d,n2d,ch3d,n3d,maxch)
      rtd=rti
      nrmem=1
   else if (hycfile==3) then
      call week_read_header(trim(fbase),rti,idm,jdm,kdm,ch2d,n2d,ch3d,n3d,maxch)
      rtd=rti
      nrmem=1
   else if (hycfile==4) then
      call read_pakheader(trim(fbase)//'.hdr',idm,jdm,kdm,rungen,plot_time,ch2d,ch3d,n2d,n3d,maxch)
      char24=plot_time
      nrmem=1


      !call get_header_info(idm,jdm,kdm,rungen,char24,n2d,n3d,ii1,jj1)
      !call get_header_info(idm,jdm,kdm,rungen,char24,ch2d,ch3d,n2d,n3d,ii1,jj1,maxch)
      print *,'File header read'
      if (char24(4:6)/= "AVE") then
         read(char24(5  :8),'(i4.4)') iyear
         read(char24(10:12),'(i3.3)') iday
         read(char24(14:15),'(i3.3)') ihour
         !rtime=iyear + min((float(iday) + float(ihour)/24)/365.,1.)
         rtd%iyy=iyear
         rtd%idd=iday
         rtd%ihh=ihour
      else
         if (char24(8:10)=="ave") then
            print '(a)','This is a multiple average file .. can not extract time ... '
            rtime=1.
         else
            read(char24( 8:11),'(i4.4)') iyear
            read(char24(13:14),'(i2.2)') imonth
            read(char24(16:16),'(i1.1)') iweek
            !rtime=iyear + float(imonth-1)/12 + float(iweek-1)/(4*12) + float(1)/(4*12*2)
            rtd%iyy=iyear
            rtd%idd=iday
            rtd%ihh=ihour
         end if
      end if
   else 
      print *,'Unknown file type ... I only accept restart/DAILY/WEEKLY or pak files'
      print *,'with their original file names...'
      print *
      print *,'Name of problematic file: '//trim(fnamein)
      print *, '(section_transport2)'
      call exit(1)
   end if
   rtime=rtd%iyy + min((float(rtd%idd) + float(rtd%ihh)/24)/365.,1.)
   !print *,rtime


   ! For safety
   if (kdm <1 .or. kdm > 50 )  then
      print *,'kdm is ',kdm
      print *,'This is probably a bug ...'
      call exit(1)
   end if

   ! Initialize IO for .ab files
   CALL XCSPMD  ! -- Requires "regional.grid.b" to be present
   CALL ZAIOST

   ! Get model grid & Depths
   print *,'Grid read'
   call get_grid()
   

   ! Read section details -- this assumes "section_intersect" is already run
   print *,'Section details read'
   call read_section_intersect()


   ! Calculate transport -- dump to files
   print *,'Reading transport details'
   call transport_specification
   print *,'Calculating transport'
   !call transport2(rtime)
   call transport2(hycfile,nrmem,kdm,rtime,trim(fbase))



   deallocate(scuy)
   deallocate(scvx)
   write(6,'(a)') 'section_transport finished -- OK'
end program section_transport2
