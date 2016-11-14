module m_construct_filename
contains

   subroutine construct_filename(hfile,ncfil)
   use mod_year_info
   use mod_hycomfile_io
   implicit none
   type (hycomfile), intent(in ) :: hfile
   character(len=*), intent(out) :: ncfil

   character(len=80) ::fbase
   integer :: find, hycfile
   type(year_info) :: rtd, rtb

   fbase=hfile%filebase
   call forecastDate(hfile,rtd)
   call startDate(hfile,rtb)

   if (trim(hfile%ftype)=='nersc_daily') then ! daily
      find=index(fbase,'DAILY_')
      find=find+5
      ncfil=fbase(1:find)//                           &
            'start'//rtb%cyy//rtb%cmm//rtb%cdm// &
            '_dump'//rtd%cyy//rtd%cmm//rtd%cdm// &
            '.nc'
   else if (trim(hfile%ftype)=='restart') then ! daily
      find=index(fbase,'restart')
      find=find+6
      ncfil=fbase(1:find)//                           &
            rtb%cyy//rtb%cmm//rtb%cdm//rtb%chh//'.nc'
   else if (trim(hfile%ftype)=='nersc_weekly') then ! daily
      !find=index(fbase,'AVE')
      !find=find+2
      !ncfil=fbase(1:find)//                           &
      !      rtb%cyy//'_'//rtb%cmm//'_'//rtb%cdm//'.nc'
      ncfil=trim(fbase)//'.nc'
   else if (trim(hfile%ftype)=='archv') then ! archv
      !TODO add dump time to filename
      find=index(fbase,'archv')
      find=find+4
      ncfil=fbase(1:find)//                &
            '_'//rtd%cyy//rtd%cmm//rtd%cdm//'_'//rtd%chh//'.nc'
   else if (trim(hfile%ftype)=='archv_wav') then ! archv_wav
      find=index(fbase,'archv_wav')
      find=find+8
      ncfil=fbase(1:find)//                &
            '_start'//rtb%cyy//rtb%cmm//rtb%cdm//'_'//hfile%start_ctime//'Z'  &
          //'_dump'//rtd%cyy//rtd%cmm//rtd%cdm//'_'//hfile%ctime//'Z.nc'
   else
      print *,'construct_filename: Unknown file type '//trim(hfile%ftype)
   end if
   end subroutine

!   subroutine construct_filename_stat(hfile,ncfil)
!   use mod_year_info
!   use mod_hycom_read
!   implicit none
!   type (hycomfile), intent(in ) :: hfile
!   character(len=*), intent(out) :: ncfil
!
!   integer :: find
!
!   if (hfile%hycfile==2) then ! daily
!      find=index(hfile%fbase,'DAILY_')
!      find=find+5
!      ncfil=hfile%fbase(1:find)//                           &
!            'start'//hfile%rtb%cyy//hfile%rtb%cmm//hfile%rtb%cdm// &
!            '_dump'//hfile%rtd%cyy//hfile%rtd%cmm//hfile%rtd%cdm// &
!            '.nc'
!   elseif (hfile%hycfile==1) then ! restart
!      find=index(hfile%fbase,'restart')
!      find=find+6
!      ncfil=hfile%fbase(1:find)//                           &
!            hfile%rtb%cyy//'_'//hfile%rtb%cmm//'_'//hfile%rtb%cdm//'.nc'
!   elseif (hfile%hycfile==3) then ! wweekly
!      find=index(hfile%fbase,'AVE')
!      find=find+2
!      ncfil=hfile%fbase(1:find)//                           &
!            hfile%rtb%cyy//'_'//hfile%rtb%cmm//'_'//hfile%rtb%cdm//'.nc'
!   end if
!   end subroutine



      

end module
