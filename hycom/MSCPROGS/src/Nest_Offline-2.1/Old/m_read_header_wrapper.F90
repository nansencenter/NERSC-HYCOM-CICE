module m_read_header_wrapper
contains
   subroutine read_header_wrapper(cfile,ftype,rt,idm,jdm,kdm)
   use mod_read_dailyab
   use mod_year_info
   implicit none
   integer, intent(in) :: idm, jdm, kdm
   character(len=*), intent(in) :: cfile, ftype

   type(year_info) :: rt, rtdump, rtinit
   integer :: idm2,jdm2,kdm2,find,nrmem

   ! Read header from file type (only daily files for now)
   if (trim(ftype)=='daily') then

      ! Get file base from filename
      find=max(index(cfile,'.a'),index(cfile,'.b'))-1

      call daily_average_read_header(cfile(1:find),rtinit,rtdump, &
         nrmem,idm2,jdm2,kdm2)

      ! Consistency check
      if (idm/=idm2 .or. jdm/=jdm2 .or. kdm/=kdm2) then
         print *,'daily file dimension mismatch'
         print *,'idm =',idm,idm2
         print *,'jdm =',jdm,jdm2
         print *,'kdm =',kdm,kdm2
         stop '(read_header_wrapper)'
      end if

      rt=rtdump 
   else
      print *,'Dont know how to handle file type '//ftype
      stop '(read_header_wrapper)'
   end if
   end subroutine

end module
