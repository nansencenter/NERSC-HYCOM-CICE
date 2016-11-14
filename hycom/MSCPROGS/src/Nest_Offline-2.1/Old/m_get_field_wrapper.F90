module m_get_field_wrapper

contains
  subroutine get_field_wrapper(fld,k,cvar,filename,filetype)
  use mod_xc
  use mod_read_dailyab
  implicit none
  real, dimension(idm,jdm), intent(out) :: fld
  integer, intent(in) :: k
  character(len=*), intent(in) :: cvar,filename,filetype

  integer :: find
  character(len=8) :: char8
  character(len=8) :: char8_2
  character(len=80) :: fbase
  real, dimension(idm,jdm) :: fldtmp


  if (trim(filetype) == 'daily') then
     find=max(index(filename,'.a')-1,index(filename,'.b')-1)
     fbase=filename(1:find)
     char8='        '
     char8(1:len_trim(cvar))=cvar(1:len_trim(cvar))


     ! Vector variables (total velocity in daily file)
     if (trim(cvar)=='u') then
        call read_field2d(trim(fbase),'utot    ',fld   ,idm,jdm,k,0.)
        call read_field2d(trim(fbase),'ubavg   ',fldtmp,idm,jdm,0,0.)
        fld=fld-fldtmp
     else if (trim(cvar)=='v') then
        call read_field2d(trim(fbase),'vtot    ',fld   ,idm,jdm,k,0.)
        call read_field2d(trim(fbase),'vbavg   ',fldtmp,idm,jdm,0,0.)
        fld=fld-fldtmp
     else
        call read_field2d(trim(fbase),char8,fld,idm,jdm,k,0.)
     end if

  else
     print *,'unknown filetype '//trim(filetype)
     stop '(get_field_wrapper)'
  end if
  end subroutine get_field_wrapper

end module m_get_field_wrapper


