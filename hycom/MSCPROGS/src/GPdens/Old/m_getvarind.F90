module m_getvarind
contains
integer function getvarind(vartofind,varname,nvars)
implicit none
integer,intent(in) :: nvars
character(len=*),intent(in) :: vartofind
character(len=12),dimension(nvars) :: varname
integer :: ivar

getvarind=-1
do ivar=1,nvars
   if (trim(vartofind)==trim(varname(ivar))) then
      getvarind=ivar
   end if
end do
end function
end module





