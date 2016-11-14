module m_filenesting
contains
function filenesting(rt,nestname)
   use mod_year_info
   implicit none
   character(len=100) filenesting

   type(year_info),  intent(in) :: rt
   character(len=*), intent(in) :: nestname
   integer l,trimlen

   !filenesting=' '
   !do l=1,40
   !   if (nestname(l:l) == ' ') exit
   !enddo
   !l=l-1

   if (len_trim(nestname)==0) then
      filenesting='nest_'//rt%cyy//'_'//rt%cdd
   else
      filenesting=nestname//'nest_'//rt%cyy//'_'//rt%cdd
   end if
end function filenesting
end module m_filenesting
