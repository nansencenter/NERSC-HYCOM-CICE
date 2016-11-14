module m_read_nest_header
contains
      subroutine read_nest_header(nestingfile,kdm,offset,cvar,klevel,year,day,hour,minv,maxv,irec)
      implicit none
         character(len=*), intent(in) :: nestingfile
         integer, intent(in) :: irec
         character(len=3), intent(out) :: cvar
         integer, intent(out) :: offset,kdm,year,day,hour,klevel
         real   , intent(out) :: minv, maxv

         integer :: ios,ir
         real    :: maxin,minin
         logical :: ex

         inquire(exist=ex,file=trim(nestingfile)//'.hdr')
         if (.not.ex) then
            print *,'Could not find nest header file '//trim(nestingfile)//'.hdr'
            stop '(get_nest_info)'
         end if

         open(16,file=trim(nestingfile)//'.hdr',form='formatted', &
                 access='direct',recl=100,status='old')

         read(16,103,rec=irec,iostat=ios) cvar,klevel,ir,offset,year,day,hour,kdm, &
                                        maxv,minv
         if (ios/=0) then
            print *,'Input error '
            stop '(read_nest_header)'
         end if
         close(16)




  103 format(a3,7x,i4,8x,i5,8x,i5,6x,i4, &
             1x,i3,1x,i2,5x,i3,9x,2e12.2)
      end subroutine


end module m_read_nest_header
