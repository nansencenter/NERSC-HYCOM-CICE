module m_save_nest_header
   contains
      subroutine save_nest_header(nestingfile,kdm,offset,cvar,klevel,year,day,hour,minv,maxv,irec)
      implicit none
         character(len=*), intent(in) :: nestingfile
         character(len=3), intent(in) :: cvar
         integer, intent(in ) :: irec
         integer, intent(in ) :: klevel,offset,year,day,hour,kdm
         real   , intent(in ) ::  minv, maxv

         logical :: ex
         character(len=80), save :: lastfile=''
         logical, save :: firstpass=.true.

         if (trim(lastfile)/=trim(nestingfile)) then
            firstpass=.true.
         else
            firstpass=.false.
         end if


         inquire(exist=ex,file=trim(nestingfile)//'.hdr')
         !print *,ex,firstpass,trim(lastfile),trim(nestingfile)
         !print *,irec
         if (firstpass) then
            if (ex) then
               open(16,file=trim(nestingfile)//'.hdr',form='formatted', &
                       access='direct',recl=100,status='replace')
            else
               open(16,file=trim(nestingfile)//'.hdr',form='formatted', &
                       access='direct',recl=100,status='new')
            end if
         else
            open(16,file=trim(nestingfile)//'.hdr',form='formatted', &
                    access='direct',recl=100,status='old')
         end if
         write(16,103,rec=irec) cvar,klevel,irec,offset,year,day,hour,kdm, &
                               maxv,minv
         close(16)

!  103 format(a3,7x,i4,8x,i5,8x,i5,6x,i4, &
!             1x,i3,1x,i2,5x,i3,9x,2e12.2)
  103 format(a3," level=",i4," record=",i5," offset=",i5," date=",i4, &
             " ",i3," ",i2," kdm=",i3," min/max:",2e12.2)

          lastfile=trim(nestingfile)

      end subroutine
end module m_save_nest_header
