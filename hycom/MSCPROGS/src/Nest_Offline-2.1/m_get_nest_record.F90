module m_get_nest_record


      contains

      subroutine get_nest_record(cvar,klevel,nestingfile,irec,frec,lrec)
      implicit none
         character(len=*), intent(in) :: cvar
         character(len=*), intent(in) :: nestingfile
         integer, intent(out) :: irec
         integer, intent(in ) :: klevel,frec,lrec

         character(len=3) :: cin
         integer ::kin,recin,offsin,yrin,dayin,hrin,kdmin,ios,ir
         real    :: maxin,minin
         logical :: found,ex

         found=.false.

         inquire(file=trim(nestingfile)//'.hdr',exist=ex)
         if (.not.ex) then
            print *,'Could not find nest header file '//trim(nestingfile)//'.hdr'
            stop '(get_nest_record)'
         end if
         open(16,file=trim(nestingfile)//'.hdr',form='formatted', &
                 access='direct',recl=100,status='old')

         ir=frec
         ios=0
         do while (.not. found .and. ios==0 .and. ir<=lrec) 
            !print *,ir
            read(16,103,rec=ir) cin,kin,recin,offsin,yrin,dayin,hrin,kdmin, &
                         maxin,minin

            found = trim(cin)==trim(cvar) .and. kin==klevel
            irec=recin
            ir=ir+1
            !print *,cin,cvar,kin,klevel,found
         end do
         close(16)

         if (ios/=0) then
            print *,'Input error '
            stop '(get_nest_record)'
         end if

         if (.not. found) then
            print '(a,a,a,i2)','Could not find variable ',cvar,' at level ',klevel
            !stop '(get_nest_record)'
            irec=-1
         end if


  103 format(a3,7x,i4,8x,i5,8x,i5,6x,i4, &
             1x,i3,1x,i2,5x,i3,9x,2e12.2)

      end subroutine

end module m_get_nest_record
