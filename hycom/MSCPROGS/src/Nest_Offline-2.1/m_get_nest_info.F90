module m_get_nest_info
contains
      subroutine get_nest_info(nestingfile,kdm2,nrec_offset,num_offset)
      implicit none
         character(len=*), intent(in) :: nestingfile
         integer, intent(out) :: nrec_offset,kdm2,num_offset  !,year,day,hour

         character(len=3) :: cin
         integer ::kin,recin,offsin,yrin,dayin,hrin,kdmin,ios,ir,maxrec
         real    :: maxin,minin
         logical :: ex

         inquire(exist=ex,file=trim(nestingfile)//'.hdr')
         if (.not.ex) then
            print *,'Could not find nest header file '//trim(nestingfile)//'.hdr'
            stop '(get_nest_info)'
         end if

         open(16,file=trim(nestingfile)//'.hdr',form='formatted', &
                 access='direct',recl=100,status='old')

         ir=1
         ios=0
         maxrec=0
         num_offset=0
         do while(ios==0)
            read(16,103,rec=ir,iostat=ios) cin,kin,recin,offsin,yrin,dayin,hrin,kdmin, &
                                           maxin,minin

            if (ios==0) then
               ir=ir+1
               num_offset=max(num_offset,offsin)
               maxrec=max(maxrec,recin)
               kdm2=kdmin
            end if
         end do
         ir=ir-1
         close(16)
         if (ios/=0.and.ir==0) then
            print *,'Input error '
            print *,'iostat is ',ios
            stop '(get_nest_info)'
         end if

         kdm2=kdmin
         nrec_offset=maxrec/num_offset
         !year=yrin
         !day=dayin
         !hour=hrin



  103 format(a3,7x,i4,8x,i5,8x,i5,6x,i4, &
             1x,i3,1x,i2,5x,i3,9x,2e12.2)
      end subroutine


end module m_get_nest_info
