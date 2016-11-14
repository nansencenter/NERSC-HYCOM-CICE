program p_tecline

! plots points which can be plotted as lines uin tecplot 

  implicit none
  integer m,i,j
  integer,parameter :: sti=1,endi=300 !NWAG
  integer,parameter :: stj=406,endj=406


  open(10,file='tecline.dat',status='unknown')
  do i=sti,endi
   do j=stj,endj
    write(10,'(2i5)') i,j   
    enddo
  enddo

!  write(10,'(30I4)')((i,i=1,ii),j=1,jj)
!  write(10,'(30I4)')((j,i=1,ii),j=1,jj)
!  write(10,900)((fld(i,j,1),i=1,ii),j=1,jj)

  close(10)
!900 format(10(1x,e12.5))

end program p_tecline


