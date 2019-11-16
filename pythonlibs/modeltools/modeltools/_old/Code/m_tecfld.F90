module m_tecfld
contains
subroutine tecfld(fname,ii,jj,lon,lat,deep,ires)
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ii,jj,ires
   real, intent(in) :: lon(ii,jj)
   real, intent(in) :: lat(ii,jj)
   real, intent(in) :: deep(ii,jj)
   integer i,j,iii,jjj

   iii=0
   do i=1,ii,ires
      iii=iii+1
   enddo

   jjj=0
   do j=1,jj,ires
      jjj=jjj+1
   enddo

   open(10,file='tec'//fname//'.dat',status='unknown')
      write(10,*)'TITLE = "',fname,'"'
      write(10,*)'VARIABLES = "i" "j" "lon" "lat" "deep"'
      write(10,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',iii,', J=',jjj,', K=1'

      write(10,'(20I6)')((i,i=1,ii,ires),j=1,jj,ires)
      write(10,'(20I6)')((j,i=1,ii,ires),j=1,jj,ires)
      write(10,900)((lon(i,j),i=1,ii,ires),j=1,jj,ires)
      write(10,900)((lat(i,j),i=1,ii,ires),j=1,jj,ires)
      write(10,900)((deep(i,j),i=1,ii,ires),j=1,jj,ires)
   close(10)
 900 format(10(1x,e12.5))

end subroutine tecfld

end module m_tecfld
