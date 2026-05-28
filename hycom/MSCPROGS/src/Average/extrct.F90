subroutine extrct(work,n,m,io,jo,array,no,mo)
 implicit none

 integer :: n,m,io,jo,no,mo
 real :: work(n,m),array(no,mo)
!
! --- array = work(io:io+no-1,jo:jo+mo-1)
!
! --- this version   c y c l i c   in i

 integer :: i,j

 do j=1,min(mo,m-jo+1)
 do i=1,no
    array(i,j)=work(mod(io+i-2,n)+1,jo+j-1)
 enddo
 enddo

 return
end subroutine extrct
