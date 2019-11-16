module m_shapiro_ini
contains
subroutine shapiro_ini(n,sh,shdim)
   implicit none
   integer n
   integer i,j,shdim
   real sh(0:shdim)
   real f2n,fj,f2nj,ff
   if(n.GT.8)then
   write(*,*)'Error in "shfact"'
   write(*,*)'n is to large (>8): n=',n
   stop
   endif
   if((n.EQ.1).OR.(n.EQ.2).OR.(n.EQ.4).OR.(n.EQ.8))then
   ff=2.0**(2*n)
   f2n=1
   do j=1,2*n
   f2n=f2n*float(j)
   enddo
   fj=1
   sh(0)=((-1.0)**(n-1))/ff
   do j=1,n

   f2nj=1
   do i=1,2*n-j
   f2nj=f2nj*float(i)
   enddo

   fj=fj*float(j)

   sh(j)=(-1.0)**(n+1-j)*f2n/(ff*fj*f2nj)
   enddo
   else
   write(*,*)'error in shfact.  n=',n
   stop
   endif
end subroutine shapiro_ini
end module m_shapiro_ini
