module m_shfilt
contains
subroutine shfilt(n,sh,dim,x,incx,y,incy,shdim)
       implicit logical (a-z)
       integer n,i,j,incx,incy,dim,ix,iy,shdim,m
       real sh(0:shdim),x(dim*incx),y(dim*incy),sum
!-------------------------------------------------------------------
!      This routine calculates the filtered vector y(i) from
!      the vector x(i), using the shapirofilter of order n,
!      and the method described in "shfact"
!
!      This routine must be called after the routine "shfact", which
!      computes the factors used in the filter.
!
!       n      The order of the filter.
!       sh     The vector of dimension sh(0:2n) which
!              returns the factors in the filter.
!       dim    The number of elements in the vectors x and y.
!       x      The vector containing the elements to be filtered.
!       incx   The increment between elements which should be filtered.
!       y      The vector containing the filtered result.
!       incy   The increment between filtered elements in y.
!
!
!                  Made by Geir Evensen juni 1991,
!                          geir@nrsc.no
!
!-------------------------------------------------------------------
       if (n.le.0) return
       
       y(1)=x(1)
       y(dim)=x(dim)

       if ((incx.eq.1).and.(incy.eq.1)) then 
!        code for both increments equal to 1
         do 111 i=2,n
           sum=0.0
           do 112 j=0,n-1
              m=1-(i-n+j)
              if (m.ge.0) then
                  sum=sum+sh(j)*(2.0*x(1)-x(1+m) + x(i+n-j))
              else
                  sum=sum+sh(j)*(x(i+n-j)+x(i-n+j))
              endif
 112       continue
           y(i)=(1.0+sh(n))*x(i)+sum
 111     continue

         do 11 i=n+1,dim-n
           sum=0.0
           do 12 j=0,n-1
              sum=sum+sh(j)*(x(i+n-j)+x(i-n+j))
 12        continue
           y(i)=(1.0+sh(n))*x(i)+sum
 11      continue
         
         do 211  i=dim-n+1,dim-1
           sum=0.0
           do 212 j=0,n-1
              m=(i+n-j)-dim
              if (m.ge.0) then
                  sum=sum+sh(j)*( 2.0*x(dim)-x(dim-m) + x(i-n+j))
              else
                  sum=sum+sh(j)*(x(i+n-j)+x(i-n+j))
              endif
 212       continue
           y(i)=(1.0+sh(n))*x(i)+sum
 211     continue


      else
!       code for unequal increments or equal increments not equal to 1
        do 13 i=n+1,dim-n
           ix=(i-1)*incx + 1
           iy=(i-1)*incy + 1

           sum=0.0
           do 14 j=0,n-1
              sum=sum+sh(j)*(x(ix+(n-j)*incx)+x(ix-(n-j)*incx))
 14        continue
           
           y(iy)=(1.0+sh(n))*x(ix)+sum
 13      continue
       endif 
       return
       end subroutine shfilt
       end module m_shfilt
