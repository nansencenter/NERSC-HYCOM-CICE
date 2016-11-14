module m_shfact
contains
subroutine shfact(n,sh)
       implicit logical (a-z)
       integer n,i,j
       real sh(0:n),f2n,fj,f2nj,ff
!-------------------------------------------------------------------
!      This routine calculates the factors in the general formula
!      for the Shapiro filter, given as
!       
!                                      (-1)^{n+j-1}(2n)!
!        y(i)=x(i) + \Sum_{j=0}^{2n} { -------------------- } 
!                                       2^{2n} j!  (2n-j)!
!      
!      which is calculated from equation (39) in 
!      (The multi-dimensional Crowley advection scheme,
!      Piotr K. Smolarkiewicz, Monthley weater review,
!      volume 110, pages 1968-1983, december 1982)
!      Note the error in the paper.
!
!       n      is the order of the filter. n=2^j  j=0,1,2, ...
!       sh     is the vector of dimension sh(0:2n) which
!              returns the factors in the filter.
!
!      Normally this routine is called first to compute the factors
!      and followed by "shfilt" for the filtering.
!
!                  Made by Geir Evensen juni 1991,
!                          geir@nrsc.no
!
!-------------------------------------------------------------------
       if (n.gt.8) then
          write(*,*)'Error in "shfact"'
          write(*,100)'n is to large (>8): n=',n
          stop
       endif
       if (amod(alog(float(n)),alog(2.0)).ne.0.0) then
          write(*,*)'Error in "shfact"'
          write(*,100)'Invalid order of filter: n=',n
 100      format(' ',a,i2)
          stop
       endif

       ff=2.0**(2*n)
       
       
       f2n=1
       do 10 j=1,2*n
          f2n=f2n*float(j)
 10    continue

       fj=1
       sh(0)=((-1.0)**(n-1))/ff
       do 11 j=1,n

!        calculates (2n-j)! for each j
         f2nj=1
         do 12 i=1,2*n-j
            f2nj=f2nj*float(i)
 12      continue

!        calculates j!
         fj=fj*float(j)
!        calculates the factors
         sh(j)=(-1.0)**(n+j-1) *real( f2n/(ff*fj*f2nj) )
 11    continue
       return
       end subroutine shfact
end module m_shfact
