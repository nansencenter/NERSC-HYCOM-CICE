module m_shapiro_filt
contains
subroutine shapiro_filt(ish,sh,x,nx,ny,shdim,af,if,il,is,jf,jl,js,ms)
   use m_sfilt
   implicit none
   integer j,i,l,yn,m,mm,ms
   integer ish,shdim
   integer max,nx,ny,af
   real sh(0:shdim)
   real x(nx,ny)    ! input field and output field
   real, allocatable ::  y(:)    ! work array
   real, allocatable ::  z(:)    ! work array

   integer if(ny,ms)
   integer il(ny,ms)
   integer jf(nx,ms)
   integer jl(nx,ms)
   integer is(ny)
   integer js(nx)

   allocate(y(max(nx,ny)))
   allocate(z(max(nx,ny)))


   if((2*ish+1.GT.nx).OR.(2*ish+1.GT.ny))then
      write(*,*)'shfilt2:  The domain is to small for ish=',ish
      return
   endif


   !Her begynner mitt program

   do i=1,af

   y=0.0

   do j=1,ny
     do l=1,is(j)
       yn=il(j,l)-if(j,l)+1
       if(yn.GT.2*ish)then
         y(1:yn)=x(if(j,l):il(j,l),j)
         call sfilt(sh,ish,shdim,y,yn,z)
         x(if(j,l):il(j,l),j)=y(1:yn)
       endif
     enddo
   enddo

   y=0.0

   do j=1,nx
     do l=1,js(j)
       yn=jl(j,l)-jf(j,l)+1
       if(yn.GT.2*ish)then
         y(1:yn)=x(jf(j,l):jl(j,l),j)
         call sfilt(sh,ish,shdim,y,yn,z)
         x(jf(j,l):jl(j,l),j)=y(1:yn)
       endif
     enddo
   enddo

   enddo

end subroutine shapiro_filt
end module m_shapiro_filt

