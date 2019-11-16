module m_sfilt
contains
subroutine sfilt(sh,ish,shdim,y,yn,z)
   implicit none
   integer ish,shdim,yn,i,j,m,mm
   real sh(0:shdim)
   real y(1:yn)
   real z(1:yn)

   z=y

   do m=2,ish
     do j=0,ish-1
       mm=1-(m-ish+j)
       if(mm.GE.0)then
         z(m)=z(m)+sh(j)*(2.0*y(1)-y(1+mm)+y(m+ish-j))
       else
         z(m)=z(m)+sh(j)*(y(m+ish-j)+y(m-ish+j))
       endif
     enddo
     z(m)=(sh(ish))*y(m)+z(m)
   enddo

   do m=ish+1,yn-ish
     do j=0,ish-1
         z(m)=z(m)+sh(j)*(y(m+ish-j)+y(m-ish+j))
     enddo
     z(m)=(sh(ish))*y(m)+z(m)
   enddo

   do m=yn-ish+1,yn-1
     do j=0,ish-1
       mm=(m+ish-j)-yn
       if(mm.GE.0)then
         z(m)=z(m)+sh(j)*(2.0*y(yn)-y(yn-mm)+y(m-ish+j))
       else
         z(m)=z(m)+sh(j)*(y(m+ish-j)+y(m-ish+j))
       endif
     enddo
     z(m)=(sh(ish))*y(m)+z(m)
   enddo

   do i=1,yn
     y(i)=z(i)
   enddo

end subroutine sfilt
end module m_sfilt


