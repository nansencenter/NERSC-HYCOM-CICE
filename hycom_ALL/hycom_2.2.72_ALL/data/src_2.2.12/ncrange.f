      subroutine ncrange(h,ii,jj,kk, fill_value, hmin,hmax)
      implicit none
c
      integer, intent(in ) ::   ii,jj,kk
      real*4,    intent(in ) :: h(ii,jj,kk)
      real,  intent(in)      :: fill_value
      real,    intent(out) :: hmin,hmax
c
c     return range of array, ignoring fill_value
c
      integer i,j,k
      real    hhmin,hhmax
c

      hhmin =  abs(fill_value)
      hhmax = -abs(fill_value)
      do k= 1,kk
        do j= 1,jj
          do i= 1,ii
            if     (h(i,j,k) .ne. fill_value) then
              hhmin = min(hhmin,h(i,j,k))
              hhmax = max(hhmax,h(i,j,k))
            endif
          enddo
        enddo
      enddo
      hmin = hhmin
      hmax = hhmax

      end subroutine ncrange
