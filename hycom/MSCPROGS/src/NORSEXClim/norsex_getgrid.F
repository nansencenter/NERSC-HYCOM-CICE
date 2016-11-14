      subroutine norsex_getgrid(lon,lat,idm,jdm)
      implicit none
      integer, intent(in) :: idm,jdm
      real, intent(out) , dimension(idm,jdm) :: lon,lat

      integer :: i,j

      do j=1,jdm
      do i=1,idm
         call norsex_locate(2,1,1,real(i),real(j),lon(i,j),lat(i,j))
         !print *,i,j,lon(i,j),lat(i,j)
      end do
      end do

      end subroutine



