module m_depth_frac

contains

   ! Routine to calculate fraction of a layer which is positioned between two
   ! depth levels
   ! * NB: Positive direction is downwards - "lower" interface > "upper interface"* 
   subroutine depth_frac(upper_z,lower_z,upper_layer,lower_layer,dfrac,idm,jdm)
      implicit none
      integer, intent(in) :: idm,jdm

      real,    intent(in) :: upper_z                       ! Upper interface for integration
      real,    intent(in) :: lower_z                       ! lower interface for integration  (> upper)
      real, dimension(idm,jdm), intent(in ) :: upper_layer ! Upper layer interfaces
      real, dimension(idm,jdm), intent(in ) :: lower_layer ! Lower layer interfaces
      real, dimension(idm,jdm), intent(out) :: dfrac       ! fraction of layer between z limits

      integer :: i,j,k
      real    :: up_int,lw_int,dphere
      integer :: itest,jtest
      real, parameter :: epsiloon=1e-3

      ! Error check
      if (lower_z < upper_z .or. upper_z<0 ) then
         print *,'Error -- inconsistent upper or lower boundary'
         print *,'upper_z:', upper_z
         print *,'lower_z:', lower_z
         print *, '(depth_frac)'
         call exit(1)
      end if



      itest=-1
      jtest=-1
      do j=1,jdm
      do i=1,idm

         ! This constrains xx_int to lie between lower_z and upper_z
         ! Note pos dir downwards
         up_int= min( max(upper_layer(i,j),upper_z) , lower_z) 
         lw_int= min( max(lower_layer(i,j),upper_z) , lower_z) 
         dphere= lower_layer(i,j)-upper_layer(i,j)
         dfrac(i,j)= (lw_int-up_int) / ( dphere + epsiloon)

         if (i==itest.and.j==jtest) then
            print '(8f10.2)',upper_layer(i,j),upper_z,up_int,lower_layer(i,j),lower_z,lw_int, &
               lw_int-up_int,dfrac(i,j)
         end if
      end do
      end do
         
   end subroutine depth_frac
end module m_depth_frac

