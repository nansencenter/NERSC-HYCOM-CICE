!=======================================================================
!BOP
!
! !IROUTINE: to_ugrid - shift from T-cell to U-cell midpoints
!
! !INTERFACE:
!
      subroutine to_ugrid(work1,work2,lmargin)
      use mod_evp
      implicit none
!
! !DESCRIPTION:
!
! Shifts quantities from the T-cell midpoint (work1) to the U-cell 
! midpoint (work2)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), intent(in)  ::
     &   work1
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), intent(out) ::
     &   work2
      integer, intent(in) :: lmargin
!
!EOP
!
      integer :: i, j

      if (lmargin>=nbdy) then
         if (mnproc==1) print *,'to_ugrid called with margin=',lmargin
         call xcstop('(to_ugrid)')
         stop '(to_ugrid)'
      end if
         

      do j=1-lmargin,jj+lmargin
      do i=1-lmargin,ii+lmargin
      if (umask(i,j)) then
       work2(i,j) = 0.25* (work1(i  ,j  )*tarea(i  ,  j) 
     &                   + work1(i+1,j  )*tarea(i+1,  j)
     &                   + work1(i  ,j+1)*tarea(i  ,j+1) 
     &                   + work1(i+1,j+1)*tarea(i+1,j+1))/uarea(i,j)
      else
         work2(i,j)=0.
      end if
      enddo
      enddo

      end subroutine to_ugrid
