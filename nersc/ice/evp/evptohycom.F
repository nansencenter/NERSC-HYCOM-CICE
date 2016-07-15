      subroutine evptohycom
      use mod_common_ice , only : iceU, iceV, tauxice, tauyice
      use mod_evp
      implicit none

      integer :: i,j
      real :: sumupoint,sumvpoint

      ! Ocean and wind velocities in EVP u-point
c --- TODO: tauxice/tauyice should perhaps be on p-grid
      imargin=min(imargin,nbdy-1)
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
C
CKAL - This can create artificial convergence when going between
CKAL   different grids. Important when HYCOM does the advection.
         iceU(i,j) = 0.5*(uvel(i-1,j-1) + uvel(i-1,j  ) )
         iceV(i,j) = 0.5*(vvel(i-1,j-1) + vvel(i  ,j-1) )
C   
CKAL - tauxice - not weighted with conc.
         tauxice(i,j) = - 0.5*( strocnx(i-1,j  ) + strocnx(i-1,j-1))
         tauyice(i,j) = - 0.5*( strocny(i-1,j-1) + strocny(i  ,j-1))
C
      end do
      end do
!$OMP END PARALLEL DO 

      call xctilr(iceU   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_uv)
      call xctilr(iceV   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_vv)
      call xctilr(tauxice( 1-nbdy,1-nbdy),1,   1, 6,6, halo_uv)
      call xctilr(tauyice( 1-nbdy,1-nbdy),1,   1, 6,6, halo_vv)

      ! 
      end subroutine evptohycom
