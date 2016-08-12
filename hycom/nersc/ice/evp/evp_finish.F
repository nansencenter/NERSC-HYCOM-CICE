!=======================================================================
!BOP
!
! !IROUTINE: evp_finish - calculates ice-ocean stress 
!
! !INTERFACE:
!
      subroutine evp_finish
!
! !DESCRIPTION:
!
! Calculation of the ice-ocean stress \\
! ...the sign will be reversed later...
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use mod_evp
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      implicit none
      integer :: i, j

      real :: vrel
      !real, dimension(itdm,jtdm) :: gtmp1,gtmp2

      imargin=6
!$OMP PARALLEL DO PRIVATE(j,i,vrel)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
      if (iceumask(i,j)) then
        ! ocean-ice stress for coupling
        vrel = dragw*sqrt((uocn(i,j) - uvel(i,j))**2 + 
     &                    (vocn(i,j) - vvel(i,j))**2)  ! m/s
        strocnx(i,j) = strocnx(i,j) 
     &               - vrel*(uvel(i,j)*cosw - vvel(i,j)*sinw)
        strocny(i,j) = strocny(i,j) 
     &               - vrel*(vvel(i,j)*cosw + uvel(i,j)*sinw)
      else
        ! set stress to zero on land and (nearly) open water
        strocnx(i,j) = 0.0
        strocny(i,j) = 0.0
      endif

      !-----------------------------------------------------------------
      ! convert strocn to T grid
      !-----------------------------------------------------------------
      !strocnxT(i,j) = strocnx(i,j)  ! prepare to shift
      !strocnyT(i,j) = strocny(i,j)

      enddo
      enddo
!$OMP END PARALLEL DO 

      call xctilr(strocnx( 1-nbdy,1-nbdy),1,   1, 6,6, halo_uv)
      call xctilr(strocny( 1-nbdy,1-nbdy),1,   1, 6,6, halo_vv)

      ! KAL - TODO - interp is done in evptohycom
      !call u2tgrid(strocnxT)        ! shift
      !call u2tgrid(strocnyT)

!      call xcaget(gtmp1,strocnx,0)
!      call xcaget(gtmp2,strocny,0)
!      if (mnproc==1) then 
!         open(10,file='test.evpfinish')
!         do j=1,jtdm
!         do i=1,itdm
!            write (10,'(2i5,2e16.3)') i,j,gtmp1(i,j),gtmp2(i,j)
!         end do
!         end do
!         close(10)
!      end if

      end subroutine evp_finish
