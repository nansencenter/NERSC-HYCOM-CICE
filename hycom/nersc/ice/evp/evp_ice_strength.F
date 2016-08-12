      subroutine evp_ice_strength(kstrngth)
      use mod_xc
      use mod_evp
#if defined(ICE_DYN_DIAG)
      use mod_common_ice, only : pice
#endif
      implicit none

      integer, intent(in) :: kstrngth

      integer :: i,j

      if (kstrngth==1) then
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
         do i = 1-imargin,ii+margin
         do j = 1-imargin,jj+margin

           ! Hibler79, Eq 17, P=Pstar*h*exp(-C*(1-A)), h=effective thick
           strength(i,j) = Pstar0*vice(i,j)*exp(-Cstar0*(1-aice(i,j)))
#if defined(ICE_DYN_DIAG)
           pice(i,j) = strength(i,j)
#endif
         enddo   
         enddo  
!$OMP END PARALLEL DO 
      else
         write(lp,*) 'Invalid ice strength flag ',kstrngth
         call xcstop('(evp_ice_strength)')
         stop '(evp_ice_strength)'
      end if


      end subroutine

