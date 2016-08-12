module m_icestate_heatsolve
   contains

   !pure function icestate_heatsolve(icem,dt,condt,condb)
   function icestate_heatsolve(icem,dt,condt,condb)
      use mod_icestate
      implicit none

      type(t_ice), intent(in)                    :: icem
      real,    intent(in)                        :: dt,condt,condb
      real, dimension(1:nlaymax) :: icestate_heatsolve

      real :: dh,rksnw,cfac,efac,rkeff,cfaceff,T_new(nlaymax)
      integer i

      ! Current vertical increment
      dh = icem%hice/icem%nlay

      ! Snow conductivity
      rksnw    = rkice * (icem%rhosnw / rhow)**1.885

      ! Get conductivity factor 
      cfac = rkice / dh

      ! Get  energy factor 
      efac = dt / (dh*rhoice*cpice)

      !if ( abs(cfac*efac)>.5) stop 'heatsolve scheme unstable'
      !print *,'cfac*efac=',cfac*efac

      ! Solve for temp in ice/snow interior
      do i = 2, icem%nlay-1
         T_new(i) = icem%vtp(i) + efac*cfac*( icem%vtp(i+1) - 2*icem%vtp(i) + icem%vtp(i-1) )
      end do

      ! Lowest layer temperature. 
      T_new(1) = icem%vtp(1) + &
             efac*(cfac*(icem%vtp(2) - icem%vtp(1)) - condb ) 

      ! Effective conductivity
      rkeff   = rksnw*rkice / (rksnw*dh/2 + rkice*icem%hsnw)
      cfaceff = rkeff / (icem%hsnw+dh/2) 

      ! Upper layer temperature
      T_new(icem%nlay) = icem%vtp(icem%nlay) + &
                        efac * ( condt - cfac*(icem%vtp(icem%nlay)-icem%vtp(icem%nlay-1)) )

      ! Set remaining layers to upper layer temp
      T_new(icem%nlay+1:nlaymax) = t_new(icem%nlay)

      ! Update
      icestate_heatsolve = T_new

   end function icestate_heatsolve

end module m_icestate_heatsolve
