module m_oldtonew
contains
subroutine oldtonew(lat_o,lon_o,lat_n,lon_n)

! this routine performes a conformal mapping of the old to the new
! coordinate system

   use mod_confmap
   implicit none

   real, intent(in)  :: lat_o,lon_o
   real, intent(out) :: lat_n,lon_n

! local variables

   real theta,phi,psi,mu
   complex z,w

! transform to spherical coordinates


   theta=mod(lon_o*rad+3.0*pi_1,2.0*pi_1)-pi_1
   phi=pi_2-lat_o*rad

! transform to the new coordinate system

   if (abs(phi-pi_1) < epsil) then
     mu=mu_s
     psi=psi_s
   elseif ((abs(phi-phi_b)<epsil).and.&
           (abs(theta-theta_b)<epsil)) then
     mu=.0
     psi=pi_1
   else
     z=tan(.5*phi)*exp(imag*theta)
     w=(z-ac)*cmnb/((z-bc)*cmna)
     mu=atan2(aimag(w),real(w))
     psi=2.*atan(abs(w))
   endif

! transform to lat/lon coordinates

   lat_n=(pi_2-psi)*deg
   lon_n=mu*deg

end subroutine oldtonew
end module m_oldtonew
