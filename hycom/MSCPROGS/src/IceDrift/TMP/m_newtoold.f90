module m_newtoold
contains
subroutine newtoold(lat_n,lon_n,lat_o,lon_o)
! this routine performes a conformal mapping of the new to the old 
! coordinate system

   use mod_confmap
   implicit none

   real lat_o,lon_o,lat_n,lon_n

! local variables

   real theta,phi,psi,mu
   complex w,z

! transform to spherical coordinates

   mu=mod(lon_n*rad+3*pi_1,2*pi_1)-pi_1
   psi=abs(pi_2-lat_n*rad)

! transform to the old coordinate system

   if (abs(psi-pi_1) < epsil) then
     theta=theta_b
     phi=phi_b
   elseif ((abs(mu-mu_s) < epsil).and.((psi-psi_s)<epsil)) then
     theta=.0
     phi=pi_1
   else
     w=tan(.5*psi)*exp(imag*mu)
     z=(ac*cmnb-w*bc*cmna)/(cmnb-w*cmna)
     theta=atan2(aimag(z),real(z))
     phi=2.*atan(abs(z))
   endif

! transform to lat/lon coordinates

   lat_o=(pi_2-phi)*deg
   lon_o=theta*deg

end subroutine newtoold
end module m_newtoold

