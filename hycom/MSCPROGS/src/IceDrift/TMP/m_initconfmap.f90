module m_initconfmap
contains
subroutine initconfmap(nx,ny)
! This routine initialize constants used in the conformal mapping
! and must be called before the routines 'oldtonew' and 'newtoold'
! are called. The arguments of this routine are the locations of
! the two poles in the old coordiante system.

   !use mod_dimensions
   use mod_confmap
   implicit none

   integer, intent(in) :: nx,ny


! local variables

   real cx,cy,cz,theta_c,phi_c
   complex c,w
   logical ass,lold

! Read info file
   open(unit=10,file='grid.info',form='formatted')
      read(10,*) lat_a,lon_a
      read(10,*) lat_b,lon_b
      read(10,*) wlim,elim,ires
      read(10,*) slim,nlim,jres
      read(10,*) ass
      read(10,*) ass
      read(10,*) ass
      read(10,*) mercator
      read(10,*) mercfac,lold
   close(10)
   if ((ires /= nx).and.(jres /= ny)) then
      print *,'initconfmap: WARNING -- the dimensions in grid.info are not'
      print *,'initconfmap: WARNING -- consistent with nx and ny'
      print *,'initconfmap: WARNING -- IGNORE IF RUNNING CURVIINT'
   endif

! some constants
   pi_1=4.*atan(1.)
   pi_2=.5*pi_1
   deg=180./pi_1
   rad=1.0/deg
   epsil=1.0E-9

   di=(elim-wlim)/float(ires-1)   ! delta lon'
   dj=(nlim-slim)/float(jres-1)   ! delta lat' for spherical grid

   if (mercator) then
      dj=di
      if (lold) then
         print *,'initconfmap: lold'
         slim=-mercfac*jres*dj
      else
         print *,'initconfmap: not lold'
         slim= mercfac
      endif
   endif

! transform to spherical coordinates

   theta_a=lon_a*rad
   phi_a=pi_2-lat_a*rad
   theta_b=lon_b*rad
   phi_b=pi_2-lat_b*rad

! find the angles of a vector pointing at a point located exactly
! between the poles

   cx=cos(theta_a)*sin(phi_a)+cos(theta_b)*sin(phi_b)
   cy=sin(theta_a)*sin(phi_a)+sin(theta_b)*sin(phi_b)
   cz=cos(phi_a)+cos(phi_b)

   theta_c=atan2(cy,cx)
   phi_c=pi_2-atan2(cz,sqrt(cx*cx+cy*cy))

! initialize constants used in the conformal mapping

   imag=(.0,1.)
   ac=tan(.5*phi_a)*exp(imag*theta_a)
   bc=tan(.5*phi_b)*exp(imag*theta_b)
   c=tan(.5*phi_c)*exp(imag*theta_c)
   cmna=c-ac
   cmnb=c-bc

   w=cmnb/cmna
   mu_s=atan2(aimag(w),real(w))
   psi_s=2.*atan(abs(w))

end subroutine initconfmap
end module m_initconfmap
