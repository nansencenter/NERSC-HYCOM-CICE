module m_initconfmap
contains
subroutine initconfmap
! This routine initialize constants used in the conformal mapping
! and must be called before the routines 'oldtonew' and 'newtoold'
! are called. The arguments of this routine are the locations of
! the two poles in the old coordiante system.
   use mod_confmap
   implicit none


! local variables

   real*8 cx,cy,cz,theta_c,phi_c
   complex*16 c,w
   character(len=80) :: reply
   logical ass

! Read info file
   open(unit=10,file='grid.info',form='formatted')
      read(10,*) lat_a,lon_a
      read(10,*) lat_b,lon_b
      read(10,*) wlim,elim,ires
      read(10,*) slim,nlim,jres
      read(10,*) dotopo
      read(10,*) dolatlon
      read(10,*) final
      read(10,*) mercator
      read(10,*) mercfac,lold
   close(10)

! some constants
   pi_1=4d0*atan(1d0)
   pi_2=.5d0*pi_1
   deg=180d0/pi_1
   rad=1.0d0/deg
   !epsil=1.0E-9
   epsil=1.0d-12

   di=(elim-wlim)/float(ires-1)   ! delta lon'
   dj=(nlim-slim)/float(jres-1)   ! delta lat' for spherical grid

   if (mercator) then
      dj=di
      dm=di
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
   theta_b=lon_b*rad
   !KAL phi_a=pi_2-lat_a*rad
   !KAL phi_b=pi_2-lat_b*rad
   phi_a=min(pi_1-real(1e-6,kind=8),pi_2-lat_a*rad) ! phi_a < pi (see below)
   phi_b=min(pi_1-real(1e-6,kind=8),pi_2-lat_b*rad) ! phi_b < pi (see below)

! find the angles of a vector pointing at a point located exactly
! between the poles

   cx=cos(theta_a)*sin(phi_a)+cos(theta_b)*sin(phi_b)
   cy=sin(theta_a)*sin(phi_a)+sin(theta_b)*sin(phi_b)
   cz=cos(phi_a)+cos(phi_b)

   theta_c=atan2(cy,cx)
   phi_c=pi_2-atan2(cz,sqrt(cx*cx+cy*cy))
#ifndef USE_BROKEN_MIDPOINT

   !KAL Fixes for "indefinite" lon/lat calculations of midpoint below...
   ! Midpoint is on z-axis. Let north pole determine theta_c
   if (sqrt(cy**2d0+cx**2d0)<1d-6) then 
      print *,'Applying longitude fix for midpoint'
      theta_c=theta_a

      ! Midpoint is actually in origo. This means the poles are opposing each
      !other on the sphere. Subtract pi/2 (90 deg) and make sure that
      !it is in the correct range
      if (abs(cz)<1d-6) then
         phi_c=cos(phi_a-pi_1/2d0)
         phi_c=acos(phi_c)
         print *,'Applying latitude fix for midpoint'
      end if

      print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
      print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
      print *,'-----------------------------------------------'
      print *,'If you are re-generating an old grid, these midpoint '
      print *,'fixes can change the grid. This is because the old way'
      print *,'of generating a grid used routines which were un-'
      print *,'predictable in some cases:'
      print *,'a) If your poles are on exact opposite sides of the earth'
      print *,'b) If your poles are symmetric about the lines '
      print *,'   joinng the real north and south pole.'
      print *
      print *,'If that is the case you can either re-generate the '
      print *,'grid (preferred solution), or you can define the flag'
      print *,'USE_BROKEN_MIDPOINT in MODEL.CPP (not preferred at all,'
      print *,'beacuse the midpoint may change from machine to machine).'
      print *,'-----------------------------------------------'
      print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
      print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   end if
#else
   print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   print *,'-----------------------------------------------'
   print *,'USE_BROKEN_MIDPOINT is set (probably in MODEL.CPP)'
   print *,'You should avoid this situation, and undef this flag'
   print *,'in MODEL.CPP. If you still want to continue'
   print *,'You should answer "Yes" below. Proceed at your'
   print *,'own risk !!!! (You may have to make changes in the HYCOM code...)'
   print *,'-----------------------------------------------'
   print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   read(*,'(a80)') reply
   if (trim(reply)/="Yes") then
      print *,reply
      print *,'Great, I will quit'
      stop '(m_initconfmap)'
   else
      print *,'Ok, have it your way...'
   end if
#endif

! initialize constants used in the conformal mapping

   imag=(.0d0,1.d0)
   ac=tan(.5d0*phi_a)*exp(imag*theta_a)
   bc=tan(.5d0*phi_b)*exp(imag*theta_b)
   c=tan(.5d0*phi_c)*exp(imag*theta_c)
   cmna=c-ac
   cmnb=c-bc

   w=cmnb/cmna
   mu_s=atan2(aimag(w),real(w))
   psi_s=2.d0*atan(abs(w))

end subroutine initconfmap
end module m_initconfmap
