module m_initconfmap
contains
subroutine initconfmap(nxg,nyg)
! This routine initialize constants used in the conformal mapping
! and must be called before the routines 'oldtonew' and 'newtoold'
! are called. The arguments of this routine are the locations of
! the two poles in the old coordiante system.

   !use mod_dimensions
   use mod_confmap
   implicit none


! local variables

   integer, intent(in) :: nxg, nyg
   real cx,cy,cz,theta_c,phi_c
   complex c,w
   logical ass,lold,ex

! Read info file
   inquire(exist=ex,file='grid.info')
   if (.not.ex) then
      !if (mnproc==1) 
      print *,'Can not find grid.info'
      !call xcstop('(initconfmap)')
      stop '(initconfmap)'
   end if
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
   if ((ires /= nxg).and.(jres /= nyg)) then
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
         !if (mnproc==1) 
         print *,'initconfmap: lold'
         slim=-mercfac*jres*dj
      else
         !if (mnproc==1)
         print *,'initconfmap: not lold'
         slim= mercfac
      endif
   endif

! transform to spherical coordinates

   theta_a=lon_a*rad
   theta_b=lon_b*rad
   !KAL phi_a=pi_2-lat_a*rad
   !KAL phi_b=pi_2-lat_b*rad
   phi_a=min(pi_1-1e-6,pi_2-lat_a*rad) ! phi_a < pi (see below)
   phi_b=min(pi_1-1e-6,pi_2-lat_b*rad) ! phi_b < pi (see below)

!   if (mnproc==1) then
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!      print *,'KAL - changes to initconfmap!'
!   end if

! find the angles of a vector pointing at a point located exactly
! between the poles

   cx=cos(theta_a)*sin(phi_a)+cos(theta_b)*sin(phi_b)
   cy=sin(theta_a)*sin(phi_a)+sin(theta_b)*sin(phi_b)
   cz=cos(phi_a)+cos(phi_b)

   ! KAL - what if cx**2+cy**2+cz**2 is very small?
   ! This is the case for the Barents Sea model.... The ifdef statement below
   ! "fixes" this
   theta_c=atan2(cy,cx)
   phi_c=pi_2-atan2(cz,sqrt(cx*cx+cy*cy))
#ifdef BARENTS
   theta_c= 1.57079632679489656
   phi_c  = 0.425580672355377887
#endif
   
   ! KAL - this should fix the problems.... But some grids will
   ! become broken - hence the ifdef..
#ifndef USE_BROKEN_MIDPOINT
   !KAL Fixes for "indefinite" lon/lat calculations of midpoint below...
   ! Midpoint is on z-axis. Let north pole determine theta_c
   if (sqrt(cy**2+cx**2)<1e-6) then 
      !if (mnproc==1) 
      print *,'Applying longitude fix for midpoint'
      theta_c=theta_a
      ! Midpoint is actually in origo. This means the poles are opposing each
      !other on the sphere. Subtract pi/2 (90 deg) and make sure that
      !it is in the correct range
      if (abs(cz)<1e-6) then
         phi_c=cos(phi_a-pi_1/2)
         phi_c=acos(phi_c)
         !if (mnproc==1) 
         print *,'Applying latitude fix for midpoint'
      end if

      !if (mnproc==1) then
         print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
         print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
         print *,'-----------------------------------------------'
         print *,'If you are using an old grid, these midpoint '
         print *,'fixes can break the grid. This is because the old way'
         print *,'of generating a grid used routines which were un-'
         print *,'predictable in some cases:'
         print *,'a) If your poles are on exact opposite sides of the earth'
         print *,'b) If your poles are symmetric about the lines '
         print *,'   joinng the real north and south pole.'
         print *
         print *,'If that is the case you can either re-generate the '
         print *,'grid (preferred solution), or you can define the flag'
         print *,'USE_CONFMAP_BROKEN_MIDPOINT in MODEL.CPP (not preferred at all,'
         print *,'beacuse the midpoint may change from machine to machine).'
         print *,'-----------------------------------------------'
         print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
         print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
      !end if

   end if
#else
   !if (mnproc==1) then
   print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   print *,'-----------------------------------------------'
   print *,'USE_CONFMAP_BROKEN_MIDPOINT is set (probably in MODEL.CPP)'
   print *,'You should avoid this situation, and undef this flag'
   print *,'in MODEL.CPP, period. If you still want to continue'
   print *,'You proceed at your own risk !!!! (You probably '
   print *,'have to make changes in HYCOM support routines)'
   print *,'-----------------------------------------------'
   print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   print *,'WARNING WARNING WARNING WARNING WARNING WARNING'
   !endif
#endif


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

!   if (mnproc==1) then
!      print *,'pi_1',pi_1
!      print *,'pi_2',pi_2
!      print *,'theta_a',theta_a
!      print *,'phi_a',phi_a
!      print *,'theta_b',theta_b
!      print *,'phi_b',phi_b
!      print *,'theta_c',theta_c
!      print *,'phi_c',phi_c
!      print *,'ac',ac
!      print *,'bc',bc
!      print *,'c',cmna+ac
!      print *,'cx',cx
!      print *,'cy',cy
!      print *,'cz',cz
!      print *,'cmna',cmna
!      print *,'cmnb',cmnb
!
!      print *,'w',w
!      !print *,'z',z
!      !print *,'psi',psi
!      !print *,'mu',mu
!   end if
!   call xcstop('initconfmap')


end subroutine initconfmap
end module m_initconfmap
