program rivers
   use mod_xc
   use mod_za
   use mod_grid
   use mod_rivers
   implicit none

   integer :: i
   real    :: amin, amax
   character(len=80) :: tmparg
#if defined (IARGC)
   integer*4, external :: iargc
#endif

   ! Get arguments
   if (iargc()==2) then
      call getarg(1,tmparg)
      read(tmparg,*) rradius
      rradius=rradius*1000.
      call getarg(2,tmparg)
      read(tmparg,*) alongshoreradius
      alongshoreradius= alongshoreradius*1000.
   elseif (iargc()==0) then
      print *,'using standard values for radii'
   else
      print *, 'Usage:'
      print *,'  rivers [radius] [alongshoreradius]'
      print *,'Units are in km'
      stop
   end if
   print *, 'River      radius [km]: ',rradius/1000.
   print *, 'Alongshore radius [km]: ',alongshoreradius/1000.


      


   ! Initialize IO
   call XCSPMD()
   call ZAIOST()
   call get_grid()
   call rivers_to_hycom(plon,plat,scpx,scpy,depths)
   print *
   print *, 'River      radius [km]: ',rradius/1000.
   print *, 'Alongshore radius [km]: ',alongshoreradius/1000.
   print *

end program
