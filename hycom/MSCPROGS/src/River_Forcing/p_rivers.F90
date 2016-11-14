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
   !elseif (iargc()==0) then
   !   print *,'using standard values for radii'
   else
      print *,'*******************************************************'
      print *,'The river forcing routine calculates monthly climatology'
      print *,'river forcing fields for use with HYCOM. The rivers are'
      print *,'specified in the file rivers.dat.'
      print *
      print *,'The input arguments are two distances. The 1st (radius) '
      print *,'indicates how far rivers should  extend away from the '
      print *,'coast. The second (alongshoreradius) indicates how far '
      print *,'along the coast the rivers should extend. These two '
      print *,'parameters are used to create a discharge area for the '
      print *,'different rivers.'
      print *,'*******************************************************'
      print *, 'Usage:'
      print *,'  rivers radius alongshoreradius'
      print *,'Units are in km'
      print *,'Example:'
      print *,'  rivers 100 200'
      print *,'Normally rivers spread more along the coast than normal to it.'
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
