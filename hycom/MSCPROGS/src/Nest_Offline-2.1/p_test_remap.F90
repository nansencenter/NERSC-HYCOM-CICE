program p_test_remap
use mod_sigma
use m_layer_remapV4
implicit none

real, parameter :: saln0=34.5

integer, parameter :: thflag=0
integer, parameter :: oldkdm=22
integer, parameter :: newkdm=26

real, dimension(oldkdm) :: olddens,oldint,oldtemp,oldsaln
real, dimension(newkdm) :: newdens,newint,newdp0

real, parameter :: dz=3.
real, parameter :: dsigold=2. , dsignew=8.5

integer :: k
logical :: fatal

newdp0=5;
oldsaln=35.
do k=1,oldkdm

   !layer thickness increases with square of k
   oldint(k)=k**2*dz

   !layer densities increase exponentially
   olddens(k)=24.0 + dsigold * exp(-real(oldkdm-k)/oldkdm)

   if (thflag==0) then
      oldtemp(k)=tofsig0(olddens(k),saln0)
   elseif (thflag==2) then
      oldtemp(k)=tofsig2(olddens(k),saln0)
   end if

end do

print *,'Old densities : ',olddens
print *
print *,'Old temp      : ',oldtemp
print *
print *,'Old interfaces: ',oldint
print *

do k=1,newkdm

   !layer densities increase exponentially
   newdens(k)=20.5 + dsignew * exp(-real(newkdm-k)/newkdm)

end do
print *
print *,'New densities : ',newdens
print *


call layer_remapV4(oldint,oldtemp,oldsaln,oldkdm,  &
                   newdens,newint,newdp0,newkdm,thflag,.true., fatal)
print *
print *,'New interfaces: ',newint

end program p_test_remap
