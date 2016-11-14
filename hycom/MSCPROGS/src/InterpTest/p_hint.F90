program hint
! A simpler version of the m2nc/m2t scripts and associated programs.
! Variables to plot are read from extract - files, and then put into
! 4D/3D netcdf variables with horizontal, vertical and time dimensions
! Everything is presented on the model grid.
! TODO: adhere to rotation/sphere etc
!
! Knut Liseter
   use mod_xc
   use mod_za
   use mod_grid
   implicit none
   character(len=80) :: c80
   real, allocatable :: twod(:,:)
   real :: tmp

   call xcspmd()
   call zaiost()
   call get_grid()
   allocate(twod(idm,jdm))
   read(*,*) twod

end program hint
