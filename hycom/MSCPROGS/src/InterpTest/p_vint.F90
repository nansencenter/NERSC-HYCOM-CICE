program vint
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
   use mod_hycomfile_io
   implicit none
   character(len=80) :: filename, ftype
   type(hycomfile) :: hfile
   real, allocatable :: twod(:,:)

   call getarg(1,filename) 
   ftype=getfiletype(filename)
   call xcspmd(10)
   call zaiost()
   allocate(twod(idm,jdm))
   call initHF(hfile,trim(filename),trim(ftype))
   call HFReadField(hfile,twod,idm,jdm,'saln    ',1,1)
   !write(6,'(10e20.10)') twod
   write(*,*) twod
end program vint
