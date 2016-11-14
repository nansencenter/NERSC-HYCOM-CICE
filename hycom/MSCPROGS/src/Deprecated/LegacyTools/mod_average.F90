module mod_average
! type defenitions and variables for  computation of averages

   use mod_dimensions

   type average_state
      character(len=3) cmm
      character(len=4) cyy
      integer counter
      real ssh(nx,ny)
      real u(nx,ny,nz)
      real v(nx,ny,nz)
      real d(nx,ny,nz)
      real t(nx,ny,nz)
      real s(nx,ny,nz)
      real w(nx,ny,nz)
   end type average_state

   type tmp_states
      real u
      real v
      real d
      real t
      real s
      real u2
      real v2
   end type tmp_states

   type average_test
      logical res
      logical ini
      logical add
      logical sav
   end type average_test

   type(average_state) ave_week
   type(average_test)  lave

   integer days_in_week(4),current_week,totit
   character(len=19) fileavew
   real t1,t2

end module mod_average
