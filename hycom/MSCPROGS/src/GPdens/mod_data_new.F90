module mod_data_new
type data_new
   real day    ! Time in days rel ref year
   real s      ! saln psu
   real t      ! temp deg C
   real dir    ! degrees 0 is North, and clockwise
   real speed  ! cm/s
   real u      ! dir east   cm/s
   real v      ! dir north  cm/s
end type data_new
real , parameter :: pi=3.1415927
end module mod_data_new
