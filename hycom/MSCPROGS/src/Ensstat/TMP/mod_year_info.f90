module mod_year_info
   type year_info
      integer iyy             ! year
      integer imm             ! month
      integer idd             ! day in year
      integer ihh             ! hours
      integer iss             ! seconds
      integer idm             ! current day in month
      character(len=4) cyy    ! year nr  (char)
      character(len=2) cmm    ! month nr (char)
      character(len=3) cdd    ! day in year  (char)
      character(len=2) chh    ! hour nr (char)
      character(len=4) css    ! second nr (char)
      character(len=3) month  ! month 'JAN' etc
      character(len=2) cdm    ! current day in month (char)
      integer totdim(12)      ! total Days In Months 
      integer daysinyear      ! total Days In year
   end type year_info
end module mod_year_info
