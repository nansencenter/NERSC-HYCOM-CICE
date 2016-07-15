module m_evp_next_step
contains
logical function EVP_next_step(rt,OGCMstep)
      use mod_evp
      use mod_year_info, only : year_info
      implicit none
      type (year_info), intent(in) :: rt
      real, intent(in) :: OGCMstep
      real :: time

      time = rt%idd*86400. + rt%ihh*3600. + rt%iss

      EVP_next_step = nint((time)/dyn_dt) /= nint((time+OGCMstep)/dyn_dt)

   end function EVP_next_step
end module m_evp_next_step

