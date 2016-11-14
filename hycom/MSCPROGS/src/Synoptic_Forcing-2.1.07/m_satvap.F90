module m_satvap
contains
real function satvap(t)
! This function calculates the saturation vapour pressure
! [Pa] from the temperature [deg K].
! Modified: Anita Jacob, June '97

! Input: t: temperature [deg K]
! Output: satvap: saturation vapour pressure at temp. t

! es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual


      implicit none
      real, intent(in):: t
      real :: aa,bb,c1,c3,c4,t00,cc
      data c1/610.78/,t00/273.16/

      if (t < t00) then
         c3 = 21.875
         c4 = 7.66
      else
         c3 = 17.269
         c4 = 35.86
      endif

      aa = c3 * (t - t00)
      bb = t - c4
      
      cc=aa/bb
      if (cc < -20.0) then
         satvap=0.0
      else
         satvap = c1 * exp(aa/bb)
      endif

end function satvap
end module m_satvap
