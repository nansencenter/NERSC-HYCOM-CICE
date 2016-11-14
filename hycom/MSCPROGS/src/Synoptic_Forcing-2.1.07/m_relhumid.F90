module m_relhumid
contains
real function relhumid(sva,svd,msl)
! This routine calculates the relative humidity by the 
! dew point temperature and the mean sea level pressure.
! Modified: Anita Jacob, June '97

! Input: sva: saturatn vapour press at air temp [K]
! svd: saturatn vapour press at dew pt temp [K]
! msl: pressure at mean sea level [Pa]
! Output: relhumid: Relative Humidity

! We use the Tetens formula:
! es(T) = C1 * exp(C3*(T - T0)/(T - C4)) from ECMWF manual
!              es(Tdew)        p - es(Tair)
! RH = 100 *  -----------   *  ------------
!             p - es(tdew)       es(Tair)
       
      implicit none

      real, intent(in):: sva
      real, intent(in):: svd
      real, intent(in):: msl
      real :: aaa,bbb
      aaa=msl - svd
      aaa = svd/aaa
      bbb = (msl - sva)/sva
      relhumid = 100. * aaa * bbb

end function relhumid
end module m_relhumid
