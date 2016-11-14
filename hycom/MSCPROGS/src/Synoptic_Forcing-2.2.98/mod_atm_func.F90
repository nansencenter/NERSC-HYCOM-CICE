module mod_atm_func

! -  A Collection of conversion formulas between humidity, vapor
! -  pressure, mixig ratio etc etc.
! -  Vapor pressure and humidity formulas (se e.g Gill and Tetens)
! -
! -  satvappi           - Saturation vapor pressure over ice
! -  satvappw           - Saturation vapor pressure over ice
! -  satvap             - Can do both of the above by checking temperature
! -  humid              - Calculates humidity
! -  relhumid           - Calculates relative humidity
! -  spechum_to_relhum  - Conversion from specific humidity to relative humidity
! -  rhtovpmix          - Converts from relative humidity to mixing ration

contains 


real function satvap(t)
! This function calculates the saturation vapour pressure
! [Pa] from the temperature [deg K].
! Modified: Anita Jacob, June '97
!
! Input: t: temperature [deg K]
! Output: satvap: saturation vapour pressure at temp. t
!
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


!#######################################################################
!Saturation vapor pressure over ice                                    
!NB: This is essentially the same equation as "satvap", but without
!    Temperature checks on coefficients "a" and "b"
!#######################################################################
real function satvappi(tx)
! This function calculates the saturation vapour pressure
! [Pa] from the temperature.
! Input is temperature in Kelvins
implicit none
real, intent(in) :: tx
satvappi=611.*10.**(9.5*(tx-273.16)/(tx-7.66)) ! tx in K
end function satvappi




!#######################################################################
!Saturation vapor pressure over water
!NB: This is essentially the same equation as "satvap", but without
!    Temperature checks on coefficients "a" and "b"
!#######################################################################
real function satvappw(tx)
! This function calculates the saturation vapour pressure
! [Pa] from the temperature.
! Input is temperature in Kelvins
implicit none
real, intent(in) :: tx 
satvappw=611.*10.**(7.5*(tx-273.16)/(tx-7.50)) ! tx in K
end function satvappw




!#######################################################################
!Specific humidity from vapor pressure
!#######################################################################
real function humid(px,es)
! Function calculates the specific humidity based on the total presure and 
! the vapor ! pressure (ie partial pressure of vapor). Nondimensional. Both
! Input is in Pa (doesnt really matter as long as they have the same unit)
! Based on formulas in Gill (1982).
implicit none
real,intent(in) :: px,es
humid=.622*es/(px-.378*es)
end function humid





real function relhumid(sva,svd,msl)
! This routine calculates the relative humidity by the 
! dew point temperature and the mean sea level pressure.
! Modified: Anita Jacob, June '97

! Input:
!    sva: saturatn vapour press at air temp [K]
!    svd: saturatn vapour press at dew pt temp [K]
!    msl: pressure at mean sea level [Pa]
! Output: 
!   relhumid: Relative Humidity

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


real function spechum_to_relhum(q,qs)
! Function converts from specific humidity to relative humidity. 
implicit none
real, intent(in) :: q   ! humidity
real, intent(in) :: qs  ! saturation humidity
real :: aa,bb

aa=q *(1-qs)
bb=qs*(1-q )
spechum_to_relhum=aa/bb
end function


real function rhtovpmix(relhum,t,slp)
! convert relative humidity to vapor mixing ratio
! KAL - Wrong! Need to fix this formula even though we dont use vapor mixing
! ratio
implicit none
real, intent(in) :: relhum ! relative humidity (range 0-1)
real, intent(in) :: t      ! air temperature [K]
real, intent(in) :: slp    ! sea level pressure [Pa]
real :: svap, sathum, rs

svap=satvap(t)         ! Saturation vapor pressure
sathum=humid(slp,svap) ! Specific humidity at saturation
rs=sathum/(1-sathum)   ! Mixing ratio at saturation
rhtovpmix = rs * relhum
end function
end module
