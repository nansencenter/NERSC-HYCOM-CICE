module m_icestate_prec
contains

! -----------------------------------------------------------------------
! -------------------------- SUBROUTINE PREC2_V2-------------------------
! -----------------------------------------------------------------------
! Subroutine which takes precipitations over ice or snow into 
! account, except the case of rain over non snow covered sea ice.  Returns
! surplus water, measured in kg/m^3
! OBS !! Salinity is not conserved here ... The best way to fix this ?
! 1) Introduce fake salinity fluxes ..
! 2) introduce a spatially dependent ice salinity
!
! To be corrected ...
subroutine icestate_prec(icem,rainfallx,tairx,rainx,snowx,rhosw, &
                         surpl_water,sal_corr)
   USE mod_icestate
#if defined(SSNOWD)
   use mod_icestate_tools
#endif
      
   IMPLICIT NONE
   type(t_ice), dimension(nthick),intent(inout) :: icem
   REAL,    INTENT(in)                       :: rainfallx,tairx,rhosw
   real,    intent(out)                      :: surpl_water
   real,    intent(out)                      :: sal_corr
   LOGICAL, INTENT(out)                      :: rainx,snowx

   REAL  :: aux(nthick), aux2(nthick), precip_temp,exc_wtr(nthick)
#if defined(SSNOWD)
   real :: snwfrac(nthick),lam(nthick),z_dm(nthick)
#endif  

   rainx  = .FALSE.
   snowx  = .FALSE.
   exc_wtr=0.
   surpl_water=0.

   ! Before any rainfall, the density of all snow layers increases towards rhosnwmax.
   where (icem%fice>epsil1)
      aux     = icem%rhosnw
      icem%rhosnw = (icem%rhosnw - rhosnwmax) * exp(- tauf * dtt / 86400.) + rhosnwmax
      icem%hsnw   = icem%hsnw * aux / icem%rhosnw
#if defined(SSNOWD)
   !2 configuration : 
   !   - if melting has not occured : snow compaction reduces accumulated snow
   !         depth
   !   - if melting has started. A part of the ice is snow free. Snow free
   !      fraction remains constant with compaction and accumulated  depth and
   !      melt depth are reduced.
      icem%hprcp = icem%hprcp * aux / icem%rhosnw
      icem%hmelt = icem%hmelt * aux / icem%rhosnw
#endif      
   endwhere

   ! Now in case there is precipitation, update rhosnwx, hsnwx and hice1x.
   aux=0.
   IF (rainfallx > 0. ) THEN
      IF (tairx < t0deg  ) THEN

         ! If air temp. < ... C, then snow is falling.
         snowx = .TRUE.

         ! On the different ice types:
         where (icem%fice>epsil1)
            icem%rhosnw = icem%hsnw*icem%rhosnw + rainfallx*rhow
            icem%hsnw   = icem%hsnw             + rainfallx*rhow/rhosnwmin
            icem%rhosnw = icem%rhosnw / icem%hsnw
#if defined (SSNOWD)
         !No snow at the previous time step
         where(icem%hprcp<epsil1)
            icem%cv    = cvsnw
            icem%hprcp = icem%hsnw
            icem%hmelt = 0.
         elsewhere
            ! 100% snow cover at previous time step. No melting has occured
            where(icem%hmelt<epsil1)
               icem%hprcp= icem%hprcp+rainfallx*rhow/rhosnwmin
            ! Less than 100% snow cover at previous time step, with new accumulation.
            ! Two options are possible (described in Liston , 2004)
            ! We use option 2 : new accumulation decreases the melt depth, 
            ! "pushing" the depletion curve back towards 100 %    
            elsewhere
             !c Two different cases :
             !     - there is enough new accumulation to push the melt depth to zero
             !     - there is not enough (the melt depth is still positive) : hsnw is updated
             ! and hmelt is reduced based on Eq 18 from Liston, 2004. (iterative solution)
               where(icem%hsnw>=icem%hprcp)
                  icem%hprcp=icem%hsnw
                  icem%hmelt=0.
               elsewhere
                  icem%hmelt = fixed_point(50,icem)
               endwhere
             endwhere 
          endwhere    
#endif            
         endwhere
      ELSE   ! It rains
         rainx = .TRUE.

         where (icem%fice>epsil1 .AND. icem%hsnw > epsil1)

            ! The superficial snow layer absorbs rain and its density increases.
            icem%rhosnw = icem%rhosnw * icem%hsnw + rhow * rainfallx
            icem%rhosnw = icem%rhosnw / icem%hsnw
            where (icem%rhosnw > rhosnwmax)
               where (icem%rhosnw < rhoice)

                  ! Snow density is more than rhomax, divide into ice and snow at rhomax
                  aux = icem%hsnw*(icem%rhosnw-rhosnwmax)/(rhoice-rhosnwmax)
                  icem%hice  = icem%hice + icem%hsnw*(icem%rhosnw-rhosnwmax)/(rhoice-rhosnwmax)
                  icem%hsnw  = icem%hsnw*(rhoice-icem%rhosnw)/(rhoice-rhosnwmax)
                  icem%rhosnw= rhosnwmax
#if defined (SSNOWD)
                  !hsnw is reduced : hmelt is updated accordingly. Formation of
                  !surimposed ice and reduction of snow cover fraction.
                  icem%hmelt = fixed_point(50,icem)
#endif                  
               elsewhere

                  ! Snow density is more than ice density, create ice and excess weater
                  aux = aux + icem%hsnw
                  icem%hice   = icem%hice + icem%hsnw
                  exc_wtr     = icem%hsnw*(icem%rhosnw-rhoice)/rhow
                  icem%hsnw   = 0.
                  icem%rhosnw = rhosnwmin
#if defined (SSNOWD)
                  !Snow is removed from the sea-ice surface. 
                  icem%hprcp  = 0.
                  icem%hmelt  = 0. 
#endif                  
               endwhere
            endwhere
         elsewhere
            exc_wtr = rainfallx
         endwhere

      ENDIF

      ! Surpl_water is the water entering the total ice fraction, hence the division
      surpl_water = rhow*sum(icem%fice*exc_wtr) / (sum(icem%fice)+epsil1)

   ENDIF

   ! If ice draft is below sea level (due to snow load) convert snow into ice to
   ! create hydrostatic balance.
   aux2=0.
   where  (icem%hsnw*icem%rhosnw>icem%hice*(rhosw-rhoice)+epsil1 .and. icem%fice>epsil1 )
     aux2      = ( icem%hsnw*icem%rhosnw - icem%hice*(rhosw-rhoice) ) / rhosw
     icem%hsnw = icem%hsnw  -aux2*rhoice/icem%rhosnw
     icem%hice = icem%hice + aux2
#if defined (SSNOWD)
     !hsnw is reduced according to volume conservation. hmelt is increased since
     !bare ice is formed and a part a ice is now snow-free
     icem%hmelt = fixed_point(50,icem)
#endif     
  endwhere

  ! Salinity correction factor (put into mixed layer)
  sal_corr = sum((aux+aux2)*icem%fice)*rhoice*sice / (sum(icem%fice)+epsil1)

  ! Snow is assumed to be limited in thickness due to wind blowing it 
  ! into leads and polynya. The following puts an upper limit on
  ! the snow thickness.
  aux=0.
#if defined (SSNOWD)
  ! No restriction is applied when SSNOWD is used. (Validation is required). 
#else  
  where (icem%hsnw>snwlim)
     aux   = (icem%hsnw-snwlim)*icem%rhosnw
     icem%hsnw = snwlim
  endwhere
#endif  

  ! Surpl_water is the water entering the total ice fraction, hence the division
  surpl_water = surpl_water + sum(aux*icem%fice)/(sum(icem%fice)+epsil1)

END SUBROUTINE icestate_prec
END module m_icestate_prec

