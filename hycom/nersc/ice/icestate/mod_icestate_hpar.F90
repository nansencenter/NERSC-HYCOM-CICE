module mod_icestate_hpar

   REAL, save      :: hice_min = .30      ! Minimum ice thickness for lateral growth [m]

contains


! =======================================================================
! The horizontal "freezing" routine. New ice created, initially with
! a thickness. hice_min. Special care is taken if that fails..
! =======================================================================
   subroutine icestate_hfreeze(icem,qrest,omficetot,fice_,hice_)
      use mod_icestate
      implicit none
      type(t_ice), intent(inout) :: icem(nthick) ! State of the ice
      real,        intent(in)    :: qrest        ! Flux available for freezing
      real,        intent(in)    :: omficetot    ! Fraction of open water
      real,        intent(out)   :: hice_        ! Thickness of newfrozen ice
      real,        intent(out)   :: fice_        ! Fraction of newfrozen ice

      real :: hice2_,fice2_,latvol,dfice_
      type(t_ice) :: icem2(nthick)

      !icem2=icem

      hice_  = hice_min

      ! Ice equiv. to the qrest energy for a unit area of open water
      dfice_ = -qrest * dtt  / ( hofusn0 * hice_ )

      ! The ice equiv. to the energy for the __lead__ area is created
      fice_  = omficetot * dfice_


      ! If total ice area exceeds fice_max, then adjust hice_
      if (fice_ + sum(icem%fice) > fice_max ) then

         ! fice2_ is the "overshoot"
         fice2_ = sum(icem%fice) + fice_ - fice_max
         fice_  = fice_ - fice2_
         !print *,fice2_

         ! The surplus volume due to the "overshoot" is divided among the classes.
         latvol    = fice2_ * hice_
         hice2_    = latvol / sum(icem%fice)
         icem%hice = icem%hice + hice2_
         !print *,hice2_


      ! If total area is lower than fice_min, adjust hice_, so that new vol has area fice_
      else if ( fice_ + sum(icem%fice) < fice_min + epsil1) then

         fice2_ = fice_min - sum(icem%fice)
         hice_  = hice_ * fice_ / fice2_
         fice_  = fice2_

      endif

      !print *,hice_,fice_,omficetot*dfice_,sum(icem2%hice-icem%hice)

   end subroutine icestate_hfreeze







! =======================================================================
! The horizontal "melting" routine. Ice is melted, taking care to not get 
! negative ice thickness in categories.
! =======================================================================
   subroutine icestate_hmelt(icem,tmlLLL,cpmlLLL,ficetot,qtotlead,hofusn)
      use mod_icestate
      implicit none
      type(t_ice), intent(inout) :: icem(nthick)   ! The initial state of the ice
      real       , intent(inout) :: hofusn(nthick) ! Heat of fusion of the ice
      real,        intent(in)    :: cpmlLLL        ! Eff. heat capacity of ml
      real,        intent(in)    :: ficetot        ! Total ice fraction [initially]
      real,        intent(inout) :: qtotlead       ! Heat flux into _LEAD_ area 
      real,        intent(inout) :: tmlLLL         ! temperature of mixed layer 
                                                   ! after qtotlead has been absorbed

      real, dimension(nthick) :: qrestLLL,dfice,icefrac,imass,smass
      real :: omficetot,tmp_ftot

!KAL: -- Lateral melt of ice disregards vtp in ice for energy budget 

      ! Scheme from Hakkinen and Mellor (1990) 
      !dficetot = min( 0.7 * omficetot *  qml_ice / sum(hice*fice + hsnw*fice*rhsnw) , qtotlead )
      !dfice    = dficetot * fice / tmp_ftot
      !tmlLLL   = tmlLLL - dficetot * dtt * sum(hice*fice+hsnw*rhosnw*fice) * latheat
      !rest as before

      

      ! temporary ficetot (total ice concentration [now])
      omficetot=1.-ficetot
      tmp_ftot = sum(icem%fice)

      ! "omficetot" fraction of energy goes to heat the lead ML. The following
      ! corrects tmlLLL calculated prior to this routine.  (Which used all the
      ! energy)
      tmlLLL = tmlLLL - ficetot * qtotlead * dtt / cpmlLLL


      ! "ficetot" fraction goes to lateral melting. If the following is not traversed, 
      ! all the energy has gone towards lead warming (remainder of this if cause could 
      ! be skipped).
      where (icem%fice>epsil1)
         icefrac = icem%fice / (tmp_ftot +1e-8)
         imass = icem%hice
         smass = icem%hsnw*icem%rhosnw/rhoice

         ! Ice equiv. to the energy "qtotlead * icem%fice * ficetot / tmp_ftot" for a unit area
         dfice       =  - ficetot * qtotlead * dtt * icefrac / &
                          ( imass*hofusn + smass*hofusn0 )

         ! The ice equiv. to the energy for the __lead__ area is subtracted
         icem%fice = icem%fice + omficetot * dfice
      endwhere


      ! Adjust if a class becomes less than zero. 
      qrestLLL=0. !! NB -- for correct totlead below
      do while (any(icem%fice<-epsil1))
         where (icem%fice<-epsil1)
            qrestLLL = - icem%fice*(imass*hofusn+smass*hofusn0)/dtt
            icem%fice = 0.
         elsewhere
            qrestLLL=0.
         endwhere

         qtotlead = sum(qrestLLL)
         tmp_ftot = sum(icem%fice)
         icefrac = icem%fice / (tmp_ftot +1e-8)

         where (icem%fice>epsil1)
            icem%fice = icem%fice - qtotlead * dtt * icefrac / &
                        (hofusn*imass +smass*hofusn0)
         endwhere

      end do


      ! We may return with qtotlead > 0  here. Also, qtotlead is now TOTAL energy
      ! over the unit area. The omficetot division gives energy going into the 
      ! lead  area again ...
      qtotlead = sum(qrestLLL) / (omficetot + 1e-8)

      ! Now, any remaining heat heats the mixed layer below the lead
      ! Adjust for energy going into _lead__ area
      !tmlLLL = tmlLLL + ficetot * qtotlead * dtt / cpmlLLL
      tmlLLL = tmlLLL + qtotlead * dtt / cpmlLLL


      ! and remaining qtotlead is zero -- ålvæis
      qtotlead=0.


   end subroutine icestate_hmelt

! =======================================================================
! The horizontal "melting" routine. Ice is melted, taking care to not get 
! negative ice thickness in categories.
! =======================================================================
   subroutine icestate_hmelt_hak(icem,smlLLL,tmlLLL,cpmlLLL,ficetot,hofusn,dhice)
      use mod_icestate
      use mod_icestate_tools
      implicit none
      type(t_ice), intent(inout) :: icem(nthick)   ! The initial state of the ice
      real       , intent(in)    :: hofusn(nthick) ! Heat of fusion of the ice
      real,        intent(in)    :: cpmlLLL        ! Eff. heat capacity of lead ml
      real,        intent(in)    :: ficetot        ! Total ice fraction [initially]
      real,        intent(inout) :: tmlLLL         ! temperature of mixed layer 
      real,        intent(in)    :: smlLLL         ! salinity    of mixed layer 
      real,        intent(in)    :: dhice          ! dh/dt from vertical melting

      real :: omficetot,tmp_ftot,hice,effhofusn,dfice,fac,avehice,dhicedt
      real :: tmlinit,ficeinit
      integer :: hk


      ! Scheme from Hakkinen and Mellor (1990) -- modified for multi--category model


      ! Locate thinnest ice category containing ice.
      hk=1
      do while (hk<nthick .and. icem(hk)%fice<1e-8)
         hk=hk+1
      end do
      hk=min(hk,nthick)
      !write(6,'(i3)',advance='no') hk
      tmlinit=tmlLLL
      ficeinit=icem(hk)%fice

      ! Ice thickness, effective heat of fusion of that category
      ! "fac" is a handy factor when converting from lead ml energy to grid cell
      ! horizontal ice melt
      avehice   = sum(icem%hice*icem%fice)/max(ficetot,1e-8)
      hice      = icem(hk)%hice
      effhofusn = icem(hk)%hice*hofusn(hk) + icem(hk)%hsnw*icem(hk)%rhosnw*hofusn0/rhoice
      omficetot = 1.-ficetot
      dhicedt   = dhice/dtt
      fac       = effhofusn / max(cpmlLLL*omficetot,1e-8)

      ! Delta fice for a unit area in the lead -- only for melt
      !dfice = max(0.,0.7*omficetot*qml_ice/max(hofusn(hk)*hice,1e-8))
      dfice = max(0.,0.7*omficetot*dhicedt/max(avehice,1e-8))
      !write(6,'(f16.8)',advance='no') dfice*dtt
      !write(6,'(f16.8)',advance='no') dhice
      !write(6,'(f16.8)',advance='no') omficetot

      ! Reduction in lead area temp due to melting 
      tmlLLL= tmlLLL - dfice * fac * dtt

      ! Correct if resulting tmlLLL < freeze_temp
      dfice  = dfice - max(0.,freeze_temp(smlLLL)-tmlLLL)/(dtt*max(fac,1e-8))
      tmlLLL = max(tmlLLL,freeze_temp(smlLLL))
      !write(6,'(f16.8)',advance='no') dfice*dtt

      ! Lead melting 
      icem(hk)%fice=icem(hk)%fice-dfice*dtt

      ! Correct if we melted too much
      dfice  = min(icem(hk)%fice/dtt,0.)
      tmlLLL = tmlLLL - dfice * fac

      !write(6,'(f16.8)',advance='no') dfice*dtt
      !write(6,'(e16.8)',advance='yes') &
      !   ((tmlinit-tmlLLL)*omficetot*cpmlLLL-(ficeinit-icem(hk)%fice)*effhofusn)/&
      !   cpmlLLL

   end subroutine icestate_hmelt_hak


end module mod_icestate_hpar
