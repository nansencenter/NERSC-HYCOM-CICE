      module mod_daily_average
c --- Public Routine:
c ---   daily_average      -- main wrapper for routines 
c ---   daily_ave_combine  -- Combining daily averages across
c ---                         an ensemble (post processing)
c ---
c --- Private Routine:
c ---   daily_ave_dump     -- Logic for handling saving daily ave
c ---                         (reading previous daily average for
c ---                         ensemble runs, ++)
c ---   daily_average_save -- This routine does the actual saving
c ---   daily_average_read -- For reading previous averages
c ---   daily_ave_add      -- For accumulating averages in this run
c ---   daily_ave_ini      -- Initializes da variable
c ---   da_div_int         -- divides da variable with integer
c ---   da_plus_da         -- adds together two daily_averages 
c ---   da_div_dp          -- Divides da variables with layer thickness
c ---   dailyfilebase      -- Creates base file names for averages
c --- -------------------------------------------------------
c ---
c --- This module accumulates daily averages. When the 
c --- day switches the average is saved and a averaging for a new 
c --- day begins.
c ---
c --- Note that only daily_average should be called by the user. 
c --- The remaining routines are private and can only be accessed 
c --- from this module.
c ---
c --- Most variables are private, except for switches which are
c --- set in m_limits.
c ---
c --- To use the routine the cpp flag "DAILY_AVERAGE" must be set. 
c --- This is done to save memory if needed.
c ---
c --- NB: Currently there are NO restart capabilities. If you 
c --- stop hycom at the middle of the day, all averaging data
c --- for that day will be lost. 
c ---
c --- TODO: Make sure everything is ok for yrflag=3
c ---
c --- Created by Knut Liseter 07.05.2003 (KAL)
c ---
c KAL Changed 20.11.2003 -- makeshift conversion for hycom 2.1.03
c KAL Changed 07.07.2004 -- Suited for MPI-runs. Global arrays are
c KAL                       gone in place of local arrays.
c KAL Changed       2005 -- Set up for hycom-type IO (.[ab] files)
c KAL Changed 27.09.2008 -- Set up for hycom 2.2. No longer needs
c KAL                       the mod_year_info module, time logic
c KAL                       is handled by hycom routines
c --- Changed 19.11.2012 -- Added feature to bypass shared daily
c ---                       average file, instead write one file
c ---                       for each ensemble member and post process
c ---                       using daily_ave_combine. To use this
c ---                       feature define METNO flag.
c ---                       Geir Arne Waagbø (met.no)
c --- -------------------------------------------------------
      use mod_xc
#if defined(NOR05)
      use mod_necessary_ecovars, only: init, ipho, isil, idet, isis,
     &                                 ifla, idia, ioxy, ised, iyel,
     &                                 icha, idetp,imeso,imicro,grosspp,
     &                                 netpp
#endif /* NOR05 */
#if defined(ECOSM)
#if defined(ECO2)
      use mod_necessary_ecovars, only: ifla, idia, imicro, imeso,
     &                                 idet, inh4, idom, init, ipho,
     &                                 isil, ioxy, iopa, ibg, ifish, 
     &                                 inpp,idiap,iflap,izlp,
     &                                 izsp,ispp,igpp,idic,ialk,
     &                                 ipco,iph,iatc,bio_diagn,ndbio
#else
      use mod_necessary_ecovars, only: ifla, idia, imicro, imeso,
     &                                 idet, inh4, idom, init, ipho,
     &                                 isil, ioxy, iopa, ibg, ifish, 
     &                                 inpp,idiap,iflap,izlp,
     &                                 izsp,ispp,bio_diagn,ndbio
#endif
#endif /* ECOSM  */
      implicit none

      private
      logical, save,public :: l_daily_average
      logical, save,public :: l_daily_accum
#ifndef DAILY_AVERAGE
c --- If the CPP flag "DAILY_AVERAGE" is undefined, this is the
c --- only variable present in this module. 
      logical, parameter,public :: l_daily_average_in_code=.false.
#else
c --- See above
      logical, parameter,public :: l_daily_average_in_code=.true.
c
c
c --- This type is used in calculations
      type daily_averages
         integer         :: nrmem
         integer         :: lmem(1000) ! To decide which members are saved
         integer         :: initday, inityear
         integer         :: dumpday, dumpyear
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: 
     &      s,t,p,u,v,w
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     &      ub, vb, ssh, heatflux, swflx, salflx, emnp, taux, tauy
#if defined (DAILY_ENSVAR)
     &      ,ssh_sq, fice_sq, hice_sq
#endif
#if defined(ICE) || defined(ICESTATE)
     &      ,fice,hice,hsnw,srfheatflux,srfalb,htndncy
#if defined(ALBSNW_EVOL)
     &      ,srfsnw   
#endif         
#endif
#if defined(EVP) 
     &      ,uice,vice,wice,tauxice,tauyice
#if defined(ICE_DYN_DIAG)
     &      ,sigp,sigm,sig12,pr,strI,strII
#endif
#if defined(TEST_ICE_AGE)
     &      ,fy_age,fy_frac,rdg_frac
#endif
#endif
#if defined (PARAM_EST)
     &      ,msshb,sstb
#endif

#if defined(NOR05)
         real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,mxtrcr) ::
     &      trcr_ave     
         real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     &      grosspp_ave, netpp_ave 
#endif /* NOR05 */
#if defined(ECOSM)
         real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,mxtrcr) ::
     &      trcr_ave     
         real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,ndbio) ::
     &      diagn_ave
#endif /* ECOSM */
      end type daily_averages
c
c --- Keep daily average here, to avoid passing it  
c --- it to subroutines. 
      type (daily_averages),save :: da
c
c --- Bookkeeping variables
      integer, save :: da_counter = 0   !Increment counter
      integer, save :: day_banned =-1   !"Banned" day  
c
      public :: daily_average,daily_ave_combine
      contains


C --- --------------------------------------------------
C --- Main routine, calls the accumulation of fields, and
C --- decides when it is appropriate to save daily files
C --- --------------------------------------------------
      subroutine daily_average(dtime,m,restart)
      implicit none

      integer,          intent(in) :: m
      logical,          intent(in) :: restart
      real*8,           intent(in) :: dtime
      integer :: iyear, iday, ihour
      include 'common_blocks.h'
c
c --- Only run if daily average is switched on
      if (.not. l_daily_average) return
c
c --- Get day info
      call forday(dtime,yrflag,iyear,iday,ihour)
c --- hycom ordinal day starts at 1
      iday=iday-1
c
c --- Special case when restarting
      if (restart) then
         if  ( ihour > 1 ) then
c ---       We do not start averaging unless we have at least 
c ---       23 hours ahead of us so this day is "banned"
            day_banned=iday
         else
            da%dumpday=iday
         end if
         call daily_ave_ini()
         da%inityear=iyear
         da%initday =iday
      end if
c
c --- Skip this day if it has been banned above
      if (day_banned==iday) then
         if (mnproc==1)  print *,'daily_average:This day is "banned"'//
     &                           '-- no daily average'
      else  
c ---    We have just switched days
         if (da%dumpday/=iday) then
            if ( day_banned/=da%dumpday) then
               if (mnproc==1)print *,'daily: dump'
               call daily_ave_dump() ! NB .. da_rtprev
            end if
            call daily_ave_ini()
         end if
c ---    Always accumulate fields - no subsampling
         call daily_ave_add(m)
      end if
c
c --- Update previous day and time info
      da%dumpday =iday
      da%dumpyear=iyear
      end subroutine daily_average



C --- --------------------------------------------------
C --- Routine adds model fields to da fields. da is shared
C --- in this module.
C --- --------------------------------------------------
      subroutine daily_ave_add(m)
      use mod_xc
      use mod_forcing_nersc
      use mod_common_ice
#if defined (EVP)
      use mod_evp , only : evp_u=> uvel , evp_v => vvel,
     & evp_strairx => strairx, evp_strairy => strairy
#endif
      implicit none
      integer             , intent(in)    :: m
c
      integer :: i,j,k,l,ktr
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     & icestrairx,icestrairy,icestrocnx,icestrocny
      real ::ttaux,ttauy,utot,vtot
      include 'common_blocks.h'
      margin=0
c --- 3D fields
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i,utot,vtot)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
            da%s(i,j,k) = da%s(i,j,k) + saln(i,j,k,m)*dp(i,j,k,m)
            da%t(i,j,k) = da%t(i,j,k) + temp(i,j,k,m)*dp(i,j,k,m)
            da%p(i,j,k) = da%p(i,j,k) + p(i,j,k+1)
c ---       NB - total velocity
            utot=u(i,j,k,m) + ubavg(i,j,m)
            vtot=v(i,j,k,m) + vbavg(i,j,m)
            da%u(i,j,k) = da%u(i,j,k) + utot*dp(i,j,k,m)
            da%v(i,j,k) = da%v(i,j,k) + vtot*dp(i,j,k,m)
c ---       Used for kinetic energy
            da%w(i,j,k) = da%w(i,j,k) + 
     &         (vtot*vtot+utot*utot)*dp(i,j,k,m)
         end do
         end do
         end do
      end do
#if defined (NOR05)
c --- Ecosystem
      do ktr=1,mxtrcr
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
            da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,k,ktr)
     &                            +tracer(i,j,k,m,ktr)*dp(i,j,k,m)
         end do
         end do
         end do
!$OMP END PARALLEL DO
      end do
      end do
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
            da%grosspp_ave(i,j,k)=da%grosspp_ave(i,j,k)
     &                            +grosspp(i,j,k)*dp(i,j,k,m)
            da%netpp_ave(i,j,k)=da%netpp_ave(i,j,k)
     &                            +netpp(i,j,k)*dp(i,j,k,m)
         end do
         end do
         end do
!$OMP END PARALLEL DO
      end do
cdiag print *,'daily_ave_add after 3D ecosyetem flds -- mnproc=', mnproc
#endif /* NOR05 */
#if defined (ECOSM)
c --- Ecosystem
      do ktr=1,mxtrcr
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
           if(ktr.eq.init.or.ktr.eq.inh4)then
            da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,k,ktr)
     &          +tracer(i,j,k,m,ktr)*1./6.625*1./12.01*dp(i,j,k,m)
           elseif(ktr.eq.ipho)then
            da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,k,ktr)
     &          +tracer(i,j,k,m,ktr)*1./106.*1./12.01*dp(i,j,k,m)
           elseif(ktr.eq.isil)then
            da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,k,ktr)
     &          +tracer(i,j,k,m,ktr)*1./6.625*1./12.01*dp(i,j,k,m)
           elseif(ktr.eq.ioxy)then
            da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,k,ktr)
     &          +tracer(i,j,k,m,ktr)*1./44.6608009*dp(i,j,k,m)
#if defined(ECO2)
           elseif(ktr.eq.idic)then
            da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,k,ktr)
     &          +tracer(i,j,k,m,ktr)*1./12.01*dp(i,j,k,m)
#endif
           else 
            da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,k,ktr)
     &                            +tracer(i,j,k,m,ktr)*dp(i,j,k,m)
           endif
         end do
         end do
         end do
!$OMP END PARALLEL DO
      end do
      end do
      do ktr=1,ndbio
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
          da%diagn_ave(i,j,k,ktr)=da%diagn_ave(i,j,k,ktr)
     &                      +bio_diagn(i,j,k,ktr)*dp(i,j,k,m)
#if defined(ECO2)
           if(ktr.eq.iatc)then
            da%diagn_ave(i,j,k,ktr)=da%diagn_ave(i,j,k,ktr)
     &          +bio_diagn(i,j,k,ktr)
           endif
#endif
         end do
         end do
         end do
!$OMP END PARALLEL DO
      end do
      end do
cdiag print *,'daily_ave_add after 3D ecosyetem flds -- mnproc=', mnproc
#endif /* ECOSM */
#if defined(EVP)
      icestrairx=evp_strairx
      icestrairy=evp_strairy
      icestrocnx=tauxice
      icestrocny=tauyice
#endif
c
c --- 2D fields
!$OMP  PARALLEL DO PRIVATE(j,i,ttaux,ttauy)
!$OMP& SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
         ttaux     =taux(i,j,l0)*w0+taux(i,j,l1)*w1
     &             +taux(i,j,l2)*w2+taux(i,j,l3)*w3
         ttauy     =tauy(i,j,l0)*w0+tauy(i,j,l1)*w1
     &             +tauy(i,j,l2)*w2+tauy(i,j,l3)*w3
#if defined(EVP)
c ---    Modify by letting taux be the net stress acting on the
c ---    combined ice/ocean surface
         ttaux     = ttaux * (1.-ficem(i,j)) + 
     &               ficem(i,j) * icestrairx(i,j)
         ttauy     = ttauy * (1.-ficem(i,j)) + 
     &               ficem(i,j) * icestrairy(i,j)
#endif
c
         da%vb      (i,j) = da%vb   (i,j)   + vbavg(i,j,m)  
         da%ub      (i,j) = da%ub   (i,j)   + ubavg(i,j,m)  
         da%ssh     (i,j) = da%ssh  (i,j)   + srfhgt (i,j)/g
         da%taux    (i,j) = da%taux(i,j)    + ttaux
         da%tauy    (i,j) = da%tauy(i,j)    + ttauy
         da%heatflux(i,j) = da%heatflux(i,j)+ surflx (i,j)   ! Heat flux into ocean
#if defined (DAILY_ENSVAR)
         da%ssh_sq  (i,j) = da%ssh_sq  (i,j)   + (srfhgt  (i,j)**2/g**2)
         da%fice_sq (i,j) = da%fice_sq (i,j)   + ficem(i,j)**2
         da%hice_sq (i,j) = da%hice_sq (i,j)   + hicem(i,j)**2
#endif

#if defined(ICE)
c ---    NB: this is downward heat flux into ice and ocean - because
c ---    of the simplicity of the standard ice module, this can be
c ---    formulated as follows:
         da%srfheatflux(i,j) = da%srfheatflux(i,j)
     &   + surflx(i,j) 
     &     - ( delta_icevol(i,j)*fusi 
     &        +delta_snwvol(i,j)*fuss ) / baclin
#elif defined(ICESTATE)
c#error " srfheatflux not def for ICESTATE"
         da%srfheatflux(i,j) = da%srfheatflux(i,j)+0.0 

#endif
#if defined(ICE) || defined(ICESTATE)
         da%swflx   (i,j) = da%swflx(i,j)+ surf_qsw_sum   (i,j)  !  Average surface sw heat flux (area average)
#else
         da%swflx   (i,j) = da%swflx(i,j)   + sswflx (i,j)  ! SW heat flux into OCEAN
#endif
         da%salflx  (i,j) = da%salflx(i,j)  + salflx (i,j)   ! Sal flux into ocean
         da%emnp    (i,j) = da%emnp  (i,j)  + pemnp  (i,j)  / onem ! sal flx interpreted as EMNP
#if defined(ICE) || defined(ICESTATE)
         da%hice(i,j)= da%hice(i,j)+ hicem(i,j)
         da%fice(i,j)= da%fice(i,j)+ ficem(i,j)
         da%hsnw(i,j)= da%hsnw(i,j)+ hsnwm(i,j)
         da%srfalb  (i,j) = da%srfalb  (i,j)+ surf_albedo_sum(i,j)  ! Average albedo (surface average)
         da%htndncy  (i,j) = 
     &      da%htndncy  (i,j)+ delta_icevol(i,j)/baclin  ! Ice growth tendency (surface average)
#if defined(ALBSNW_EVOL)
         da%srfsnw(i,j) = da%srfsnw(i,j) + surf_fscov(i,j) ! Snow cover  fraction
#endif         
#endif
#if defined(EVP)
         da%uice(i,j)= da%uice(i,j)+ iceu(i,j)
         da%vice(i,j)= da%vice(i,j)+ icev(i,j)
         da%wice(i,j)= da%wice(i,j)+ sqrt(icev(i,j)**2 + iceu(i,j)**2)
         da%tauxice(i,j)= da%tauxice(i,j)+icestrocnx(i,j)
         da%tauyice(i,j)= da%tauyice(i,j)+icestrocny(i,j)
#if defined(ICE_DYN_DIAG)
         da%sigp(i,j)= da%sigp(i,j)+stressp(i,j)
         da%sigm(i,j)= da%sigm(i,j)+stressm(i,j)
         da%sig12(i,j)= da%sig12(i,j)+stress12(i,j)
         da%pr(i,j)= da%pr(i,j)+pice(i,j)
         da%strI(i,j)= da%strI(i,j)+strainI(i,j)
         da%strII(i,j)= da%strII(i,j)+strainII(i,j)
#endif
#if defined(TEST_ICE_AGE)
         da%fy_frac (i,j) = da%fy_frac (i,j) + fy_frac (i,j)
         da%fy_age  (i,j) = da%fy_age  (i,j) + fy_age  (i,j)
         da%rdg_frac(i,j) = da%rdg_frac(i,j) + rdg_frac(i,j)
#endif
#if defined(PARAM_EST)
         da%sstb (i,j) = da%sstb (i,j) + sstb (i,j)
         da%msshb(i,j) = da%msshb(i,j) + msshb(i,j)
#endif


#endif
      end do
      end do
      end do
c
c --- Augment counter
      da_counter=da_counter+1
      end subroutine daily_ave_add


C --- --------------------------------------------------
C --- Routine initializes da fields. da variable is shared 
C --- in this module
C --- --------------------------------------------------
      subroutine daily_ave_ini()
      implicit none
      integer :: i,j,k,ktr
c Fanf:Initialisation so we do the whole table
      margin=0
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP& SCHEDULE(STATIC,jblk)
         do j=1-margin,jdm+margin
         do i=1-margin,idm+margin
            da%s(i,j,k) = 0.
            da%t(i,j,k) = 0.
            da%p(i,j,k) = 0.
            da%u(i,j,k) = 0.
            da%v(i,j,k) = 0.
            da%w(i,j,k) = 0.
         end do
         end do
      end do
#if defined (NOR05)
c Fanf:Initialisation so we do the whole table
      do ktr=1,mxtrcr
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP& SCHEDULE(STATIC,jblk)
         do j=1-margin,jdm+margin
         do i=1-margin,idm+margin
            da%trcr_ave(i,j,k,ktr) = 0.
         end do
         end do
      end do
      end do
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP& SCHEDULE(STATIC,jblk)
         do j=1-margin,jdm+margin
         do i=1-margin,idm+margin
            da%grosspp_ave(i,j,k) = 0.
            da%netpp_ave(i,j,k) = 0.
         end do
         end do
!OMP END PARALLEL DO
      end do
#endif /* NOR05 */
#if defined (ECOSM)
c Fanf:Initialisation so we do the whole table
      do ktr=1,mxtrcr
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP& SCHEDULE(STATIC,jblk)
         do j=1-margin,jdm+margin
         do i=1-margin,idm+margin
            da%trcr_ave(i,j,k,ktr) = 0.
         end do
         end do
      end do
      end do
      do ktr=1,ndbio
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP& SCHEDULE(STATIC,jblk)
         do j=1-margin,jdm+margin
         do i=1-margin,idm+margin
            da%diagn_ave(i,j,k,ktr) = 0.
         end do
         end do
         end do
!OMP END PARALLEL DO
      end do
#endif /* ECOSM */
c
c Fanf:Initialisation so we do the whole table
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP& SCHEDULE(STATIC,jblk)
      do j=1-margin,jdm+margin
      do i=1-margin,idm+margin
         da%vb      (i,j)  = 0.
         da%ub      (i,j)  = 0.
         da%ssh     (i,j)  = 0.
#if defined (DAILY_ENSVAR)
         da%ssh_sq   (i,j)  = 0.
         da%fice_sq  (i,j)  = 0.
         da%hice_sq  (i,j)  = 0.
#endif
         da%taux    (i,j)  = 0.
         da%tauy    (i,j)  = 0.
         da%heatflux(i,j)  = 0.
         da%swflx   (i,j)  = 0.
         da%salflx  (i,j)  = 0.
         da%emnp    (i,j)  = 0.
#if defined(ICE) || defined(ICESTATE)
         da%hice(i,j)= 0.
         da%fice(i,j)= 0.
         da%hsnw(i,j)= 0.
         da%srfalb(i,j)= 0.
         da%htndncy(i,j)= 0.
         da%srfheatflux(i,j)= 0.
#if defined(ALBSNW_EVOL)        
         da%srfsnw(i,j)= 0.
#endif
#endif
#if defined(EVP)
         da%uice(i,j)= 0.
         da%vice(i,j)= 0.
         da%wice(i,j)= 0.
         da%tauxice(i,j)= 0.
         da%tauyice(i,j)= 0.
#if defined(ICE_DYN_DIAG)
         da%sigp(i,j) = 0.
         da%sigm(i,j) = 0.
         da%sig12(i,j) = 0.
         da%pr(i,j) = 0.
         da%strI(i,j) = 0.
         da%strII(i,j) = 0.
#endif
#if defined(TEST_ICE_AGE)
         da%fy_age(i,j)= 0.
         da%fy_frac(i,j)= 0.
         da%rdg_frac(i,j)= 0.
#endif
#if defined(PARAM_EST)
         da%sstb (i,j)= 0.
         da%msshb(i,j)= 0.
#endif
#endif
      end do
      end do
c
      da_counter=0
      end subroutine daily_ave_ini




C --- --------------------------------------------------
C --- Dumping of average fields -- with a crude implementation 
C --- of a file lock
C --- --------------------------------------------------
      subroutine daily_ave_dump()
      use mod_hycom_nersc
      implicit none
c
      integer ::j,iostat,lmem(1000),numtries, ios
      character(len=80) :: lockname
      type (daily_averages) :: da_tmp
      logical :: ex,ex2,islocked
      character(len=60) :: filenamea,filenameb, filebase
      integer*4, parameter :: sleeptime=30
      real :: rios(1)
      integer :: irecl
#if defined (EVP_DRIFT_ASSIMILATION)
      character(len=60) :: ficedrft
      real*4, dimension(itdm,jtdm) :: dmpu,dmpv
      real  , dimension(itdm,jtdm) :: tmpuv
#endif
      include 'common_blocks.h'
#if defined (TRACK_ASSIM)
      character(len=60) :: trackfile
      real*4, dimension(itdm,jtdm) :: dmpssh
      real  , dimension(itdm,jtdm) :: tmpssh
#endif
#if defined (PROBACAST)
      character(len=60) :: ufile
      character(len=60) :: ficefile
      real*4, dimension(itdm,jtdm) :: dmput,dmpvt
      real*4, dimension(itdm,jtdm) :: dmpfice
      real  , dimension(itdm,jtdm) :: tmpuvt
      real  , dimension(itdm,jtdm) :: tmpfice
#endif
c --- Divide by counter
      call da_div_int(da_counter)
c 
c --- Divide by dp
      call da_div_dp()
c
c --- Get average
      da%lmem    = 0 ; da%lmem(imem) = 1
      lmem=da%lmem
      da%nrmem   = 1  ! By default, we have one member pr record
c
c --- Generate file lock and file names
      filebase=dailyfilebase(da%inityear, da%initday, 
     &                       da%dumpyear, da%dumpday)
      lockname =trim(filebase)//'.lock'
      filenamea=trim(filebase)//'.a'
      filenameb=trim(filebase)//'.b'
c
#if defined (EVP_DRIFT_ASSIMILATION)
c --- This is needed to produce drift trajectories
c --- from the EVP model
      ficedrft=trim(filebase)//'_ICEDRIFT.uf'
      call xcaget(tmpuv,da%uice,0) ; dmpu=tmpuv
      call xcaget(tmpuv,da%vice,0) ; dmpv=tmpuv
      if (mnproc==1) then
         print *,'Dumping EVP drift for assimilation'
         inquire(iolength=irecl)dmpu,dmpv
         open(20,file=trim(ficedrft),status='unknown',
     &           form='unformatted',recl=irecl,access='direct')
c ---    Always replace record
         write(20,rec=imem) dmpu, dmpv
         close(20)
      end if
#endif
#if defined (PROBACAST)
      ufile=trim(dailyfilebase(da%inityear, da%initday,
     &                         da%dumpyear, da%dumpday))//
     &   '_UV01.uf'
      call xcaget(tmpuvt,da%u,0) ; dmput=tmpuvt
      call xcaget(tmpuvt,da%v,0) ; dmpvt=tmpuvt
      if (mnproc==1) then
         print *,'Dumping Ut and Vt for probacast'
         inquire(iolength=irecl)dmput,dmpvt
         open(20,file=trim(ufile),status='unknown',
     &           form='unformatted',recl=irecl,access='direct')
         ! Always overwrite
         write(20,rec=imem) dmput, dmpvt
         close(20)
      end if
       ficefile=trim(dailyfilebase(da%inityear, da%initday,
     &                             da%dumpyear, da%dumpday))//
     &   '_FICE.uf'
      call xcaget(tmpfice,da%fice,0) ; dmpfice=tmpfice
      if (mnproc==1) then
         print *,'Dumping fice for probacast'
         inquire(iolength=irecl)dmpfice
         open(20,file=trim(ficefile),status='unknown',
     &           form='unformatted',recl=irecl,access='direct')
         ! Always overwrite
         write(20,rec=imem) dmpfice
         close(20)
      end if
#endif
#if defined (TRACK_ASSIM)
      ! We store daily the model output (here SSH)
      trackfile=trim(filebase)//'_TSSH.uf'
      call xcaget(tmpssh,da%ssh,0) 
      dmpssh=tmpssh
      if (mnproc==1) then
         print *,'Dumping SSH for Track assimilation'
         inquire(iolength=irecl)dmpssh
         open(33,file=trim(trackfile),status='unknown',
     &           form='unformatted',recl=irecl,access='direct')
         ! Always overwrite
         write(33,rec=imem) dmpssh
         close(33)
      end if
#endif
#if defined (METNO)
c --- With the METNO flag: bypass file locking mechanism
c --- and shared daily mean files, it creates sporadic
c --- disturbed fields for F07.
c ---
c --- Instead: save individual files for each ensemble member and
c --- postprocess using the subroutine daily_ave_combine
      call daily_ave_save(imem)
#else
c
c --- Crude file locking mechanism -- holds until lock file is removed 
c --- A safer (and less verbose) approach may be to use system-level 
c --- file locks (eg NFS)
      islocked=.true.
      numtries=0
      if (mnproc==1) then
         print *,filenamea
         print *,filenameb
         write (lp,'(a)')'daily_ave_dump: Waiting for file lock...'
         write (lp,'(a)')'daily_ave_dump: lockfile is '//trim(lockname)
      end if
c
c --- Loop while file is locked
      do while (islocked)
         numtries=numtries+1
c ---    Abort on too many tries
         if (numtries>20) then
            if (mnproc==1) then
            write (lp,'(a,i5,a)') 'daily_ave_dump:'//
     &      ' Tried ',numtries,' times but no success'//
     &      ' in opening lockfile. Write for this member is skipped '
            end if
            return
         end if
c
c ---    If file lock exists, file is in use by other executable
         call xcsync(flush_lp)
         inquire(exist=islocked,file=trim(lockname))
         call xcsync(flush_lp)
         if (islocked) then
            if (mnproc==1) then
               write (lp,'(a,i4)')'daily_ave_dump:'//
     &         ' Waiting for file lock.. attempt number ', numtries
            end if
c ---       Sleep for a while before trying again
#if defined(AIX)
            call sleep_(sleeptime)
#else
            call sleep(sleeptime)
#endif
            cycle
         end if
c
c ---    If file doesnt exist, prepare for writing. Communicate 
c ---    the result to all cpus.
         ios=0
         if (mnproc==1) 
     &      open (11,file=trim(lockname),status='new',iostat=ios)
         rios(1)=ios
         call xcastr(rios,1)
         ios=nint(rios(1))
c
c ---    Opening lock file failed - try again
         if (ios/=0) then
            islocked =.true. ! Stay in loop
            if (mnproc==1) then
               write (lp,'(a,i4)')'daily_ave_dump:'//
     &         ' Rechecking file lock.. attempt number ', numtries
            end if
#if defined(AIX)
            call sleep_(sleeptime)
#else
            call sleep(sleeptime)
#endif
            cycle
c ---    Open succeeded - unset lock flag (exit loop)
         else
            islocked=.false.
         end if
      end do ! if islocked
c
c --- At this stage the file should be free of locks, and
c --- we are ready for reading/writing
      if (mnproc==1) then
         write(11,'(i4)') imem
      end if
c
c --- Find unformatted file(s) -- this finds the .b-file
      inquire(exist=ex2,file=trim(filenameb))
c
c --- Must be in sync brefore proceeding
      call xcsync(flush_lp)
c
c --- If old file, we must read data from it. Override if daily_accum is switched off
      if (ex2.and.l_daily_accum) then
         if(mnproc==1)write(lp,'(a)') 'daily_ave_dump: File exists'
c
c ---    Read old daily average from file
         call daily_ave_read(da_tmp,da%inityear, da%initday,
     &                              da%dumpyear, da%dumpday)
c
c ---    NB -- lmem holds members we have dumped. If 
c ---    lmem is not equal to 0, we have already dumped it.
c ---    In that case skip writing. 
c ---    (This comes in handy for a restart etc...)
         if (da_tmp%lmem(imem) == 0 ) then
            if (mnproc==1)
     &         write(*,'(a)')'daily_ave_dump: Adding to existing member'
            call da_plus_da(da,da_tmp)
            call daily_ave_save()
         else
            if (mnproc==1)
     &         write(lp,'(a)') 'daily_ave_dump: WARNING this '//
     &              ' member already dumped -- item skipped'
            lmem=da_tmp%lmem
         end if
c
c --- New file  first write
      else
         if (mnproc==1 .and. .not. ex2) 
     &      write(lp,'(a)') 'daily_ave_dump: File does not exist'
         if (mnproc==1 .and. .not. l_daily_accum) 
     &      write(lp,'(a)') 'daily_ave_dump: Accumulation is off'
         lmem=da%lmem
         call daily_ave_save()
      end if
c
c --- Close lock file
      if (mnproc==1) then
         close(11)
         write(lp,'(a,i4)')'daily_ave_dump:Total number of members ',
     &         sum(lmem)
      end if
c
c --- When finished remove lockfile
      if (mnproc==1) then
        call unlink_wrap(trim(lockname)//char(0))
        write(lp,'(a,21i2)') 
     &    'daily_ave_dump: lmem around this member:',
     &    lmem( (max(1,imem-10)):min(imem+10,1000))
      end if
#endif /* METNO */
      end subroutine daily_ave_dump


c --- ----------------------------------------------
c --- Divide da variables with counter. da is global 
c --- in this module
c --- ----------------------------------------------
      subroutine da_div_int(divisor)
      implicit none
      integer             , intent(in   ) :: divisor
      integer :: i,j,k,ktr
c --- 3D variables
      margin=0
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         da%s(i,j,k)       = da%s(i,j,k) /divisor
         da%t(i,j,k)       = da%t(i,j,k) /divisor
         da%p(i,j,k)       = da%p(i,j,k) /divisor
         da%u(i,j,k)       = da%u(i,j,k) /divisor 
         da%v(i,j,k)       = da%v(i,j,k) /divisor 
         da%w(i,j,k)       = da%w(i,j,k) /divisor 
      end do
      end do
      end do
#if defined (NOR05)
      do ktr=1,mxtrcr
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
            da%trcr_ave(i,j,k,ktr)=
     &      da%trcr_ave(i,j,k,ktr) / divisor
      end do
      end do
      end do
      end do
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
            da%grosspp_ave(i,j,k)=
     &      da%grosspp_ave(i,j,k) / divisor
            da%netpp_ave(i,j,k)=
     &      da%netpp_ave(i,j,k) / divisor
      end do
      end do
      end do
#endif /* NOR05 */
#if defined (ECOSM)
      do ktr=1,mxtrcr
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
            da%trcr_ave(i,j,k,ktr)=
     &      da%trcr_ave(i,j,k,ktr) / divisor
      end do
      end do
      end do
      end do
      do ktr=1,ndbio
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
            da%diagn_ave(i,j,k,ktr)=
     &      da%diagn_ave(i,j,k,ktr) / divisor
      end do
      end do
      end do
      end do
#endif /* ECOSM */
c --- 2D variables
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         da%ub      (i,j)= da%ub      (i,j)/divisor
         da%vb      (i,j)= da%vb      (i,j)/divisor
         da%ssh     (i,j)= da%ssh     (i,j)/divisor
#if defined (DAILY_ENSVAR)
         da%ssh_sq  (i,j)= da%ssh_sq  (i,j)/divisor
         da%fice_sq (i,j)= da%fice_sq (i,j)/divisor
         da%hice_sq (i,j)= da%hice_sq (i,j)/divisor
#endif
         da%heatflux(i,j)= da%heatflux(i,j)/divisor
         da%swflx   (i,j)= da%swflx   (i,j)/divisor
         da%salflx  (i,j)= da%salflx  (i,j)/divisor
         da%emnp    (i,j)= da%emnp    (i,j)/divisor
         da%taux    (i,j)= da%taux    (i,j)/divisor
         da%tauy    (i,j)= da%tauy    (i,j)/divisor
#if defined(ICE) || defined(ICESTATE)
         da%fice       (i,j)= da%fice       (i,j)/divisor
         da%hice       (i,j)= da%hice       (i,j)/divisor
         da%hsnw       (i,j)= da%hsnw       (i,j)/divisor
         da%srfalb     (i,j)= da%srfalb     (i,j)/divisor
         da%htndncy    (i,j)= da%htndncy    (i,j)/divisor
         da%srfheatflux(i,j)= da%srfheatflux(i,j)/divisor
#if defined(ALBSNW_EVOL)   
         da%srfsnw     (i,j)= da%srfsnw     (i,j)/divisor
#endif
#endif
#if defined(EVP)
         da%uice      (i,j)= da%uice   (i,j)  /divisor
         da%vice      (i,j)= da%vice   (i,j)  /divisor
         da%wice      (i,j)= da%wice   (i,j)  /divisor
         da%tauxice   (i,j)= da%tauxice(i,j)  /divisor
         da%tauyice   (i,j)= da%tauyice(i,j)  /divisor
#if defined(ICE_DYN_DIAG)
         da%sigp   (i,j)   = da%sigp   (i,j)  /divisor
         da%sigm   (i,j)   = da%sigm   (i,j)  /divisor
         da%sig12  (i,j)   = da%sig12  (i,j)  /divisor
         da%pr     (i,j)   = da%pr     (i,j)  /divisor
         da%strI   (i,j)   = da%strI   (i,j)  /divisor
         da%strII  (i,j)   = da%strII  (i,j)  /divisor
#endif
#if defined(TEST_ICE_AGE)
         da%fy_age      (i,j)= da%fy_age   (i,j)  /divisor
         da%fy_frac     (i,j)= da%fy_frac  (i,j)  /divisor
         da%rdg_frac    (i,j)= da%rdg_frac (i,j)  /divisor
#endif
#if defined(PARAM_EST)
         da%msshb     (i,j)= da%msshb  (i,j)  /divisor
         da%sstb      (i,j)= da%sstb   (i,j)  /divisor
#endif
#endif
      end do
      end do
      end subroutine da_div_int


c --- ----------------------------------------------
c --- Divide 3D variables with accumulated layer  thickness.
c --- da is shared variable in this module.
c --- ----------------------------------------------
      subroutine da_div_dp()
      implicit none
      real :: tdp
      integer :: lastmass(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer i,j,k,l,klast,ktr
      include 'common_blocks.h'
c
      margin=0
      lastmass=0
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i,tdp,klast)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
            klast=lastmass(i,j)
            if (k==1) then
               tdp= da%p(i,j,k) 
            else
               tdp= da%p(i,j,k)-da%p(i,j,k-1) 
            end if
            da%s(i,j,k)= da%s(i,j,k)/(tdp+onemm)
            da%t(i,j,k)= da%t(i,j,k)/(tdp+onemm)
            da%u(i,j,k)= da%u(i,j,k)/(tdp+onemm) 
            da%v(i,j,k)= da%v(i,j,k)/(tdp+onemm) 
            da%w(i,j,k)= da%w(i,j,k)/(tdp+onemm) 
            if (tdp<tencm .and. klast>0) then
               ! Fill with sensible values
               da%s(i,j,k)=da%s(i,j,klast)
               da%t(i,j,k)=da%t(i,j,klast)
               da%u(i,j,k)=da%u(i,j,klast)
               da%v(i,j,k)=da%v(i,j,klast)
               da%w(i,j,k)=da%w(i,j,klast)
            elseif (tdp<tencm .and. klast==0) then
               ! Fill with zerooo values
               da%s(i,j,k)=0.
               da%t(i,j,k)=0.
               da%u(i,j,k)=0.
               da%v(i,j,k)=0.
               da%w(i,j,k)=0.
            else
               ! Values already ok, and this is a mass-
               ! containing layer
               lastmass(i,j)=k
            end if
         end do
         end do
         end do
      end do
c
c --- ecosys 3D variables
#if defined (NOR05)
      do ktr=1,ntracr
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i,tdp)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
            klast=lastmass(i,j)
            if (k==1) then
               tdp= da%p(i,j,k) 
            else
               tdp= da%p(i,j,k)-da%p(i,j,k-1) 
            end if
            da%trcr_ave(i,j,k,ktr)=
     &      da%trcr_ave(i,j,k,ktr) / (tdp+onemm)
            if (tdp<tencm .and. klast>0) then
               ! Fill with sensible values
               da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,klast,ktr)
            elseif (tdp<tencm .and. klast==0) then
               ! Fill with zerooo values
               da%trcr_ave(i,j,k,ktr)=0
            else
               ! Values already ok, and this is a mass-
               ! containing layer
               lastmass(i,j)=k
            end if
         end do
         end do
         end do
      end do
      end do
!!primary production
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i,tdp)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
            klast=lastmass(i,j)
            if (k==1) then
               tdp= da%p(i,j,k) 
            else
               tdp= da%p(i,j,k)-da%p(i,j,k-1) 
            end if
            da%grosspp_ave(i,j,k)=
     &      da%grosspp_ave(i,j,k) / (tdp+onemm)
            da%netpp_ave(i,j,k)=
     &      da%netpp_ave(i,j,k) / (tdp+onemm)
            if (tdp<tencm .and. klast>0) then
               ! Fill with sensible values
               da%grosspp_ave(i,j,k)=da%grosspp_ave(i,j,klast)
               da%netpp_ave(i,j,k)=da%netpp_ave(i,j,klast)
            elseif (tdp<tencm .and. klast==0) then
               ! Fill with zerooo values
               da%grosspp_ave(i,j,k)=0
               da%netpp_ave(i,j,k)=0
            else
               ! Values already ok, and this is a mass-
               ! containing layer
               lastmass(i,j)=k
            end if
         end do
         end do
         end do
!$OMP END PARALLEL DO
      end do
#endif /* NOR05 */
#if defined (ECOSM)
      do ktr=1,ntracr
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i,tdp)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
            klast=lastmass(i,j)
            if (k==1) then
               tdp= da%p(i,j,k) 
            else
               tdp= da%p(i,j,k)-da%p(i,j,k-1) 
            end if
            da%trcr_ave(i,j,k,ktr)=
     &      da%trcr_ave(i,j,k,ktr) / (tdp+onemm)
            if (tdp<tencm .and. klast>0) then
               ! Fill with sensible values
               da%trcr_ave(i,j,k,ktr)=da%trcr_ave(i,j,klast,ktr)
            elseif (tdp<tencm .and. klast==0) then
               ! Fill with zerooo values
               da%trcr_ave(i,j,k,ktr)=0
            else
               ! Values already ok, and this is a mass-
               ! containing layer
               lastmass(i,j)=k
            end if
         end do
         end do
         end do
      end do
      end do
!!primary production
      do ktr=1,ndbio
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(j,l,i,tdp)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
            klast=lastmass(i,j)
            if (k==1) then
               tdp= da%p(i,j,k) 
            else
               tdp= da%p(i,j,k)-da%p(i,j,k-1) 
            end if
            da%diagn_ave(i,j,k,ktr)=
     &da%diagn_ave(i,j,k,ktr) / (tdp+onemm)
            if (tdp<tencm .and. klast>0) then
               ! Fill with sensible values
            da%diagn_ave(i,j,k,ktr)=da%diagn_ave(i,j,klast,ktr)
            elseif (tdp<tencm .and. klast==0) then
               ! Fill with zerooo values
            da%diagn_ave(i,j,k,ktr)=0
            else
               ! Values already ok, and this is a mass-
               ! containing layer
               lastmass(i,j)=k
            end if
         end do
         end do
         end do
!$OMP END PARALLEL DO
      end do
      end do
#endif /* ECOSM */
      end subroutine da_div_dp
 
 
c --- ----------------------------------------------
c --- Add two daily averages. Here we do not use the
c --- shared variable da but pass as arguments
c --- ----------------------------------------------
      subroutine da_plus_da(da1,da2)
      implicit none
      type(daily_averages), intent(inout) :: da1
      type(daily_averages), intent(in   ) :: da2
      integer :: i,j,k, ktr
c --- Consistency checks on dates, they must match
      if (da1%inityear == da2%inityear  .and. 
     &    da1%initday  == da2%initday ) then
         continue
      else
         if (mnproc==1) then
            print *,'(Error on da adding rtinit)'
            print *,da1%inityear,da1%initday
            print *,da2%inityear,da2%initday
         end if
         call xcstop('(mod_daily_average:da_plus_da)')
         stop '(mod_daily_average:da_plus_da)'
      end if
      if (da1%dumpyear   == da2%dumpyear    .and. 
     &    da1%dumpday    == da2%dumpday   ) then
         continue
      else
         if (mnproc==1) then 
            print *,'(Error on da adding rtdump)'
            print *,da1%dumpyear,da1%dumpday
            print *,da2%dumpyear,da2%dumpday
         end if
         call xcstop('(mod_daily_average:da_plus_da)')
         stop '(mod_daily_average:da_plus_da)'
      end if
c --- Update member indicators and total # of members
      da1%nrmem   = da1%nrmem   + da2%nrmem
      da1%lmem    = da1%lmem    + da2%lmem 
c
      margin=0
c --- Sum 3D ocean fields
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         da1%s(i,j,k)       = da1%s(i,j,k) + da2%s(i,j,k)
         da1%t(i,j,k)       = da1%t(i,j,k) + da2%t(i,j,k)
         da1%p(i,j,k)       = da1%p(i,j,k) + da2%p(i,j,k)
         da1%u(i,j,k)       = da1%u(i,j,k) + da2%u(i,j,k)
         da1%v(i,j,k)       = da1%v(i,j,k) + da2%v(i,j,k)
         da1%w(i,j,k)       = da1%w(i,j,k) + da2%w(i,j,k)
      end do
      end do
      end do
c
#if defined (NOR05)
c --- Sum 3D ecosys fields
      do ktr=1,mxtrcr
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
            da1%trcr_ave(i,j,k,ktr)=
     &      da1%trcr_ave(i,j,k,ktr)+
     &      da2%trcr_ave(i,j,k,ktr)
      end do
      end do
      end do
      end do
!primary production
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
            da1%grosspp_ave(i,j,k)=
     &      da1%grosspp_ave(i,j,k)+
     &      da2%grosspp_ave(i,j,k)
            da1%netpp_ave(i,j,k)=
     &      da1%netpp_ave(i,j,k)+
     &      da2%netpp_ave(i,j,k)
      end do
      end do
      end do
#endif /* NOR05 */
#if defined (ECOSM)
c --- Sum 3D ecosys fields
      do ktr=1,mxtrcr
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
            da1%trcr_ave(i,j,k,ktr)=
     &      da1%trcr_ave(i,j,k,ktr)+
     &      da2%trcr_ave(i,j,k,ktr)
      end do
      end do
      end do
      end do
!primary production
      do ktr=1,ndbio
      do k=1,kdm
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
            da1%diagn_ave(i,j,k,ktr)=
     &      da1%diagn_ave(i,j,k,ktr)+
     &      da2%diagn_ave(i,j,k,ktr)
      end do
      end do
      end do
      end do
#endif /* ECOSM */
c --- Sum various 2D fields
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         da1%ub      (i,j)= da1%ub      (i,j)+ da2%ub      (i,j)
         da1%vb      (i,j)= da1%vb      (i,j)+ da2%vb      (i,j)
         da1%ssh     (i,j)= da1%ssh     (i,j)+ da2%ssh     (i,j)
#if defined (DAILY_ENSVAR)
         da1%ssh_sq   (i,j)= da1%ssh_sq   (i,j)+ da2%ssh_sq   (i,j) 
         da1%fice_sq  (i,j)= da1%fice_sq  (i,j)+ da2%fice_sq  (i,j) 
         da1%hice_sq  (i,j)= da1%hice_sq  (i,j)+ da2%hice_sq  (i,j) 
#endif
         da1%heatflux(i,j)= da1%heatflux(i,j)+ da2%heatflux(i,j)
         da1%swflx   (i,j)= da1%swflx   (i,j)+ da2%swflx   (i,j)
         da1%salflx  (i,j)= da1%salflx  (i,j)+ da2%salflx  (i,j)
         da1%emnp    (i,j)= da1%emnp    (i,j)+ da2%emnp    (i,j)
         da1%taux    (i,j)= da1%taux    (i,j)+ da2%taux    (i,j)
         da1%tauy    (i,j)= da1%tauy    (i,j)+ da2%tauy    (i,j)
#if defined(ICE) || defined(ICESTATE)
         da1%fice       (i,j)=da1%fice (i,j)      +da2%fice(i,j)
         da1%hice       (i,j)=da1%hice (i,j)      +da2%hice(i,j)
         da1%hsnw       (i,j)=da1%hsnw (i,j)      +da2%hsnw(i,j)
         da1%srfalb     (i,j)=da1%srfalb(i,j)     +da2%srfalb(i,j)
         da1%htndncy    (i,j)=da1%htndncy(i,j)    +da2%htndncy(i,j)
         da1%srfheatflux(i,j)=da1%srfheatflux(i,j)+da2%srfheatflux(i,j)
#if defined(ALBSNW_EVOL)  
         da1%srfsnw     (i,j)=da1%srfsnw(i,j)     +da2%srfsnw(i,j)
#endif         
#endif
#if defined(EVP)
         da1%uice   (i,j) = da1%uice   (i,j) + da2%uice(i,j)
         da1%vice   (i,j) = da1%vice   (i,j) + da2%vice(i,j)
         da1%wice   (i,j) = da1%wice   (i,j) + da2%wice(i,j)
         da1%tauxice(i,j) = da1%tauxice(i,j) + da2%tauxice(i,j)
         da1%tauyice(i,j) = da1%tauyice(i,j) + da2%tauyice(i,j)
#if defined(ICE_DYN_DIAG)
         da1%sigp   (i,j) = da1%sigp   (i,j) + da2%sigp(i,j)
         da1%sigm   (i,j) = da1%sigm   (i,j) + da2%sigm(i,j)
         da1%sig12  (i,j) = da1%sig12  (i,j) + da2%sig12(i,j)
         da1%pr     (i,j) = da1%pr     (i,j) + da2%pr(i,j)
         da1%strI   (i,j) = da1%strI   (i,j) + da2%strI(i,j)
         da1%strII  (i,j) = da1%strII  (i,j) + da2%strII(i,j)
#endif
#if defined(TEST_ICE_AGE)
         da1%fy_frac (i,j) = da1%fy_frac (i,j) + da2%fy_frac (i,j)
         da1%fy_age  (i,j) = da1%fy_age  (i,j) + da2%fy_age  (i,j)
         da1%rdg_frac(i,j) = da1%rdg_frac(i,j) + da2%rdg_frac(i,j)
#endif
#if defined(PARAM_EST)
         da1%sstb (i,j) = da1%sstb (i,j) + da2%sstb (i,j)
         da1%msshb(i,j) = da1%msshb(i,j) + da2%msshb(i,j)
#endif
#endif
      end do
      end do
      end subroutine da_plus_da




c --- ----------------------------------------------
c --- Save averaged daily state. Uses .ab-file IO
c --- If rmem provided and rmem>0 add '_rmem' to the
c --- daily average file name.
c --- ----------------------------------------------
      subroutine daily_ave_save(rmem,rstep,rtime)
      use mod_za
      implicit none
c
      integer, optional, intent(in) :: rmem
      integer, optional, intent(in) :: rstep
      real, optional, intent(in) :: rtime
      real :: xmin,xmax,xmin2,xmax2,coord,loctime
      integer :: ifld,n,k,nop,ktr,locstep
      character(len=60) :: filenamea,filenameb, filebase
      character(len=8)  :: char8
      integer :: n2dfld,n3dfld
      logical :: written
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmpfld
      include 'common_blocks.h'
c
      if(present(rstep))then
        locstep = rstep
      else
c ---   Use global value
        locstep = nstep
      endif
      if(present(rtime))then
        loctime = rtime
      else
c ---   Use global value
        loctime = time
      endif
c
      n2dfld=23 
#if defined(TEST_ICE_AGE)
      n2dfld=26
#endif
#if defined(ALBSNW_EVOL)
      n2dfld=27
#endif
#if defined(ICE_DYN_DIAG)
      n2dfld=33
#endif
#if defined(PARAM_EST)
      n2dfld=35
#endif
      n3dfld=6+ntracr
#if defined (NOR05)
cnumber of fields totally written out ntracr number of tracers in bio 2 -- rates Ute
      if (ntracr.ne.0) then
        n3dfld=6+ntracr+2
      end if
#endif
#if defined (ECOSM)
cnumber of fields totally written out ntracr number of tracers in bio(ntracr)plus potential rates e.g. production 2 -- rates Ute
      if (ntracr.ne.0) then
        n3dfld=6+ntracr+ndbio
      end if
#endif
c --- Filenames
      if (present(rmem)) then
        if (rmem>0) then
         filebase=dailyfilebase(da%inityear, da%initday,
     &                         da%dumpyear, da%dumpday, rmem)
        else
         filebase=dailyfilebase(da%inityear, da%initday,
     &                         da%dumpyear, da%dumpday)
        end if
      else
        filebase=dailyfilebase(da%inityear, da%initday,
     &                         da%dumpyear, da%dumpday)
      end if
      filenamea=trim(filebase)//'.a'
      filenameb=trim(filebase)//'.b'
c
c --- Open files .. 
      nop=711
      if (mnproc==1) 
     &   open (unit=nop,file=trim(filenameb),
     &         status='replace')
      call zaiopf(trim(filenamea),'replace',nop)
c
c --- Write .b Header
      if (mnproc==1) 
     &   write(nop,116) ctitle,iversn,iexpt,yrflag,itdm,jtdm, 
     &      kdm,da%inityear,da%initday,da%dumpyear,da%dumpday,
     &      da%nrmem
c
c --- Write the local ensemble members into ".b" header
      if (mnproc==1) then
         print *,filenamea
         print *,filenameb
         print *,da%inityear,da%initday
         print *,da%dumpyear,da%dumpday
         do n=1,1000
             write(nop,118) n,da%lmem(n)
          end do
          write(nop,119)
       end if
c
c ---  Dump all 2D fields. If you add more fields 
c ---  you may want to augment n2dfld
       coord = 0
       ifld=1
       do ifld=1,n2dfld
c --     initialize write flag
         written = .false.
c
         if     (ifld==1) then
            tmpfld=da%ub
            char8='ubavg   '
            written=.true.
         elseif (ifld==2) then
            tmpfld=da%vb
            char8='vbavg   '
            written=.true.
         elseif (ifld==3) then
            tmpfld=da%ssh
            char8='ssh     '
            written=.true.
         elseif (ifld==4) then
            tmpfld=da%heatflux
            char8='surflx  '
            written=.true.
         elseif (ifld==5) then
            tmpfld=da%swflx
            char8='swflx   '
            written=.true.
         elseif (ifld==6) then
            tmpfld=da%salflx
            char8='salflx  '
            written=.true.
         elseif (ifld==7) then
            tmpfld=da%emnp
            char8='emnp    '
            written=.true.
         elseif (ifld==8) then
            tmpfld=da%taux
            char8='taux    '
            written=.true.
         elseif (ifld==9) then
            tmpfld=da%tauy
            char8='tauy    '
            written=.true.
#if defined(ICE) || defined(ICESTATE)
         elseif (ifld==10) then
            tmpfld=da%fice
            char8='fice    '
            written=.true.
         elseif (ifld==11) then
            tmpfld=da%hice
            char8='hice    '
            written=.true.
         elseif (ifld==12) then
            tmpfld=da%hsnw
            char8='hsnw    '
            written=.true.
         elseif (ifld==13) then
            tmpfld=da%srfalb
            char8='albedo  '
            written=.true.
         elseif (ifld==14) then
            tmpfld=da%srfheatflux
            char8='surflxs '
            written=.true.
         elseif (ifld==15) then
            tmpfld=da%htndncy
            char8='htndncy '
            written=.true.
#endif
#if defined (DAILY_ENSVAR)
         elseif (ifld==21) then
            tmpfld=da%ssh_sq  
            char8='ssh_sq  '
            written=.true.
         elseif (ifld==22) then
            tmpfld=da%fice_sq  
            char8='fice_sq '
            written=.true.
         elseif (ifld==23) then
            tmpfld=da%hice_sq  
            char8='hice_sq '
            written=.true.
#endif
#if defined(EVP)
         elseif (ifld==16) then
            tmpfld=da%uice
            char8='uice    '
            written=.true.
         elseif (ifld==17) then
            tmpfld=da%vice
            char8='vice    '
            written=.true.
         elseif (ifld==18) then
            tmpfld=da%wice
            char8='spdice  '
            written=.true.
         elseif (ifld==19) then
            tmpfld=da%tauxice
            char8='tauxi   '
            written=.true.
         elseif (ifld==20) then
            tmpfld=da%tauyice
            char8='tauyi   '
            written=.true.
#if defined(TEST_ICE_AGE)
         elseif (ifld==24) then
            tmpfld=da%fy_age
            char8='fy_age  '
            written=.true.
         elseif (ifld==25) then
            tmpfld=da%fy_frac
            char8='fy_frac '
            written=.true.
         elseif (ifld==26) then
            tmpfld=da%rdg_frac
            char8='rdg_frac'
#endif
#if defined(ALBSNW_EVOL)  
         elseif (ifld==27) then
            tmpfld=da%srfsnw
            char8='fsnw_cov'
            written=.true.
#endif             
#if defined(ICE_DYN_DIAG)
         elseif (ifld==28) then
            tmpfld=da%sigp
            char8='sigp    '
            written=.true.
         elseif (ifld==29) then
            tmpfld=da%sigm
            char8='sigm    '
            written=.true.
         elseif (ifld==30) then
            tmpfld=da%sig12
            char8='sig12   '
            written=.true.
         elseif (ifld==31) then
            tmpfld=da%pr
            char8='pr      '
            written=.true.
         elseif (ifld==32) then
            tmpfld=da%strI
            char8='strI    '
            written=.true.
         elseif (ifld==33) then
            tmpfld=da%strII
            char8='strII   '
            written=.true.
#endif
#endif
#if defined(PARAM_EST)
         elseif (ifld==34) then
            tmpfld=da%sstb
            char8='sstb    '
            written=.true.
         elseif (ifld==35) then
            tmpfld=da%msshb
            char8='msshb   '
           written=.true.
#endif
         end if
         if     (written) then
            call zaiowr(tmpfld(1-nbdy,1-nbdy),ip,.false.,
     &                        xmin,xmax, nop, .false.)
            if     (mnproc.eq.1) then
               write (nop,117) char8,locstep,loctime,0,coord,xmin,xmax
               call flush(nop)
            end if
         endif 
      end do
c
c ---  Dump all 3D fields. If you add more fields 
c ---  you may want to augment n3dfld
       do k=1,kdm
       coord=sigma(k)
       do ifld=1,n3dfld
         if (ifld==1) then
            call zaiowr(da%u(1-nbdy,1-nbdy,k),iu,.false.,
     &                        xmin,xmax, nop, .false.)
            char8='utot    '
         elseif (ifld==2) then
            call zaiowr(da%v(1-nbdy,1-nbdy,k),iv,.false.,
     &                        xmin,xmax, nop, .false.)
            char8='vtot    '
         elseif (ifld==3) then
            call zaiowr(da%s(1-nbdy,1-nbdy,k),ip,.false.,
     &                        xmin,xmax, nop, .false.)
            char8='saln    '
         elseif (ifld==4) then
            call zaiowr(da%t(1-nbdy,1-nbdy,k),ip,.false.,
     &                        xmin,xmax, nop, .false.)
            char8='temp    '
         elseif (ifld==5) then
            call zaiowr(da%p(1-nbdy,1-nbdy,k),ip,.false.,
     &                        xmin,xmax, nop, .false.)
            char8='pres    '
         elseif (ifld==6) then
            call zaiowr(da%w(1-nbdy,1-nbdy,k),ip,.false.,
     &                        xmin,xmax, nop, .false.)
            char8='kinetic '
#if defined (NOR05)
c   UTe set names according to your setup
         elseif (ifld>6 .and. ifld <= 6+ntracr) then
            ktr=ifld-6
            !KAL -- I snagged this from prtsol
            select case (ktr)
                case(init)  ; char8='NIT     '   
                case(ipho)  ; char8='PHO     '   
                case(isil)  ; char8='SIL     '   
                case(idet)  ; char8='DET     '   
                case(isis)  ; char8='SIS     '   
                case(ifla)  ; char8='FLA     '   
                case(idia)  ; char8='DIA     '   
                case(ioxy)  ; char8='OXY     '   
                case(ised)  ; char8='SED     '   
                case(iyel)  ; char8='YEL     '   
                case(idetp) ; char8='DETP    '   
                case(imicro); char8='MICRO   '   
                case(imeso) ; char8='MESO    '   
          end select
            call zaiowr(da%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,.false.,
     &                        xmin,xmax, nop, .false.)
         elseif (ifld == 7+ntracr) then
             char8='primprod'
             call zaiowr(da%grosspp_ave(1-nbdy,1-nbdy,k),ip,.false.,
     &                        xmin,xmax, nop, .false.)
         elseif (ifld == 8+ntracr) then
             char8='netpp'
             call zaiowr(da%netpp_ave(1-nbdy,1-nbdy,k),ip,.false.,
     &                        xmin,xmax, nop, .false.)
#endif /* NOR05 */
#if defined (ECOSM)
c   UTe set names according to your setup
         elseif (ifld>6 .and. ifld <= 6+ntracr) then
            ktr=ifld-6
            !KAL -- I snagged this from prtsol
            select case (ktr)
                case(ifla)  ; char8='FLA     '   
                case(idia)  ; char8='DIA     '   
                case(imicro); char8='MICRO   '   
                case(imeso) ; char8='MESO    '   
                case(idet)  ; char8='DET     '   
                case(inh4)  ; char8='NH4     '   
                case(idom)  ; char8='DOM     '   
                case(init)  ; char8='NIT     '   
                case(ipho)  ; char8='PHO     '   
                case(isil)  ; char8='SIL     '   
                case(ioxy)  ; char8='OXY     '   
                case(iopa)  ; char8='OPAL    '   
                case(ibg)   ; char8='CYAN    '   
                case(ifish) ; char8='FISH    '   
c                case(ised1) ; char8='SED1    '   
c                case(ised2) ; char8='SED_Pho '   
c                case(ised3) ; char8='SED_Sil '   
#if defined(ECO2)
                case(idic) ; char8='DIC    '   
                case(ialk) ; char8='ALK    '   
#endif
            end select
            call zaiowr(da%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,.false.,
     &                        xmin,xmax, nop, .false.)
         elseif (ifld>6+ntracr .and. ifld <= 6+ntracr+ndbio) then
            ktr=ifld-6-ntracr
            select case (ktr)
                case(iflap)  ; char8='flapp     '   
                case(idiap)  ; char8='diapp     '   
                case(izsp); char8='zspp   '   
                case(izlp) ; char8='zlpp    '   
                case(ispp)  ; char8='netsp     '   
                case(inpp)  ; char8='netpp     '   
#if defined(ECO2)
                case(ipco) ; char8='PCO2    '   
                case(iph) ; char8='PH    '   
                case(iatc) ; char8='CO2flx    '   
#endif
            end select

            call zaiowr(da%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,.false.,
     &                        xmin,xmax, nop, .false.)
#endif /* ECOSM */
         end if
c
         if     (mnproc.eq.1) then
            write (nop,117) char8,locstep,loctime,k,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
      end do
      end do
      call zaiocl(nop)
      if (mnproc==1) close(nop)
c
 116  format (a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & i5,4x,'''kdm   '' = Vertical     array size'/
     & i5,4x,'''syear '' = Year of integration start '/
     & i5,4x,'''sday  '' = Day of integration start'/
     & i5,4x,'''dyear '' = Year of this dump      '/
     & i5,4x,'''dday  '' = Day of this dump     '/
     & i5,4x,'''count '' = Ensemble counter       ')
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
 118  format ('member ',i5.5,' = ',i1,' Ensemble member flag')
 119  format('field         time step  model day',
     & '  k  dens        min              max')
      end subroutine daily_ave_save


      
c --- ----------------------------------------------
c --- Read a previous averaged daily state. Needed when 
c --- Running ensembles. (And may be used if restart
c --- capability is added)
c --- ----------------------------------------------
      subroutine daily_ave_read(datmp, inityear, initday,
     &                          dumpyear, dumpday,
     &                          rmem,rstep,rtime)
      use mod_za
      implicit none
      type(daily_averages), intent(inout) :: datmp
      integer, intent(in) :: inityear, initday,dumpyear,dumpday
      integer, optional, intent(in) :: rmem
      integer, optional, intent(out) :: rstep
      real, optional, intent(out) :: rtime
      real :: xmin,xmax,xmin2,xmax2,coord,loctime
      integer :: idummy, l_itdm,l_jtdm,l_kdm,k,nop, ios,n,
     &           locnstep, lociversn, lociexpt, locyrflag,
     &           ktr
      character(len=80) :: locctitle(4)
      character(len=60) :: filenamea,filenameb, filebase
      character(len=8) :: char8
      logical :: no2dfld,no3dfld
      include 'common_blocks.h'
c
c --- Set filenames
      filebase=dailyfilebase(inityear, initday,
     &                       dumpyear, dumpday, rmem)
      filenamea=trim(filebase)//'.a'
      filenameb=trim(filebase)//'.b'
      !print *,filenamea
      !print *,filenameb
      if (mnproc==1) then
         write(lp,'(a)') 'Reading '//
     &      trim(filebase)//'.[ab]'
         call flush(lp)
      end if
c
c --- Open files 
      nop=710
      open (unit=nop,file=trim(filenameb),
     &      status='old')
      call zaiopf(trim(filenamea),'old',nop)
c
c --- Read .b Header
      read(nop,116) locctitle,lociversn,lociexpt,locyrflag,
     &   l_itdm,l_jtdm, 
     &   l_kdm,datmp%inityear,datmp%initday,
     &   datmp%dumpyear,datmp%dumpday,datmp%nrmem
c
c --  Security Check
      if (l_itdm/=itdm.or.l_jtdm/=jtdm.or.l_kdm/=kdm) then
         Print *,'Grid dimensions are wrong!!!'
         call xcstop ('(mod_daily_average:daily_ave_read)')
         stop         '(mod_daily_average:daily_ave_read)'
      end if
c
c --- Read the local ensemble members into ".b" header
      do n=1,1000
          read(nop,118) idummy,datmp%lmem(n)
       end do
       read(nop,119)
c
c ---  Start reading fields Read until we fail
       ios=0
       do while (ios==0)
          read (nop,117,iostat=ios) char8,locnstep,loctime,
     &                              k,coord,xmin2,xmax2
          if (present(rstep)) then
            rstep = locnstep
          end if
          if (present(rtime)) then
            rtime = loctime
          end if
          no2dfld=.false.
          no3dfld=.false.
          ktr=0
          if (ios==0) then
c ---     Read 2D fields
          if (k==0) then
             if      (adjustl(trim(char8))=='ubavg') then
               call zaiord(datmp%ub(1-nbdy,1-nbdy),iu,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='vbavg') then
               call zaiord(datmp%vb(1-nbdy,1-nbdy),iv,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='ssh'  ) then
               call zaiord(datmp%ssh(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
#if defined (DAILY_ENSVAR)
             else if (adjustl(trim(char8))=='ssh_sq'  ) then
               call zaiord(datmp%ssh_sq(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='fice_sq'  ) then
               call zaiord(datmp%fice_sq(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='hice_sq'  ) then
               call zaiord(datmp%hice_sq(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
#endif
             else if (adjustl(trim(char8))=='surflx') then
               call zaiord(datmp%heatflux(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='swflx') then
               call zaiord(datmp%swflx(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='salflx') then
               call zaiord(datmp%salflx(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='emnp') then
               call zaiord(datmp%emnp(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='taux') then
               call zaiord(datmp%taux(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='tauy') then
               call zaiord(datmp%tauy(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
#if defined (ICE) || defined(ICESTATE)
             else if (adjustl(trim(char8))=='fice') then
               call zaiord(datmp%fice(1-nbdy,1-nbdy),ip,.false.,
     &                            xmin,xmax, nop)
             else if (adjustl(trim(char8))=='hice') then
               call zaiord(datmp%hice(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='hsnw') then
               call zaiord(datmp%hsnw(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='surflxs') then
               call zaiord(datmp%srfheatflux(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='albedo') then
               call zaiord(datmp%srfalb(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='htndncy') then
               call zaiord(datmp%htndncy(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
#if defined(ALBSNW_EVOL)
             else if (adjustl(trim(char8))=='fsnw_cov') then
               call zaiord(datmp%srfsnw(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
#endif               
#endif
#if defined(EVP)
             else if (adjustl(trim(char8))=='uice') then
               call zaiord(datmp%uice(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='vice') then
               call zaiord(datmp%vice(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='spdice') then
               call zaiord(datmp%wice(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='tauxi  ') then
               call zaiord(datmp%tauxice(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='tauyi  ') then
               call zaiord(datmp%tauyice(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
#if defined(ICE_DYN_DIAG)
             else if (adjustl(trim(char8))=='sigp   ') then
               call zaiord(datmp%sigp(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='sigm   ') then
               call zaiord(datmp%sigm(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='sig12  ') then
               call zaiord(datmp%sig12(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='pr     ') then
               call zaiord(datmp%pr(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='strI   ') then
               call zaiord(datmp%strI(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='strII  ') then
               call zaiord(datmp%strII(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
#endif
#if defined(TEST_ICE_AGE)
             else if (adjustl(trim(char8))=='fy_frac') then
               call zaiord(datmp%fy_frac(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='fy_age ') then
               call zaiord(datmp%fy_age(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='rdg_frac') then
               call zaiord(datmp%rdg_frac(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)

#endif
#if defined(PARAM_EST)
             else if (adjustl(trim(char8))=='sstb') then
               call zaiord(datmp%sstb(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='msshb ') then
               call zaiord(datmp%msshb(1-nbdy,1-nbdy),ip,.false.,
     &                           xmin,xmax, nop)
#endif
#endif
             else
                no2dfld=.true.
             end if
c ---     Read 3D fields
          else if (k<=kdm .and. k>0) then

             if      (adjustl(trim(char8))=='utot') then
               call zaiord(datmp%u(1-nbdy,1-nbdy,k),iu,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='vtot') then
               call zaiord(datmp%v(1-nbdy,1-nbdy,k),iv,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='saln'  ) then
               call zaiord(datmp%s(1-nbdy,1-nbdy,k),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='temp') then
               call zaiord(datmp%t(1-nbdy,1-nbdy,k),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='pres') then
               call zaiord(datmp%p(1-nbdy,1-nbdy,k),ip,.false.,
     &                           xmin,xmax, nop)
             else if (adjustl(trim(char8))=='kinetic') then
               call zaiord(datmp%w(1-nbdy,1-nbdy,k),ip,.false.,
     &                           xmin,xmax, nop)
#if defined (NOR05)
            elseif (adjustl(trim(char8))=='NIT') then
               ktr=init
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='PHO') then
               ktr=ipho
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='SIL') then
               ktr=isil
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='DET') then
               ktr=idet
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='SIS') then
               ktr=isis
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='FLA') then
               ktr=ifla
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='DIA') then
               ktr=idia
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='OXY') then
               ktr=ioxy
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='SED') then
               ktr=ised
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='YEL') then
               ktr=iyel
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='DETP') then
               ktr=idetp
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='MICRO') then
               ktr=imicro
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='MESO') then
               ktr=imeso
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='primprod') then
               call zaiord(datmp%grosspp_ave(1-nbdy,1-nbdy,k),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='netpp') then
               call zaiord(datmp%netpp_ave(1-nbdy,1-nbdy,k),ip,
     &                      .false.,xmin,xmax, nop)
#endif /* NOR05 */
#if defined (ECOSM)
            elseif (adjustl(trim(char8))=='FLA') then
               ktr=ifla
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='DIA') then
               ktr=idia
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='MICRO') then
               ktr=imicro
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='MESO') then
               ktr=imeso
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='DET') then
               ktr=idet
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='NH4') then
               ktr=inh4
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='DOM') then
               ktr=idom
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='NIT') then
               ktr=init
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='PHO') then
               ktr=ipho
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='SIL') then
               ktr=isil
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='OXY') then
               ktr=ioxy
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='OPAL') then
               ktr=iopa
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='CYAN') then
               ktr=ibg
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='FISH') then
               ktr=ifish
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
#if defined(ECO2)
            elseif (adjustl(trim(char8))=='DIC') then
               ktr=idic
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='ALK') then
               ktr=ialk
               call zaiord(datmp%trcr_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='PCO2') then
               ktr=ipco
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='PH') then
               ktr=iph
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='CO2flx') then
               ktr=iatc
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
#endif
            elseif (adjustl(trim(char8))=='diapp') then
               ktr=idiap
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='flapp') then
               ktr=iflap
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='netpp') then
               ktr=inpp
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='netzs') then
               ktr=izsp
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='netzl') then
               ktr=izlp
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
            elseif (adjustl(trim(char8))=='netsp') then
               ktr=ispp
               call zaiord(datmp%diagn_ave(1-nbdy,1-nbdy,k,ktr),ip,
     &                      .false.,xmin,xmax, nop)
#endif /*ECOSM  */
             else
                no3dfld=.true.
            end if
c ---    We just read a field with k> kdm or k<0. No good
         else
            if (mnproc==1) then
               write(lp,'(a)') 'Illegal specification of k'
               write(lp,'(a,i3)') 'Offending field '//char8,k
            end if
            call xcstop('(mod_daily_average:daily_ave_read)')
            stop        '(mod_daily_average:daily_ave_read)'
         end if
c
c ---    Sanity check:
         if (no2dfld.and.no3dfld) then
            if (mnproc==1) then
               write(lp,'(a)') 'No flds corresponding to ',
     &            adjustl(trim(char8)), ' found'
            end if
            call xcstop('(mod_daily_average:daily_ave_read)')
            stop        '(mod_daily_average:daily_ave_read)'
         end if
c
c ---    Check for coherency between min and maxes from a and b files
         if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.
     &       abs(xmax-xmax2)>abs(xmax)*1e-4) then
            if     (mnproc.eq.1) then
               write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           'error - .a and .b files not consistent:',
     &           '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,
     &           '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
               print *,'Variable :',char8,k
            endif
            call xcstop('(mod_daily_average:daily_ave_read)')
            stop        '(mod_daily_average:daily_ave_read)'
         end if
         end if  !ios==0
       enddo
c
c ---  Finished - close files
       call zaiocl(nop)
       if (mnproc==1) close(nop)
c
 116  format (a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & i5,4x,'''kdm   '' = Vertical     array size'/
     & i5,4x,'''syear '' = Year of integration start '/
     & i5,4x,'''sday  '' = Day of integration start'/
     & i5,4x,'''dyear '' = Year of this dump      '/
     & i5,4x,'''dday  '' = Day of this dump     '/
     & i5,4x,'''count '' = Ensemble counter       ')
C    & 'field       time step  model day',
C    & '  k  dens        min              max')
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
 118  format ('member ',i5.5,' = ',i1,' Ensemble member flag')
 119  format('field         time step  model day',
     & '  k  dens        min              max')
      end subroutine daily_ave_read

c --- ----------------------------------------------
c --- Combine daily average files for a number of
c --- ensemble members.
c --- ----------------------------------------------
      subroutine daily_ave_combine(byear,bjuld,fyear,fjuld,enssize)
      implicit none
      integer, intent(in) :: byear,bjuld,fyear,fjuld,enssize
      integer :: locmem, locstep
      real :: loctime
      type(daily_averages) :: datmp

c --- Initialize the daily average structure
      call daily_ave_ini()
c --- Read the average for the first ensemble member
      call daily_ave_read(da,byear,bjuld,fyear,fjuld,1,locstep,loctime)
c --- Iterate through and add the other ensemble members
      do locmem=2,enssize
        call daily_ave_read(datmp,byear,bjuld,fyear,fjuld,locmem)
        call da_plus_da(da,datmp)
      end do
c --- Save the combined average to a new file
      call daily_ave_save(0,locstep,loctime)
      end subroutine daily_ave_combine


c --- ----------------------------------------------
c --- Produce a file name for the daily average, based 
c --- on start of run and the day it is valid for
c --- ----------------------------------------------
      function dailyfilebase(inityear, initday,dumpyear, dumpday, rm)
      use mod_hycom_nersc
      implicit none
      integer, intent(in ) :: inityear, initday,
     &                        dumpyear, dumpday
      integer, optional, intent(in ) :: rm
      character(len=60) :: dailyfilebase

      if (present(rm)) then
        if (rm>0) then
           write(dailyfilebase,121) rungen, inityear, initday,
     &                           dumpyear, dumpday, rm
        else
           write(dailyfilebase,120) rungen, inityear, initday,
     &                           dumpyear, dumpday
        end if
      else
        write(dailyfilebase,120) rungen, inityear, initday,
     &                           dumpyear, dumpday
      end if
 120  format (a3,'DAILY_',i4.4,'_',i3.3,'_',i4.4,'_',i3.3)
 121  format (a3,'DAILY_',i4.4,'_',i3.3,'_',i4.4,'_',i3.3,'_mem',i3.3)
      end function dailyfilebase
#endif
      end module mod_daily_average
