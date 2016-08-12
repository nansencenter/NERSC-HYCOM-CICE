      module mod_average
c --- --------------------------------------------------------
c --- Module for processing weekly averages. 
c ---
c --- Public routines:
c ---   ave_process      - Processes the averaging logic each 
c ---                      time step. Sets when to accumulate/read/save
c ---   read_restart_ave - Read weekly average restart into week_ave
c ---   save_restart_ave - Save weekly average restart from week_ave
c ---
c --- Private routines:
c ---    ave_ini         - Initializes week_ave and counters
c ---    aveweek_setup   - Sets up week info
c ---    average_dump    - Saves a weeks worth of averages
c ---    save_ave2       - Saves ave_week state
c ---                      (Used by save_restart_ave, average_dump)
c ---    read_ave2       - Reads ave_week state
c ---                      (Used by read_restart_ave)
c ---    fileaveres      - Returns file name used by restart files
c --- --------------------------------------------------------
c --- Flag (laverage) and averaging interval (hours) are public 
c --- and set in m_limits.
c ---
c --- --------------------------------------------------------
c --- Personnel: Geir Evensen (GE) and Knut Liseter (KAL)
c ---
c GE  Some time in the 1990s : Initial version
c KAL 2003(?)    : Modified for MPI runs (uses tiled arrays)
c KAL 2008.09.28 : Modified for HYCOM 2.2. Tighter wrapping 
c ...              of module. One-line calls from hycom. Moved
c ...              restart logic in here from mod_restart
c KAL 2008.09.28 : Will work without restart files (incomplete averages)
c ...              Previously this stopped hycom.
c --- --------------------------------------------------------
      use mod_xc
#if defined(ICESTATE)
      use mod_icestate
      use mod_icestate_fluxes
#endif
      implicit none
c --- Shield everything not explicitly declared public
      private
c
      logical, save, public :: laverage
      integer, save, public :: average_dt
      integer, save :: average_frq
#ifndef WEEKLY_AVERAGE
      logical,parameter, public :: l_weekly_average_in_code=.false.
#else
      ! type defenitions and variables for  computation of averages
      logical,parameter, public :: l_weekly_average_in_code=.true.



      ! This state is used for intermediate calculations.
      ! Introduced to do the work on tiled regions in place
      ! of the global region (MPI).
      type average_state_tile
         character(len=3) cmm
         character(len=4) cyy
         integer counter
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &      ssh, relsalflx, relsurflx,salflx,surflx, dpmix,ub,vb
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     &      u,v,d,t,s,w,difs
         real, 
     &   dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,mxtrcr) ::
     &   tracer
#if defined(EVP) 
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &      uice,vice,wice
#if defined(ICE_DYN_DIAG)
     &      ,sigp,sigm,sig12,pr,strI,strII
#endif
#endif
#if defined(TEST_ICE_AGE)
     &      ,fy_age,fy_frac,rdg_frac
#endif
#if defined (PARAM_EST)
     &  ,msshb,sstb
#endif      
#if defined(ICE)||defined(ICESTATE)
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &      hicem,ficem,hsnwm
#if defined(ALB_SNW)
     &    salb,fscov
#endif  
#if defined(SSNOWD_ICE)
     &      ,hprcp,hmelt
#endif         
#endif
#if defined(ICESTATE)
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nthick) ::
     &      hice,fice,hsnw,rsnw,ice_tsrf,ice_alb,ice_qstr, 
     &      ice_swfl, ice_swtr, ice_trb, ice_lw, ice_brfl, 
     &      ice_mlfl,ice_ctop,ice_cbot,ice_nbot,ice_ntop,
     &      ice_grw, ice_lgrw
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)        ::
     &      lead_tot, lead_sw, lead_lw, lead_trb
#if defined (SSNOWD)
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nthick) ::
     &      hprcp,hmelt
#endif  
#if defined (ICEAGE)
         real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nthick) ::
     &      age
#endif         
#endif    
      end type average_state_tile


      ! Temp state
      type average_test
         logical res
         logical ini
         logical add
         logical sav
      end type average_test

c --- Type containing weekly average. Shared among routines
      type(average_state_tile), save :: ave_week

c --- Type denoting when to add, save, read restart etc
      type(average_test),       save ::  lave

c --- Total number of records in an average. set by aveweek_setup
      integer, save :: totit
      integer, save :: estit

c --- File name of average file. Set by aveweek_setup
      character(len=16) fileavew 

      public :: ave_process, read_restart_ave, save_restart_ave
                

      contains

C --- -----------------------------------------------------------
C --- Moved some averaging logic from hycom in here
C --- -----------------------------------------------------------
      subroutine ave_process(rt,m,n)
      use mod_year_info, only: year_info
      implicit none
      type(year_info),  intent(in) :: rt
      integer,          intent(in) :: m,n
      include 'common_blocks.h'


c --- For calculating/dumping to average 
      if (laverage) then
         call aveweek_setup(rt,baclin)
         if (lave%ini) then
            call ave_ini
            if (mnproc==1) write(lp,'(a)') 
     &         'AVERAGE: setting all ave fields to zero'
         endif
         if (lave%add) then
            if (mnproc==1) write(lp,'(a,3i5,a)')  'AVERAGE: ',
     &         ave_week%counter,totit,estit,'---Accumulating averages'
            call ave_add(m,n)
         endif
         if (lave%sav) then
            if (ave_week%counter==totit) then
               if (mnproc==1) write(lp,'(a,3i5,a)')  'AVERAGE: ',
     &            ave_week%counter,totit,estit,'---Dumping average'
            else
               if (mnproc==1) write(lp,'(a,3i5,a)')  'AVERAGE: ',
     &            ave_week%counter,totit,estit,
     &            '---Dumping partial average'
            end if
            call average_dump
         end if
      endif
      end subroutine

C --- -----------------------------------------------------------
C --- Restart logic moved from hycom and mod_restart in here
C --- -----------------------------------------------------------
      subroutine read_restart_ave(rt,baclin)
      use mod_year_info, only: year_info
      implicit none
      type(year_info),  intent(in) :: rt
      real,             intent(in) :: baclin
      character(len=80) :: tmpchar
      integer :: find
      logical :: exa, exb
      if (laverage) then
         call ave_ini
         call aveweek_setup(rt,baclin)
         if (lave%res) then
            tmpchar=fileaveres(rt)

            inquire(exist=exa,file=trim(tmpchar)//'.a')
            inquire(exist=exb,file=trim(tmpchar)//'.b')
            if (.not.exa.and..not.exb) then
               if (mnproc==1) write(lp,'(a)') 'no weekly restart --'//
     &            'continuing with ave=0 (incomplete average)'
               return
            end if
c
            if (mnproc==1) write(lp,'(a)') 'reading:'//trim(tmpchar)//
     &                                      '.[ab]'
            call read_ave2(trim(tmpchar))
         end if
      endif
      end subroutine


C --- -----------------------------------------------------------
C --- Restart logic moved from hycom and mod_restart in here
C --- -----------------------------------------------------------
      subroutine save_restart_ave(rt)
      use mod_year_info, only: year_info
      implicit none
      type(year_info),  intent(in) :: rt
      character(len=80) :: tmpchar
      integer :: find
      if (laverage) then
         if (mnproc==1) print *,'Dumping ave fields'
         tmpchar=fileaveres(rt)
         if (mnproc==1) 
     &      write(lp,'(a)',advance='no') 'save_restart_ave: saving:'//
     &      trim(tmpchar)//'.[ab]'
         call save_ave2(trim(tmpchar))
         if (mnproc==1)
     &    write(lp,'(a,i5,a)') '   counter= ',ave_week%counter, 
     &      '   AVE done'
      endif
      end subroutine






C --- -----------------------------------------------------------
C --- Ave_add routine adds current state to the average week_ave.
c --- week_ave is shared in this module
C --- -----------------------------------------------------------
      subroutine ave_add(m,n)
      use mod_forcing_nersc
C     use m_sshcalc
#if defined(EVP)
      use mod_evp , only :  evp_u=>uvel , evp_v=>vvel
#if defined(ICE_DYN_DIAG)
      use mod_common_ice , only : stressp,stressm,stress12,pice
     &  ,strainI,strainII
#endif
#endif
#if defined(ICE)||defined(ICESTATE)
      use mod_common_ice , only : ficem, hicem,hsnwm
#if defined(ALBSNW_EVOL)
     &  ,surf_albedo_sum,surf_fscov
#endif
#if defined(TEST_ICE_AGE)
     &  ,fy_age,fy_frac,rdg_frac
#endif
#if defined (PARAM_EST)
     &  ,msshb,sstb
#endif      
#if defined (SSNOWD_ICE)
     &      ,hprcp,hmelt
#endif      
#endif

      implicit none
      integer, intent(in) :: m,n
c
      integer i,j,k,l,it,ip1,jp1
      integer ktr
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   tmpfld1,tmpfld2
      real th,tmpu,tmpv,tmpd
      include 'common_blocks.h'
      include 'stmt_fns.h'

      ! Tile before averaging below
      call xctilr(u(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
      call xctilr(v(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
      call xctilr(ubavg(  1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_uv)
      margin=1

      do k=1,kk
!$OMP PARALLEL DO PRIVATE(j,l,i,ip1,jp1,tmpu,tmpv,tmpd)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         jp1=j+1 ! ok for margin < nbdy
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
         ip1=i+1 ! ok for margin < nbdy


            tmpu = 0.5*(u  (i,j,k,n) + u  (ip1,j,k,n))
     &           + 0.5*(ubavg(i,j,n) + ubavg(ip1,j,n))
            tmpv = 0.5*(v  (i,j,k,n) + v  (i,jp1,k,n))
     &           + 0.5*(vbavg(i,j,n) + vbavg(i,jp1,n))
            tmpd =  dp(i,j,k,n)/onem

            ave_week%u(i,j,k)=ave_week%u(i,j,k) + tmpu*tmpd
            ave_week%v(i,j,k)=ave_week%v(i,j,k) + tmpv*tmpd
            ave_week%t(i,j,k)=ave_week%t(i,j,k) + temp(i,j,k,n)*tmpd
            ave_week%s(i,j,k)=ave_week%s(i,j,k) + saln(i,j,k,n)*tmpd
            ave_week%d(i,j,k)=ave_week%d(i,j,k) + tmpd
            ave_week%w(i,j,k)=ave_week%w(i,j,k) 
     &                       +(tmpu**2 + tmpv**2)*tmpd
            ave_week%difs(i,j,k)=ave_week%difs(i,j,k) + difs(i,j,k) ! NB - has kk+1 interfaces
            if (ntracr.ne.0) then
              do ktr=1,ntracr
                ave_week%tracer(i,j,k,ktr)=ave_week%tracer(i,j,k,ktr) + 
     &                                    tracer(i,j,k,n,ktr)*tmpd
              enddo
           end if

         enddo
         enddo
         enddo
!$OMP END PARALLEL DO 
      enddo


!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)

         ave_week%ssh  (i,j) =ave_week%ssh  (i,j) + srfhgt(i,j)/g
         ave_week%dpmix(i,j) =ave_week%dpmix(i,j) + dpmixl(i,j,n)/onem

         ! KAL -- get relaxation flux
         ave_week%relsalflx(i,j) =ave_week%relsalflx(i,j) + sssflx(i,j)
         ave_week%relsurflx(i,j) =ave_week%relsurflx(i,j) + sstflx(i,j)
         ave_week%salflx(i,j) =ave_week%salflx(i,j) + salflx(i,j)
         ave_week%surflx(i,j) =ave_week%surflx(i,j) + surflx(i,j)
         ave_week%vb    (i,j) = ave_week%vb(i,j)   + vbavg(i,j,n)  
         ave_week%ub    (i,j) = ave_week%ub(i,j)   + ubavg(i,j,n)  

#if defined(EVP)
         ave_week%uice(i,j)=ave_week%uice(i,j)+evp_u(i,j)
         ave_week%vice(i,j)=ave_week%vice(i,j)+evp_v(i,j)
         ave_week%wice(i,j)=ave_week%wice(i,j)+sqrt(evp_u(i,j)**2 
     &                     + evp_v(i,j)**2)
#if defined(ICE_DYN_DIAG)
         ave_week%sigp(i,j)= ave_week%sigp(i,j)+stressp(i,j)
         ave_week%sigm(i,j)= ave_week%sigm(i,j)+stressm(i,j)
         ave_week%sig12(i,j)= ave_week%sig12(i,j)+stress12(i,j)
         ave_week%pr(i,j)= ave_week%pr(i,j)+pice(i,j)
         ave_week%strI(i,j)= ave_week%strI(i,j)+strainI(i,j)
         ave_week%strII(i,j)= ave_week%strII(i,j)+strainII(i,j)
#endif
#if defined(TEST_ICE_AGE)
         ave_week%fy_age(i,j) =ave_week%fy_age(i,j) +fy_age(i,j)
         ave_week%fy_frac(i,j)=ave_week%fy_frac(i,j)+fy_frac(i,j)
         ave_week%rdg_frac(i,j)=ave_week%rdg_frac(i,j)+rdg_frac(i,j)
#endif
#if defined (PARAM_EST)
         ave_week%sstb(i,j) =ave_week%sstb(i,j) +sstb(i,j)
         ave_week%msshb(i,j)=ave_week%msshb(i,j)+msshb(i,j)
#endif      
#endif
#if defined(ICE)||defined(ICESTATE)
         ave_week%hicem(i,j)=ave_week%hicem(i,j)+hicem(i,j)
         ave_week%hsnwm(i,j)=ave_week%hsnwm(i,j)+hsnwm(i,j)
         ave_week%ficem(i,j)=ave_week%ficem(i,j)+ficem(i,j)
#if defined(ALBSNW_EVOL)      
         ave_week%salb(i,j)=ave_week%salb(i,j)+surf_albedo_sum(i,j)
         ave_week%fscov(i,j)=ave_week%fscov(i,j)+surf_fscov(i,j)
#endif
#if defined(SSNOWD_ICE)
         ave_week%hmelt(i,j)=ave_week%hmelt(i,j)+hmelt(i,j)
         ave_week%hprcp(i,j)=ave_week%hprcp(i,j)+hprcp(i,j)
#endif         
#endif
#if defined(ICESTATE)
         !Icestate 2d vars
         ave_week%lead_tot(i,j)=ave_week%lead_tot(i,j)+
     &      icestate_lead_tot(i,j)
         ave_week%lead_sw (i,j)=ave_week%lead_sw (i,j)+
     &      icestate_lead_sw (i,j)
         ave_week%lead_lw (i,j)=ave_week%lead_lw (i,j)+
     &      icestate_lead_lw (i,j)
         ave_week%lead_trb(i,j)=ave_week%lead_trb(i,j)+
     &      icestate_lead_trb(i,j)
#endif
      end do
      end do
      end do


#if defined(ICESTATE)
      do k=1,nthick
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
         ave_week%hice(i,j,k)=ave_week%hice(i,j,k)
     &                       +icestate(i,j)%ice(k)%hice
         ave_week%fice(i,j,k)=ave_week%fice(i,j,k)
     &                       +icestate(i,j)%ice(k)%fice
         ave_week%hsnw(i,j,k)=ave_week%hsnw(i,j,k)
     &                       +icestate(i,j)%ice(k)%hsnw
         ave_week%rsnw(i,j,k)=ave_week%rsnw(i,j,k)
     &                       +icestate(i,j)%ice(k)%rhosnw
         ave_week%ice_tsrf(i,j,k)=ave_week%ice_tsrf(i,j,k)
     &                           +icestate(i,j)%ice(k)%tsrf
         ave_week%ice_alb (i,j,k)=ave_week%ice_alb (i,j,k)
     &                           +icestate(i,j)%ice(k)%albs
         ave_week%ice_qstr(i,j,k)=ave_week%ice_qstr(i,j,k)
     &                           +icestate(i,j)%ice(k)%qstore 
         ave_week%ice_swfl(i,j,k)=ave_week%ice_swfl(i,j,k)
     &                           +icestate_swfl(i,j,k)
         ave_week%ice_swtr(i,j,k)=ave_week%ice_swtr(i,j,k)
     &                           +icestate_swtr(i,j,k)
         ave_week%ice_trb (i,j,k)=ave_week%ice_trb (i,j,k)
     &                           +icestate_trb (i,j,k)
         ave_week%ice_lw  (i,j,k)=ave_week%ice_lw  (i,j,k)
     &                           +icestate_lw  (i,j,k)
         ave_week%ice_brfl(i,j,k)=ave_week%ice_brfl(i,j,k)
     &                           +icestate_brfl(i,j,k)
         ave_week%ice_mlfl(i,j,k)=ave_week%ice_mlfl(i,j,k)
     &                           +icestate_mlfl(i,j,k)
         ave_week%ice_ctop(i,j,k)=ave_week%ice_ctop(i,j,k)
     &                           +icestate_ctop(i,j,k)
         ave_week%ice_cbot(i,j,k)=ave_week%ice_cbot(i,j,k)
     &                           +icestate_cbot(i,j,k)
         ave_week%ice_ntop(i,j,k)=ave_week%ice_ntop(i,j,k)
     &                           +icestate_ntop(i,j,k)
         ave_week%ice_nbot(i,j,k)=ave_week%ice_nbot(i,j,k)
     &                           +icestate_nbot(i,j,k)
         ave_week%ice_grw (i,j,k)=ave_week%ice_grw (i,j,k)
     &                           +icestate_grw (i,j,k)
         ave_week%ice_lgrw(i,j,k)=ave_week%ice_lgrw(i,j,k)
     &                           +icestate_lgrw(i,j,k)
#if defined (SSNOWD)
         ave_week%hprcp(i,j,k)=ave_week%hprcp(i,j,k)
     &                            +icestate(i,j)%ice(k)%hprcp
         ave_week%hmelt(i,j,k)=ave_week%hmelt(i,j,k)
     &                            +icestate(i,j)%ice(k)%hmelt
#endif 
#if defined (ICEAGE)
         ave_week%age(i,j,k)=ave_week%age(i,j,k)
     &                            +icestate(i,j)%ice(k)%age
#endif          
      end do
      end do
      end do
      end do
#endif
      !print *,'mod_average ssh :',minval(ssh),maxval(ssh)
      !print *,'mod_average ave ssh
      !:',minval(ave_week%ssh),maxval(ave_week%ssh)
      !print *,'mod_average tslp:',minval(tslp),maxval(tslp)
      ! Increment counter 
      ave_week%counter=ave_week%counter+1
      end subroutine ave_add



C --- -----------------------------------------------------------
C --- Ave_ini routine initializes the state week_ave. week_ave
c --- is shared in this module
C --- -----------------------------------------------------------
      subroutine ave_ini
      implicit none
      include 'common_blocks.h'
      integer :: i,j,k, ktr
      average_frq=24/average_dt
      do k=1,kk
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
         do i=1-margin,ii+margin
            ave_week%u(i,j,k)=0.0
            ave_week%v(i,j,k)=0.0
            ave_week%d(i,j,k)=0.0
            ave_week%t(i,j,k)=0.0
            ave_week%s(i,j,k)=0.0
            ave_week%difs(i,j,k)=0.0
            ave_week%w(i,j,k)=0.0
            if (ntracr.ne.0) then
               do ktr=1,ntracr
                  ave_week%tracer(i,j,k,ktr)=0.0
               end do
            end if
         end do
         end do
!$OMP END PARALLEL DO 
      end do
      ave_week%counter=0
      ave_week%cyy='0'
      ave_week%cmm='0'
      ave_week%ssh=0.
      ave_week%ub=0.
      ave_week%vb=0.
      ave_week%dpmix=0.
      ave_week%relsalflx=0.
      ave_week%relsurflx=0.
      ave_week%salflx=0.
      ave_week%surflx=0.
#if defined(EVP) 
      ave_week%uice=0.
      ave_week%vice=0.
      ave_week%wice=0.
#if defined(ICE_DYN_DIAG)
      ave_week%sigp =0.
      ave_week%sigm =0.
      ave_week%sig12=0.
      ave_week%pr   =0.
      ave_week%strI =0.
      ave_week%strII=0.
#endif
#endif
#if defined (TEST_ICE_AGE)
      ave_week%fy_age=0.
      ave_week%fy_frac=0.
      ave_week%rdg_frac=0.
#endif
#if defined (PARAM_EST)
      ave_week%sstb=0.
      ave_week%msshb=0.
#endif
#if defined(ICE)||defined(ICESTATE)
      ave_week%ficem=0.
      ave_week%hicem=0.
      ave_week%hsnwm=0.
#if defined(ALBSNW_EVOL) 
      ave_week%salb=0.
      ave_week%fscov=0.
#endif 
#if defined (SSNOWD_ICE)
      ave_week%hprcp=0.
      ave_week%hmelt=0.
#endif       
#endif
#if defined(ICESTATE)
      !Icestate 2d vars
      ave_week%lead_tot=0.
      ave_week%lead_sw =0.
      ave_week%lead_lw =0.
      ave_week%lead_trb=0.

      ave_week%hice=0.
      ave_week%fice=0.
      ave_week%hsnw=0.
      ave_week%rsnw=0.
      ave_week%ice_tsrf=0.
      ave_week%ice_alb =0.
      ave_week%ice_qstr=0.
      ave_week%ice_swfl=0.
      ave_week%ice_swtr=0.
      ave_week%ice_trb =0.
      ave_week%ice_lw  =0.
      ave_week%ice_brfl=0.
      ave_week%ice_mlfl=0.
      ave_week%ice_ctop=0.
      ave_week%ice_cbot=0.
      ave_week%ice_ntop=0.
      ave_week%ice_nbot=0.
      ave_week%ice_grw =0.
      ave_week%ice_lgrw=0.
#if defined (SSNOWD)
      ave_week%hprcp=0.
      ave_week%hmelt=0.
#endif      
#if defined (ICEAGE)
      ave_week%age   =0.
#endif 
#endif
      end subroutine ave_ini




C --- -----------------------------------------------------------
C --  Aveweek_setup routine Sets up for averaging
C --- -----------------------------------------------------------
      subroutine aveweek_setup(rt,dt)
      use mod_hycom_nersc
      use mod_year_info, only: year_info
      implicit none
      type(year_info), intent(in) :: rt
      real, intent(in) :: dt
c
      integer i
      integer it
      integer aveit
      integer current_day, current_week
      character(len=1) tag1
      logical, save :: lfirst=.true.
      integer  week(4)

      aveit=ave_week%counter

c --- Checking if time step is an integer multiple of the average_dt
      if (mod(float(average_dt*3600),dt) /= 0) then
         if (mnproc==1) then
            print *,'aveweek_setup: time step is not an '//
     &              'integer multiple of average_dt'
            print *,'aveweek_setup:         dt= ',dt
            print *,'aveweek_setup: average_dt= ',average_dt*3600
         end if
         call xcstop ('(aveweek_setup)')
         stop '(aveweek_setup)'
      endif
c
c --- Checking if average_dt integer multiple of 24
      if (mod(24,average_dt) /= 0) then
         if (mnproc==1) then
            print *,'aveweek_setup: average_dt should be a factor in 24'
            print *,'aveweek_setup: average_dt=1, 2, 3, 4, 6, 8, 12, 24'
            print *,'aveweek_setup: average_dt<= 4 is recomended to '//
     &              'resolve inertial oscillations'
         end if
         call xcstop ('(aveweek_setup)')
         stop '(aveweek_setup)'
      endif
c
c --- Setting number of days in each week.   Gives always four weeks in each 
c --- month independent of number of days in month.
      select case (rt%totdim(rt%imm))
      case (28)
         week(1)=7; week(2)=7; week(3)=7; week(4)=7
      case (29)
         week(1)=7; week(2)=7; week(3)=7; week(4)=8
      case (30)
         week(1)=7; week(2)=8; week(3)=7; week(4)=8
      case (31) 
         week(1)=7; week(2)=8; week(3)=8; week(4)=8
      end select
c
c --- Returns the number of the current week in this month (1-4)
      if (rt%idm+1 <= sum(week(1:1)) ) then
            current_week=1
      elseif (rt%idm+1 <= sum(week(1:2)) ) then
            current_week=2
      elseif (rt%idm+1 <= sum(week(1:3)) ) then
            current_week=3
      else
            current_week=4
      endif
c
c --- Returns current day in this week (1-7)
      current_day=rt%idm-sum(week(1:current_week-1))
      current_day=current_day+1
c
c --- Resets all variables to zero before the accumulation of fields starts.
      if ((current_day == 1).and.(rt%ihh == average_dt).and.
     &    (rt%iss < nint(dt))) then
         lave%ini=.true.
      else
         lave%ini=.false.
      endif
c
c --- If restart after beginning of week a restart file for the 
c --- average data must be read.
      lave%res=.false.
      if (lfirst) then
      if (.not.((current_day == 1).and.(rt%ihh == 0).and.
     &    (rt%iss < nint(dt)))) then
            lave%res=.true.
            lave%ini=.false.
      endif
      endif
c
c --- Set filename and total counter
      if (lave%ini.or.lave%res) then
         write(tag1,'(i1.1)')current_week
         fileavew(1:3)=rungen(1:3)
         fileavew(4:7)='AVE_'
         fileavew(8:11)=rt%cyy(1:4)
         fileavew(12:12)='_'
         fileavew(13:14)=rt%cmm(1:2)
         fileavew(15:15)='_'
         fileavew(16:16)=tag1

c ---    The total number of fields to be accumulated this week
         totit=week(current_week)*average_frq
      endif
c
c --- estit should be the number of the average counter. On a week
c --- switch it will be -1. use this to trigger save of an imcomplete 
c --- average
      estit=(current_day-1)*average_frq+rt%ihh/average_dt-1
c
c --- Accumulates fields every average_dt and dumps results after
c --- final field has been accumulated for the current week.
      if ( (.not.lfirst).and.(rt%iss < nint(dt)).and.
     &     (mod(rt%ihh,average_dt)==0) ) then
         lave%add=.true.
         if (aveit == totit-1) then
            lave%sav=.true.
         else
            lave%sav=.false.
         endif
c ---    TODO: Fix averaging logic also for incomplete averages
c ---          (when missing restart)
         if (estit == -1) then
            lave%sav=.true.
         else
            lave%sav=.false.
         endif
      else
         lave%add=.false.
         lave%sav=.false.
      endif
      lfirst=.false.
      end subroutine aveweek_setup



C --- -----------------------------------------------------------
C --- Average_dump Routine dumps average to a file
C --- -----------------------------------------------------------
      subroutine average_dump()
      implicit none
      integer :: find
      character(len=80) :: tmpchar
      tmpchar=fileavew
      call save_ave2(trim(fileavew))
      end subroutine

C --- -----------------------------------------------------------
C --- -----------------------------------------------------------
C --- save_ave2 -- Save daily average and collect it from nodes
C --- -----------------------------------------------------------
C --- -----------------------------------------------------------
      ! Gather daily average and save it to file
      ! NB -- file saving mechanism is now similar
      ! to the one used in standard hycom...
      subroutine save_ave2(filebase)
         use mod_za
         implicit none
         character(len=*), intent(in) :: filebase
c
         real :: coord,xmin,xmax
         integer :: k,iyy,imm,nop,ifld,hk
         character(len=8) :: char8
         character(len=2) :: ccat

         include 'common_blocks.h'

         ! Integer representation of week month and counter
         !print *,ave_week%cmm
         !print *,ave_week%cyy
         read(ave_week%cmm,'(i3)') imm
         read(ave_week%cyy,'(i4)') iyy


         ! Open files .. 
         nop=710
         if (mnproc==1) 
     &      open (unit=nop,file=filebase//'.b',
     &            status='replace')
         call zaiopf(filebase//'.a','replace',nop)
         
         ! Header
         if (mnproc==1) 
     &      write(nop,116) ctitle,iversn,iexpt,yrflag,itdm,jtdm, 
     &         kdm,imm,iyy,ave_week%counter

! --     Surface fields
         coord=0.0
         call zaiowr(ave_week%ssh(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'ssh     ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%dpmix(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'dpmix   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%relsalflx(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'rsalflx ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%relsurflx(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'rsurflx ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%salflx(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'salflx  ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%surflx(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'surflx  ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
         call zaiowr(ave_week%ub(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'ubavg   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
         call zaiowr(ave_week%vb(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'vbavg   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
#if defined(EVP) 
         call zaiowr(ave_week%uice(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'uice    ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%vice(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'vice    ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%wice(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'spdice  ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
#if defined(ICE_DYN_DIAG)
         call zaiowr(ave_week%sigp(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'sigp    ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%sigm(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'sigm    ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%sig12(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'sig12   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%pr(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'pr      ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%strI(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'strI    ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%strII(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'strII    ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
#endif
#endif
#if defined(TEST_ICE_AGE) 
         call zaiowr(ave_week%fy_age(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'fy_age  ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%fy_frac(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'fy_frac ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%rdg_frac(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'rdg_frac',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
#endif
#if defined (PARAM_EST)
         call zaiowr(ave_week%sstb(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'sstb    ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
         call zaiowr(ave_week%msshb(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'msshb   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
#endif

#if defined(ICE)||defined(ICESTATE)
         call zaiowr(ave_week%hicem(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'hicem   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%ficem(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'ficem   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%hsnwm(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'hsnwm   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
#if defined(ALBSNW_EVOL)
         call zaiowr(ave_week%salb(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'surf_alb',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif

         call zaiowr(ave_week%fscov(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'frac_snw',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif
#endif 
#if defined (SSNOWD_ICE)
         call zaiowr(ave_week%hprcp(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'hprcp   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif

         call zaiowr(ave_week%hmelt(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'hmelt   ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif
#endif         
#endif
#if defined(ICESTATE)
         call zaiowr(ave_week%lead_tot(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'lead_tot',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%lead_sw(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'lead_sw ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%lead_lw(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'lead_lw ',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile

         call zaiowr(ave_week%lead_trb(1-nbdy,1-nbdy),ip,.false.,
     &       xmin,xmax, nop, .false.)
         if     (mnproc.eq.1) then
            write (nop,117) 'lead_trb',nstep,time,0,coord,xmin,xmax
            call flush(nop)
         endif !1st tile
#endif




! --     3D fields
         do k=1,kdm
         do ifld=1,7+ntracr
            coord = sigma(k)
            ! U - velocity
            if (ifld==1) then
               call zaiowr(ave_week%u(1-nbdy,1-nbdy,k),iu,.false.,
     &               xmin,xmax, nop, .false.)
               char8='utot    '
            elseif (ifld==2) then
               call zaiowr(ave_week%v(1-nbdy,1-nbdy,k),iv,.false.,
     &               xmin,xmax, nop, .false.)
               char8='vtot    '
            elseif (ifld==3) then
               call zaiowr(ave_week%d(1-nbdy,1-nbdy,k),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='pres    '
            elseif (ifld==4) then
               call zaiowr(ave_week%t(1-nbdy,1-nbdy,k),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='temp    '
            elseif (ifld==5) then
               call zaiowr(ave_week%s(1-nbdy,1-nbdy,k),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='saln    '
            elseif (ifld==6) then
               call zaiowr(ave_week%w(1-nbdy,1-nbdy,k),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='kinetic '
            elseif (ifld==7) then
               call zaiowr(ave_week%difs(1-nbdy,1-nbdy,k),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='diffs   '
            elseif ( ifld==8) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,1),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr1   '
            elseif ( ifld==9) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,2),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr2   '
            elseif ( ifld==10) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,3),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr3   '
            elseif ( ifld==11) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,4),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr4   '
            elseif ( ifld==12) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,5),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr5   '
            elseif ( ifld==13) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,6),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr6   '
            elseif ( ifld==14) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,7),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr7   '
            elseif ( ifld==15) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,8),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr8   '
            elseif ( ifld==16) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,9),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr9   '
            elseif ( ifld==17) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,10),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr10  '
            elseif ( ifld==18) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,11),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr11  '
            elseif ( ifld==19) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,12),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr12  '
            elseif ( ifld==20) then
               call zaiowr(ave_week%tracer(1-nbdy,1-nbdy,k,13),ip,
     &              .false.,xmin,xmax, nop, .false.)
               char8='trcr13  '
            end if
            if     (mnproc.eq.1) then
               write (nop,117) char8,nstep,time,k,coord,xmin,xmax
               call flush(nop)
            endif !1st tile
         end do
         end do

#if defined(ICESTATE)
         do hk=1,nthick
#if defined(SSNOWD) && defined(ICEAGE)
         do ifld=1,22
#elif defined(SSNOWD) 
         do ifld=1,21
#elif defined(ICEAGE)
          do ifld=1,20
#else
         do ifld=1,19
#endif          
            coord = hk
            write(ccat,'(i2.2)') hk

            if (   ifld==1) then
               call zaiowr(ave_week%hice(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_hice'//ccat
            elseif (ifld==2) then
               call zaiowr(ave_week%fice(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_fice'//ccat
            elseif (ifld==3) then
               call zaiowr(ave_week%hsnw(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_hsnw'//ccat
            elseif (ifld==4) then
               call zaiowr(ave_week%rsnw(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_rsnw'//ccat
            elseif (ifld==5) then
               call zaiowr(ave_week%ice_tsrf(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_tsrf'//ccat
            elseif (ifld==6) then
               call zaiowr(ave_week%ice_alb (1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_albs '//ccat
            elseif (ifld==7) then
               call zaiowr(ave_week%ice_qstr(1-nbdy,1-nbdy,hk),
     &            ip,.false., xmin,xmax, nop, .false.)
               char8='i_qbrn'//ccat
            elseif (ifld==8) then
               call zaiowr(ave_week%ice_swfl(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_swfl'//ccat
            elseif (ifld==9) then
               call zaiowr(ave_week%ice_swtr(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_swtr'//ccat
            elseif (ifld==10) then
               call zaiowr(ave_week%ice_trb (1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_turb'//ccat
            elseif (ifld==11) then
               call zaiowr(ave_week%ice_lw  (1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_lw'//ccat
            elseif (ifld==12) then
               call zaiowr(ave_week%ice_brfl(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_brfl'//ccat
            elseif (ifld==13) then
               call zaiowr(ave_week%ice_mlfl(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_mlfl'//ccat
            elseif (ifld==14) then
               call zaiowr(ave_week%ice_ctop(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_ctop'//ccat
            elseif (ifld==15) then
               call zaiowr(ave_week%ice_cbot(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_cbot'//ccat
            elseif (ifld==16) then
               call zaiowr(ave_week%ice_nbot(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_nbot'//ccat
            elseif (ifld==17) then
               call zaiowr(ave_week%ice_ntop(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_ntop'//ccat
            elseif (ifld==18) then
               call zaiowr(ave_week%ice_grw (1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_grw'//ccat
            elseif (ifld==19) then
               call zaiowr(ave_week%ice_lgrw(1-nbdy,1-nbdy,hk),
     &            ip,.false.,xmin,xmax, nop, .false.)
               char8='i_lgrw'//ccat
#if defined(SSNOWD) && defined(ICEAGE)
            elseif (ifld==20) then
               call zaiowr(ave_week%hprcp(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_hpcp'//ccat
            elseif (ifld==21) then
               call zaiowr(ave_week%hmelt(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_hmlt'//ccat
            elseif (ifld==22) then
               call zaiowr(ave_week%age(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_iage'//ccat
#elif defined(SSNOWD) 
            elseif (ifld==20) then
               call zaiowr(ave_week%hprcp(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_hpcp'//ccat
            elseif (ifld==21) then
               call zaiowr(ave_week%hmelt(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_hmlt'//ccat               
#elif defined(ICEAGE)
            elseif (ifld==20) then
               call zaiowr(ave_week%age(1-nbdy,1-nbdy,hk),ip,.false.,
     &               xmin,xmax, nop, .false.)
               char8='i_iage'//ccat
#endif               
            else
               if (mnproc ==1) write(lp,*) 'You shouldnt be here ...'
               call xcstop ('(save_ave2)')
               stop '(save_ave2)'
            end if
            if     (mnproc.eq.1) then
c ---          write (nop,117) char8,nstep,time,hk,coord,xmin,xmax
               write (nop,117) char8,nstep,time,0,coord,xmin,xmax
               call flush(nop)
            endif !1st tile
         end do
         end do
#endif


         ! Close unit 
         close(unit=nop)
         call zaiocl(nop)

 116  format (a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & i5,4x,'''kdm   '' = Vertical     array size'/
     & i5,4x,'''month '' = Month of this dump     '/
     & i5,4x,'''year  '' = Year of this dump      '/
     & i5,4x,'''count '' = Averaging counter      '/
     & 'field       time step  model day',
     & '  k  dens        min              max')
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)

      end subroutine


C --- -----------------------------------------------------------
C --- -----------------------------------------------------------
C --- read_ave2 -- Read weekly average and distribute it to nodes
C --- Result is put in ave_week, private to this module
C --- -----------------------------------------------------------
C --- -----------------------------------------------------------
      subroutine read_ave2(filebase)
         use mod_za
         implicit none
         character(len=*), intent(in) :: filebase
c
         real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmpfld
         real :: xmin,xmax,xmin2,xmax2,coord
         integer :: iyy, imm, nop,dummyint,tmpcount
         integer :: l_kdm, l_itdm, l_jtdm,ifld,ios,klevel
         character(len=8) a8
         character(len=80) a80
         logical :: exa,exb

         include 'common_blocks.h'

         ! Open files .. 
         inquire(exist=exa,file=filebase//'.a')
         inquire(exist=exb,file=filebase//'.b')
         if (.not.(exa.and.exb)) then
            if (mnproc==1) 
     &      write(lp,*) 'Can not find '//trim(filebase)//'.b'
            call xcstop('(read_ave2)')
            stop '(read_ave2)'
         end if


         nop=99
         open (unit=nop,file=filebase//'.b',
     &         status='old')
         call zaiopf(filebase//'.a','old',nop)
         
         ! Header
!         read(nop,116) ctitle,iversn,iexpt,yrflag,l_itdm,l_jtdm, 
!     &      l_kdm,imm,iyy,tmpcount
         read(nop,'(a80)') a80;! print *,a80 !Header
         read(nop,'(a80)') a80;! print *,a80 !Header
         read(nop,'(a80)') a80;! print *,a80!Header
         read(nop,'(a80)') a80;! print *,a80!Header
         read(nop,'(a80)') a80;! print *,a80!iversn
         read(nop,'(a80)') a80;! print *,a80!iexpt
         read(nop,'(a80)') a80;! print *,a80!yrflag
         call aveini(l_itdm   ,'idm   ')
         call aveini(l_jtdm   ,'jdm   ')
         call aveini(l_kdm    ,'kdm   ')
         call aveini(imm      ,'month ')
         call aveini(iyy      ,'year  ')
         call aveini(tmpcount ,'count ')
         read(nop,*) a80!header

         ! Check dims....
         if (l_kdm/=kdm .or. l_itdm/=itdm .or. l_jtdm /= jtdm) then
            if (mnproc==1) then
               write(lp,'(a)') 'read_ave2:dimension mismatch....'
               print *,itdm,l_itdm
               print *,jtdm,l_jtdm
               print *,kdm ,l_kdm
            end if
            call xcstop('(read_ave2)')
            stop '(read_ave2)'
         end if

         ! Put into ave_week
         write(ave_week%cmm,'(i3.3)') imm
         write(ave_week%cyy,'(i4.4)') iyy
         ave_week%counter=tmpcount


         ! Start reading -- for each field check layer and var name
         ios=0
         do while (ios==0)
            read (nop,117,iostat=ios) a8,nstep,time,klevel,coord,
     &                                xmin2,xmax2
Cdiag       print *,a8,time,klevel

            if (ios/=0) exit

            call zaiord(tmpfld(1-nbdy,1-nbdy),ip,.false.,
     &          xmin,xmax, nop)

            ! Check xmin/xmax
            if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.
     &          abs(xmax-xmax2)>abs(xmax)*1e-4) then
               if     (mnproc.eq.1) then
               write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           'error - .a and .b files not consistent:',
     &           '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,
     &           '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
                  print *,'Variable :',a8
               endif
               call xcstop('(read_ave2)')
               stop '(read_ave2)'
            end if


            ! Put into various fields...
            if (trim(a8)=='ssh'.and.klevel==0) then
               ave_week%ssh=tmpfld
            elseif (trim(a8)=='dpmix'.and.klevel==0) then
               ave_week%dpmix=tmpfld
            elseif ( trim(a8) == 'rsalflx'.and. klevel==0) then
               ave_week%relsalflx=tmpfld
            elseif ( trim(a8) == 'salflx' .and. klevel==0) then
               ave_week%salflx=tmpfld
            elseif ( trim(a8) == 'rsurflx'.and. klevel==0) then
               ave_week%relsurflx=tmpfld
            elseif ( trim(a8) == 'surflx' .and. klevel==0) then
               ave_week%surflx=tmpfld
            elseif ( trim(a8) == 'ubavg' .and. klevel==0) then
               ave_week%ub=tmpfld
            elseif ( trim(a8) == 'vbavg' .and. klevel==0) then
               ave_week%vb=tmpfld
#if defined(EVP) 
            elseif ( trim(a8) == 'uice' .and. klevel==0) then
               ave_week%uice=tmpfld
            elseif ( trim(a8) == 'vice' .and. klevel==0) then
               ave_week%vice=tmpfld
            elseif ( trim(a8) == 'spdice' .and. klevel==0) then
               ave_week%wice=tmpfld
#if defined(ICE_DYN_DIAG)
            elseif ( trim(a8) == 'sigp' .and. klevel==0) then
               ave_week%sigp=tmpfld
            elseif ( trim(a8) == 'sigm' .and. klevel==0) then
               ave_week%sigm=tmpfld
            elseif ( trim(a8) == 'sig12' .and. klevel==0) then
               ave_week%sig12=tmpfld
            elseif ( trim(a8) == 'pr' .and. klevel==0) then
               ave_week%pr=tmpfld
            elseif ( trim(a8) == 'strI' .and. klevel==0) then
               ave_week%strI=tmpfld
            elseif ( trim(a8) == 'strII' .and. klevel==0) then
               ave_week%strII=tmpfld
#endif

#endif
#if defined(TEST_ICE_AGE) 
            elseif ( trim(a8) == 'fy_age' .and. klevel==0) then
               ave_week%fy_age=tmpfld
            elseif ( trim(a8) == 'fy_frac' .and. klevel==0) then
               ave_week%fy_frac=tmpfld
            elseif ( trim(a8) == 'rdg_frac' .and. klevel==0) then
               ave_week%rdg_frac=tmpfld
#endif
#if defined (PARAM_EST)
            elseif ( trim(a8) == 'sstb' .and. klevel==0) then
               ave_week%sstb=tmpfld
            elseif ( trim(a8) == 'msshb' .and. klevel==0) then
               ave_week%msshb=tmpfld
#endif
#if defined(ICE)||defined(ICESTATE)
            elseif ( trim(a8) == 'hicem' .and. klevel==0) then
               ave_week%hicem=tmpfld
            elseif ( trim(a8) == 'ficem' .and. klevel==0) then
               ave_week%ficem=tmpfld
            elseif ( trim(a8) == 'hsnwm' .and. klevel==0) then
               ave_week%hsnwm=tmpfld 
#if defined(ALBSNW_EVOL)
            elseif ( trim(a8) == 'salb' .and. klevel==0) then
               ave_week%salb=tmpfld
            elseif ( trim(a8) == 'fscov' .and. klevel==0) then
               ave_week%fscov=tmpfld    
#endif
#if defined (SSNOWD_ICE)
            elseif ( trim(a8) == 'hprcp' .and. klevel==0) then
               ave_week%hprcp=tmpfld
            elseif ( trim(a8) == 'hmelt' .and. klevel==0) then
               ave_week%hmelt=tmpfld    
#endif                
#endif
#if defined(ICESTATE)
            elseif ( trim(a8) == 'lead_tot' .and. klevel==0) then
               ave_week%lead_tot=tmpfld
            elseif ( trim(a8) == 'lead_sw'  .and. klevel==0) then
               ave_week%lead_sw =tmpfld
            elseif ( trim(a8) == 'lead_lw'  .and. klevel==0) then
               ave_week%lead_lw =tmpfld
            elseif ( trim(a8) == 'lead_trb' .and. klevel==0) then
               ave_week%lead_trb =tmpfld
            elseif ( trim(a8) == 'hice' .and. klevel/=0) then
               ave_week%hice (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'fice' .and. klevel/=0) then
               ave_week%fice (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'hsnw' .and. klevel/=0) then
               ave_week%hsnw (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'rsnw' .and. klevel/=0) then
               ave_week%rsnw (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'ice_tsrf' .and. klevel/=0) then
               ave_week%ice_tsrf (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'ice_alb' .and. klevel/=0) then
               ave_week%ice_alb  (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'ice_qstr' .and. klevel/=0) then
               ave_week%ice_qstr (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'ice_swfl' .and. klevel/=0) then
               ave_week%ice_swfl (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'ice_swtr' .and. klevel/=0) then
               ave_week%ice_swtr (:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'ice_trb' .and. klevel/=0) then
               ave_week%ice_trb  (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'ice_lw' .and. klevel/=0) then
               ave_week%ice_lw   (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'ice_brfl' .and. klevel/=0) then
               ave_week%ice_brfl (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'ice_mlfl' .and. klevel/=0) then
               ave_week%ice_mlfl (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'ice_ctop' .and. klevel/=0) then
               ave_week%ice_ctop (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'ice_cbot' .and. klevel/=0) then
               ave_week%ice_cbot (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'ice_nbot' .and. klevel/=0) then
               ave_week%ice_nbot (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'ice_ntop' .and. klevel/=0) then
               ave_week%ice_ntop (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'ice_grw ' .and. klevel/=0) then
               ave_week%ice_grw  (:,:,klevel)  =tmpfld
            elseif ( trim(a8) == 'ice_lgrw' .and. klevel/=0) then
               ave_week%ice_lgrw (:,:,klevel) =tmpfld
#if defined (SSNOWD)
            elseif ( trim(a8) == 'hprcp' .and. klevel/=0) then
               ave_week%hprcp (:,:,klevel) =tmpfld
            elseif ( trim(a8) == 'hmelt' .and. klevel/=0) then
               ave_week%hmelt (:,:,klevel) =tmpfld      
#endif   
#if defined (ICEAGE)
            elseif ( trim(a8) == 'age' .and. klevel/=0) then
               ave_week%age (:,:,klevel) =tmpfld
#endif               
#endif
            elseif ( trim(a8) == 'utot' .and. klevel/=0) then
               ave_week%u(:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'vtot' .and. klevel/=0) then
               ave_week%v(:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'saln' .and. klevel/=0) then
               ave_week%s(:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'diffs' .and. klevel/=0) then
               ave_week%difs(:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'temp' .and. klevel/=0) then
               ave_week%t(:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'pres' .and. klevel/=0) then
               ave_week%d(:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'kinetic' .and. klevel/=0) then
               ave_week%w(:,:,klevel)=tmpfld
            elseif ( trim(a8) == 'trcr1' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,1)=tmpfld
            elseif ( trim(a8) == 'trcr2' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,2)=tmpfld
            elseif ( trim(a8) == 'trcr3' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,3)=tmpfld
            elseif ( trim(a8) == 'trcr4' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,4)=tmpfld
            elseif ( trim(a8) == 'trcr5' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,5)=tmpfld
            elseif ( trim(a8) == 'trcr6' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,6)=tmpfld
            elseif ( trim(a8) == 'trcr7' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,7)=tmpfld
            elseif ( trim(a8) == 'trcr8' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,8)=tmpfld
            elseif ( trim(a8) == 'trcr9' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,9)=tmpfld
            elseif ( trim(a8) == 'trcr10' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,10)=tmpfld
            elseif ( trim(a8) == 'trcr11' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,11)=tmpfld 
            elseif ( trim(a8) == 'trcr12' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,12)=tmpfld
            elseif ( trim(a8) == 'trcr13' .and. klevel/=0) then
               ave_week%tracer(:,:,klevel,13)=tmpfld
           else
               if (mnproc==1) 
     &            write(lp,*) 'read_ave2: Unknown field '//trim(a8)
            end if

         end do

         ! Close unit 
         close(unit=nop)
         call zaiocl(nop)




 116  format (a80/a80/a80/a80/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & 'field       time step  model day',
     & '  k  dens        min              max')
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
      end subroutine


      function fileaveres(rt)
         use mod_year_info, only: year_info
         use mod_hycom_nersc
         implicit none
         character(len=24) fileaveres
         type(year_info) rt
         fileaveres(1:3)=rungen(1:3)
         fileaveres(4:10)='restart'
         fileaveres(11:14)=rt%cyy(1:4)
         fileaveres(15:15)='_'
         fileaveres(16:18)=rt%cdd(1:3)
         fileaveres(19:19)='_'
         fileaveres(20:21)=rt%chh(1:2)
         fileaveres(22:24)='AVE'
      end function fileaveres


      subroutine aveini(ivar,cvar)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,ivar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in aveini - input ',cvarin,
     +                       ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(aveini)')
               stop
      endif
      return
 6000 format(a6,' =',i10)
      end subroutine aveini

#endif
      end module mod_average
