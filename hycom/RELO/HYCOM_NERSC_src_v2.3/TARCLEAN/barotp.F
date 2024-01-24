#if defined(ROW_LAND)
#define SEA_P .true.
#define SEA_U .true.
#define SEA_V .true.
#elif defined(ROW_ALLSEA)
#define SEA_P allip(j).or.ip(i,j).ne.0
#define SEA_U alliu(j).or.iu(i,j).ne.0
#define SEA_V alliv(j).or.iv(i,j).ne.0
#else
#define SEA_P ip(i,j).ne.0
#define SEA_U iu(i,j).ne.0
#define SEA_V iv(i,j).ne.0
#endif

#if defined(ICE_INFLUENCE_BAROTROPIC)
      module mod_seaice_draft
      !TODO: Handle situautions where ice thickness exceeds draft...
      use mod_xc
      implicit none
      real,external :: pbudel_seaice, pbvdel_seaice
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::
     &   draft_u,draft_u_thk,draft_v,draft_v_thk
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &        draft_u,draft_v,draft_u_thk,draft_v_thk
#endif

      contains
      subroutine init_seaice_draft()
      implicit none
#if defined(RELO)
        allocate(
     &          draft_u    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          draft_v    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          draft_u_thk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          draft_v_thk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
        call mem_stat_add( 4*(idm+2*nbdy)*(jdm+2*nbdy) )
        draft_u=r_init
        draft_v=r_init
        draft_u_thk=r_init
        draft_v_thk=r_init
#endif
      end subroutine
C
      subroutine baclin_velocity_over_draft(n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
      integer, intent(in) :: n
      integer :: i,j,k
      real    :: usur1,dpu1,vsur1,dpv1, iceth,draft,psur1

      do j=1,jj
      do i=1,ii
         draft_u    (i,j)=0.0
         draft_u_thk(i,j)=0.0
         draft_v    (i,j)=0.0
         draft_v_thk(i,j)=0.0
C        
         if (SEA_U) then
            psur1=0.0
            usur1=0.0
            iceth=.5*(thkice(i,j)+thkice(i-1,j))
            draft=0.9*iceth*onem
            do k=1,kk
              dpu1 = min( dpu(i,j,k,n), max( 0.0, draft-psur1))
              usur1 = usur1 + dpu1*u(i,j,k,n)
              psur1 = psur1 + dpu1
              if     (psur1.ge.draft) then
                exit
              endif
            end do
            draft_u_thk(i,j)=psur1/onem
            draft_u    (i,j)=usur1/max(onemm,psur1)
CKAL  if (mnproc==1) print '(a,2i5,2f12.3)',"draft u:",i,j,
CKAL &   draft_u_thk(i,j),draft_u(i,j)
         end if
C
         if (SEA_V) then
            psur1=0.0
            vsur1=0.0
            iceth=.5*(thkice(i,j)+thkice(i,j-1))
            draft=0.9*iceth*onem
            do k=1,kk
              dpv1 = min( dpv(i,j,k,n), max( 0.0, draft-psur1))
              vsur1 = vsur1 + dpv1*v(i,j,k,n)
              psur1 = psur1 + dpv1
              if     (psur1.ge.draft) then
                exit
              endif
            end do
            draft_v_thk(i,j)=psur1/onem
            draft_v    (i,j)=vsur1/max(onemm,psur1)
CKAL  if (mnproc==1) print '(a,2i5,2f12.3)',"draft v:",i,j,
CKAL &   draft_v_thk(i,j),draft_v(i,j)
         end if
C
      end do
      end do
      call xctilr( draft_u    ,1,1, nbdy,nbdy, halo_uv)
      call xctilr( draft_u_thk,1,1, nbdy,nbdy, halo_uv)
      call xctilr( draft_v    ,1,1, nbdy,nbdy, halo_vv)
      call xctilr( draft_v_thk,1,1, nbdy,nbdy, halo_vv)
      end subroutine
c---- Include effect of moving ice. Ice divergence which differs from
c --- ocean divergence averaged over ice draft will now contribute to
c --- the evolution of the barotropic pressure.
c
      real function pbudel_seaice(i,j,n)
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
      integer, intent(in) :: i,j,n
C 
C --- Sea ice drift in u-point
      pbudel_seaice = .5* (
     &  + (si_u(i+1,j)+si_u(i,j))*scuy(i+1,j)*draft_u_thk(i+1,j)
     &  - (si_u(i-1,j)+si_u(i,j))*scuy(i  ,j)*draft_u_thk(i  ,j)
     &     )
C
C --- Subtract Total velocity over draft
      pbudel_seaice  = pbudel_seaice   - (
     &  + (draft_u(i+1,j)+ubavg(i+1,j,n))*scuy(i+1,j)*draft_u_thk(i+1,j)
     &  - (draft_u(i  ,j)+ubavg(i  ,j,n))*scuy(i  ,j)*draft_u_thk(i  ,j)
     &   )
CKAL
CKAL  if (mnproc==1) print '(a,2i5,3i5)',"pbudel_si 0",i,j,
CKAL &   ip(i-1,j),ip(i,j),ip(i+1,j)
CKAL  if (mnproc==1) print '(a,2i5,3f12.3)',"pbudel_si 1",i,j,
CKAL &   si_u(i-1,j),si_u(i,j),si_u(i+1,j)
CKAL  if (mnproc==1) print '(a,2i5,2f12.3)',"pbudel_si 2",i,j,
CKAL &   scuy(i+1,j),scuy(i,j)
CKAL  if (mnproc==1) print '(a,2i5,2f12.3)',"pbudel_si 3",i,j,
CKAL &   draft_u_thk(i+1,j),draft_u_thk(i,j)
CKAL  if (mnproc==1) print '(a,2i5,f12.3)', "pbudel_si 4",i,j,
CKAL &   pbudel_seaice
      end function
c
c
      real function pbvdel_seaice(i,j,n)
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
      integer, intent(in) :: i,j,n
C 
C --- Sea ice drift in v-point
      pbvdel_seaice = .5* (
     &  + (si_v(i,j+1)+si_v(i,j))*scvx(i,j+1)*draft_v_thk(i,j+1)
     &  - (si_v(i,j-1)+si_v(i,j))*scvx(i  ,j)*draft_v_thk(i  ,j)
     &     )
C
C --- Subtract Total velocity over draft
      pbvdel_seaice  = pbvdel_seaice - (
     &  +(draft_v(i,j+1)+vbavg(i,j+1,n))*scvx(i,j+1)*draft_v_thk(i,j+1)
     &  -(draft_v(i,j  )+vbavg(i,j  ,n))*scvx(i,j  )*draft_v_thk(i  ,j)
     &  )
      end function
c 
      end module mod_seaice_draft
#endif /*defined(ICE_INFLUENCE_BAROTROPIC)*/

      subroutine barotp(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_pipe       ! HYCOM debugging interface
      use mod_tides      ! HYCOM tides
#if defined(STOKES)
      use mod_stokes     !    HYCOM Stokes Drift
#endif
#if defined(ICE_INFLUENCE_BAROTROPIC)
      use mod_seaice_draft
#endif
      implicit none
c
      integer m,n
c
c --- ------------------------------------------------------------------------
c --- advance barotropic equations.
c ---   on entry: -n- is time t-dt, -m- is time t
c ---   on exit:                    -m- is time t, -n- is time t+dt
c ---   time level 3 is only used internally (n and m are always 1 or 2).
c
c --- LeapFrog version based on:
c ---   Y. Morel, Baraille, R., Pichon A. (2008) "Time splitting and
c ---   linear stability of the slow part of the barotropic component", 
c ---   Ocean Modeling, 23, pp 73-81.
c --- ------------------------------------------------------------------------
c
      logical    lpipe_barotp
      parameter (lpipe_barotp=.false.)
      logical    ldebug_barotp
      parameter (ldebug_barotp=.false.)
c
      real    q,pbudel,pbvdel,utndcy,vtndcy,wblpf
      real    d11,d12,d21,d22,ubp,vbp
      real*8  sump
      integer i,j,l,lll,ml,nl,mn,lstep1,margin,mbdy
c	 & ,iffstep
      logical ldrag
c	  data iffstep/0/
c	  save iffstep
c
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::
     &        pbavo,ubavo,vbavo,displd
c
      if     (.not.allocated(pbavo)) then
        allocate(
     &          pbavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          ubavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          vbavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         displd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( 4*(idm+2*nbdy)*(jdm+2*nbdy) )
                pbavo = r_init
                ubavo = r_init
                vbavo = r_init
               displd = r_init
#if defined(ICE_INFLUENCE_BAROTROPIC)
        call init_seaice_draft
#endif
      endif
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &        pbavo,ubavo,vbavo,displd
#endif
c
      mbdy = 6
c

C$RMY c   >>>>>>>>>>>>>>>>wave induced pressure      
C$RMY #if defined(STOKES_3D)
C$RMY       if (stpres) then
C$RMY          margin = 0
C$RMY !$OMP PARALLEL DO PRIVATE(j,i)
C$RMY !$OMP&         SCHEDULE(STATIC,jblk)
C$RMY          do j = 1-margin, jj+margin
C$RMY          do i = 1-margin, ii+margin
C$RMY             if (SEA_U) then
C$RMY                utotn(i,j) = utotn(i,j)
C$RMY      .                 - (ws0*(sj(i,j,ls0)-sj(i-1,j,ls0)) +
C$RMY      .                    ws1*(sj(i,j,ls1)-sj(i-1,j,ls1)))*scuxi(i,j)
C$RMY             endif !iu
C$RMY             if (SEA_V) then
C$RMY                vtotn(i,j) = vtotn(i,j)
C$RMY      .                 - (ws0*(sj(i,j,ls0)-sj(i,j-1,ls0)) +
C$RMY      .                    ws1*(sj(i,j,ls1)-sj(i,j-1,ls1)))*scvyi(i,j)
C$RMY             endif !iv
C$RMY          enddo !i
C$RMY          enddo !j
C$RMY       endif
C$RMY #endif
C$RMY c
C$RMY c   >>>>>>>>>>>>>>>>wave induced radiation stress
C$RMY #if defined(STOKES_2D)
C$RMY       margin = 0
C$RMY !$OMP PARALLEL DO PRIVATE(j,i)
C$RMY !$OMP&         SCHEDULE(STATIC,jblk)
C$RMY       do j = 1-margin, jj+margin
C$RMY       do i = 1-margin, ii+margin
C$RMY          if (SEA_U) then
C$RMY             utotn(i,j) = utotn(i,j)
C$RMY      .              - ( ws0*d_sx(i,j,ls0) + ws1*d_sx(i,j,ls1) )
C$RMY      .               / max(onem,
C$RMY      .                 0.5*(pbot(i-1,j)*onetai(i-1,j)+pbavg(i-1,j,m) + 
C$RMY      .                      pbot(i,j)*onetai(i,j)    +pbavg(i,j,m)))
C$RMY          endif !iu
C$RMY          if (SEA_V) then
C$RMY             vtotn(i,j) = vtotn(i,j)
C$RMY      .              - ( ws0*d_sy(i,j,ls0) + ws1*d_sy(i,j,ls1) )
C$RMY      .               / max(onem,
C$RMY      .                 0.5*(pbot(i-1,j)*onetai(i,j-1)+pbavg(i,j-1,m) + 
C$RMY      .                      pbot(i,j)*onetai(i,j)    +pbavg(i,j,m)))
C$RMY          endif !iv
C$RMY       enddo !i
C$RMY       enddo !j
C$RMY #endif

c
c --- utotn,vtotn from momtum is time step t-1 to t+1 barotropic tendency
      call xctilr(utotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_uv)
      call xctilr(vtotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_vv)
#if defined(ICE_INFLUENCE_BAROTROPIC)
c
c --- pre-compute velocities over ice draft
      call baclin_velocity_over_draft(m)
#endif
c
      if     (lpipe .and. lpipe_barotp) then
c ---   compare two model runs.
        call pipe_compare_sym2(utotn, iu,'barotp:utotn',
     &                         vtotn, iv,'barotp:vtotn')
        call pipe_compare_sym1(pvtrop,iq,'barotp:pvtrp')
      endif
c
c --- explicit time integration of barotropic flow (forward-backward scheme)
c --- in order to combine forward-backward scheme with leapfrog treatment of
c --- coriolis term, v-eqn must be solved before u-eqn every other time step
c
      if     (btrlfr) then
        if     (delt1.ne.baclin) then  !not on very 1st time step
c ---     start at time level t-dt and go to t+dt.
          lstep1 = lstep + lstep  !more stable, but also more expensive
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              pbavo(i,j)   = pbavg(i,j,n)  !save t-1 for RA filter
              ubavo(i,j)   = ubavg(i,j,n)  !save t-1 for RA filter
              vbavo(i,j)   = vbavg(i,j,n)  !save t-1 for RA filter
c
              pbavg(i,j,3) = pbavg(i,j,n)
              ubavg(i,j,3) = ubavg(i,j,n)
              vbavg(i,j,3) = vbavg(i,j,n)
            enddo !i
          enddo !j
        else !1st time step
c ---     start at time level t and go to t+dt.
          lstep1 = lstep
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              pbavo(i,j)   = 0.0 !makes correct mean height safe
              pbavg(i,j,n) = pbavg(i,j,m)
              ubavg(i,j,n) = ubavg(i,j,m)
              vbavg(i,j,n) = vbavg(i,j,m)
              pbavg(i,j,3) = pbavg(i,j,m)
              ubavg(i,j,3) = ubavg(i,j,m)
              vbavg(i,j,3) = vbavg(i,j,m)
            enddo !i
          enddo !j
        endif !usual:1st time step
      else
c ---   start at time level t    and go to t+dt.
        lstep1 = lstep          !original, less stable, method
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            pbavo(i,j)   = 0.0 !makes correct mean height safe
            pbavg(i,j,n) = pbavg(i,j,m)
            ubavg(i,j,n) = ubavg(i,j,m)
            vbavg(i,j,n) = vbavg(i,j,m)
          enddo !i
        enddo !j
      endif !btrlfr
c
      ldrag = tidflg.gt.0 .and. drgscl.ne.0.0 .and. thkdrg.eq.0.0
c
      if     (ldrag) then
        displd(:,:) = 0.0
      endif
c
c --- time step loop
c
      if     (btrlfr) then
        wblpf = 0.0   !1st minor time step, lll=1, only
      else
        wblpf = wbaro
      endif
c
      do 840 lll=1,lstep1,2
c
      call xctilr(pbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_ps)
      call xctilr(ubavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_vv)
c
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,nl),ip,'barot+:pbavn')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,nl),iu,'barot+:ubavn',
     &    vbavg(1-nbdy,1-nbdy,nl),iv,'barot+:vbavn')
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,ml),ip,'barot+:pbavm')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,ml),iu,'barot+:ubavm',
     &    vbavg(1-nbdy,1-nbdy,ml),iv,'barot+:vbavm')
      endif
c
c --- odd minor time step.
c
      ml=n
      nl=3
c
c --- continuity equation, and tidal drag on p-grid
c
c --- rhs: pbavg, ubavg+, vbavg+
c --- lhs: pbavg
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,pbudel,pbvdel,
!$OMP&                    ubp,vbp,d11,d12,d21,d22,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
#if defined(STOKES)
c
c   Barotropic Stokes flow included here
c
            pbudel = (ubavg(i+1,j,ml)+usdbavg(i+1,j))*
     &                    (depthu(i+1,j)*scuy(i+1,j))
     &              -(ubavg(i,  j,ml)+usdbavg(i,  j))*
     &                    (depthu(i,  j)*scuy(i,  j))
            pbvdel = (vbavg(i,j+1,ml)+vsdbavg(i,j+1))*
     &                    (depthv(i,j+1)*scvx(i,j+1))
     &              -(vbavg(i,j,  ml)+vsdbavg(i,j  ))*
     &                    (depthv(i,j  )*scvx(i,j  ))
#else
            pbudel =  ubavg(i+1,j,ml)*(depthu(i+1,j)*scuy(i+1,j))
     &               -ubavg(i  ,j,ml)*(depthu(i  ,j)*scuy(i  ,j))
            pbvdel =  vbavg(i,j+1,ml)*(depthv(i,j+1)*scvx(i,j+1))
     &               -vbavg(i,j  ,ml)*(depthv(i,j  )*scvx(i,j  ))
#endif
#if defined(ICE_INFLUENCE_BAROTROPIC)
c --------- Include effect of moving ice.
            pbudel=pbudel+pbudel_seaice(i,j,ml)
            pbvdel=pbvdel+pbvdel_seaice(i,j,ml)
CKAL  if (mnproc==1) print '(i5,a3,i5,4f12.3)' ,i,"-1 ",j,
CKAL &   pbudel_seaice(i-1,j,ml),covice(i-1,j),thkice(i-1,j),si_u(i-1,j)
CKAL  if (mnproc==1) print '(i5,a3,i5,4f12.3)' ,i,"   ",j,
CKAL &   pbudel_seaice(i,j,ml),covice(i,j),thkice(i,j),si_u(i,j)
CKAL  if (mnproc==1) print '(i5,a3,i5,4f12.3)' ,i,"+1 ",j,
CKAL &   pbudel_seaice(i+1,j,ml),covice(i+1,j),thkice(i+1,j),si_u(i+1,j)
#endif
c
            pbavg(i,j,nl)=
     &        ((1.-wblpf)*pbavg(i,j,ml)+
     &             wblpf *pbavg(i,j,nl) )-
     &         (1.+wblpf)*dlt*(pbudel + pbvdel)*scp2i(i,j)
c
            if     (ldrag) then
c
c ---         tidal drag tensor on p-grid:
c ---           ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
c ---           vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
c ---         solve implicitly by inverting the matrix:
c ---            1+(dlt/H)*t.11    (dlt/H)*t.12
c ---              (dlt/H)*t.21  1+(dlt/H)*t.22
c ---         use depths (H) rather than onem*pbavg (h) for stability.
c
              ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
              vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
              d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
              d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
              d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
              d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
              q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
c ---         set util5,util6 to the ubavg,vbavg drag increment
              util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
              util6(i,j) = q*(ubp*d21+vbp*(1.0-d11)) - vbp
c ---         add an explicit antidrag correction
*             util5(i,j) = util5(i,j) - (d11*untide(i,j)+
*    &                                   d12*vntide(i,j) )
*             util6(i,j) = util6(i,j) - (d21*untide(i,j)+
*    &                                   d22*vntide(i,j) )
c ---         dissipation per m^2
              displd(i,j) = displd(i,j) +
     &                      (ubp*util5(i,j) + vbp*util6(i,j))*
     &                      depths(i,j)*qthref/dlt
c
*             if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*               write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
*    &            nstep,i+i0,j+j0,lll,
*    &            'ubp,new,vbp,new =',
*    &          ubp,ubp+util5(i,j),
*    &          vbp,vbp+util6(i,j)
*             endif !debug
            else
              util5(i,j) = 0.0
              util6(i,j) = 0.0
            endif !ldrag
          endif !ip
        enddo !i
      enddo !j
c
      mn=ml
c
c --- u momentum equation, 1st
c
c --- rhs: pbavg+, vbavg+, pvtrop+
c --- lhs: ubavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,utndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            utndcy=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)+
     &       ((vbavg(i  ,j,  mn)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,mn)*depthv(i  ,j+1))+
     &        (vbavg(i-1,j,  mn)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i,j+1)))
c
            ubavg(i,j,nl)=
     &        ((1.-wblpf)*ubavg(i,j,ml)+
     &             wblpf *ubavg(i,j,nl) )+
     &         (1.+wblpf)*dlt*(utndcy+utotn(i,j))+
     &                0.5*(util5(i,j)+util5(i-1,j))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,6g15.6)')
*    &          nstep,i+i0,j+j0,lll,
*    &          'u_old,u_new,p_grad,corio,u_star,drag =',
*    &          ubavg(i,j,ml),ubavg(i,j,nl),
*    &           -thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
*    &          (vbavg(i  ,j,  mn)*depthv(i  ,j)
*    &          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
*    &          +vbavg(i-1,j,  mn)*depthv(i-1,j)
*    &          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
*    &          *(pvtrop(i,j)+pvtrop(i,j+1))
*    &          *.125 * dlt,utotn(i,j) * dlt,
*    &          0.5*(util5(i,j)+util5(i-1,j))
*           endif !debug
          endif !iu
        enddo !i
      enddo !j
c
      mn = nl
c
c --- v momentum equation, 2nd
c --- rhs: pbavg+, ubavg+, pvtrop+
c --- lhs: vbavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,vtndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_V) then
            vtndcy=-thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)-
     &       ((ubavg(i,  j  ,mn)*depthu(i,  j  )
     &        +ubavg(i+1,j  ,mn)*depthu(i+1,j  ))+
     &        (ubavg(i,  j-1,mn)*depthu(i,  j-1)
     &        +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i+1,j)))
c
            vbavg(i,j,nl)=
     &        ((1.-wblpf)*vbavg(i,j,ml)+
     &             wblpf *vbavg(i,j,nl) )+
     &         (1.+wblpf)*dlt*(vtndcy+vtotn(i,j))+
     &                0.5*(util6(i,j)+util6(i,j-1))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,6g15.6)')
*    &          nstep,i+i0,j+j0,lll,
*    &          'v_old,v_new,p_grad,corio,v_star,drag =',
*    &          vbavg(i,j,ml),vbavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
*    &          -(ubavg(i,  j  ,mn)*depthu(i,j  )
*    &           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
*    &           +ubavg(i,  j-1,mn)*depthu(i,j-1)
*    &           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
*    &          *(pvtrop(i,j)+pvtrop(i+1,j))
*    &          *.125 * dlt, vtotn(i,j) * dlt,
*    &          0.5*(util6(i,j)+util6(i,j-1))
*           endif !debug
          endif !iv
        enddo !i
      enddo !j
c
      if     (ldebug_barotp) then
        call xcsync(flush_lp)
      endif
c
#if ! defined(RELO)
      if     (lbflag.eq.1) then
        call latbdp( nl)
      elseif (lbflag.eq.3) then
        call latbdf( nl,lll)
      endif
#endif
      if     (lbflag.eq.2) then
        call latbdt( nl,lll)
      elseif (lbflag.eq.4) then
        call latbdtf(nl,lll)
      endif
c
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,nl),ip,'barot+:pbavn')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,nl),iu,'barot+:ubavn',
     &    vbavg(1-nbdy,1-nbdy,nl),iv,'barot+:vbavn')
        call pipe_compare_sym1(
     &    pbavg(1-nbdy,1-nbdy,ml),ip,'barot+:pbavm')
        call pipe_compare_sym2(
     &    ubavg(1-nbdy,1-nbdy,ml),iu,'barot+:ubavm',
     &    vbavg(1-nbdy,1-nbdy,ml),iv,'barot+:vbavm')
      endif
c
c --- even minor time step.
c
      ml=3
      nl=n
      wblpf = wbaro  !used for all subsequent time steps: lll=2,lstep1
c
c --- continuity equation
c
c --- rhs: pbavg, ubavg+, vbavg+
c --- lhs: pbavg
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,pbudel,pbvdel,
!$OMP&                    ubp,vbp,d11,d12,d21,d22,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
#if defined(STOKES)
c
c           Barotropic Stokes flow included here
c
            pbudel = (ubavg(i+1,j,ml)+usdbavg(i+1,j))*
     &                    (depthu(i+1,j)*scuy(i+1,j))
     &              -(ubavg(i,  j,ml)+usdbavg(i,  j))*
     &                     (depthu(i ,j)*scuy(i,  j))
            pbvdel = (vbavg(i,j+1,ml)+vsdbavg(i,j+1))*
     &                    (depthv(i,j+1)*scvx(i,j+1))
     &              -(vbavg(i,j,  ml)+vsdbavg(i,j  ))*
     &                    (depthv(i,j  )*scvx(i,j  ))
#else
            pbudel =  ubavg(i+1,j,ml)*(depthu(i+1,j)*scuy(i+1,j))
     &               -ubavg(i  ,j,ml)*(depthu(i  ,j)*scuy(i  ,j))
            pbvdel =  vbavg(i,j+1,ml)*(depthv(i,j+1)*scvx(i,j+1))
     &               -vbavg(i,j  ,ml)*(depthv(i,j  )*scvx(i,j  ))

#endif
#if defined(ICE_INFLUENCE_BAROTROPIC)
c --------- Include effect of moving ice.
            pbudel=pbudel+pbudel_seaice(i,j,ml)
            pbvdel=pbvdel+pbvdel_seaice(i,j,ml)
#endif
            pbavg(i,j,nl)=
     &        ((1.-wblpf)*pbavg(i,j,ml)+
     &             wblpf *pbavg(i,j,nl) )-
     &         (1.+wblpf)*dlt*(pbudel + pbvdel)*scp2i(i,j)
c
            if     (ldrag) then
c ---         tidal drag tensor on p-grid:
c ---           ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
c ---           vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
c ---         solve implicitly by inverting the matrix:
c ---            1+(dlt/H)*t.11    (dlt/H)*t.12
c ---              (dlt/H)*t.21  1+(dlt/H)*t.22
c ---         use depths (H) rather than onem*pbavg (h) for stability.
c
              ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
              vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
              d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
              d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
              d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
              d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
              q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
c ---         set util5,util6 to the ubavg,vbavg drag increment
              util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
              util6(i,j) = q*(ubp*d21+vbp*(1.0-d11)) - vbp
c ---         add an explicit antidrag correction
*             util5(i,j) = util5(i,j) - (d11*untide(i,j)+
*    &                                   d12*vntide(i,j) )
*             util6(i,j) = util6(i,j) - (d21*untide(i,j)+
*    &                                   d22*vntide(i,j) )
c ---         dissipation per m^2
              displd(i,j) = displd(i,j) +
     &                      (ubp*util5(i,j) + vbp*util6(i,j))*
     &                      depths(i,j)*qthref/dlt
c
*             if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*               write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
*    &            nstep,i+i0,j+j0,lll+1,
*    &            'ubp,new,vbp,new =',
*    &          ubp,ubp+util5(i,j),
*    &          vbp,vbp+util6(i,j)
*             endif !debug
            else
              util5(i,j) = 0.0
              util6(i,j) = 0.0
            endif !ldrag
          endif !ip
        enddo !i
      enddo !j
c
      mn=ml
c
c --- v momentum equation, 1st
c
c --- rhs: pbavg+, ubavg+, pvtrop+
c --- lhs: vbavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,vtndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_V) then
            vtndcy=-thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)-
     &       ((ubavg(i,  j  ,mn)*depthu(i,  j  )
     &        +ubavg(i+1,j  ,mn)*depthu(i+1,j  ))+
     &        (ubavg(i,  j-1,mn)*depthu(i,  j-1)
     &        +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i+1,j)))
c
            vbavg(i,j,nl)=
     &        ((1.-wblpf)*vbavg(i,j,ml)+
     &             wblpf *vbavg(i,j,nl))+
     &         (1.+wblpf)*dlt*(vtndcy+vtotn(i,j))+
     &                0.5*(util6(i,j)+util6(i,j-1))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,6g15.6)')
*    &          nstep,i+i0,j+j0,lll+1,
*    &          'v_old,v_new,p_grad,corio,v_star,drag =',
*    &          vbavg(i,j,ml),vbavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
*    &          -(ubavg(i,  j  ,mn)*depthu(i,j  )
*    &           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
*    &           +ubavg(i,  j-1,mn)*depthu(i,j-1)
*    &           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
*    &          *(pvtrop(i,j)+pvtrop(i+1,j))
*    &          *.125 * dlt, vtotn(i,j) * dlt,
*    &          0.5*(util6(i,j)+util6(i,j-1))
*           endif !debug
          endif !iv
        enddo !i
      enddo !j
c
      mn=nl
c
c --- u momentum equation, 2nd
c
c --- rhs: pbavg+, vbavg+, pvtrop+
c --- lhs: ubavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,i,utndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            utndcy=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)+
     &       ((vbavg(i  ,j,  mn)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,mn)*depthv(i  ,j+1))+
     &        (vbavg(i-1,j,  mn)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))*
     &       (0.125*(pvtrop(i,j)+pvtrop(i,j+1)))
c
            ubavg(i,j,nl)=
     &        ((1.-wblpf)*ubavg(i,j,ml)+
     &             wblpf *ubavg(i,j,nl) )+
     &         (1.+wblpf)*dlt*(utndcy+utotn(i,j))+
     &                0.5*(util5(i,j)+util5(i-1,j))
c
*           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
*             write (lp,'(i9,2i5,i3,3x,a,6g15.6)')
*    &          nstep,i+i0,j+j0,lll+1,
*    &          'u_old,u_new,p_grad,corio,u_star,drag =',
*    &          ubavg(i,j,ml),ubavg(i,j,nl),
*    &          -thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
*    &          (vbavg(i  ,j,  mn)*depthv(i  ,j)
*    &          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
*    &          +vbavg(i-1,j,  mn)*depthv(i-1,j)
*    &          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
*    &          *(pvtrop(i,j)+pvtrop(i,j+1))
*    &          *.125 * dlt,utotn(i,j) * dlt,
*    &          0.5*(util5(i,j)+util5(i-1,j))
*           endif !debug
          endif !iu
        enddo !i
      enddo !j
c
      if     (ldebug_barotp) then
        call xcsync(flush_lp)
      endif
c
#if ! defined(RELO)
      if     (lbflag.eq.1) then
        call latbdp( nl)
      elseif (lbflag.eq.3) then
        call latbdf( nl,lll+1)
      endif
#endif
      if     (lbflag.eq.2) then
        call latbdt( nl,lll+1)
      elseif (lbflag.eq.4) then
        call latbdtf(nl,lll+1)
      endif
c
 840  continue  ! lll=1,lstep1,2
c
      if     (ldrag) then  !disp_count updated in momtum
        displd_mn(:,:) = displd_mn(:,:) + displd(:,:)/real(lstep1)
      endif
c
      if     (lbflag.eq.1) then
c
c ---   correct mean height.
c ---   this should not be required - so there may be a bug in the bc.
c
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1-margin,ii+margin
            if (SEA_P) then
              util1(i,j) = pbavg(i,j,n)*scp2(i,j)
            endif !ip
          enddo !i
        enddo !j
        call xcsum(sump, util1,ipa)
        q = sump/area
c
c ---   rhs: pbavg
c ---   lhs: pbavg
c
        margin = 0
c
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              pbavo(i,j)   = pbavo(i,j)   - q
              pbavg(i,j,1) = pbavg(i,j,1) - q
              pbavg(i,j,2) = pbavg(i,j,2) - q
              pbavg(i,j,3) = pbavg(i,j,3) - q
            endif !ip
          enddo !i
        enddo !j
      endif
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare(pbavg(1-nbdy,1-nbdy,1), ip,'barotp:pbav1')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,2), ip,'barotp:pbav2')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,3), ip,'barotp:pbav3')
      endif
c
      if     (btrlfr .and. delt1.ne.baclin) then  !not on very 1st time step
c ---   Robert-Asselin time filter 
!$OMP   PARALLEL DO PRIVATE(j,i,q)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            q = 0.5*ra2fac*(    pbavo(i,j)   +   !t-1
     &                          pbavg(i,j,n) -   !t+1
     &                      2.0*pbavg(i,j,m)  )  !t
            pbavg(i,j,m) = pbavg(i,j,m) + q
            q = 0.5*ra2fac*(    ubavo(i,j)   +   !t-1
     &                          ubavg(i,j,n) -   !t+1
     &                      2.0*ubavg(i,j,m)  )  !t
            ubavg(i,j,m) = ubavg(i,j,m) + q
            q = 0.5*ra2fac*(    vbavo(i,j)   +   !t-1
     &                          vbavg(i,j,n) -   !t+1
     &                      2.0*vbavg(i,j,m)  )  !t
            vbavg(i,j,m) = vbavg(i,j,m) + q
          enddo !i
        enddo !j
      endif !btrlfr & not on very 1st time step
c
      return
      end subroutine barotp
c
c> Revision history:
c>
c> Mar. 1995 - changed vertical velocity averaging interval from 10 cm to 1 m
c>             (loops 33,35)
c> Mar. 1995 - changed order of loop nesting in loop 842
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c> Aug. 1997 - transferred loops preceding loop 840 to momeq2.f
c> Jan. 2000 - added latbdp for lateral boundary ports
c> Aug. 2001 - two barotropic time steps per loop, for halo efficiency
c> Nov. 2006 - added lbflag==3 (latbdf) and thref_bt (mod_tides)
c> Nov. 2006 - removed thref_bt (and mod_tides)
c> Apr. 2007 - added btrlfr: leapfrog time step; see also momtum
c> Apr. 2010 - bugfixes for 1st time step and 1st miner time step
c> Apr  2011 - added    Robert-Asselin filtering for btrlfr
c> Aug  2011 - reworked Robert-Asselin filtering for btrlfr
c> Mar. 2012 - added latbdtf for nesting with Flather b.c.'s.
c> Jan. 2013 - added tidal drag tensor
c> June 2013 - added   lbflag==6 for latbdtc
c> Apr. 2014 - replace ip with ipa for mass sums
c> May  2014 - use land/sea masks (e.g. ip) to skip land
c> May  2014 - removed lbflag==6 for latbdtc
