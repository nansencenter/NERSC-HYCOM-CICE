module mod_mean
  implicit none

  ! --- HYCOM mean: array allocation and calculation interface.

  ! --- ii    = 1st dimension of array (==idm)
  ! --- jj    = 2nd dimension of array (==jdm)
  ! --- kk    = number of layers (typically 1)
  ! --- nmean = number of archive records in the mean
  ! --- ntracr= number of tracers to form the mean of

  integer, save :: ii,ii1,ii2,iorign,jj,jj1,jj2,jorign,kk,ntracr
  integer, save :: nmean,nstep,nstep_m

  ! --- lsteric = steric in input (and output) means
  ! --- lwtrflx = wtrflx in input (and output) means

  logical, save :: loneta, lsteric, lwtrflx
  logical, save :: meansq, single, trcout, icegln, hisurf, bgcout

  ! --- archive header

  character(len=80), save, dimension(4) :: ctitle

  ! --- arrays:

  real, save, allocatable, dimension (:,:,:,:) :: &
         tracer, &
         tracer_m

  real, save, allocatable, dimension (:,:,:) :: &
         u,  v,  ke,  temp,  saln,  th3d,  dp  ,visc,  tdff,  sdff,  dw,p,  &
         u_m,v_m,ke_m,temp_m,saln_m,th3d_m,dp_m,visc_m,tdff_m,sdff_m,dpu_m,dpv_m,dw_m

  real, save, allocatable, dimension (:,:) :: & 
         depths,depthu,depthv, &
         ubaro,vbaro,pbaro,kebaro, &
         montg,srfht,steric,oneta,onetaw,dpbl,dpmixl, &
         tmix,smix,thmix,umix,vmix,kemix, &
         surflx,salflx,wtrflx,covice,thkice,temice, &
         ubaro_m,vbaro_m,pbaro_m,kebaro_m, &
         montg_m,srfht_m,steric_m,dpbl_m,dpmixl_m, &
         tmix_m,smix_m,thmix_m,umix_m,vmix_m,kemix_m, &
         surflx_m,salflx_m,wtrflx_m,covice_m,thkice_m,temice_m, &
         oneta_m,onetaw_m,onetaw_u,onetaw_v, &
         si_u, si_um, si_v, si_vm, &
         surtx, surtxm, surty, surtym
  ! --- add two sea-ice variables for drift
  ! --- add x-/y surface stress

  real, save, allocatable, dimension (:) :: & 
         theta

  integer, save, allocatable, dimension (:,:) :: &
         ip,iq,iu,iv

  ! --- ECOSMO

  integer, parameter :: ntracr_bgc_2d = 10, ntracr_bgc_3d = 45
  
  character(len=8), save, allocatable, dimension(:) :: &
         nvar_bgc_2d

  character(len=8), save, allocatable, dimension(:) :: &
         nvar_bgc_3d

  real, save, allocatable, dimension (:,:,:,:) :: &
         tracer_bgc_3d, &
         tracer_bgc_3dm
  
  real, save, allocatable, dimension (:,:,:) :: &
         tracer_bgc_2d, &
         tracer_bgc_2dm

  ! --- module subroutines

  contains

  subroutine mean_alloc
  implicit none
  real, parameter :: spval = 2.0**100

  ! --- initialize allocatable arrays.

  ii1 = ii - 1
  ii2 = ii - 2
  jj1 = jj - 1
  jj2 = jj - 2

  nmean = 0

  lsteric = .false. !default
  lwtrflx = .false. !default
  loneta  = .false. !default

  ! ---
  
  allocate(      p(  ii,  jj,kk+1) )
  allocate(  dpu_m(  ii,  jj,kk  ) ); dpu_m = 0.0
  allocate(  dpv_m(  ii,  jj,kk  ) ); dpv_m = 0.0
  allocate( depthu(  ii,  jj     ) )
  allocate( depthv(  ii,  jj     ) )
  allocate( depths(0:ii,0:jj     ) )
  
  allocate( onetaw_u(ii,jj) );   onetaw_u = 0.0
  allocate( onetaw_v(ii,jj) );   onetaw_v = 0.0

  allocate( ip(ii,jj) )
  allocate( iq(ii,jj) )
  allocate( iu(ii,jj) )
  allocate( iv(ii,jj) )

  allocate( theta(kk) )

  ! --- 2D
  
  allocate(  montg(ii,jj),  montg_m(ii,jj) );  montg_m = 0.0
  allocate(  srfht(ii,jj),  srfht_m(ii,jj) );  srfht_m = 0.0
  allocate( steric(ii,jj), steric_m(ii,jj) ); steric_m = 0.0
  allocate(  oneta(ii,jj),  oneta_m(ii,jj) );  oneta_m = 0.0
  allocate( onetaw(ii,jj), onetaw_m(ii,jj) ); onetaw_m = 0.0
  allocate( surflx(ii,jj), surflx_m(ii,jj) ); surflx_m = 0.0
  allocate( wtrflx(ii,jj), wtrflx_m(ii,jj) ); wtrflx_m = 0.0
  allocate( salflx(ii,jj), salflx_m(ii,jj) ); salflx_m = 0.0
  allocate(   dpbl(ii,jj),   dpbl_m(ii,jj) );   dpbl_m = 0.0
  allocate( dpmixl(ii,jj), dpmixl_m(ii,jj) ); dpmixl_m = 0.0
  allocate(   tmix(ii,jj),   tmix_m(ii,jj) );   tmix_m = 0.0
  allocate(   smix(ii,jj),   smix_m(ii,jj) );   smix_m = 0.0
  allocate(  thmix(ii,jj),  thmix_m(ii,jj) );  thmix_m = 0.0
  allocate(   umix(ii,jj),   umix_m(ii,jj) );   umix_m = 0.0
  allocate(   vmix(ii,jj),   vmix_m(ii,jj) );   vmix_m = 0.0
  allocate(  kemix(ii,jj),  kemix_m(ii,jj) );  kemix_m = 0.0
  
  allocate( covice(ii,jj), covice_m(ii,jj) ); covice_m = 0.0
  allocate( thkice(ii,jj), thkice_m(ii,jj) ); thkice_m = 0.0
  allocate( temice(ii,jj), temice_m(ii,jj) ); temice_m = 0.0
  allocate(   si_u(ii,jj),    si_um(ii,jj) );    si_um = 0.0
  allocate(   si_v(ii,jj),    si_vm(ii,jj) );    si_vm = 0.0
  allocate(  surtx(ii,jj),   surtxm(ii,jj) );   surtxm = 0.0
  allocate(  surty(ii,jj),   surtym(ii,jj) );   surtym = 0.0

  allocate(  ubaro(ii,jj),  ubaro_m(ii,jj) );  ubaro_m = 0.0
  allocate(  vbaro(ii,jj),  vbaro_m(ii,jj) );  vbaro_m = 0.0
  allocate(  pbaro(ii,jj),  pbaro_m(ii,jj) );  pbaro_m = 0.0
  allocate( kebaro(ii,jj), kebaro_m(ii,jj) ); kebaro_m = 0.0

  ! 3D

  allocate(    u(ii,jj,kk),    u_m(ii,jj,kk) );    u_m = 0.0
  allocate(    v(ii,jj,kk),    v_m(ii,jj,kk) );    v_m = 0.0
  allocate(   ke(ii,jj,kk),   ke_m(ii,jj,kk) );   ke_m = 0.0
  allocate(   dp(ii,jj,kk),   dp_m(ii,jj,kk) );   dp_m = 0.0
  allocate(   dw(ii,jj,kk),   dw_m(ii,jj,kk) );   dw_m = 0.0
  allocate( temp(ii,jj,kk), temp_m(ii,jj,kk) ); temp_m = 0.0
  allocate( saln(ii,jj,kk), saln_m(ii,jj,kk) ); saln_m = 0.0
  allocate( th3d(ii,jj,kk), th3d_m(ii,jj,kk) ); th3d_m = 0.0
  allocate( visc(ii,jj,kk), visc_m(ii,jj,kk) ); visc_m = 0.0
  allocate( tdff(ii,jj,kk), tdff_m(ii,jj,kk) ); tdff_m = 0.0
  allocate( sdff(ii,jj,kk), sdff_m(ii,jj,kk) ); sdff_m = 0.0

  if (trcout) then
     allocate( tracer(  ii,jj,kk,ntracr) )
     allocate( tracer_m(ii,jj,kk,ntracr) ); tracer_m = 0.0
  endif

  ! --- ECOSMO

  if (bgcout) then
    allocate( nvar_bgc_2d(ntracr_bgc_2d) )
    allocate( nvar_bgc_3d(ntracr_bgc_3d) )
    
    allocate( tracer_bgc_2d(ii,jj,   ntracr_bgc_2d)) 
    allocate( tracer_bgc_2dm(ii,jj,   ntracr_bgc_2d)); tracer_bgc_2dm(:,:,:) = 0.0
    allocate( tracer_bgc_3d(ii,jj,kk,ntracr_bgc_3d)) 
    allocate( tracer_bgc_3dm(ii,jj,kk,ntracr_bgc_3d)); tracer_bgc_3dm(:,:,:,:) = 0.0

    ! register variable names (TODO: register in namelist)

    nvar_bgc_3d( 1) = 'CO2_c   ' 
    nvar_bgc_3d( 2) = 'CO2_TA  ' 
    nvar_bgc_3d( 3) = 'ECO_no3 '  
    nvar_bgc_3d( 4) = 'ECO_nh4 ' 
    nvar_bgc_3d( 5) = 'ECO_pho ' 
    nvar_bgc_3d( 6) = 'ECO_sil ' 
    nvar_bgc_3d( 7) = 'ECO_oxy ' 
    nvar_bgc_3d( 8) = 'ECO_fla ' 
    nvar_bgc_3d( 9) = 'ECO_dia ' 
    nvar_bgc_3d(10) = 'ECO_ccl ' 
    nvar_bgc_3d(11) = 'ECO_cclc' 
    nvar_bgc_3d(12) = 'ECO_caco' 
    nvar_bgc_3d(13) = 'ECO_diac' 
    nvar_bgc_3d(14) = 'ECO_flac' 
    nvar_bgc_3d(15) = 'ECO_micr' 
    nvar_bgc_3d(16) = 'ECO_meso' 
    nvar_bgc_3d(17) = 'ECO_det ' 
    nvar_bgc_3d(18) = 'ECO_opa ' 
    nvar_bgc_3d(19) = 'ECO_dom ' 
    nvar_bgc_3d(20) = 'ECO_dsnk' 
    nvar_bgc_3d(21) = 'CO2_pH  ' 
    nvar_bgc_3d(22) = 'CO2_pCO2' 
    nvar_bgc_3d(23) = 'CO2_Carb' 
    nvar_bgc_3d(24) = 'CO2_BiCa' 
    nvar_bgc_3d(25) = 'CO2_Carb' 
    nvar_bgc_3d(26) = 'CO2_Om_c' 
    nvar_bgc_3d(27) = 'CO2_Om_a' 
    nvar_bgc_3d(28) = 'ECO_prim' 
    nvar_bgc_3d(29) = 'ECO_secp' 
    nvar_bgc_3d(30) = 'ECO_netp' 
    nvar_bgc_3d(31) = 'ECO_parm' 
    nvar_bgc_3d(32) = 'ECO_Nlim' 
    nvar_bgc_3d(33) = 'ECO_Plim' 
    nvar_bgc_3d(34) = 'ECO_Slim' 
    nvar_bgc_3d(35) = 'ECO_Llim' 
    nvar_bgc_3d(36) = 'ECO_deni' 
    nvar_bgc_3d(37) = 'ECO_snks' 
    nvar_bgc_3d(38) = 'ECO_c2ch' 
    nvar_bgc_3d(39) = 'ECO_c2ch' 
    nvar_bgc_3d(40) = 'ECO_c2ch' 
    nvar_bgc_3d(41) = 'light_sw' 
    nvar_bgc_3d(42) = 'light_pa' 
    nvar_bgc_3d(43) = 'attenuat' 
    nvar_bgc_3d(44) = 'total_ch' 
    nvar_bgc_3d(45) = 'total_ca'
        
    nvar_bgc_2d( 1) = 'ECO_sed4'
    nvar_bgc_2d( 2) = 'ECO_sed1'
    nvar_bgc_2d( 3) = 'ECO_sed2'
    nvar_bgc_2d( 4) = 'ECO_sed3'
    nvar_bgc_2d( 5) = 'CO2_fair'
    nvar_bgc_2d( 6) = 'CO2_wind'
    nvar_bgc_2d( 7) = 'ECO_bots'
    nvar_bgc_2d( 8) = 'light_pa'
    nvar_bgc_2d( 9) = 'surface_'
    nvar_bgc_2d(10) = 'surface_'
  endif
  
  end subroutine mean_alloc

  subroutine mean_add(iweight)
  implicit none

  integer, intent(in) :: iweight

  ! --- add an archive to the mean.
  ! --- layer quantities weighted by layer thickness (i.e. by dw).

  integer             :: i,im,j,jm,k,ktr
  real                :: s,swk
  real, dimension(kk) :: sw

  nmean = nmean + iweight

  s = iweight

  p(:,:,:) = 0.0
  do j= 1,jj
  do i= 1,ii
     if (ip(i,j).eq.1) then
        p(i,j,1) = 0.0
        do k= 1,kk
          p(i,j,k+1) = p(i,j,k) + dw(i,j,k)
        enddo
  !     else
  !     p(i,j,:) = 0.0
     endif
  enddo
  enddo

  do j= 1,jj
  do i= 1,ii
     if (iu(i,j).eq.1) then
        ubaro_m(i,j) =  ubaro_m(i,j) + ubaro(i,j) * s
         umix_m(i,j) =   umix_m(i,j) +  umix(i,j) * s

        if (i.ne.1) then
           im = i-1
        else
           im = ii
        endif
        ! ---       depthu is either depths(i,j) or depths(i-1,j)
        if (depths(i,j).eq.depths(im,j)) then
           onetaw_u(i,j) = 0.5*(onetaw(i,j)+onetaw(im,j))
        elseif (depths(i,j).eq.depthu(i,j)) then
           onetaw_u(i,j) =      onetaw(i,j)
        else
           onetaw_u(i,j) =                  onetaw(im,j)
        endif
        do k= 1,kk
           swk = s*max(0.0, min(depthu(i,j),0.5*(p(i,j,k+1)+p(im,j,k+1))) &
                          - min(depthu(i,j),0.5*(p(i,j,k  )+p(im,j,k  ))) )*onetaw_u(i,j)
           dpu_m(i,j,k) = dpu_m(i,j,k) +            swk
             u_m(i,j,k) =   u_m(i,j,k) + u(i,j,k) * swk
        enddo
     endif !iu

     if (iv(i,j).eq.1) then
        vbaro_m(i,j) = vbaro_m(i,j) + vbaro(i,j) * s
         vmix_m(i,j) =  vmix_m(i,j) +  vmix(i,j) * s

        if (j.ne.1) then
           jm = j-1
        else
           jm = jj
        endif
        ! --- depthv is either depths(i,j) or depths(i,j-1)
        if     (depths(i,j).eq.depths(i,jm)) then
              onetaw_v(i,j) = 0.5*(onetaw(i,j)+onetaw(i,jm))
        elseif (depths(i,j).eq.depthv(i,j)) then
              onetaw_v(i,j) =      onetaw(i,j)
        else
              onetaw_v(i,j) =                  onetaw(i,jm)
        endif
        do k= 1,kk
          swk = s*max(0.0, min(depthv(i,j),0.5*(p(i,j,k+1)+p(i,jm,k+1))) &
                         - min(depthv(i,j),0.5*(p(i,j,k  )+p(i,jm,k  ))) )*onetaw_v(i,j)
          dpv_m(i,j,k) = dpv_m(i,j,k) +            swk
            v_m(i,j,k) =   v_m(i,j,k) + v(i,j,k) * swk
        enddo
     endif !iv

     if (ip(i,j).eq.1) then
         pbaro_m(i,j) =  pbaro_m(i,j)   +  pbaro(i,j)   * s
        kebaro_m(i,j) = kebaro_m(i,j)   + kebaro(i,j)   * s
         montg_m(i,j) =  montg_m(i,j)   +  montg(i,j)   * s
         srfht_m(i,j) =  srfht_m(i,j)   +  srfht(i,j)   * s
        steric_m(i,j) = steric_m(i,j)   + steric(i,j)   * s
          dpbl_m(i,j) =   dpbl_m(i,j)   +   dpbl(i,j)   * s
        dpmixl_m(i,j) = dpmixl_m(i,j)   + dpmixl(i,j)   * s
          tmix_m(i,j) =   tmix_m(i,j)   +   tmix(i,j)   * s
          smix_m(i,j) =   smix_m(i,j)   +   smix(i,j)   * s
         thmix_m(i,j) =  thmix_m(i,j)   +  thmix(i,j)   * s
         kemix_m(i,j) =  kemix_m(i,j)   +  kemix(i,j)   * s
        surflx_m(i,j) = surflx_m(i,j)   + surflx(i,j)   * s
        salflx_m(i,j) = salflx_m(i,j)   + salflx(i,j)   * s
        wtrflx_m(i,j) = wtrflx_m(i,j)   + wtrflx(i,j)   * s
        if (icegln) then
           covice_m(i,j) = covice_m(i,j) + covice(i,j) * s
           thkice_m(i,j) = thkice_m(i,j) + thkice(i,j) * s
           temice_m(i,j) = temice_m(i,j) + temice(i,j) * s
              si_um(i,j) =    si_um(i,j) +   si_u(i,j) * s
              si_vm(i,j) =    si_vm(i,j) +   si_v(i,j) * s
             surtxm(i,j) =   surtxm(i,j) +  surtx(i,j) * s
             surtym(i,j) =   surtym(i,j) +  surty(i,j) * s
        endif
        if (bgcout) then
          do ktr= 1,ntracr_bgc_2d
             tracer_bgc_2dm(i,j,ktr) = tracer_bgc_2dm(i,j,ktr) + tracer_bgc_2d(i,j,ktr) * s
          enddo !ktr
        endif

         oneta_m(i,j) =  oneta_m(i,j)   +  oneta(i,j)   * s
        onetaw_m(i,j) = onetaw_m(i,j)   + onetaw(i,j)   * s

                sw(:) = onetaw(i,j) * dw(i,j,:) * s

          dw_m(i,j,:) =     dw_m(i,j,:) +                 sw(:)
          dp_m(i,j,:) =     dp_m(i,j,:) +                 sw(:)
        temp_m(i,j,:) =   temp_m(i,j,:) +   temp(i,j,:) * sw(:)
        saln_m(i,j,:) =   saln_m(i,j,:) +   saln(i,j,:) * sw(:)
        th3d_m(i,j,:) =   th3d_m(i,j,:) +   th3d(i,j,:) * sw(:)
          ke_m(i,j,:) =     ke_m(i,j,:) +     ke(i,j,:) * sw(:)
        visc_m(i,j,:) =   visc_m(i,j,:) +   visc(i,j,:) * sw(:)
        tdff_m(i,j,:) =   tdff_m(i,j,:) +   tdff(i,j,:) * sw(:)
        sdff_m(i,j,:) =   sdff_m(i,j,:) +   sdff(i,j,:) * sw(:)

       if (trcout) then
          do ktr= 1,ntracr
             tracer_m(i,j,:,ktr) = tracer_m(i,j,:,ktr) + tracer(i,j,:,ktr) * sw(:)
          enddo !ktr
       endif
       
       if (bgcout) then
          do ktr= 1,ntracr_bgc_3d
             tracer_bgc_3dm(i,j,:,ktr) = tracer_bgc_3dm(i,j,:,ktr) + tracer_bgc_3d(i,j,:,ktr) * sw(:)
          enddo !ktr
       endif
     endif !ip
  enddo !i
  enddo !j

  !     write(6,*) 'mean_add   -    dp_m = ',   dp_m(54, 1,1),
  !    &                                        dp(  54, 1,1)

  end subroutine mean_add

  subroutine mean_addsq(iweight)
  implicit none

  integer, intent(in) :: iweight

  ! --- add an archive squared to the mean.
  ! --- layer quantities weighted by layer thickness (i.e. by dw).

  integer             :: i,im,j,jm,k,ktr
  real                :: s,swk
  real, dimension(kk) :: sw

  nmean = nmean + iweight

  s = iweight

  do j= 1,jj
  do i= 1,ii
     if (ip(i,j).eq.1) then
        p(i,j,1) = 0.0
        do k= 1,kk
           p(i,j,k+1) = p(i,j,k) + dw(i,j,k)
        enddo
  !  else
  !     p(i,j,:) = 0.0
     endif
  enddo
  enddo

  do j= 1,jj
  do i= 1,ii
     if (iu(i,j).eq.1) then
        ubaro_m(i,j) =  ubaro_m(i,j) + ubaro(i,j)**2 * s
         umix_m(i,j) =   umix_m(i,j) +  umix(i,j)**2 * s

        if (i.ne.1) then
           im = i-1
        else
           im = ii
        endif
        ! ---       depthu is either depths(i,j) or depths(i-1,j)
        if (depths(i,j).eq.depths(im,j)) then
           onetaw_u(i,j) = 0.5*(onetaw(i,j)+onetaw(im,j))
        elseif (depths(i,j).eq.depthu(i,j)) then
           onetaw_u(i,j) =      onetaw(i,j)
        else
           onetaw_u(i,j) =                  onetaw(im,j)
        endif
        do k= 1,kk
           swk = s*max(0.0, min(depthu(i,j), 0.5*(p(i,j,k+1)+p(im,j,k+1))) &
                          - min(depthu(i,j), 0.5*(p(i,j,k  )+p(im,j,k  ))) )*onetaw_u(i,j)
           dpu_m(i,j,k) = dpu_m(i,j,k) +               swk
             u_m(i,j,k) =   u_m(i,j,k) + u(i,j,k)**2 * swk
        enddo
     endif !iu

     if (iv(i,j).eq.1) then
        vbaro_m(i,j) = vbaro_m(i,j) + vbaro(i,j)**2 * s
         vmix_m(i,j) =  vmix_m(i,j) +  vmix(i,j)**2 * s

        if (j.ne.1) then
           jm = j-1
        else
           jm = jj
        endif
        ! ---       depthv is either depths(i,j) or depths(i,j-1)
        if     (depths(i,j).eq.depths(i,jm)) then
           onetaw_v(i,j) = 0.5*(onetaw(i,j)+onetaw(i,jm))
        elseif (depths(i,j).eq.depthv(i,j)) then
           onetaw_v(i,j) =      onetaw(i,j)
        else
           onetaw_v(i,j) =                  onetaw(i,jm)
        endif
        do k= 1,kk
           swk = s*max( 0.0, min(depthv(i,j),0.5*(p(i,j,k+1)+p(i,jm,k+1))) &
                           - min(depthv(i,j),0.5*(p(i,j,k  )+p(i,jm,k  ))) )*onetaw_v(i,j)
           dpv_m(i,j,k) = dpv_m(i,j,k) +               swk
             v_m(i,j,k) =   v_m(i,j,k) + v(i,j,k)**2 * swk
        enddo
     endif !iv

     if (ip(i,j).eq.1) then
         pbaro_m(i,j)   =  pbaro_m(i,j)   +  pbaro(i,j)**2   * s
        kebaro_m(i,j)   = kebaro_m(i,j)   + kebaro(i,j)**2   * s
         montg_m(i,j)   =  montg_m(i,j)   +  montg(i,j)**2   * s
         srfht_m(i,j)   =  srfht_m(i,j)   +  srfht(i,j)**2   * s
        steric_m(i,j)   = steric_m(i,j)   + steric(i,j)**2   * s
          dpbl_m(i,j)   =   dpbl_m(i,j)   +   dpbl(i,j)**2   * s
        dpmixl_m(i,j)   = dpmixl_m(i,j)   + dpmixl(i,j)**2   * s
          tmix_m(i,j)   =   tmix_m(i,j)   +   tmix(i,j)**2   * s
          smix_m(i,j)   =   smix_m(i,j)   +   smix(i,j)**2   * s
         thmix_m(i,j)   =  thmix_m(i,j)   +  thmix(i,j)**2   * s
         kemix_m(i,j)   =  kemix_m(i,j)   +  kemix(i,j)**2   * s
        surflx_m(i,j)   = surflx_m(i,j)   + surflx(i,j)**2   * s
        salflx_m(i,j)   = salflx_m(i,j)   + salflx(i,j)**2   * s
        wtrflx_m(i,j)   = wtrflx_m(i,j)   + wtrflx(i,j)**2   * s
        if (icegln) then
           covice_m(i,j) = covice_m(i,j) + covice(i,j)**2 * s
           thkice_m(i,j) = thkice_m(i,j) + thkice(i,j)**2 * s
           temice_m(i,j) = temice_m(i,j) + temice(i,j)**2 * s
              si_um(i,j) =    si_um(i,j) +   si_u(i,j)**2 * s
              si_vm(i,j) =    si_vm(i,j) +   si_v(i,j)**2 * s
             surtxm(i,j) =   surtxm(i,j) +  surtx(i,j)**2 * s
             surtym(i,j) =   surtym(i,j) +  surty(i,j)**2 * s
        endif
        if (bgcout) then
           do ktr= 1,ntracr_bgc_2d
              tracer_bgc_2dm(i,j,ktr) = tracer_bgc_2dm(i,j,ktr) + tracer_bgc_2d(i,j,ktr)**2 *s
           enddo !ktr
        endif
         oneta_m(i,j)   =  oneta_m(i,j)   +  oneta(i,j)**2 * s
        onetaw_m(i,j)   = onetaw_m(i,j)   + onetaw(i,j)    * s
            dp_m(i,j,:) =     dp_m(i,j,:) + dp(i,j,:)**2   * s

                  sw(:) =   onetaw(i,j) * dw(i,j,:) * s
            dw_m(i,j,:) =   dw_m(i,j,:) +                  sw(:)
            ke_m(i,j,:) =   ke_m(i,j,:) +   ke(i,j,:)**2 * sw(:)
          temp_m(i,j,:) = temp_m(i,j,:) + temp(i,j,:)**2 * sw(:)
          saln_m(i,j,:) = saln_m(i,j,:) + saln(i,j,:)**2 * sw(:)
          th3d_m(i,j,:) = th3d_m(i,j,:) + th3d(i,j,:)**2 * sw(:)
          visc_m(i,j,:) = visc_m(i,j,:) + visc(i,j,:)**2 * sw(:)
          tdff_m(i,j,:) = tdff_m(i,j,:) + tdff(i,j,:)**2 * sw(:)
          sdff_m(i,j,:) = sdff_m(i,j,:) + sdff(i,j,:)**2 * sw(:)
        if (trcout) then
           do ktr= 1,ntracr
              tracer_m(i,j,:,ktr) = tracer_m(i,j,:,ktr) + tracer(i,j,:,ktr)**2 * sw(:)
           enddo !ktr
        endif
        if (bgcout) then
           do ktr= 1,ntracr_bgc_3d
              tracer_bgc_3dm(i,j,:,ktr) = tracer_bgc_3dm(i,j,:,ktr) + tracer_bgc_3d(i,j,:,ktr)**2 * sw(:)
           enddo !ktr
        endif
     endif !ip
  enddo
  enddo

  ! write(6,*) 'mean_addsq -    dp_m = ',   dp_m(54, 1,1),
  !    &                                        dp(  54, 1,1)**2

  end subroutine mean_addsq

  subroutine mean_box(a,li,lj,lk,nb,larctic,lperiod)
  implicit none

  logical :: larctic,lperiod
  integer :: li,lj,lk,nb
  real, dimension(li,lj,lk) :: a

  ! --- form the 2*nb+1 square running average of a(:,:,k)

  real, parameter   :: spval = 2.0**100
  integer           :: i,iq,j,jq,k
  real              :: rs,qc
  real, allocatable :: b(:,:),s(:),q(:)

  allocate( b(1-nb:li+nb,1-nb:lj+nb),s(1-nb:li+nb),q(1-nb:li+nb) )

  do k= 1,lk

     do j= 1,lj
     do i= 1,li
        b(i,j) = a(i,j,k)
     enddo
     enddo

     if (.not.lperiod) then !closed domain
        do j = 1,lj
        do iq = 1,nb
           b( 1-iq,j) = spval
           b(li+iq,j) = spval
        enddo !iq
        enddo !j
        do i = 1-nb,li+nb
        do jq = 1,nb
           b(i, 1-jq) = spval
           b(i,lj+jq) = spval
        enddo !jq
        enddo !i
     elseif (.not.larctic) then !lperiod only, i.e. near-global
        do jq = 1,nb
        do i = 1,li
           b(i, 1-jq) = spval          !closed bottom boundary
           b(i,lj+jq) = spval          !closed top    boundary
        enddo !i
        enddo !jq
        do j = 1-nb,lj+nb
        do iq = 1,nb
           b(-nb+iq,j) = b(li-nb+iq,j)  !periodic in longitude
           b( li+iq,j) = b(      iq,j)  !periodic in longitude
        enddo !iq
        enddo !j
     else !global with arctic patch
        do jq = 1,nb
        do i = 1,li
           b(i,1-jq) = spval !closed bottom boundary
        enddo !i
        enddo !jq
        do j= lj+1,lj+nb
           jq = lj-1-(j-lj)
           do i= 1,li
              iq = li-mod(i-1,li)
              b(i,j) = b(iq,jq)  !arctic patch across top boundary
           enddo !i
        enddo !j
        do j= 1-nb,lj+nb
        do iq= 1,nb
           b(-nb+iq,j) = b(li-nb+iq,j)  !periodic in longitude
           b( li+iq,j) = b(      iq,j)  !periodic in longitude
        enddo !iq
        enddo !j
     endif !domain type

     do j= 1,lj
        do i= 1-nb,li+nb
        rs = 0.0
        qc = 0.0
        do jq= -nb,nb
           if (b(i,j+jq).ne.spval) then
              rs = rs + b(i,j+jq)
              qc = qc + 1.0
           endif
        enddo !jq
        s(i) = rs
        q(i) = qc
        enddo !i
        do i= 1,li
           if (b(i,j) .ne. spval) then
              rs = 0.0
              qc = 0.0
              do iq= -nb,nb
                 rs = rs + s(i+iq)
                 qc = qc + q(i+iq)
              enddo !iq
              a(i,j,k) = rs/qc  !qc can't be zero, since b(i,j).ne.spval
           else
              a(i,j,k) = spval
           endif
        enddo !i
     enddo !j

  enddo !k

  deallocate( b,s,q )
  end subroutine mean_box

  subroutine mean_copy
  implicit none

  ! --- copy archive to mean archive

     nmean = nstep

       u_m =      u
       v_m =      v
      ke_m =     ke
    temp_m =   temp
    saln_m =   saln
    th3d_m =   th3d
      dp_m =     dp
      dw_m =     dw
    visc_m =   visc
    tdff_m =   tdff
    sdff_m =   sdff

  if (trcout) then
     tracer_m = tracer
  endif   

  if (bgcout) then
     tracer_bgc_3dm = tracer_bgc_3d
  endif   

  ubaro_m =  ubaro
   vbaro_m =  vbaro
   pbaro_m =  pbaro
  kebaro_m = kebaro
   montg_m =  montg
   srfht_m =  srfht
  steric_m = steric
    dpbl_m =   dpbl
  dpmixl_m = dpmixl
    tmix_m =   tmix
    smix_m =   smix
   thmix_m =  thmix
    umix_m =   umix
    vmix_m =   vmix
   kemix_m =  kemix
  surflx_m = surflx
  salflx_m = salflx
  wtrflx_m = wtrflx
  if (icegln) then
     covice_m = covice
     thkice_m = thkice
     temice_m = temice
      oneta_m =  oneta
     onetaw_m = onetaw
  endif
  if (bgcout) then
     tracer_bgc_2dm = tracer_bgc_2d
  endif   
  if (icegln) then
     si_um =   si_u
     si_vm =   si_v 
    surtxm =  surtx
    surtym =  surty
  endif

  !     write(6,*) 'mean_copy  -    dp_m = ',   dp_m(54, 1,1),
  !    &                                        dp(  54, 1,1)

  end subroutine mean_copy

  subroutine mean_depths
  implicit none

  ! --- calculate depthu and depthv

  integer :: i,im,j,jm

  depths(:,:) = 9806.0 * depths(:,:)  ! convert to pressure units

  do j= 1,jj
  do i= 1,ii
     if (i.ne.1) then
        im = i-1
     else
        im = ii
     endif
     if (min(ip(i,j),ip(im,j)).eq.1) then
        depthu(i,j) = min(depths(i,j),depths(im,j))
     elseif (ip(i ,j).eq.1) then
        depthu(i,j) = depths(i ,j)
     elseif (ip(im,j).eq.1) then
        depthu(i,j) = depths(im,j)
     else
        depthu(i,j) = 0.0
     endif

     if (j.ne.1) then
        jm = j-1
     else
        jm = jj
     endif
     if (min(ip(i,j),ip(i,jm)).eq.1) then
        depthv(i,j) = min(depths(i,j),depths(i,jm))
     elseif (ip(i,j) .eq.1) then
        depthv(i,j) = depths(i,j)
     elseif (ip(i,jm).eq.1) then
        depthv(i,j) = depths(i,jm)
     else
        depthv(i,j) = 0.0
     endif
  enddo
  enddo

  end subroutine mean_depths

  subroutine mean_diff(nscale,nbox)
  implicit none

  integer nscale,nbox

  ! --- form the difference of two archives, 1st already in _m

  real, parameter :: zero = 0.0

  real :: q
  logical :: larctic,lperiod
  integer :: i,j,k,ktr

  nmean = 2  !always 2 for diff archives
  q = 1.0/real(nscale)

  do j= 1,jj
  do i= 1,ii
     do k= 1,kk
        if (iu(i,j).eq.1) then
           u_m(i,j,k)    = q*(   u_m(i,j,k) -    u(i,j,k))
        endif
        if (iv(i,j).eq.1) then
           v_m(i,j,k)    = q*(   v_m(i,j,k) -    v(i,j,k))
        endif
        if (ip(i,j).eq.1) then
             dw_m(i,j,k) = dp(i,j,k)
             dp_m(i,j,k) = q*(  dp_m(i,j,k) -   dp(i,j,k))
             ke_m(i,j,k) = q*(  ke_m(i,j,k) -   ke(i,j,k))
           temp_m(i,j,k) = q*(temp_m(i,j,k) - temp(i,j,k))
           saln_m(i,j,k) = q*(saln_m(i,j,k) - saln(i,j,k))
           th3d_m(i,j,k) = q*(th3d_m(i,j,k) - th3d(i,j,k))
           visc_m(i,j,k) = q*(visc_m(i,j,k) - visc(i,j,k))
           tdff_m(i,j,k) = q*(tdff_m(i,j,k) - tdff(i,j,k))
           sdff_m(i,j,k) = q*(sdff_m(i,j,k) - sdff(i,j,k))
           if (trcout) then
              do ktr= 1,ntracr
                 tracer_m(i,j,k,ktr) = q*(tracer_m(i,j,k,ktr) - tracer(i,j,k,ktr))
              enddo !ktr
           endif
           if (bgcout) then
              do ktr= 1,ntracr_bgc_3d
                 tracer_bgc_3dm(i,j,k,ktr) = q*(tracer_bgc_3dm(i,j,k,ktr) - tracer_bgc_3d(i,j,k,ktr))
              enddo !ktr
           endif
        endif
     enddo !k

     if (iu(i,j).eq.1) then
         ubaro_m(i,j) = q*( ubaro_m(i,j) -  ubaro(i,j))
          umix_m(i,j) = q*(  umix_m(i,j) -   umix(i,j))
     endif
     if (iv(i,j).eq.1) then
         vbaro_m(i,j) = q*( vbaro_m(i,j) -  vbaro(i,j))
          vmix_m(i,j) = q*(  vmix_m(i,j) -   vmix(i,j))
     endif
     if (ip(i,j).eq.1) then
         pbaro_m(i,j) = q*( pbaro_m(i,j) -  pbaro(i,j))
        kebaro_m(i,j) = q*(kebaro_m(i,j) - kebaro(i,j))
         montg_m(i,j) = q*( montg_m(i,j) -  montg(i,j))
         srfht_m(i,j) = q*( srfht_m(i,j) -  srfht(i,j))
        steric_m(i,j) = q*(steric_m(i,j) - steric(i,j))
          dpbl_m(i,j) = q*(  dpbl_m(i,j) -   dpbl(i,j))
        dpmixl_m(i,j) = q*(dpmixl_m(i,j) - dpmixl(i,j))
          tmix_m(i,j) = q*(  tmix_m(i,j) -   tmix(i,j))
          smix_m(i,j) = q*(  smix_m(i,j) -   smix(i,j))
         thmix_m(i,j) = q*( thmix_m(i,j) -  thmix(i,j))
         kemix_m(i,j) = q*( kemix_m(i,j) -  kemix(i,j))
        surflx_m(i,j) = q*(surflx_m(i,j) - surflx(i,j))
        salflx_m(i,j) = q*(salflx_m(i,j) - salflx(i,j))
        wtrflx_m(i,j) = q*(wtrflx_m(i,j) - wtrflx(i,j))
        if (icegln) then
           covice_m(i,j) = q*(covice_m(i,j) - covice(i,j))
           thkice_m(i,j) = q*(thkice_m(i,j) - thkice(i,j))
           temice_m(i,j) = q*(temice_m(i,j) - temice(i,j))
        endif
        if (bgcout) then
           do ktr= 1,ntracr_bgc_2d
              tracer_bgc_2dm(i,j,ktr) = q*(tracer_bgc_2dm(i,j,ktr) - tracer_bgc_2d(i,j,ktr))
           enddo !ktr
        endif
        if (icegln) then
           si_um(i,j) = q*(   si_um(i,j) -   si_u(i,j))
           si_vm(i,j) = q*(   si_vm(i,j) -   si_v(i,j))
          surtxm(i,j) = q*(  surtxm(i,j) -  surtx(i,j))
          surtym(i,j) = q*(  surtym(i,j) -  surty(i,j))
        endif
     endif
  enddo !i
  enddo !j

  if (nbox.ne.0) then
     larctic = any( ip(:,jj).eq.1 )  !land on last row
     lperiod = any( ip(ii,:).eq.1 )  !land on last column

     call mean_box(     u_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(     v_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(  temp_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(  saln_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(  th3d_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(    dp_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(    ke_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(  visc_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(  tdff_m,ii,jj,kk,nbox,larctic,lperiod)
     call mean_box(  sdff_m,ii,jj,kk,nbox,larctic,lperiod)

     if (trcout) then
        do ktr= 1,ntracr
           call mean_box(tracer_m(1,1,1,ktr),ii,jj,kk,nbox,larctic,lperiod)
        enddo !ktr
     endif
     if (bgcout) then
        do ktr= 1,ntracr_bgc_3d
           call mean_box(tracer_bgc_3dm(1,1,1,ktr),ii,jj,kk,nbox,larctic,lperiod)
        enddo !ktr
     endif

     call mean_box( ubaro_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(  umix_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box( vbaro_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(  vmix_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box( pbaro_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(kebaro_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box( montg_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box( srfht_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(steric_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(  dpbl_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(dpmixl_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(  tmix_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(  smix_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box( thmix_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box( kemix_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(surflx_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(salflx_m,ii,jj, 1,nbox,larctic,lperiod)
     call mean_box(wtrflx_m,ii,jj, 1,nbox,larctic,lperiod)
     if (icegln) then
        call mean_box(covice_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(thkice_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(temice_m,ii,jj, 1,nbox,larctic,lperiod)
     endif
     if (bgcout) then
        do ktr= 1,ntracr_bgc_2d
           call mean_box(tracer_bgc_2dm(1,1,ktr),ii,jj,1,nbox,larctic,lperiod)
        enddo !ktr
     endif
     if (icegln) then
        call mean_box(   si_um,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(   si_vm,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(  surtxm,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(  surtym,ii,jj, 1,nbox,larctic,lperiod)
     endif
  endif !nbox

  !     write(6,*) 'mean_diff  -    dp_m = ',   dp_m(54, 1,1),
  !    &                                        dw_m(54, 1,1),
  !    &                                        dp(  54, 1,1)

  end subroutine mean_diff

  subroutine mean_end
  implicit none

! --- reduce sum of archives to their mean.

  real, parameter :: spval = 2.0**100

  integer :: i,im,j,jm,k,ktr
  real :: s,swk
  real, dimension(kk) :: sw(kk)

  s = 1.0/nmean

  do j= 1,jj
  do i= 1,ii
     if (iu(i,j).eq.1) then
        if (i.ne.1) then
           im = i-1
        else
           im = ii
        endif
        do k= 1,kk
           swk = dpu_m(i,j,k) * s
           if (swk.ge.0.000001) then
              swk = s/swk
              u_m(i,j,k) = u_m(i,j,k) * swk
           else  ! project into zero thickness layers
              u_m(i,j,k) = u_m(i,j,k-1)
           endif
        enddo
        ubaro_m(i,j)   = ubaro_m(i,j) * s
         umix_m(i,j)   =  umix_m(i,j) * s
     else
            u_m(i,j,:) = spval
        ubaro_m(i,j)   = spval
         umix_m(i,j)   = spval
     endif !iu

     if (iv(i,j).eq.1) then
        if (j.ne.1) then
           jm = j-1
        else
           jm = jj
        endif
        do k= 1,kk
           swk = dpv_m(i,j,k) * s
           if (swk.ge.0.000001) then
              swk = s/swk
              v_m(i,j,k) = v_m(i,j,k) * swk
           else  ! project into zero thickness layers
              v_m(i,j,k) = v_m(i,j,k-1)
           endif
        enddo
           vbaro_m(i,j)   = vbaro_m(i,j) * s
            vmix_m(i,j)   =  vmix_m(i,j) * s
     else
               v_m(i,j,:) = spval
           vbaro_m(i,j)   = spval
            vmix_m(i,j)   = spval
     endif

     if (ip(i,j).eq.1) then
         pbaro_m(i,j)  =  pbaro_m(i,j) * s
        kebaro_m(i,j)  = kebaro_m(i,j) * s
         montg_m(i,j)  =  montg_m(i,j) * s
         srfht_m(i,j)  =  srfht_m(i,j) * s
        steric_m(i,j)  = steric_m(i,j) * s
          dpbl_m(i,j)  =   dpbl_m(i,j) * s
        dpmixl_m(i,j)  = dpmixl_m(i,j) * s
          tmix_m(i,j)  =   tmix_m(i,j) * s
          smix_m(i,j)  =   smix_m(i,j) * s
         thmix_m(i,j)  =  thmix_m(i,j) * s
         kemix_m(i,j)  =  kemix_m(i,j) * s
        surflx_m(i,j)  = surflx_m(i,j) * s
        salflx_m(i,j)  = salflx_m(i,j) * s
        wtrflx_m(i,j)  = wtrflx_m(i,j) * s
        if (icegln) then
           covice_m(i,j)  = covice_m(i,j) * s
           thkice_m(i,j)  = thkice_m(i,j) * s
           temice_m(i,j)  = temice_m(i,j) * s
        endif
         oneta_m(i,j)  =  oneta_m(i,j) * s
        onetaw_m(i,j)  = onetaw_m(i,j) * s
        if (bgcout) then
           do ktr= 1,ntracr_bgc_2d
              tracer_bgc_2dm(i,j,ktr) = tracer_bgc_2dm(i,j,ktr) * s
           enddo !ktr
        endif
        if (icegln) then
           si_um(i,j)  =    si_um(i,j) * s
           si_vm(i,j)  =    si_vm(i,j) * s
          surtxm(i,j)  =   surtxm(i,j) * s
          surtym(i,j)  =   surtym(i,j) * s
        endif

        do k= 1,kk
           dw_m(i,j,k) = dw_m(i,j,k) * s
           dp_m(i,j,k) = dp_m(i,j,k) * s
           if (dw_m(i,j,k).ge.0.000001) then
                        swk = s/dw_m(i,j,k)
                ke_m(i,j,k) =     ke_m(i,j,k) * swk
              temp_m(i,j,k) =   temp_m(i,j,k) * swk
              saln_m(i,j,k) =   saln_m(i,j,k) * swk
              th3d_m(i,j,k) =   th3d_m(i,j,k) * swk
              visc_m(i,j,k) =   visc_m(i,j,k) * swk
              tdff_m(i,j,k) =   tdff_m(i,j,k) * swk
              sdff_m(i,j,k) =   sdff_m(i,j,k) * swk
              if (trcout) then
                 do ktr= 1,ntracr
                    tracer_m(i,j,k,ktr) = tracer_m(i,j,k,ktr) * swk
                 enddo !ktr
              endif
              if (bgcout) then
                 do ktr= 1,ntracr_bgc_3d
                    tracer_bgc_3dm(i,j,k,ktr) = tracer_bgc_3dm(i,j,k,ktr) * swk
                 enddo !ktr
              endif
           else  ! project into zero thickness layers
                ke_m(i,j,k) =     ke_m(i,j,k-1)
              temp_m(i,j,k) =   temp_m(i,j,k-1)
              saln_m(i,j,k) =   saln_m(i,j,k-1)
              th3d_m(i,j,k) =   th3d_m(i,j,k-1)
              visc_m(i,j,k) =   visc_m(i,j,k-1)
              tdff_m(i,j,k) =   tdff_m(i,j,k-1)
              sdff_m(i,j,k) =   sdff_m(i,j,k-1)
              if (trcout) then
                 do ktr= 1,ntracr
                    tracer_m(i,j,k,ktr) = tracer_m(i,j,k-1,ktr)
                 enddo !ktr
              endif
              if (bgcout) then
                 do ktr= 1,ntracr_bgc_3d
                    tracer_bgc_3dm(i,j,k,ktr) = tracer_bgc_3dm(i,j,k-1,ktr)
                 enddo !ktr
              endif
           endif
           ! ---         archived dp_m is based on dp' (dp_m/oneta_m)
           dp_m(i,j,k)   = dp_m(i,j,k)/onetaw_m(i,j)
           dw_m(i,j,k)   = dw_m(i,j,k)/onetaw_m(i,j)
              p(i,j,k+1) = dw_m(i,j,k) + p(i,j,k)
        enddo
     else
         pbaro_m(i,j)   = spval
        kebaro_m(i,j)   = spval
         montg_m(i,j)   = spval
         srfht_m(i,j)   = spval
        steric_m(i,j)   = spval
          dpbl_m(i,j)   = spval
        dpmixl_m(i,j)   = spval
          tmix_m(i,j)   = spval
          smix_m(i,j)   = spval
         thmix_m(i,j)   = spval
         kemix_m(i,j)   = spval
        surflx_m(i,j)   = spval
        salflx_m(i,j)   = spval
        wtrflx_m(i,j)   = spval
        if (icegln) then
           covice_m(i,j)   = spval
           thkice_m(i,j)   = spval
           temice_m(i,j)   = spval
        endif
        if (bgcout) then
           do ktr= 1,ntracr_bgc_2d
              tracer_bgc_2dm(i,j,ktr) = spval
           enddo !ktr
        endif 
        if (icegln) then
           si_um(i,j)   = spval
           si_vm(i,j)   = spval
          surtxm(i,j)   = spval
          surtym(i,j)   = spval
        endif
         oneta_m(i,j)   = spval
        onetaw_m(i,j)   = spval

            dw_m(i,j,:) = spval
            dp_m(i,j,:) = spval
            ke_m(i,j,:) = spval
          temp_m(i,j,:) = spval
          saln_m(i,j,:) = spval
          th3d_m(i,j,:) = spval
          visc_m(i,j,:) = spval
          tdff_m(i,j,:) = spval
          sdff_m(i,j,:) = spval
        if (trcout) then
           do ktr= 1,ntracr
              tracer_m(i,j,:,ktr) = spval
           enddo !ktr
        endif 
        if (bgcout) then
           do ktr= 1,ntracr_bgc_3d
              tracer_bgc_3dm(i,j,:,ktr) = spval
           enddo !ktr
        endif 
     endif
  enddo
  enddo

  ! write(6,*) 'mean_end   -    dp_m = ',   dp_m(54, 1,1)

  end subroutine mean_end

  subroutine mean_std
  implicit none

  ! --- form the std.dev = sqrt(mnsq-mean**2)

  real, parameter :: zero = 0.0
  integer :: i,j,k,ktr
  real :: std,x
  
  std(x) = sqrt(max(zero,x))

  do j= 1,jj
  do i= 1,ii
     do k= 1,kk
        if (iu(i,j).eq.1) then
              u_m(i,j,k) = std(   u(i,j,k) -    u_m(i,j,k)**2)
        endif
        if (iv(i,j).eq.1) then
              v_m(i,j,k) = std(   v(i,j,k) -    v_m(i,j,k)**2)
        endif
        if (ip(i,j).eq.1) then
             dw_m(i,j,k) =                     dp_m(i,j,k)
             dp_m(i,j,k) = std(  dp(i,j,k) -   dp_m(i,j,k)**2)
             ke_m(i,j,k) = std(  ke(i,j,k) -   ke_m(i,j,k)**2)
           temp_m(i,j,k) = std(temp(i,j,k) - temp_m(i,j,k)**2)
           saln_m(i,j,k) = std(saln(i,j,k) - saln_m(i,j,k)**2)
           th3d_m(i,j,k) = std(th3d(i,j,k) - th3d_m(i,j,k)**2)
           visc_m(i,j,k) = std(visc(i,j,k) - visc_m(i,j,k)**2)
           tdff_m(i,j,k) = std(tdff(i,j,k) - tdff_m(i,j,k)**2)
           sdff_m(i,j,k) = std(sdff(i,j,k) - sdff_m(i,j,k)**2)
           if (trcout) then
              do ktr= 1,ntracr
                 tracer_m(i,j,k,ktr) = std(tracer(i,j,k,ktr) - tracer_m(i,j,k,ktr)**2)
              enddo !ktr
           endif
           if (bgcout) then
              do ktr= 1,ntracr_bgc_3d
                 tracer_bgc_3dm(i,j,k,ktr) = std(tracer_bgc_3d(i,j,k,ktr) - tracer_bgc_3dm(i,j,k,ktr)**2)
              enddo !ktr
           endif
        endif
     enddo

     if (iu(i,j).eq.1) then
        ubaro_m(i,j) = std(  ubaro(i,j) -  ubaro_m(i,j)**2)
         umix_m(i,j) = std(   umix(i,j) -   umix_m(i,j)**2)
     endif
     if (iv(i,j).eq.1) then
        vbaro_m(i,j) = std(  vbaro(i,j) -  vbaro_m(i,j)**2)
         vmix_m(i,j) = std(   vmix(i,j) -   vmix_m(i,j)**2)
     endif
     if (ip(i,j).eq.1) then
         pbaro_m(i,j) = std( pbaro(i,j) -  pbaro_m(i,j)**2)
        kebaro_m(i,j) = std(kebaro(i,j) - kebaro_m(i,j)**2)
         montg_m(i,j) = std( montg(i,j) -  montg_m(i,j)**2)
         srfht_m(i,j) = std( srfht(i,j) -  srfht_m(i,j)**2)
        steric_m(i,j) = std(steric(i,j) - steric_m(i,j)**2)
          dpbl_m(i,j) = std(  dpbl(i,j) -   dpbl_m(i,j)**2)
        dpmixl_m(i,j) = std(dpmixl(i,j) - dpmixl_m(i,j)**2)
          tmix_m(i,j) = std(  tmix(i,j) -   tmix_m(i,j)**2)
          smix_m(i,j) = std(  smix(i,j) -   smix_m(i,j)**2)
         thmix_m(i,j) = std( thmix(i,j) -  thmix_m(i,j)**2)
         kemix_m(i,j) = std( kemix(i,j) -  kemix_m(i,j)**2)
        surflx_m(i,j) = std(surflx(i,j) - surflx_m(i,j)**2)
        salflx_m(i,j) = std(salflx(i,j) - salflx_m(i,j)**2)
        wtrflx_m(i,j) = std(wtrflx(i,j) - wtrflx_m(i,j)**2)
        if (icegln) then
           covice_m(i,j) = std(covice(i,j) - covice_m(i,j)**2)
           thkice_m(i,j) = std(thkice(i,j) - thkice_m(i,j)**2)
           temice_m(i,j) = std(temice(i,j) - temice_m(i,j)**2)
        endif
        onetaw_m(i,j) = onetaw_m(i,j)
         oneta_m(i,j) = std( oneta(i,j) -  oneta_m(i,j)**2)
        if (bgcout) then
           do ktr= 1,ntracr_bgc_2d
              tracer_bgc_2dm(i,j,ktr) = std(tracer_bgc_2d(i,j,ktr) - tracer_bgc_2dm(i,j,ktr)**2)
           enddo !ktr
        endif
        if (icegln) then
           si_um(i,j) = std(  si_u(i,j) -    si_um(i,j)**2)
           si_vm(i,j) = std(  si_v(i,j) -    si_vm(i,j)**2)
          surtxm(i,j) = std( surtx(i,j) -   surtxm(i,j)**2)
          surtym(i,j) = std( surty(i,j) -   surtym(i,j)**2)
        endif
     endif
  enddo
  enddo

  !     write(6,*) 'mean_std   -    dp_m = ',   dp_m(54, 1,1),
  !    &                                        dw_m(54, 1,1),
  !    &                                        dp(  54, 1,1)

  end subroutine mean_std

  subroutine mean_velocity
  implicit none

  ! --- update velocity to include depth averaged component, and
  ! --- calculate kinetic energy.
  ! --- only called for standard archive fields.

  integer :: i,ia,ip1,j

  do j= 1,jj
  do i= 1,ii
     if (iu(i,j).eq.1) then
           u(i,j,:) =    u(i,j,:) + ubaro(i,j)
        umix(i,j)   = umix(i,j)   + ubaro(i,j)
     endif
     if (iv(i,j).eq.1) then
           v(i,j,:) =    v(i,j,:) + vbaro(i,j)
        vmix(i,j)   = vmix(i,j)   + vbaro(i,j)
     endif
  enddo
  enddo

  do j= 1,jj-1
  do i= 1,ii
     if (i.ne.ii) then
        ip1 = i+1
     else
        ip1 = 1  !global periodic region, also works for closed domains since ip(ii,:)=0
     endif
     if (ip(i,j).eq.1) then 
     ! --- kinetic energy / mass (m**2/s**2)
            ke(i,j,:) = 0.5*((0.5*(    u(i,j,:) +     u(ip1,j,:)))**2 + &
                             (0.5*(    v(i,j,:) +     v(i,j+1,:)))**2  )
         kemix(i,j)   = 0.5*((0.5*( umix(i,j)   +  umix(ip1,j)  ))**2 + &
                             (0.5*( vmix(i,j)   +  vmix(i,j+1)  ))**2  )
        kebaro(i,j)   = 0.5*((0.5*(ubaro(i,j)   + ubaro(ip1,j)  ))**2 + &
                             (0.5*(vbaro(i,j)   + vbaro(i,j+1)  ))**2  )
     endif
  enddo
  enddo
! --- arctic patch, also works for closed domains since ip(:,jj)=0
  do i= 1,ii
     ia = ii-mod(i-1,ii)
     if (ip(i,jj).eq.1) then
            ke(i,jj,:) =     ke(ia,jj-1,:)
         kemix(i,jj)   =  kemix(ia,jj-1)
        kebaro(i,jj)   = kebaro(ia,jj-1)
     endif
  enddo !i

  ! write(6,*) 'mean_velocity -   ke = ',     ke(54, 1,1)

  end subroutine mean_velocity

end module mod_mean

