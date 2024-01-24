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
      module mod_seaice_draft
      !TODO: Handle situautions where ice thickness exceeds draft...
      use mod_xc
      implicit none
      real,external :: pbudel_seaice, pbvdel_seaice
#if defined(RELO)
      real, save, allocatable, dimension(:,:) :: &
         draft_u,draft_u_thk,draft_v,draft_v_thk
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
              draft_u,draft_v,draft_u_thk,draft_v_thk
#endif

      contains
      subroutine init_seaice_draft()
      implicit none
#if defined(RELO)
        allocate( &
                draft_u    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                draft_v    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                draft_u_thk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                draft_v_thk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
        call mem_stat_add( 4*(idm+2*nbdy)*(jdm+2*nbdy) )
        draft_u=r_init
        draft_v=r_init
        draft_u_thk=r_init
        draft_v_thk=r_init
#endif
      end subroutine
!
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
!        
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
         end if
!
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
         end if
!
      end do
      end do
      call xctilr( draft_u    ,1,1, nbdy,nbdy, halo_uv)
      call xctilr( draft_u_thk,1,1, nbdy,nbdy, halo_uv)
      call xctilr( draft_v    ,1,1, nbdy,nbdy, halo_vv)
      call xctilr( draft_v_thk,1,1, nbdy,nbdy, halo_vv)
      end subroutine
!---- Include effect of moving ice. Ice divergence which differs from
! --- ocean divergence averaged over ice draft will now contribute to
! --- the evolution of the barotropic pressure.
!
      real function pbudel_seaice(i,j,n)
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
      integer, intent(in) :: i,j,n
! 
! --- Sea ice drift in u-point
      pbudel_seaice = .5* ( &
        + (si_u(i+1,j)+si_u(i,j))*scuy(i+1,j)*draft_u_thk(i+1,j) &
        - (si_u(i-1,j)+si_u(i,j))*scuy(i  ,j)*draft_u_thk(i  ,j) &
           )
!
! --- Subtract Total velocity over draft
      pbudel_seaice  = pbudel_seaice   - ( &
        + (draft_u(i+1,j)+ubavg(i+1,j,n))*scuy(i+1,j)*draft_u_thk(i+1,j) &
        - (draft_u(i  ,j)+ubavg(i  ,j,n))*scuy(i  ,j)*draft_u_thk(i  ,j) &
         )
      end function
!
!
      real function pbvdel_seaice(i,j,n)
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
      integer, intent(in) :: i,j,n
! 
! --- Sea ice drift in v-point
      pbvdel_seaice = .5* ( &
        + (si_v(i,j+1)+si_v(i,j))*scvx(i,j+1)*draft_v_thk(i,j+1) &
        - (si_v(i,j-1)+si_v(i,j))*scvx(i  ,j)*draft_v_thk(i  ,j) &
           )
!
! --- Subtract Total velocity over draft
      pbvdel_seaice  = pbvdel_seaice - ( &
        +(draft_v(i,j+1)+vbavg(i,j+1,n))*scvx(i,j+1)*draft_v_thk(i,j+1) &
        -(draft_v(i,j  )+vbavg(i,j  ,n))*scvx(i,j  )*draft_v_thk(i  ,j) &
        )
      end function
! 
      end module mod_seaice_draft

