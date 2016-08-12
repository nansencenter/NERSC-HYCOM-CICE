!!****f* HYCOM_2.1.03/EVP_MPI/hycomtoevp
!!
!! NAME
!!    hycomtoevp - Transfer needed hycom variables from hycom to the evp model
!!
!! SYNOPSIS
!!    subroutine hycomtoevp(AVEu,AVEv,m,n)
!!
!!
!! DESCRIPTION
!!    This routine transfers variables from HYCOM to the EVP model.
!!    Variables include atmospheric winds, ocean currens, ssh gradient,
!!    upper ocean layer thickness, ice concentration and average ice
!!    thickness.
!!
!! INPUTS   
!!    AVEu, AVEv      Time averaged surface velocities of ocean
!!    m,n             Old and new time step indices in HYCOM
!!
!! SIDE EFFECTS
!!    No changes are made to ocean variables, except a tiling step.
!!
!!    Routine sets up the following EVP variables (from mod_evp):
!!      uair
!!      vair
!!      uocn
!!      vocn
!!      ss_tltx
!!      ss_tlty
!!    Routine sets up the following EVP variables (from mod_ice_common), 
!!    if the CPP flag ICE is set:
!!      aice
!!      vice
!!
!! WARNINGS
!!    EVP_MPI routines are still in testing
!!
!!
!! PARAMETERS
!!    None
!!
!! AUTHOR
!!    Knut Arild Liseter
!!
!! CREATION DATE
!!    Jan 16th 2006
!!
!! HISTORY
!!
!!    Moved into CVS - Feb 2007
!!    Parallel testing/debugging, Mar 2007
!!    12.03.2007 - Added OMP directives
!! 
!! TODO
!!    Check transfer to/from hycom and mask setup
!!
!!
!! SOURCE
!!


      subroutine hycomtoevp(AVEu,AVEv,m,n)
      use mod_xc
      use mod_za
      use mod_evp   , only : uocn, vocn, uair, vair,
     &   umask, dxu, dyu, imargin, ss_tltx, ss_tlty,
     &   aice, vice
#if defined (ICESTATE) && defined (ICE)
#error CPP error - both ICE and ICESTATE is defined
#elif defined (ICESTATE)
      use mod_icestate  , only : icestate, nthick
#elif defined (ICE)
      use mod_common_ice  , only : ficem, hicem
#endif
      use mod_forcing_nersc
      implicit none
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), intent(in) :: 
     &   AVEu,AVEv
      integer, intent(in) :: m,n
c
      integer :: i,j,hk
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     &   ssh, dp1, tuwind, tvwind
      real :: xmin,xmax
      include 'common_blocks.h'

c --- Tile arrays for evp input
      call xctilr(montg1( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_ps)
      call xctilr(dp    ( 1-nbdy,1-nbdy,1,m),1,   1, 6,6, halo_ps)
      call xctilr(pbavg ( 1-nbdy,1-nbdy,  m),1,   1, 6,6, halo_ps)
c
c --- Get tiled hycom data
      imargin=nbdy
c$OMP PARALLEL DO PRIVATE(j,i)
c$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
         if (iu(i,j)==1) then
            tuwind(i,j)= uwind(i,j,l0)*w0+uwind(i,j,l1)*w1
     &                  +uwind(i,j,l2)*w2+uwind(i,j,l3)*w3
         else
            tuwind(i,j)=0.
         end if
         if (iv(i,j)==1) then
            tvwind(i,j)= vwind(i,j,l0)*w0+vwind(i,j,l1)*w1
     &                  +vwind(i,j,l2)*w2+vwind(i,j,l3)*w3
         else
            tvwind(i,j)=0.
         end if

         ssh(i,j) = (montg1(i,j)+thref*pbavg(i,j,m))/onem
         dp1(i,j) = (dp(i,j,1,m)/onem)
      end do
      end do
c$OMP END PARALLEL DO 
c
c --- Ocean and wind velocities in EVP u-point
      uocn=0.
      uair=0.
      vocn=0.
      vair=0.
      imargin=nbdy-1
c$OMP PARALLEL DO PRIVATE(j,i)
c$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
         if (umask(i,j)) then
            uocn(i,j) = 0.5*(  AVEu(i+1,j  ) + AVEu(i+1,j+1))
            vocn(i,j) = 0.5*(  AVEv(i  ,j+1) + AVEv(i+1,j+1))

            if (windf) then
               uair(i,j) = (tuwind(i+1,j  )*iu(i+1,j  ) +
     &                      tuwind(i+1,j+1)*iu(i+1,j+1)) /
     &                      max(1,iu(i+1,j  )+iu(i+1,j+1))
               vair(i,j) = (tvwind(i  ,j+1)*iu(i  ,j+1) +
     &                      tvwind(i+1,j+1)*iu(i+1,j+1)) /
     &                      max(1,iv(i  ,j+1)+iv(i+1,j+1))
            else
               uair(i,j) = 0.0
               vair(i,j) = 0.0
            end if
         else
            uocn(i,j)=0.
            vocn(i,j)=0.
            uair(i,j)=0.
            vair(i,j)=0.
         end if
      end do
      end do
!$OMP END PARALLEL DO 
      call xctilr(uocn( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_qv)
      call xctilr(vocn( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_qv)
      call xctilr(uair( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_qv)
      call xctilr(vair( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_qv)

cdiag call zaiopf('tuwnd.a','replace',117)
cdiag call zaiowr(tuwind,ip,.false.,xmin,xmax,117,.true.)
cdiag call zaiowr(tvwind,ip,.false.,xmin,xmax,117,.true.)
cdiag call zaiowr(uair  ,ip,.false.,xmin,xmax,117,.true.)
cdiag call zaiowr(vair  ,ip,.false.,xmin,xmax,117,.true.)
cdiag call zaiowr(tmpmsk,ip,.false.,xmin,xmax,117,.true.)
cdiag call zaiocl(117)


c --- sea surface height in u-points
      imargin=nbdy-1
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin

c ---   dxu, dyu not defined outside of physical dom
c ---   along model boundaries (which are not periodic)
        if (iu(i,j)==1 .and. iv(i,j)==1) then
          ss_tltx(i,j)=  0.5*(
     .          (ssh(i+1,j  ) - ssh(i  ,j  )) + 
     .          (ssh(i+1,j+1) - ssh(i  ,j+1)) 
     .       ) / dxu(i,j)
c
          ss_tlty(i,j)=  0.5*(
     .           (ssh(i+1,j+1) - ssh(i+1,j  )) + 
     .           (ssh(i  ,j+1) - ssh(i  ,j  )) 
     .       ) / dyu(i,j)
        else
           ss_tltx(i,j)=0.
           ss_tlty(i,j)=0.
        end if
      end do
      end do
!$OMP END PARALLEL DO 
      call xctilr(ss_tltx( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_qs)
      call xctilr(ss_tlty( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_qs)


#if defined (ICE)
c --- Ice concentration and volume from standard 1-category ice model
      call xctilr(ficem ( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_ps)
      call xctilr(hicem ( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_ps)
      imargin=nbdy
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
         aice(i,j) = ficem(i,j)
         vice(i,j) = ficem(i,j)*hicem(i,j)
      end do
      end do
!$OMP END PARALLEL DO 
      call xctilr(aice( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_ps)
      call xctilr(vice( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_ps)
c
#elif defined (ICESTATE)
c --- Ice concentration and volume from icestate multi-category ice model
      imargin=nbdy
      aice=0.
      vice=0.
      do hk=1,nthick
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
         aice(i,j) = aice(i,j)+icestate(i,j)%ice(hk)%fice
         vice(i,j) = vice(i,j)+icestate(i,j)%ice(hk)%fice*
     &                         icestate(i,j)%ice(hk)%hice
      end do
      end do
!$OMP END PARALLEL DO 
      end do
      call xctilr(aice( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_ps)
      call xctilr(vice( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_ps)
#endif
       end subroutine hycomtoevp
!!******
