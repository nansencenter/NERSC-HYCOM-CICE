      module m_layer_remapV2
      contains

      ! KAL - This remapping routine is based finding layers, then
      ! umixing layers to get target densities

      subroutine layer_remapV2(oldint,oldtem,oldsal,kdm, 
     &                       targetdens,newint,newdp0,newkdm)
      implicit none

      integer, intent(in)                     :: kdm, newkdm
      real,    intent(in) , dimension(kdm   ) :: oldint,oldtem,oldsal
      real,    intent(in) , dimension(newkdm) :: targetdens,newdp0
      real,    intent(out), dimension(newkdm) :: newint

      integer :: k,kold
      integer :: klist(kdm,2)
      real, dimension(kdm) :: olddens

      real,parameter :: qthref=1e-3
      include '../stmt_fns_SIGMA0.h'


      ! Densities on old grid
      do kold=1,kdm
         olddens(kold)=sig(oldtem(kold),oldsal(kold))

         ! safety check
         if (kold>1) then
            if (olddens(kold-1)>=olddens(kold)) then
               stop '(layer_remapV2 - dens error !'
            end if
         end if
      end do

      ! Density limits - new grid
      do k=1,newkdm+1

         if(k==1) then
            denslim(k)=targetdensity(k)
         else
            denslim(k)=targetdensity(k)+targetdensity(k-1)
         end if

         if(k==newkdm+1) then
            denslim(k)=targetdensity(k)
         else
            denslim(k)=targetdensity(k)+targetdensity(k+1)
         end if

      end do

      ! Cycle new layers
      kolast=1
      do k=1,newkdm

         ! find old density interval within density limits of this one
         klist(:)=(/oldkdm+1,0/)
         do kold=1,oldkdm
            if (olddens(kold) > denslim(k  ) .and.
     &          olddens(kold) < denslim(k+1)  ) then

               klist(1) = min(kold,klist(1))
               klist(2) = max(kold,klist(2))
            end if
         end do


            
