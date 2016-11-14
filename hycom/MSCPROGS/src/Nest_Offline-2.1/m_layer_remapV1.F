      module m_layer_remapV1
      contains
      subroutine layer_remapV1(oldint,oldtem,oldsal,kdm, 
     &                       targetdens,newint,newdp0,newkdm)
      implicit none

      integer, intent(in)                     :: kdm, newkdm
      real,    intent(in) , dimension(kdm   ) :: oldint,oldtem,oldsal
      real,    intent(in) , dimension(newkdm) :: targetdens,newdp0
      real,    intent(out), dimension(newkdm) :: newint

      integer :: k,kold,klist
      real, dimension(kdm) :: olddens

      real,parameter :: qthref=1e-3
      include 'stmt_fns_SIGMA0.h'


      do kold=1,kdm
         olddens(kold)=sig(oldtem(kold),oldsal(kold))
      end do

      do k=1,newkdm


         ! Go through old layers, find densities which are lower than 
         ! target density
         klist=0
         do kold=1,kdm
            if (olddens(kold)<=targetdens(k)) then
               klist=kold
            end if
         end do

         ! Every layer above klist is lighter than targetdensity. 
         ! This becomes the new interface 
         if (klist>0) then
            newint(k)=oldint(klist)

         ! Couldnt find a lighter layer. The new interface must still obey the 
         ! depth restrictions, though
         else 

            if (k==1) then
               newint(k)= newdp0(k)
            else
               newint(k)= newint(k-1)+newdp0(k)
            end if

         end if

         ! Check that new layer interface is below last layer
         if (k>1) then
            newint(k)=max(newint(k-1)+newdp0(k),newint(k))
         else
            newint(k)=max(newdp0(k),newint(k))
         end if

         ! Check that we are not below sea floor
         newint(k)=min(newint(k),oldint(kdm))
      end do


      end subroutine layer_remapV1
      end module m_layer_remapV1


            

