module m_layer_mixV1

contains

   subroutine  layer_mixV1(oldint,oldvar,oldkdm,newint,newvar,newkdm)
   implicit none

   integer, intent(in) :: oldkdm, newkdm
   real, dimension(oldkdm), intent(in ) :: oldvar, oldint
   real, dimension(newkdm), intent(in ) :: newint
   real, dimension(newkdm), intent(out) :: newvar

   integer :: ko,kn,lastfilled
   real, dimension(oldkdm) :: oldfrac
   real :: olddn,newdn,oldup,newup,upper,lower,newdp

   lastfilled=0
   do kn=1,newkdm

      if (kn==1) then
         newup=0.
      else
         newup=newint(kn-1)
      end if
      newdn=newint(kn)
      newdp=newdn-newup

      ! Find fraction of old layer k within interfaces
      ! of new layer kn
      newvar(kn)=0.
      do ko=1,oldkdm

         if (ko==1) then
            oldup=0.
         else
            oldup=oldint(ko-1)
         end if
         olddn=oldint(ko)


         upper=max(newup,oldup)
         lower=min(newdn,olddn)


         ! Fraction of old layer within new layer
         oldfrac(ko)=max(0.,lower-upper)/max(.001,newdp)
         !print *,upper,lower,oldfrac(ko),newdp

         newvar(kn)= newvar(kn)+oldfrac(ko)*oldvar(ko)
      end do

      if (newdp>1.) lastfilled=kn

      ! Consistency check
      if (newdp>1.) then
         if (abs(1.-sum(oldfrac))>.01) then
            print *,'Layers do not add up to 1...'
            print *,'kn         =',kn
            print *,'sum oldfrac=',sum(oldfrac)
            print *,'newdp      =',newdp
            print *,'oldint     =',oldint
            print *,'newint     =',newint
            stop '(layer_mixV1)'
         end if

      ! Fill with last mass-filled layer values
      else if (lastfilled>0) then
         newvar(kn)=newvar(lastfilled)
      end if
   end do

   end subroutine layer_mixV1
end module m_layer_mixV1
