      module m_layer_remapV2
      contains

      ! KAL - This remapping routine is based on mixing over several
      ! KAL - layers to get target densities. Staircase profiles are
      ! KAL - assumed. Mixing and nonlinearities of the state function
      ! KAL - are not taken into account ....

      subroutine layer_remapV2(oldint,oldtem,oldsal,oldkdm, 
     &                       newdens,newint,newdp0,newkdm,thflag,
     &                       tstpoint)
      use mod_sigma
      implicit none

      integer, intent(in)                     :: oldkdm, newkdm,thflag
      real,    intent(in) , dimension(oldkdm) :: oldint,oldtem,oldsal
      real,    intent(in) , dimension(newkdm) :: newdens,newdp0
      real,    intent(out), dimension(newkdm) :: newint
      logical, intent(in) :: tstpoint

      integer :: k,kold, kolast,match
      integer :: klist(newkdm,2)
      real, dimension(oldkdm) :: olddens

      real :: newdp, oldup, sumdp, initdens, tmpdens, newup, oldpfrac
      integer :: iderr
      integer,save :: callcount=1
      integer, parameter :: modcallcount=1000
      character(len=5) :: ccall

      integer, parameter :: nop=234

      real,parameter :: tol=1e-2
      real,parameter :: tol2=1e-4


      ! Densities on old grid
      iderr=0
      do kold=1,oldkdm
         if (thflag==0) then
            olddens(kold)=sig0(oldtem(kold),oldsal(kold))
         elseif (thflag==2) then
            olddens(kold)=sig2(oldtem(kold),oldsal(kold))
         else
            print *,'Unknown thflag ',thflag
            call exit(1)
         end if
         !print *,thflag,olddens(kold)

         ! safety check
         if (kold>1) then
C           print *,kold,olddens(kold-1),olddens(kold),oldtem(kold),
C    &              oldsal(kold)
            if (olddens(kold-1)>olddens(kold)+tol) then
               if (iderr==0) then
                  iderr=kold
                  print *,kold,oldkdm
                  print *,kold+1,olddens(kold+1),
     &                    oldtem(kold+1), oldsal(kold+1), oldint(kold+1)
                  print *,kold,olddens(kold),
     &                    oldtem(kold), oldsal(kold), oldint(kold)
                  print *,kold-1,olddens(kold-1),
     &                    oldtem(kold-1), oldsal(kold-1), oldint(kold-1)
               end if
            else if (olddens(kold-1)>olddens(kold)) then
               olddens(kold)=olddens(kold-1)

            end if
         end if
      end do
      if (iderr>0) then
         print *,'dens error:',iderr
         stop 'layer_remapV2 - dens error !'
      end if

      if (tstpoint) then
         print *,'olddens   :',olddens
         print *,'targetdens:',newdens
      end if


      ! Set up some useful quantities
      klist=0
      do k=1,newkdm

         ! -- Go through old layers, find densities which are lower than 
         ! -- target density (1:klist)
         do kold=1,oldkdm
            if (olddens(kold)<=newdens(k)) then
               klist(k,1)=kold
            end if
         end do

         ! -- Go through old layers, find densities which are higher than 
         ! -- target density.
         do kold=oldkdm,1,-1
            if (olddens(kold)>=newdens(k)) then
               klist(k,2)=kold
            end if
         end do
      end do

      if (tstpoint) then
         print *,'klist1    :',klist(:,1)
         print *,'klist2    :',klist(:,2)
      end if

      ! Cycle new layers
      kolast=1
      do k=1,newkdm

         ! Last new layer upper interface
         if (k==1) then
            newup=0.;
         else
            newup=newint(k-1);
         end if
         !print *,'newup:',newup


         if (klist(k,1)==0)  then
            ! klist is zero - could not find layers lighter 
            ! than the target density -> layer outcrops at surface, dp
            ! should be zero in an isopycnic model, here we must consider
            ! layer thickness restrictions
            !newint(k)=newint(k)+newdp0(k);
            newint(k)=newup+newdp0(k);
            if (tstpoint) print *,'klist(k,1)=0 in :',k

         elseif (klist(k,1)<kolast)  then
            ! There are water masses lighter than this one, but they are
            ! above the current interface. We can not reach the
            ! target densty. The next iteration may succeed, however, so
            !newint(k)=newint(k)+newdp0(k);
            newint(k)=newup+newdp0(k);
            if (tstpoint) print *,'klist(k,1)<kolast in :',k


         elseif (klist(k,2)==0)  then
            ! klist2 is zero - no water masses are heavier than the target 
            ! density. There is no hope of reaching our target density.... This
            ! layer fills  the water column
            newint(k)=oldint(oldkdm);
            if (tstpoint) print *,'klist(k,2)=0 in :',k

         else
            ! There are elements heavier and lighter than the target density 
            ! in the currently "available" old water column - further processing needed...

            if (tstpoint) print *,'klist match in :',k

            tmpdens=0.
            initdens=0
            sumdp=0.
            match=0
            do kold=1,oldkdm

               ! Last old layer upper interface
               if (kold==1)  then
                  oldup=0.
               else
                  oldup=oldint(kold-1)
               end if
               !print *,'   oldup/dnw:',oldup,' ',oldint(kold)

         
               ! fraction of this old layer available for mixing
               oldpfrac=max(0.,oldint(kold) - max(newup,oldup));

               if (sumdp==0.)  then
                  sumdp=oldpfrac;
                  tmpdens=olddens(kold);
                  initdens=tmpdens;
               else
                  ! mixing a "newdp" layer of density olddens(kold) will result
                  ! in target density

                  ! Various checks
                  if (abs(newdens(k)-tmpdens)<tol)  then
                     !print *,'   match';
                     newdp=0.;
                     !if (match==0)
                     !   kolast=kold;
                     !end
                     match=1;

                  elseif (newdens(k)<tmpdens)  then
                     newdp=oldpfrac;
                  elseif (newdens(k)>olddens(kold))  then
                     newdp=oldpfrac;
                  else
                     newdp=sumdp*(newdens(k)-tmpdens)/
     &                     (olddens(kold)-newdens(k));
                     newdp=max(0.,min(newdp,oldint(kold)-oldup));
                  end  if

                  tmpdens=(olddens(kold)*newdp + tmpdens*sumdp)/
     &                    (sumdp+newdp);
                  sumdp=sumdp+newdp;
               end if
            end do ! kold

            newint(k)=newup+sumdp;
            if (match==0 .and. sumdp>1e-4) then  ! No match -- use dp0
               newint(k)=newint(k)-sumdp;
               newint(k)=newint(k)+newdp0(k);
            end if



         end if


         newint(k)=min(newint(k),oldint(oldkdm))
         if (k==newkdm) newint(k)=oldint(oldkdm)
      end  do
      !figure(3) ; clf; 
      !stairs(olddens,-oldint) ; hold on
      !stairs(newdens,-newint,'Color','r') ; hold on


      if (mod(callcount,modcallcount)==0.and.oldint(oldkdm)>0) then
         write(ccall,'(i5.5)') callcount/modcallcount
         open(nop,file='oldprof'//ccall//'.dat',form='formatted',   
     &         status='replace')
         do k=1,oldkdm
            write(nop,*) oldint(k),olddens(k)
         end do
         close(nop)
         open(nop,file='newprof'//ccall//'.dat',form='formatted', 
     &        status='replace')
         do k=1,newkdm
            write(nop,*) newint(k),newdens(k)
         end do
         close(nop)
      end if
      callcount=callcount+1




      end subroutine layer_remapV2
      end module m_layer_remapV2
