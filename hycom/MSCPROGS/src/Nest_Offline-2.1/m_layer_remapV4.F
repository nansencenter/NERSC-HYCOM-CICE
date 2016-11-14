      module m_layer_remapV4
      contains

      ! KAL - This remapping routine is based on mixing over several
      ! KAL - layers to get target densities. Linear profiles are
      ! KAL - set up. Mixing and nonlinearities of the state function
      ! KAL - are not taken into account ....

      subroutine layer_remapV4(oldint,olddens,oldkdm, 
     &                       targetdens,newint,newdp0,newkdm,thflag,
     &                       depth,tstpoint,lfatal)
      use mod_sigma
      implicit none

      integer, intent(in)                     :: oldkdm, newkdm,thflag
      real   , intent(in)                     :: depth
      real,    intent(in) , dimension(oldkdm) :: oldint,olddens
      real,    intent(in) , dimension(newkdm) :: targetdens,newdp0
      real,    intent(out), dimension(newkdm) :: newint
      logical, intent(in) :: tstpoint
      logical, intent(out):: lfatal

      integer :: k,kold, kolast,match
      integer :: klist(newkdm,2)

      real :: newdp, oldup, sumdp, initdens, tmpdens, newup, oldpfrac
      integer :: iderr
      integer,save :: callcount=1
      integer, parameter :: modcallcount=10000
      character(len=5) :: ccall

      integer, parameter :: nop=234
      real, parameter :: onem=9806.


      real, dimension(1:oldkdm+1) :: tmpoldint
      real, dimension(1:oldkdm+1) :: idnso     ! old interface dens
      real, dimension(1:newkdm+1) :: tmpnewint
      real, dimension(1:newkdm+1) :: idns      ! new interface dens

      real,parameter :: tol=1e-3, ptol=0.01
      real,parameter :: tol2=1e-4

      integer :: lastk
      real :: pseek, frac_olddp, olddensdiff, oldpresdiff

      lfatal=.false.




! --- Set up target  densities  at layer interfaces
      idns(1)=targetdens(1)
      idns(newkdm+1)=targetdens(newkdm)
      do k=2,newkdm
         ! Linear 
         idns(k) = 0.5*(targetdens(k-1)+targetdens(k))
      end do

! --- Set up old  densities and work pressures at layer interfaces
      idnso(1)=olddens(1)
      idnso(oldkdm+1)=olddens(oldkdm)
      do k=2,oldkdm
         ! Linear 
         idnso(k) = 0.5*(olddens(k-1)+olddens(k))
      end do
      tmpoldint=0.
      do k=1,oldkdm
         ! Linear 
         tmpoldint(k+1) = oldint(k)
      end do

      !if (tstpoint) then
      !   print *
      !   print *,'newdp0         :',newdp0 
      !   print *
      !   print *,'olddens        :',olddens
      !   print *,'olddens   (tmp):',idnso
      !   print *,'oldint    (tmp):',tmpoldint
      !   print *
      !   print *,'targetdens(tmp):',idns
      !   print *
      !end if

! --- Find layers

      tmpnewint=0.
      do k=1,newkdm

         ! Initalize sought after interface (upper interface of current layer)
         pseek=tmpnewint(k)

         ! Loop old layers
         do kold=1,oldkdm

            
            ! Pressure and density difference betwen lower and upper interface of
            ! this layer
            olddensdiff =  (idnso(kold+1) - idnso(kold))
            oldpresdiff =  (tmpoldint(kold+1) - tmpoldint(kold))

            ! Fraction of old layer available at this point
            frac_olddp  = tmpoldint(kold+1)-tmpnewint(k)
!            print '(2i5,4f12.2)',k,kold,frac_olddp, tmpoldint(kold+1),
!     &              tmpnewint(k),newdp0(k)

            ! Proceed to next layer if fraction is non-negative
            if (frac_olddp>0.) then

               ! - Find location where old piecewise linear density  
               ! - is equal to lower interface target density

               ! use linear weighting relationship
               if (olddensdiff>tol .and. oldpresdiff > ptol) then  ! real density difference

                  ! Sought after pressure of lower target interface
                  pseek = (oldpresdiff)*idns(k+1)           - 
     &                    tmpoldint(kold+1)*idnso(kold  )  + 
     &                    tmpoldint(kold  )*idnso(kold+1)      
                  pseek = pseek / olddensdiff

                  !print *,'lin:',idns(k+1),pseek,olddensdiff,oldpresdiff


                  ! pseek is higher in w ccolumn than upper target interface - cant find layers
                  ! lighter than this target density -- exit loop
                  if (pseek < tmpnewint(k)) then
                     pseek=tmpnewint(k)
                     !exit kold loop at this point
                     goto 100

                  ! pseek is lower than lower interface of the old layer --
                  ! must continue seek in next layer
                  else if (pseek > tmpoldint(kold+1)) then
                     pseek=tmpoldint(kold+1)
                     !proceed to next kold layer
                     cycle

                  ! pseek is located within this layer 
                  else
                     !exit kold loop at this point 
                     goto 100

                  end if


               ! non-zero layer, but density is near-constant in layer
               else if (oldpresdiff>ptol .and.  olddensdiff<=tol ) then


                  !print *,'cns:',idns(k+1),pseek,olddensdiff,oldpresdiff

                  ! density of old layer higher than lower interface density
                  ! exit kold loop and adhere to minimum layer thickness
                  if (idnso(kold+1) > idns(k+1) ) then
                     !exit loop at this point, dont touch pseek
                     pseek=tmpnewint(k)

                  ! density of old layer lower than lower interface density
                  ! use entire layer, continue seek in next layer
                  else 
                     pseek = tmpoldint(kold+1)
                  end if

               ! old layer is zero thickness - proceed to next layer
               else

               end if


            else  ! frac_dp<=0 
               ! proceed to next layer
               cycle

            end if !frac_dp 

         end do ! kold
  100    continue

         ! Adhere to minimum thickness here
         !print *,'final pseek:',tmpnewint(k),pseek,newdp0(k)

         pseek=max(newdp0(k)+tmpnewint(k),pseek)

         ! Adhere to max depth
         pseek=min(pseek,tmpoldint(oldkdm+1))

         tmpnewint(k+1)=pseek

      enddo ! k


      ! Ensure lowest mass filled layer stretches to the bottom
      lastk=newkdm+4
      do k=1,newkdm
         if (tmpnewint(k+1)-tmpnewint(k)>ptol) then
            lastk = k
         end if
      end do
      lastk=min(newkdm,k)
      tmpnewint(lastk:newkdm+1)=tmpoldint(oldkdm+1)


      ! return this array
      newint=tmpnewint(2:newkdm+1)



      if ((mod(callcount,modcallcount)==0.and.oldint(oldkdm)>0).or.
     &    tstpoint) then
        write(ccall,'(i5.5)') callcount/modcallcount
        open(nop,file='oldprof'//ccall//'.dat',form='formatted',   
     &        status='replace')
        do k=1,oldkdm
           write(nop,*) oldint(k),olddens(k)
        end do
        close(nop)
        open(nop,file='newprof'//ccall//'.dat',form='formatted', 
     &       status='replace')
        do k=1,newkdm
           write(nop,*) newint(k),targetdens(k)
        end do
        close(nop)
      end if
      callcount=callcount+1

      end subroutine layer_remapV4
      end module m_layer_remapV4
