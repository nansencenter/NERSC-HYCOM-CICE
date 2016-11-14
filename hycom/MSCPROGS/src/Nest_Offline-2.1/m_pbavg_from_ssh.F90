module m_pbavg_from_ssh

contains


   subroutine pbavg_from_ssh(saln,temp,dp,ssh,psikk,thkk,depths,pbavg,idm,jdm,kdm)
   use mod_sigma
   use mod_za, only : zaiopf, zaiocl, zaiowr
   implicit none

   integer, intent(in) :: idm,jdm,kdm
   real, intent(out), dimension(idm,jdm) ::  pbavg
   real, intent( in), dimension(idm,jdm) ::  depths,ssh, psikk, thkk
   real, intent( in), dimension(idm,jdm,kdm) ::  saln,temp, dp


   real, parameter :: onem  = 9806.
   real, parameter :: undef = -1e14

   real, dimension(idm,jdm,kdm) :: thstar
   real, dimension(idm,jdm,kdm+1) :: p

   integer :: i,j,k,ia,ja
   real, dimension(idm,jdm)     :: montg,oneta,th3d, &
      montgkk,ssh2, tmp
   logical, dimension(idm,jdm) :: gmsk
   integer, dimension(idm,jdm) :: idummy

   real :: rvar, thbase, amin, amax
   integer :: kapflg, thflag
   logical :: tbaric
   character(len=6)  :: cvarin

   logical :: ex
   integer, parameter :: lp=6

   !include 'stmt_funcs.h'

   inquire(file='blkdat.input',exist=ex)
   if (.not.ex) then
      write(lp,*) 'ERROR: blkdat.input does not exist'
      stop '(pbavg_from_state:blkdat.input parse)'
   endif
   open(10,file='blkdat.input')
   cvarin=''
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)

   do while (cvarin/='kapflg')
      read(10,*) rvar,cvarin
   end do
   if (cvarin/='kapflg') then
      write(lp,*) 'ERROR: Could not get kapflg'
      stop '(pbavg_from_state:blkdat.input parse)'
   else
      kapflg=nint(rvar)
   end if

   do while (cvarin/='thflag')
      read(10,*) rvar,cvarin
   end do
   if (cvarin/='thflag') then
      write(lp,*) 'ERROR: Could not get thbase'
      stop '(pbavg_from_state:blkdat.input parse)'
   else
      thflag=nint(rvar)
   end if

   do while (cvarin/='thbase')
      read(10,*) thbase,cvarin
   end do
   if (cvarin/='thbase') then
      write(lp,*) 'ERROR: Could not get thbase'
      stop '(pbavg_from_state:blkdat.input parse)'
   end if

   tbaric=kapflg==thflag
   print *,'blkdat.input -- kapflg:',kapflg
   print *,'blkdat.input -- thflag:',thflag
   print *,'blkdat.input -- thbase:',thbase
   print *,'blkdat.input -- tbaric:',tbaric


   gmsk=(depths>.1 .and. depths<1e29) 

   ! Get restart fields psikk and thkk
   !print *,rstbase

   ! cumulative pressure
   p(:,:,:)=0.
   pbavg=0.
   do k=1,kdm

      ! use upper interface pressure in converting sigma to sigma-star.
      ! this is to avoid density variations in layers intersected by bottom
      ! NB - KAL - thstar/th3d are densities with thbase = 0.
      do j=1,jdm
      do i=1,idm

         ia=mod(idm+i-2,idm)+1
         ja=max(1,j-1)
         if (gmsk(i,j)) then
            if (thflag==0) then
               th3d(i,j) = sig0(temp(i,j,k),saln(i,j,k))
            else if (thflag==2) then
               th3d(i,j) = sig2(temp(i,j,k),saln(i,j,k))
            else if (thflag==4) then
               th3d(i,j) = sig4(temp(i,j,k),saln(i,j,k))
            else
               print *,'Unknown thflag : ',thflag
               stop '(m_bavg_flds_from_ave)'
            end if
            if (tbaric) then
               if (thflag==0) then
                  thstar(i,j,k)=th3d(i,j)+kappaf0(temp(i,j,k), &
                                         saln(i,j,k),p(i,j,k))
               elseif (thflag==2) then
                  thstar(i,j,k)=th3d(i,j)+kappaf2(temp(i,j,k), &
                                         saln(i,j,k),p(i,j,k))
               elseif (thflag==4) then
                  thstar(i,j,k)=th3d(i,j)+kappaf4(temp(i,j,k), &
                                         saln(i,j,k),p(i,j,k))
               else
                  print *,'Unknown thflag : ',thflag
                  stop '(m_bavg_flds_from_ave)'
               end if
            else
               thstar(i,j,k)=th3d(i,j)
            end if
            thstar(i,j,k)=thstar(i,j,k)-thbase
            !p(i,j,k+1)=p(i,j,k)+dp(i,j,k)*onem
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
            !if (i==200.and.j==300) print *,k,thstar(i,j,k), p(i,j,k+1)/onem
         end if

      enddo
      enddo
         !print *,'dp   max/min:',maxval(dp  ),minval(dp  )
         !print *,'p    max/min:',maxval(p   ),minval(p   )
         !print *,'temp max/min:',maxval(temp),minval(temp)
         !print *,'saln max/min:',maxval(saln),minval(saln)
         !print *,'utot max/min:',maxval(utot),minval(utot)
         !print *,'vtot max/min:',maxval(vtot),minval(vtot)
         !print *,'ssh  max/min:',maxval(ssh ),minval(ssh )
   end do



   do j=1,jdm
   do i=1,idm
      if (gmsk(i,j)) then
      ! store (1+eta) (= p_total/p_prime) in -oneta-
      ! KAL -- oneta is UNKNOWN !
      !oneta(i,j)=1.+pbavg(i,j)/p(i,j,kdm+1)


      ! m_prime in lowest layer:
      !montgkk(i,j)=psikk(i,j)+(p(i,j,kdm+1)*(thkk(i,j)+thbase-thstar(i,j,kdm)) &
      !            -pbavg(i,j)*(thstar(i,j,kdm)))*thref**2
      montgkk(i,j)=psikk(i,j)+( &
                                p(i,j,kdm+1)*(thkk(i,j)-thstar(i,j,kdm)) &
                               )*thref**2
   end if
   enddo
   enddo


   ! m_prime in remaining layers:
   montg=0.
   do k=kdm-1,1,-1
      do j=1,jdm
      do i=1,idm
      if (gmsk(i,j)) then
         !montg(i,j,k)=montg(i,j,k+1)+p(i,j,k+1)*oneta(i,j) &
         !        *(thstar(i,j,k+1)-thstar(i,j,k))*thref**2

         ! KAL - oneta is unknown -- "divide" by thref and an unknown oneta
         !montg(i,j,k)=montg(i,j,k+1)+p(i,j,k+1)*oneta(i,j) &
         !        *(thstar(i,j,k+1)-thstar(i,j,k))*thref**2
         ! Keeps oneta out of the game
         montg(i,j)=montg(i,j)+p(i,j,k+1)*(thstar(i,j,k+1)-thstar(i,j,k))
      end if
      enddo
      enddo
   enddo


!   ! SSH
!   do  j=1,jdm
!   do  i=1,idm
!      if (gmsk(i,j))then
!      ssh(i,j)=(montg(i,j,1)/thref+pbavg(i,j))/onem 
!   end if
!   end do
!   end do

   ! pbavg
   ssh2=0.
   do  j=1,jdm
   do  i=1,idm
      if (gmsk(i,j))then
         !ssh(i,j)=(montg(i,j,1)/thref+pbavg(i,j))/onem 
         pbavg(i,j)=ssh(i,j)*onem - montgkk(i,j)/thref - montg(i,j)*thref
         pbavg(i,j)=pbavg(i,j)/(1. + thref * montg(i,j)/p(i,j,kdm+1) - thstar(i,j,kdm)*thref )

         ! oneta
         oneta(i,j)=1.+pbavg(i,j)/p(i,j,kdm+1)

         ! ssh2 - should be same as ssh
         ssh2(i,j)=( montgkk(i,j) - pbavg(i,j)*thstar(i,j,kdm)*thref**2  &
                    +oneta(i,j)*montg(i,j)*thref**2 )  / thref
         ssh2(i,j)=(ssh2(i,j)+pbavg(i,j))/onem
         
      end if
   end do
   end do

   !print *,'M1/p     ',maxval(montg/p(:,:,kdm+1),mask=p(:,:,kdm+1)>1.)
   !print *,'Mkk-M1   ',maxval(abs(montgkk - montg)*thref**2)
   !print *,'Mkk/thref',maxval(abs(montgkk/thref))
   !print *,'M1 *thref',maxval(abs(montg*thref))
   !print *,'ssh*onem ',maxval(abs(ssh*onem))
   !print *,'ssh2*onem ',maxval(abs(ssh2*onem))
   !print *,'thkk     ',maxval(abs(thkk))
   !print *,'thref*montg/p ',maxval(abs(thref*montg/p(:,:,kdm+1)),mask=p(:,:,kdm+1)>1.)

   ! Open .a file and dump new fields
   idummy=0
   call zaiopf('tstbavg2.a','replace',456)
   call zaiowr(pbavg  ,idummy,.false.,amin,amax,456,.false.)
   tmp=ssh
   call zaiowr(tmp    ,idummy,.false.,amin,amax,456,.false.)
   !call zaiowr(ssh2   ,idummy,.false.,amin,amax,456,.false.)
   tmp=psikk
   call zaiowr(tmp    ,idummy,.false.,amin,amax,456,.false.)
   tmp=thkk
   print *,'thkk range : ', minval(thkk),maxval(thkk)
   call zaiowr(tmp    ,idummy,.false.,amin,amax,456,.false.)

   call zaiocl(456)


   print *,'************************************************'
   print *,'************************************************'
   print *,'Barotropic calculations : '
   print *,'ssh(known) max/min:',maxval(ssh  ),minval(ssh  )
   print *,'pbavg/onem max/min:',maxval(pbavg/9806.),minval(pbavg/9806.)
   print *,'Diagnostics dumped to tstbavg2.a'
   print *,'************************************************'
   print *,'************************************************'
   print *
   print *


   end subroutine pbavg_from_ssh
end module m_pbavg_from_ssh



