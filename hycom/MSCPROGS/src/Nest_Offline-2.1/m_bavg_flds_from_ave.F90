module m_bavg_flds_from_ave

! TODO - meant to be called for nersc_daily, nersc_weekly

contains


   subroutine bavg_flds_from_ave(ubavg,vbavg,pbavg,depths,rstbase,avebase, &
      idm,jdm,kdm)
   use mod_sigma
   use mod_za, only : zaiopf, zaiocl, zaiowr
   use mod_hycomfile_io
   implicit none

   integer, intent(in) :: idm,jdm,kdm
   real, intent(out), dimension(idm,jdm) ::  ubavg, vbavg, pbavg
   real, intent( in), dimension(idm,jdm) ::  depths
   character(len=*), intent(in) :: rstbase,avebase

   real, dimension(idm,jdm,kdm) :: thstar
   real, dimension(idm,jdm,kdm+1) :: p

   integer :: i,j,k,ia,ja
   real, dimension(idm,jdm)     :: montg,oneta,temp,saln,th3d,utot,vtot, &
      montgkk,psikk,thkk, dp, ssh, ssh2,prsold,prs
   logical, dimension(idm,jdm) :: gmsk
   integer, dimension(idm,jdm) :: idummy

   real :: rvar, thbase, amin, amax
   integer :: kapflg, thflag
   logical :: tbaric
   character(len=6)  :: cvarin

   logical :: ex
   integer, parameter :: lp=6
   type(hycomfile) :: hfile
   character(len=80) :: ftype

   !include 'stmt_funcs.h'

   print *,'Bottom properties retrieved from: ',trim(rstbase)
   print *,'SSH/temp/saln/etc retrieved from: ',trim(avebase)

   

   inquire(file='blkdat.input',exist=ex)
   if (.not.ex) then
      write(lp,*) 'ERROR: blkdat.input does not exist'
      stop '(ssh_from_state:blkdat.input parse)'
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
      write(lp,*) 'ERROR: Could not get jdm'
      stop '(ssh_from_state:blkdat.input parse)'
   else
      kapflg=nint(rvar)
   end if

   do while (cvarin/='thflag')
      read(10,*) rvar,cvarin
   end do
   if (cvarin/='thflag') then
      write(lp,*) 'ERROR: Could not get thbase'
      stop '(ssh_from_state:blkdat.input parse)'
   else
      thflag=nint(rvar)
   end if

   do while (cvarin/='thbase')
      read(10,*) thbase,cvarin
   end do
   if (cvarin/='thbase') then
      write(lp,*) 'ERROR: Could not get thbase'
      stop '(ssh_from_state:blkdat.input parse)'
   end if

   tbaric=kapflg==thflag
   print *,'blkdat.input -- kapflg:',kapflg
   print *,'blkdat.input -- thflag:',thflag
   print *,'blkdat.input -- thbase:',thbase
   print *,'blkdat.input -- tbaric:',tbaric


   gmsk=(depths>.1 .and. depths<0.5*huge) 

   ! Get restart fields psikk and thkk
   call initHF(hfile,trim(rstbase)//'.a','restart')
   call HFReadField(hfile,psikk,idm,jdm,'psikk   ',0,1)
   call HFReadField(hfile,thkk ,idm,jdm,'thkk    ',0,1)

   ! Get daily (weekly average?) field ssh
   ftype=getfiletype(trim(avebase)//'.a')
   call initHF(hfile,trim(avebase)//'.a',trim(ftype))
   call HFReadField(hfile,ssh,idm,jdm,'ssh     ',0,1)
   where (.not. gmsk) ssh=0.




   ! cumulative pressure
   p(:,:,:)=0.
   ubavg=0.
   vbavg=0.
   pbavg=0.
   do k=1,kdm

      ! Read layer thickness
      !call HFReadDPField(hfile,dp,idm,jdm,k,1)
      !dp = dp/onem
      call HFReadDPField_m(hfile,dp,idm,jdm,k,1) ! Returns in meters

      ! Read remaining variables (this wont work for archv !)
      call HFReadField(hfile,saln,idm,jdm,'saln    ',k,1)
      call HFReadField(hfile,temp,idm,jdm,'temp    ',k,1)
      call HFReadField(hfile,utot,idm,jdm,'utot    ',k,1)
      call HFReadField(hfile,vtot,idm,jdm,'vtot    ',k,1)


      ! use upper interface pressure in converting sigma to sigma-star.
      ! this is to avoid density variations in layers intersected by bottom
      ! NB - KAL - thstar/th3d are densities with thbase = 0.
      do j=1,jdm
      do i=1,idm

         ia=mod(idm+i-2,idm)+1
         ja=max(1,j-1)
         if (gmsk(i,j)) then
            if (thflag==0) then
               th3d(i,j) = sig0(temp(i,j),saln(i,j))
            else if (thflag==2) then
               th3d(i,j) = sig2(temp(i,j),saln(i,j))
            else if (thflag==4) then
               th3d(i,j) = sig4(temp(i,j),saln(i,j))
            else
               print *,'Unknown thflag : ',thflag
               stop '(m_bavg_flds_from_ave)'
            end if
            if (tbaric) then
               if (thflag==0) then
                  thstar(i,j,k)=th3d(i,j)+kappaf0(temp(i,j), &
                                         saln(i,j),p(i,j,k))
               elseif (thflag==2) then
                  thstar(i,j,k)=th3d(i,j)+kappaf2(temp(i,j), &
                                         saln(i,j),p(i,j,k))
               elseif (thflag==4) then
                  thstar(i,j,k)=th3d(i,j)+kappaf4(temp(i,j), &
                                         saln(i,j),p(i,j,k))
               else
                  print *,'Unknown thflag : ',thflag
                  stop '(m_bavg_flds_from_ave)'
               end if
            else
               thstar(i,j,k)=th3d(i,j)
            end if
            thstar(i,j,k)=thstar(i,j,k)-thbase
            p(i,j,k+1)=p(i,j,k)+dp(i,j)*onem
            !if (i==200.and.j==300) print *,k,thstar(i,j,k), p(i,j,k+1)/onem
         end if

      ! TODO -- fix these
      ! nb ave files interp u to p points ?
      !if (gmsk(i,j).and.gmsk(ia,j)) then
      if (gmsk(i,j)) then
         ubavg(i,j)=ubavg(i,j)+utot(i,j)*dp(i,j)*onem
      end if
      ! nb ave files interp u to p points ? 
      !if (gmsk(i,j).and.gmsk(i,ja)) then
      if (gmsk(i,j)) then
         vbavg(i,j)=vbavg(i,j)+vtot(i,j)*dp(i,j)*onem
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

   ! Ubavg/vbavg
   do  j=1,jdm
   do  i=1,idm
   if (gmsk(i,j))then
      ubavg(i,j)=ubavg(i,j)/p(i,j,kdm+1)
      vbavg(i,j)=vbavg(i,j)/p(i,j,kdm+1)
   end if
   end do
   end do

   ! Open .a file and dump new fields
   idummy=0
   call zaiopf('tstbavg.a','replace',456)
   call zaiowr(ubavg  ,idummy,.false.,amin,amax,456,.false.)
   call zaiowr(vbavg  ,idummy,.false.,amin,amax,456,.false.)
   call zaiowr(pbavg  ,idummy,.false.,amin,amax,456,.false.)
   call zaiowr(ssh    ,idummy,.false.,amin,amax,456,.false.)
   call zaiowr(ssh2   ,idummy,.false.,amin,amax,456,.false.)
   call zaiocl(456)


   print *,'************************************************'
   print *,'************************************************'
   print *,'Barotropic calculations : '
   print *,'ssh(known) max/min:',maxval(ssh  ),minval(ssh  )
   print *,'ubavg      max/min:',maxval(ubavg),minval(ubavg)
   print *,'vbavg      max/min:',maxval(vbavg),minval(vbavg)
   print *,'pbavg/onem max/min:',maxval(pbavg/9806.),minval(pbavg/9806.)
   print *,'Diagnostics dumped to tstbavg.a'
   print *,'************************************************'
   print *,'************************************************'
   print *
   print *


   end subroutine bavg_flds_from_ave
end module m_bavg_flds_from_ave



