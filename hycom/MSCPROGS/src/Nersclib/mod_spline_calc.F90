module mod_spline_calc
! Module is used to do vertical interpolation of a variable onto
! specified depth levels. There are routines for vertical interpolation
! of 1D fields amd 2D fields
! and then 


use mod_parameters
   integer, save :: ndeep
   real, save, allocatable :: deeps(:)
   character(len=20) ,save :: vimethod

contains

subroutine spline_calc_ini_fromfile(cmethod)
   implicit none
   character(len=*), intent(in) :: cmethod

   call get_depth_levels()

   if (trim(cmethod)=='spline' .or. &
       trim(cmethod)=='linear' .or. &
       trim(cmethod)=='staircase') then
      vimethod=trim(cmethod)
   else
      print *,'Vertical interpolation method '//trim(cmethod)//'unknown'
      print *,'(mod_spline_calc:spline_calc_ini_fromfile)'
      call exit(1)
   end if


end subroutine

subroutine spline_calc_ini_frominput(cmethod,deepin,ndeepin)
   implicit none
   character(len=*), intent(in) :: cmethod
   integer, intent(in) :: ndeepin
   real   , intent(in) :: deepin(ndeepin)
   allocate(deeps(ndeepin))
   deeps=deepin
   ndeep=ndeepin
   if (trim(cmethod)=='spline' .or. &
       trim(cmethod)=='linear' .or. &
       trim(cmethod)=='staircase') then
      vimethod=trim(cmethod)
   else
      print *,'Vertical interpolation method '//trim(cmethod)//'unknown'
      print *,'(mod_spline_calc:spline_calc_ini_frominput)'
      call exit(1)
   end if
end subroutine



subroutine spline_calc(fldin,depthin,nlons,nlats,mask,fldout,nzlev,kdm,deepsin)
implicit none
integer, intent(in) :: nlons,nlats,nzlev,kdm
real, dimension(nlons,nlats,kdm)  ,intent(in) :: fldin,depthin
real, dimension(nlons,nlats,nzlev),intent(out) :: fldout !,depthout

!KAL20151110 Allow interpolation depths to be set based on input
real, dimension(nzlev),intent(in),optional :: deepsin
!KAL20151110
logical, intent(in) :: mask(nlons,nlats)

integer :: ilon,jlat,sub_last,k,klevel,kzlevel,klevel2,km1
real, dimension(kdm) :: tmpfld,tmpd
real, dimension(kdm) :: subfld,subd
real, dimension(nzlev) :: sfld
!KAL20151110 Allow interpolation depths to be set based on input. used_deeps is
!KAL20151110 depth array used in routine (can be either module array or input
!KAL20151110 array). actual_deeps is used to account for negative depths (indicates height
!KAL20151110 from bottomn)
real, dimension(nzlev) :: actual_deeps, used_deeps
!KAL20151110
real :: maxdin,stepk,dp,w1,w2, lwint, upint, fac1, fac2
real, dimension(nlons,nlats,kdm) :: depthmid

integer, parameter :: ilontest=-1,jlattest=-1

!print *,'enter spline_calc'

!KAL20151110 Set used deeps (can be argument deepsin , or use module array deeps)
if (present(deepsin)) then
   used_deeps = deepsin
elseif (nzlev/=ndeep) then
   print *,'Dimension mismatch in'
   print *,nzlev,ndeep
   stop '(spline_calc)'
else
   used_deeps=deeps
end if

if  (trim(vimethod)=='spline') then
   write(6,'(a)',advance='no') 'Spline ticker (10 dots):'

   ! NB - splines can be slow and is point-by-point. To speed it up, group
   ! more operations together (slowness mainly due to LAPACK/ESSL call). 
   do jlat=1,nlats


      ! Ticker for the impatient ....
      if (mod(jlat,nlats/10)==0) then
         write(6,'(a1)',advance='no') '.'
         call flush(6)
      end if

      do ilon=1,nlons
      if(mask(ilon,jlat) .and. depthin(ilon,jlat,kdm)>0.1 ) then

         tmpfld=fldin  (ilon,jlat,:)
         tmpd  =depthin(ilon,jlat,:)
         !maxdin=maxval(depthin(ilon,jlat,:))
         call ini_spline(tmpd,tmpfld,kdm,subd,subfld,sub_last)

!KAL20151110
         do k=1,nzlev
            actual_deeps(k) = actual_depth(used_deeps(k),depthin(ilon,jlat,kdm))
         end do
!KAL20151110


!KAL20151110
         !call compute_spline(subd,subfld,kdm,sub_last,sfld,used_deeps,nzlev)
         call compute_spline(subd,subfld,kdm,sub_last,sfld,actual_deeps,nzlev)
!KAL20151110
         fldout(ilon,jlat,:)=sfld
      else
         fldout(ilon,jlat,:)=undef
      endif
   enddo
   enddo


else if  (trim(vimethod)=='staircase') then

   ! Fast option - use layer value, not optimized
   fldout=undef
   do kzlevel=1,nzlev
   do klevel=1,kdm
   do jlat=1,nlats
   do ilon=1,nlons
   if (depthin(ilon,jlat,kdm)>.1 .and. mask(ilon,jlat)) then
!KAL20151110
      actual_deeps(kzlevel) = actual_depth(used_deeps(kzlevel),depthin(ilon,jlat,kdm))
!KAL20151110


!      ! If shallower than 1st layer lower interface use 1st layer values
!      if (klevel==1) then
!         if (depthin(ilon,jlat,klevel)>used_deeps(kzlevel)) then
!            fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
!         end if
!
!      ! If shallower than X layer lower interface and deeper than
!      ! X layer upper interface use layer values
!      else 
!         if (depthin(ilon,jlat,klevel  )> used_deeps(kzlevel) .and. &
!             depthin(ilon,jlat,klevel-1)<=used_deeps(kzlevel)) then
!            fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
!         end if
!      end if

      
      ! Does the same as above but without branching
      lwint=depthin(ilon,jlat,klevel)            ! deepest level from model
      upint=depthin(ilon,jlat,max(klevel-1,1))   ! shallowest level from model
      upint=min(klevel-1,1)*upint                ! = upint unless klevel=1, when its 0
!KAL20151110
!     fac1=(1.+sign(1.,used_deeps(kzlevel)-upint))*.5 ! 1 if upint<used_deeps(kzlevl), 0 otherwise
!     fac2=(1.+sign(1.,lwint-used_deeps(kzlevel)))*.5 ! 1 if lwint>used_deeps(kzlevl), 0 otherwise
      fac1=(1.+sign(1.,actual_deeps(kzlevel)-upint))*.5 ! 1 if upint<actual_deeps(kzlevl), 0 otherwise
      fac2=(1.+sign(1.,lwint-actual_deeps(kzlevel)))*.5 ! 1 if lwint>actual_deeps(kzlevl), 0 otherwise
!KAL20151110
      fac1=fac1*fac2

      ! If fac1 = 1, use fldin. Otherwise use old value of fldout
      fldout(ilon,jlat,kzlevel) =     fac1 *fldin (ilon,jlat,klevel) + &
                                  (1.-fac1)*fldout(ilon,jlat,kzlevel)


   end if
   end do
   end do
   end do
   end do

elseif  (trim(vimethod)=='linear') then

   ! Fast option - use weighted value, not optimized

   ! Layer mid points
   do klevel=1,kdm
   do jlat=1,nlats
   do ilon=1,nlons
      if (klevel==1) then
         depthmid(ilon,jlat,klevel)=depthin(ilon,jlat,klevel)*.5
      else
         depthmid(ilon,jlat,klevel)= &
         (depthin(ilon,jlat,klevel-1)+depthin(ilon,jlat,klevel))*.5
      end if
   end do
   end do
   end do


   fldout=undef
   do kzlevel=1,nzlev
   do klevel=1,kdm
   do jlat=1,nlats
   do ilon=1,nlons
   if (depthin(ilon,jlat,kdm)>.1 .and. mask(ilon,jlat)) then
!KAL20151110
      actual_deeps(kzlevel) = actual_depth(used_deeps(kzlevel),depthin(ilon,jlat,kdm))
!KAL20151110

!      ! If shallower than 1st layer midpoint use 1st layer values
!      if (klevel==1) then
!         if (depthmid(ilon,jlat,klevel)>used_deeps(kzlevel)) then
!            fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
!         end if
!
!      ! If deeper than last layer midpoint use last layer values
!      ! (if model is deep enough!!)
!      else if (klevel==kdm) then
!         if (depthin (ilon,jlat,klevel)>dused_eeps(kzlevel) .and. &
!             depthmid(ilon,jlat,klevel)<dused_eeps(kzlevel)) then
!            fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
!         end if
!
!      ! Otherwise we are straddled between two interior layers
!      else 
!         if (depthmid(ilon,jlat,klevel  )> dused_eeps(kzlevel) .and. &
!             depthmid(ilon,jlat,klevel-1)<=dused_eeps(kzlevel)) then
!
!            dp=1./(depthmid(ilon,jlat,klevel) - depthmid(ilon,jlat,klevel-1))
!            w1= (dused_eeps(kzlevel) - depthmid(ilon,jlat,klevel-1))
!            w2=-(dused_eeps(kzlevel) - depthmid(ilon,jlat,klevel))
!            w1=w1*dp
!            w2=w2*dp
!            fldout(ilon,jlat,kzlevel)= w1*fldin(ilon,jlat,klevel  ) + &
!                                       w2*fldin(ilon,jlat,klevel-1)
!         end if
!      end if

      ! Does the same as above but with less branching
      km1=max(klevel-1,1)
      upint=depthmid(ilon,jlat,km1)
      upint=min(klevel-1,1)*upint        ! = upint unless klevel=1, when its 0
      fac1 = real(klevel/kdm)            ! = 1 when klevel==kdm
      lwint=depthmid(ilon,jlat,klevel)   ! deepest level from model
      lwint=(1.-fac1)*lwint+fac1*depthin(ilon,jlat,klevel)

!KAL20151110
!     w1=used_deeps(kzlevel)-upint
!     w2=lwint-used_deeps(kzlevel)
      w1=actual_deeps(kzlevel)-upint
      w2=lwint-actual_deeps(kzlevel)
!KAL20151110
      dp=1./max(lwint-upint,1e-4)
      fac1=(1.+sign(1.,w1))*.5 ! 1 if upint<actual_deeps(kzlevl), 0 otherwise
      fac2=(1.+sign(1.,w2))*.5 ! 1 if lwint>actual_deeps(kzlevl), 0 otherwise
      w1=w1*dp
      w2=w2*dp
      fac1=fac1*fac2           ! 1 if in correct depth range

      ! If fac1 = 1, use fldin weighted average. Otherwise use old value of fldout
      fldout(ilon,jlat,kzlevel) =      &
              fac1* (w1*fldin(ilon,jlat,klevel) + w2*fldin (ilon,jlat,km1)) + &
          (1.-fac1)*fldout(ilon,jlat,kzlevel)

   end if
   end do
   end do
   end do
   end do
else
   print *,'Unknown vertical Interpolation method '//trim(vimethod)
end if

end subroutine spline_calc

subroutine spline_calc_1d(fldin,depthin,fldout,nzlev,kdm)
implicit none
integer, intent(in) :: nzlev,kdm
real, dimension(kdm)  ,intent(in) :: fldin,depthin
real, dimension(nzlev),intent(out) :: fldout !,depthout

integer :: ilon,jlat,sub_last,k
real, dimension(kdm) :: tmpfld,tmpd
real, dimension(kdm) :: subfld,subd
real, dimension(ndeep) :: sfld
real :: maxdin


if (nzlev/=ndeep) then
   print *,'Dimension mismatch in'
   print *,nzlev,ndeep
   stop '(mod_spline_calc:spline_calc_1d)'
end if


if  (trim(vimethod)/='spline') then
   print *,'1D vertical interpolation can only use spline interpolation method'
   stop '(mod_spline_calc:spline_calc_1d)'
end if


tmpfld=fldin  
tmpd  =depthin
maxdin=maxval(depthin)
call ini_spline(tmpd,tmpfld,kdm,subd,subfld,sub_last)
call compute_spline(subd,subfld,kdm,sub_last,sfld,deeps,ndeep)
fldout=sfld
end subroutine spline_calc_1d


subroutine ini_spline(tmpd,tmpfld,nz,subd,subfld,sub_last)
! setting up vectors to interpolate from.
! neglecting thin layers
! For spline version
   implicit none
   integer, intent(in)  :: nz
   real, intent(in) , dimension(nz) :: tmpd, tmpfld
   real, intent(out), dimension(nz) :: subd, subfld
   integer, intent(out) :: sub_last

   integer ip(nz)
   integer l,k,nip, tmpip

   ! nip is lowest mass-containing layer
   ip(1)=1
   l=1
   do k=2,nz
      !print *,'diff',(tmpd(k)-tmpd(k-1))

      !KAL -- added check for "valid" horizontal points here
      !    -- A depth is always defined, but neighbouring points may be below
      !    -- sea floor, in which case they may be undefined
      if ((tmpd(k)-tmpd(k-1)) > 0.1 .and. tmpfld(k)/=undef) then
         l=l+1
         ip(l)=k
      endif
   enddo
   nip=l
   !print *,'nip,nz',nip,nz

   ! Put into array
   l=0
   do k=1,nip
      l=l+1
      tmpip=ip(k) 
      subfld(l)=tmpfld(tmpip)
      subd  (l)=tmpd  (tmpip)
   enddo
   sub_last=l
end subroutine ini_spline

subroutine compute_spline(subd,subfld,nz,sub_last,sfld,deeps2,ndepths)
   !use m_calcconst
   implicit none
   integer, intent(in) :: sub_last,ndepths,nz
   real, intent(in), dimension(nz) :: subd, subfld
   real, intent(out),  dimension(ndepths) :: sfld
   real, intent(in),  dimension(ndepths) :: deeps2

   real, dimension(sub_last):: a, b
   real, dimension(0:sub_last):: c
   integer ic
   integer k

   call calcconst(sub_last,subfld,subd,a,b,c,'A','a')
   sfld=undef
   do k=1,ndepths
      sfld(k)=calcu(sub_last,a,b,c,subd,deeps2(k),ic)
   enddo
end subroutine compute_spline


subroutine calcconst(kk,U,deep,a,b,c,cflag,bflag)
   implicit none 
   integer, intent(in) :: kk
   character(len=1), intent(in) :: cflag
   character(len=1), intent(in) :: bflag
   real, intent(in), dimension(kk)::U,deep
   real, intent(out), dimension(kk)::a,b
   real, intent(out), dimension(0:kk)::c

   integer i,j,k, dimen, eq, info
   real*8, dimension(3*kk+1)::x,rhs,Z
   real*8, dimension(4*(3*kk+1))::lapwork
   integer, dimension(3*kk+1)::IPVT
   real RCOND
   real*8, dimension(3*kk+1,3*kk+1)::Amatr
     
   dimen = 3*kk+1

! Generate Amatr and rhs

   rhs = 0.
   Amatr = 0.

! First layer
   Amatr(1,1) = -1.                ! c_0
   Amatr(1,2) = 1.                 ! c_1
   Amatr(1,3) = .5*deep(1)         ! b_1
   Amatr(1,4) = (.5*deep(1))**2    ! a_1

   Amatr(2,3) = 1.                 ! b_1
   Amatr(2,4) = deep(1)              ! a_1

   Amatr(3,1) = .5                   ! c_0
   Amatr(3,2) = .5                   ! c_1
   Amatr(3,3) = 3./8.*deep(1)      ! b_1
   Amatr(3,4) = 7./24.*deep(1)**2  ! a_1
   rhs(3) = U(1)

   eq = 3
! Interior layers
! 
   do i=2,kk 
! Cont i function 

      eq = eq + 1
      Amatr(eq,3*(i-1)-1) =-1.            ! C(i-1)
      Amatr(eq,3*(i-1)  ) =-deep(i-1)       ! b(i-1)
      Amatr(eq,3*(i-1)+1) =-deep(i-1)**2    ! a(i-1)

      Amatr(eq,3*(i  )-1) = 1.            ! C(i)
      Amatr(eq,3*(i  )  ) = deep(i-1)       ! b(i)
      Amatr(eq,3*(i  )+1) = deep(i-1)**2    ! a(i)

! Cont derivative

      eq = eq + 1
      if (cflag == 'A') then
         ! Continuity in derivative
         Amatr(eq,3*(i-1)  ) =-1.            ! b(i-1)
         Amatr(eq,3*(i-1)+1) =-deep(i-1)*2.  ! a(i-1)
         Amatr(eq,3*(i  )  ) = 1.            ! b(i)
         Amatr(eq,3*(i  )+1) = deep(i-1)*2.  ! a(i)
      elseif (cflag == 'B') then
         ! interpolates average interface velocity (break continuity in derivative)
         Amatr(eq,3*(i  )-1) = 1.            ! C(i)
         Amatr(eq,3*(i  )  ) = deep(i-1)       ! b(i)
         Amatr(eq,3*(i  )+1) = deep(i-1)**2    ! a(i)
         rhs(eq)=0.5*(U(i-1) + U(i))
      elseif (cflag == 'C') then
         ! interpolates weighted average interface velocity (break continuity in derivative)
         Amatr(eq,3*(i  )-1) = 1.            ! C(i)
         Amatr(eq,3*(i  )  ) = deep(i-1)       ! b(i)
         Amatr(eq,3*(i  )+1) = deep(i-1)**2    ! a(i)
         if (i == 2) then
            rhs(eq)=0.5*(deep(i-1)*U(i-1) +  (deep(i)-deep(i-1))*U(i) )/deep(i)
         else
            rhs(eq)=0.5*( (deep(i-1)-deep(i-2))*U(i-1) + (deep(i)-deep(i-1))* U(i))  / &
                                  (deep(i)-deep(i-2))
         endif
      else
         stop 'invalid cflag in calcconst' 
      endif
 

! Mean 
      eq = eq + 1
      Amatr(eq,3*(i  )-1) = 1.                                                 ! C(i)
      Amatr(eq,3*(i  )  ) = 0.5*(deep(i-1)+deep(i))                            ! b(i)
      Amatr(eq,3*(i  )+1) = (deep(i-1)**2 + deep(i)*deep(i-1)+deep(i)**2)/3.0  ! a(i)
      rhs(eq) = U(i)
   enddo

! Closure

   eq = eq + 1 
   i = kk  
   if (bflag == 'a') then
      Amatr(eq,3*(i  )-1) = 1.          ! C(kk)
      Amatr(eq,3*(i  )  ) = deep(i)       ! b(kk)
      Amatr(eq,3*(i  )+1) = deep(i)**2    ! a(kk)
      rhs(eq)=U(i)
   elseif (bflag == 'b') then
      Amatr(eq,3*(i  )  ) = 1.            ! b(i)
      Amatr(eq,3*(i  )+1) = deep(i)*2.    ! a(i)
   else
       stop 'invalid bflag in calcconst'
   endif
          

! Solve

!   write(*,'(13g10.2)')transpose(Amatr)

   ! Arch definition time ...
#if defined (SGI)
   call DGECO(Amatr,dimen,dimen,IPVT,RCOND,Z)
   call DGESL(Amatr,dimen,dimen,IPVT,rhs,0)
#elif defined (ESSL)
   call dgef(Amatr,dimen,dimen,IPVT)
   call dges(Amatr,dimen,dimen,IPVT,rhs,0)
#elif defined (LAPACK)
   !call DGECON('1',dimen,Amatr,dimen,'I',RCOND,lapwork,ipvt,info)
   call DGESV(dimen,1,Amatr,dimen,ipvt,rhs,dimen,info)
#else
#error calcconst has not a library routine defined for this machine
   print *,'calcconst has not a library routine'
   print *,'defined for this machine'
   stop '(m_calcconst.F90)'
#endif


   c(0) = rhs(1)
   c(1:kk) = rhs(1+1:3*kk+1-2:3)
   b(1:kk) = rhs(1+2:3*kk+1-1:3)
   a(1:kk) = rhs(1+3:3*kk+1-0:3)
end subroutine  calcconst

   function calcu(k,a,b,c,d,z,ic)
   implicit none
   integer, intent(in)::k
   real, intent(in), dimension(k)::a,b,d
   real, intent(in), dimension(0:k)::c
   real, intent(in)::z
   real calcu
   integer, intent(out):: ic
   integer i

   if (z .le. d(1)) then
      if (z.le.0.5*d(1)) then
         calcu = c(0)
         ic = 1
         return
      else
         calcu = a(1)*z**2+b(1)*z+c(1)
         ic = 1
         return
      endif
   endif

   do i = 2,k
      if (z .le. d(i)) then
         calcu = a(i)*z**2+b(i)*z+c(i)
         ic = i
         return
      endif
   enddo
   calcu=undef
!   stop 'missed in calcu'
end function calcu




subroutine  get_depth_levels()
   implicit none
   character(len=*), parameter  :: infile_depths   = 'depthlevels.in'
   character(len=20) :: tmpchar

   logical :: ex
   integer :: ios,counter
   real :: v1,v2

   ! Get depth levels to use
   inquire(exist=ex,file=infile_depths)
   if (.not. ex) then
      print '(a)','You must specify depth levels to use in'
      print '(a)','the file called  '//infile_depths
      stop '(get_depth_levels)'
   endif

   ! Start reading
   ios=0
   counter=1
   open(10,file=infile_depths,form='formatted')
   read(10,*,iostat=ios) ndeep,tmpchar
   print '(i3,a)',ndeep,' Depth levels:'
   allocate(deeps(ndeep))
   do while (ios==0 .and. counter<=ndeep)
      read(10,*,iostat=ios) deeps(counter)
      if (deeps(counter)<0.) then
         print '(a,i3,f10.2,"(=height from bottom)")','Depth level :',counter,deeps(counter)
      else 
         print '(a,i3,f10.2)','Depth level :',counter,deeps(counter)
      end if
      counter=counter+1
   end do
   close(10)
   if (counter/=ndeep+1) then
      print *,'ERROR: Not all depth levels set'
      stop '(get_depth_levels)'
   endif
   if (ios/=0) then
      print *,'ERROR: occured when reading '//infile_depths
      stop '(get_depth_levels)'
   endif
end subroutine get_depth_levels


subroutine  clear_spline_calc()
implicit none
   if (allocated(deeps)) then
      deallocate(deeps)
   end if
   ndeep=0
end subroutine


real function actual_depth(mydepth,bottom_depth) 
implicit none
real, intent(in) :: mydepth,bottom_depth
actual_depth = mydepth
! Negative values - interpret as height from bottom
if (mydepth< 0.0 ) then
   actual_depth = bottom_depth + mydepth
   ! Above sea level - reset to be below sea floor (masked)
   if (actual_depth< 0.0 ) then
      actual_depth = bottom_depth + 10.
   end if
end if
end function

end module mod_spline_calc
