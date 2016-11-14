module mod_spline_calc
use mod_parameters
   integer, save :: ndeep
   real, save, allocatable :: deeps(:)
   integer, dimension(:,:,:), allocatable, save :: kindex
   logical,save :: splineflag

contains

subroutine spline_calc_ini(lsplinein,deepin,ndeepin)
   implicit none
   logical, intent(in) :: lsplinein
   integer, intent(in), optional :: ndeepin
   real   , intent(in), optional :: deepin(ndeepin)

   ! If present in arg call set directly
   if (present(deepin) .and. present(ndeepin)) then
      allocate(deeps(ndeepin))
      deeps=deepin
      ndeep=ndeepin

   ! Otherwise read depthlevels.in
   else
      call get_depth_levels()
      splineflag=lsplinein
   end if

end subroutine


subroutine spline_calc(fldin,depthin,nlons,nlats,mask,fldout,nzlev,kdm)
implicit none
integer, intent(in) :: nlons,nlats,nzlev,kdm
real, dimension(nlons,nlats,kdm)  ,intent(in) :: fldin,depthin
real, dimension(nlons,nlats,nzlev),intent(out) :: fldout !,depthout
logical, intent(in) :: mask(nlons,nlats)

integer :: ilon,jlat,sub_last,k,klevel,kzlevel,klevel2
real, dimension(kdm) :: tmpfld,tmpd
real, dimension(kdm) :: subfld,subd
real, dimension(ndeep) :: sfld
real :: maxdin,stepk

integer, parameter :: ilontest=-1,jlattest=-1

!print *,'enter spline_calc'

if (nzlev/=ndeep) then
   print *,'Dimension mismatch in'
   print *,nzlev,ndeep
   stop '(spline_calc)'
end if

if  (splineflag) then
!print *,'loop begins'
do jlat=1,nlats
!if (mod(jlat,nlats/10)==0) print *,jlat,nlats
if (mod(jlat,nlats/10)==0) write(6,'(a1)',advance='no') '.'
do ilon=1,nlons
!print *,'beg',jlat,ilon,nlats,nlons,mask(ilon,jlat)
if(mask(ilon,jlat) .and. depthin(ilon,jlat,2)>0 ) then

      tmpfld=fldin  (ilon,jlat,:)
      tmpd  =depthin(ilon,jlat,:)
      maxdin=maxval(depthin(ilon,jlat,:))
      if (ilon==ilontest .and. jlat==jlattest) then
         print *
         print *,'tmpd  ',tmpd
         print *,'tmpfld',tmpfld
         print *
      end if

      call ini_spline(tmpd,tmpfld,kdm,subd,subfld,sub_last)
      if (ilon==ilontest .and. jlat==jlattest) then
         print *,'subd  ',subd(1:sub_last)
         print *,'subfld',subfld(1:sub_last)
         print *,'sub_last',sub_last
      end if

      call compute_spline(subd,subfld,kdm,sub_last,sfld,deeps,ndeep)
      if (ilon==ilontest .and. jlat==jlattest) then
         print *,'deep  ',deeps
         print *,'sfld  ',sfld
      end if

      fldout(ilon,jlat,:)=sfld
      !do k=1,ndeep
      !   depthout(ilon,jlat,k)=deeps(k)
      !   depthout(ilon,jlat,k)=min(depthout(ilon,jlat,k),maxdin)
      !end do

   else
      fldout(ilon,jlat,:)=undef
   endif
   !print *,'end',jlat,ilon,nlats,nlons
enddo
enddo

else

   !print *,'splineflag',splineflag
   fldout=undef

   ! Fast option - use layer value, can be improved a lot further
   if (.not.allocated(kindex)) then
      allocate(kindex(nlons,nlats,ndeep))
      kindex=-1
      !print *,'splineflag',splineflag

      do kzlevel=1,ndeep
      do klevel=1,kdm
      do jlat=1,nlats
      do ilon=1,nlons
         if (klevel==1) then
            if (depthin(ilon,jlat,klevel)>deeps(kzlevel)) then
               fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
               kindex(ilon,jlat,kzlevel)=klevel
            end if
         else 
            if (depthin(ilon,jlat,klevel  )> deeps(kzlevel) .and. &
                depthin(ilon,jlat,klevel-1)<=deeps(kzlevel)) then
               fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
               kindex(ilon,jlat,kzlevel)=klevel
            end if
         end if
      end do
      end do
      end do
      end do
      !print *,'End First kindex loop'

   else

      !print *,'2nd kindex loop'

      do kzlevel=1,ndeep
      do jlat=1,nlats
      do ilon=1,nlons
         klevel=kindex(ilon,jlat,kzlevel)
         klevel2=max(klevel,1)
         stepk=(1.+sign(1.,float(klevel)))*.5


         fldout(ilon,jlat,kzlevel)= &
            stepk*fldin(ilon,jlat,klevel2) + (1.-stepk)*undef
         !if (klevel/=-1) then
         !   fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
         !else
         !   fldout(ilon,jlat,kzlevel)=undef
         !end if
      end do
      end do
      end do
   end if
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
   stop '(spline_calc)'
end if


tmpfld=fldin  
tmpd  =depthin
maxdin=maxval(depthin)

call ini_spline(tmpd,tmpfld,kdm,subd,subfld,sub_last)
!print *,'subd  ',subd
!print *,'subfld',subfld
!print *,'sub_last',sub_last

call compute_spline(subd,subfld,kdm,sub_last,sfld,deeps,ndeep)
!print *,'deep  ',deeps
!print *,'sfld  ',sfld

fldout=sfld
!do k=1,ndeep
!   depthout(ilon,jlat,k)=deeps(k)
!   depthout(ilon,jlat,k)=min(depthout(ilon,jlat,k),maxdin)
!end do
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
         if ((tmpd(k)-tmpd(k-1)) > 1.0 .and. tmpfld(k)/=undef) then
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
#elif defined (IBM)
   call dgef(Amatr,dimen,dimen,IPVT)
   call dges(Amatr,dimen,dimen,IPVT,rhs,0)
#elif defined (LAPACK)
   !call DGECON('1',dimen,Amatr,dimen,'I',RCOND,lapwork,ipvt,info)
   call DGESV(dimen,1,Amatr,dimen,ipvt,rhs,dimen,info)
#else
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
      print '(a,i3,f10.2)','Depth level :',counter,deeps(counter)
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

end module mod_spline_calc
