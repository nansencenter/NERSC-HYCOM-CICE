module m_density
contains

   subroutine density(new,nrm,cdata,theta,deep,maxA,cid)
   use mod_data_new
   use mod_netcdf_ops
   implicit none
   integer, intent(in) :: nrm
   type(data_new), intent(in) :: new(nrm)
   character(len=2), intent(in) :: cdata
   integer , intent(in) :: deep
   real, intent(out) :: theta
   real, intent(inout) :: maxA
   character(len=1) cid
   character(len=80) cflt

   integer, parameter :: nx=21    ! 31
   integer, parameter :: ny=21    ! 31
   real   , parameter :: dv=0.5     !  2
   integer ia(-nx:nx,-ny:ny)
   real     A(-nx:nx,-ny:ny)
   real    as(-nx:nx,-ny:ny)

   real*8 AP(3)

   real*8 p(2,2),lambda(2),zz(2,2),w(10),t1,t2
   integer info
   integer i,j,m
   real tmp
   real time,u,v
   real scale,scale1
   real varu,varv,covuv,aveu,avev
   real torad,todeg
   real x(2)
   integer :: ncid,vid,vertind, nxid, nyid, tid
   integer, allocatable :: tmp2(:)
   torad=pi/180.0
   todeg=1.0/torad


   !print '(3a)','Statistics for ',cdata,' data'
   !read(chdeep,*) deep


   ! Calculate variance, average and covariance between u and v
   ! Also calculate histogram (?) ia
   aveu=0.0; avev=0.0; varu=0.0; varv=0.0; covuv=0.0
   ia=0
   do m=1,nrm
      u=new(m)%u
      v=new(m)%v


      aveu=aveu+u
      varu=varu+u*u
      avev=avev+v
      varv=varv+v*v
      covuv=covuv+u*v

      ! KAL - histogram?
      if (u >= 0.0) then
         i=int(u/dv)
      else
         i=int(u/dv)-1
      endif
      if (v >= 0.0) then
         j=int(v/dv)
      else
         j=int(v/dv)-1
      endif
      i=min(max(-nx,i),nx)
      j=min(max(-ny,j),ny)
      !print *,i,j
      ia(i,j)=ia(i,j)+1

   enddo
101 close(10)

   ! Finalize histogram values
   A=float(ia)/float(nrm)
   maxA=max(maxval(A),maxA)

   ! Finalize variance ++ values
   scale=1.0/float(nrm)
   scale1=1.0/float(nrm-1)
   varu=scale1*(varu-scale*aveu*aveu)
   varv=scale1*(varv-scale*avev*avev)
   covuv=scale1*(covuv-scale*aveu*avev)
   aveu=aveu/float(nrm)
   avev=avev/float(nrm)
!   print '(a,f10.4)','   Mean(u) is: ',aveu
!   print '(a,f10.4)','   Mean(v) is: ',avev
!   print *
   p(1,1)=varu; p(2,2)=varv; p(1,2)=covuv; p(2,1)=covuv
!   print '(a,2f10.4)','   Covar matrix is: ',p(1,1),p(1,2)
!   print '(a,2f10.4)','                    ',p(2,1),p(2,2)
!   print *

   ! Eigenvectors - dominant directions in u/v space
   AP(1)=p(1,1)
   AP(2)=p(2,1)
   AP(3)=p(2,2)
!On fimm no DSPEV, try with DSYEV!!!
#if defined(AIX)
   call DSPEV(1,AP,lambda,p,2,2,w,9)
#else
   call DSYEV('V','U',2,p,2,lambda,w,9,info)
#endif
!   if (info /= 0) print *,'info=',info
!   print '(a,2f10.4)','   eigenvalues:',lambda(1:2)
!   print '(a)','   rsm eignevectors:'
!   print '(tr3,f10.4,a,f10.4)',p(1,1),'   ',p(1,2)
!   print '(tr3,f10.4,a,f10.4)',p(2,1),'   ',p(2,2)
!   print *


   !tmp=0.5*sqrt((varu+varv)**2-4.0*(varu*varv-covuv**2))
   !print *,'lambda1=',0.5*(varu+varv)+tmp
   !print *,'lambda2=',0.5*(varu+varv)-tmp


!   ! Write to average file  
!   open(10,file=trim(fname)//'.ave',status='replace')
!   write(10,'(2f10.4)')aveu,avev
!   close(10)

   theta=atan2(p(2,2),p(1,2))*todeg
   if (theta > 90.0) theta=theta-180.0
   if (theta < -90.0) theta=theta+180.0

!   ! Append to table 
!   open(10,file='table.txt',position='append')
!   write(10,'(a,a,8(f6.1,a),f6.0)')fname(14:),'&',aveu,'&',avev,'&', &
!      varu,'&',varv,'&',covuv,'&',lambda(1),'&',lambda(2),'&',theta
!   close(10)


!   open(10,file=trim(fname)//'_eqn.mcr',status='unknown')
!   write(10,'(a)')'#!MC 700'
!   write(10,'(a)')'$!MACROFUNCTION'
!   select case (cid)
!   case('U')
!      write(10,'(5a)')'  NAME = "station',chstat,'.',chdeep,'_eqn"'
!   case('F')
!      write(10,'(5a)')'  NAME = "station',chstat,'.',chdeep,'_F_eqn"'
!   case('T')
!      write(10,'(5a)')'  NAME = "station',chstat,'.',chdeep,'_T_eqn"'
!   end select
!   write(10,'(a)')'$!ALTERDATA'
!   write(10,'(2(a,f10.4),a)',advance='no')'EQUATION = ''V6=(',varv,'*(V3-',aveu,')**2'
!   write(10,'(2(a,f10.4),a)',advance='no')              '+',varu,'*(V4-',avev,')**2'
!   write(10,'(6(a,f10.4),a)')              '-2.0*',covuv,'*(V3-',aveu,')*(V4-',avev,'))/(',varu,'*',varv,'-',covuv,'**2)'''
!   write(10,'(a)')'$!ENDMACROFUNCTION'
!   close(10)
!
!   open(10,file=trim(fname)//'.dens',status='unknown')
!      write(10,*)'TITLE = "',trim(fname),'"'
!      write(10,*)'VARIABLES = "i"  "j"  "x" "y" "field"'
!      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',2*nx+1,', J=',2*ny+1,', K=1'
!      write(10,'(30I4)')((i,i=-nx,nx),j=-ny,ny)
!      write(10,'(30I4)')((j,i=-nx,nx),j=-ny,ny)
!      write(10,900)((float(i*dv),i=-nx,nx),j=-ny,ny)
!      write(10,900)((float(j*dv),i=-nx,nx),j=-ny,ny)
!      write(10,900)((A(i,j),i=-nx,nx),j=-ny,ny)
!   close(10)
!   900  format(10(1x,e12.5))


! KAL - added netcdf option - keeps all info in one file

   call ncopencreate(ncfile,ncid)
   call ncinqdefdim (ncid,'vlevel'  ,NF90_UNLIMITED,vid)
   call ncinqdefvertvar(ncid,vid,deep,vertind)
   call ncinqputvar(ncid,'uave_'//cdata//'_'//cid,(/vid/),(/aveu/),start=(/vertind/))
   call ncinqputvar(ncid,'vave_'//cdata//'_'//cid,(/vid/),(/avev/),start=(/vertind/))
   call ncinqputvar(ncid,'uvar_'//cdata//'_'//cid,(/vid/),(/varu/),start=(/vertind/))
   call ncinqputvar(ncid,'vvar_'//cdata//'_'//cid,(/vid/),(/varv/),start=(/vertind/))
   call ncinqputvar(ncid,'uvcov_'//cdata//'_'//cid,(/vid/),(/covuv/),start=(/vertind/))
   call ncinqputvar(ncid,'theta_'//cdata//'_'//cid,(/vid/),(/theta/),start=(/vertind/))

   call ncinqdefdim (ncid,'nx'  ,2*nx+1,nxid)
   call ncinqdefdim (ncid,'ny'  ,2*ny+1,nyid)
   allocate(tmp2(-nx:nx))
   do i=-nx,nx
      tmp2(i)=i
   end do
   call ncinqputvar(ncid,'nx',(/nxid/),real(tmp2)*dv,(/1/))
   deallocate(tmp2)
   allocate(tmp2(-ny:ny))
   do i=-ny,ny
      tmp2(i)=i
   end do
   call ncinqputvar(ncid,'ny',(/nyid/),real(tmp2)*dv,(/1/))
   call ncinqputvarslice(ncid,'A_'//cdata//'_'//cid,(/nxid,nyid,vid/),A,(/1,1,vertind/))

!  Also add time series here - not the time is not record dimension, vertical
!  level is
   call ncinqdefdim (ncid,'time'  ,nrm,tid)
   call ncinqputvarslice(ncid,'u_'//cdata//'_'//cid,(/tid,vid/),new(:)%u,start=(/1,vertind/))
   call ncinqputvarslice(ncid,'v_'//cdata//'_'//cid,(/tid,vid/),new(:)%v,start=(/1,vertind/))
   call ncinqputvarslice(ncid,'spd_'//cdata//'_'//cid,(/tid,vid/),new(:)%speed,start=(/1,vertind/))

   cflt=''
   if (cid=='U') then
      cflt=' unfiltered'
   else if (cid=='T') then
      cflt=' tidal component'
   else if (cid=='F') then
      cflt=' residual component'
   end if
   call ncputatt(ncid,'uave_'//cdata//'_'//cid,'comment', &
     'Mean u for '//cdata//trim(cflt))
   call ncputatt(ncid,'vave_'//cdata//'_'//cid,'comment', &
     'Mean v for '//cdata//trim(cflt))
   call ncputatt(ncid,'uvar_'//cdata//'_'//cid,'comment', &
     'Variance of u for '//cdata//trim(cflt))
   call ncputatt(ncid,'vvar_'//cdata//'_'//cid,'comment', &
     'Variance of v for '//cdata//trim(cflt))
   call ncputatt(ncid,'uvcov_'//cdata//'_'//cid,'comment',&
     'covariance of u and v components '//cdata//trim(cflt))
   call ncputatt(ncid,'theta_'//cdata//'_'//cid,'comment', &
      'principal direction  for '//cdata//trim(cflt))
   call ncputatt(ncid,'u_'//cdata//'_'//cid,'comment', &
      'u timeseries  for '//cdata//trim(cflt))
   call ncputatt(ncid,'v_'//cdata//'_'//cid,'comment', &
      'v  timeseries  for '//cdata//trim(cflt))
   call ncputatt(ncid,'spd_'//cdata//'_'//cid,'comment', &
      'speed timeseries  for '//cdata//trim(cflt))
   call ncerr(nf90_close(ncid))

end subroutine density
end module m_density
