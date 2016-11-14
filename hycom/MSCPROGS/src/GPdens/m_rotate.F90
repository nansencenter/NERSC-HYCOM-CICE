module m_rotate
contains


subroutine rotate(cmfilt,gpfilt,nrm,thetacm,thetagp,deep)
   use mod_data_new
   use mod_netcdf_ops
   implicit none
   integer, intent(in) :: nrm
   integer, intent(in) :: deep
   type(data_new), intent(in)  :: cmfilt(nrm),gpfilt(nrm)
   type(data_new), dimension(nrm) :: tmpcm,tmpgp
   real, intent(in) :: thetacm,thetagp

   integer m
   real torad,todeg,theta
   integer ncid, vid, tid,vertind

   torad=pi/180.0
   todeg=1.0/torad

   
!  open(10,file='rotated.out',form='formatted',status='unknown',access='sequential')
   do m=1,nrm
      tmpcm(m)=cmfilt(m)
      theta=thetacm*torad
      tmpcm(m)%u=+cos(theta)*cmfilt(m)%u+sin(theta)*cmfilt(m)%v
      tmpcm(m)%v=-sin(theta)*cmfilt(m)%u+cos(theta)*cmfilt(m)%v
      tmpcm(m)%t=thetacm

      tmpgp(m)=gpfilt(m)
      theta=thetagp*torad
      tmpgp(m)%u=+cos(theta)*gpfilt(m)%u+sin(theta)*gpfilt(m)%v
      tmpgp(m)%v=-sin(theta)*gpfilt(m)%u+cos(theta)*gpfilt(m)%v
      tmpgp(m)%t=thetagp

      !write(10,'(i6,2(f9.4,2f9.3,4f8.2))')deep,tmpgp,tmpcm
      enddo
!   close(10)

   call ncopencreate(ncfile,ncid)
   call ncinqdefdim (ncid,'vlevel'  ,NF90_UNLIMITED,vid)
   call ncinqdefdim (ncid,'time'  ,nrm,tid)
   call ncinqdefvertvar(ncid,vid,deep,vertind)
   call ncinqputvarslice(ncid,'u_CM_princ_F',(/tid,vid/),tmpcm(:)%u  ,start=(/1,vertind/))
   call ncinqputvarslice(ncid,'v_CM_princ_F',(/tid,vid/),tmpcm(:)%v  ,start=(/1,vertind/))
   call ncinqputvar(ncid,'v_CM_princ_theta',(/vid/),(/thetacm/),start=(/vertind/))
   call ncinqputvarslice(ncid,'u_GP_princ_F',(/tid,vid/),tmpgp%u  ,start=(/1,vertind/))
   call ncinqputvarslice(ncid,'v_GP_princ_F',(/tid,vid/),tmpgp%v  ,start=(/1,vertind/))
   call ncinqputvar(ncid,'v_GP_princ_theta',(/vid/),(/thetagp/),start=(/vertind/))
   call ncerr(NF90_close(ncid))

end subroutine rotate
end module m_rotate
