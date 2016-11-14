module m_correlation
contains
subroutine correlation(cm2,gp2,cmfilt,gpfilt,cmtide,gptide,nrm,deep)
   use mod_data_new
   use mod_netcdf_ops
   implicit none
   integer, intent(in) ::  nrm
   integer, intent(in) ::  deep
   type(data_new), intent(in) :: cm2(nrm),gp2(nrm)
   type(data_new), intent(in) :: cmtide(nrm),gptide(nrm)
   type(data_new), intent(in) :: cmfilt(nrm),gpfilt(nrm)
   integer m
   complex w1(nrm),w2(nrm)
   complex w1filt(nrm),w2filt(nrm)
   complex w1tide(nrm),w2tide(nrm)
   complex w12,rho
   complex w12filt,rhofilt
   complex w12tide,rhotide
   real w11,w22,absrho,argrho
   real w11filt,w22filt,absrhofilt,argrhofilt
   real w11tide,w22tide,absrhotide,argrhotide
   integer  :: ncid, vid, vertind

   do m=1,nrm
      w1(m)=cmplx(cm2(m)%u,cm2(m)%v)
      w2(m)=cmplx(gp2(m)%u,gp2(m)%v)
      w1filt(m)=cmplx(cmfilt(m)%u,cmfilt(m)%v)
      w2filt(m)=cmplx(gpfilt(m)%u,gpfilt(m)%v)
      w1tide(m)=cmplx(cmtide(m)%u,cmtide(m)%v)
      w2tide(m)=cmplx(gptide(m)%u,gptide(m)%v)
   enddo

   w12=dot_product(w1,w2)
   w11=dot_product(w1,w1)
   w22=dot_product(w2,w2)
   rho=w12/(sqrt(w11)*sqrt(w22))
   absrho=abs(rho)
   argrho=atan(aimag(rho)/real(rho))

   w12filt=dot_product(w1filt,w2filt)
   w11filt=dot_product(w1filt,w1filt)
   w22filt=dot_product(w2filt,w2filt)
   rhofilt=w12filt/(sqrt(w11filt)*sqrt(w22filt))
   absrhofilt=abs(rhofilt)
   argrhofilt=atan(aimag(rhofilt)/real(rhofilt))

   w12tide=dot_product(w1tide,w2tide)
   w11tide=dot_product(w1tide,w1tide)
   w22tide=dot_product(w2tide,w2tide)
   rhotide=w12tide/(sqrt(w11tide)*sqrt(w22tide))
   absrhotide=abs(rhotide)
   argrhotide=atan(aimag(rhotide)/real(rhotide))

!   open(10,file='compcorr.dat')
!      write(10,'(i5,6f10.2)')deep,absrho,argrho*180.0/3.14159, &
!                                  absrhofilt,argrhofilt*180.0/3.14159, &
!                                  absrhotide,argrhotide*180.0/3.14159
!   close(10)


   call ncopencreate(ncfile,ncid)
   call ncinqdefdim (ncid,'vlevel'  ,NF90_UNLIMITED,vid)
   call ncinqdefvertvar(ncid,vid,deep,vertind)
   call ncinqputvar(ncid,'absrho_U',(/vid/),(/absrho            /),start=(/vertind/))
   call ncinqputvar(ncid,'argrho_U',(/vid/),(/argrho*180/3.14159/),start=(/vertind/))
   call ncinqputvar(ncid,'absrho_F',(/vid/),(/absrhofilt        /),start=(/vertind/))
   call ncinqputvar(ncid,'argrho_F',(/vid/),(/argrhofilt*180/3.14159/),start=(/vertind/))
   call ncinqputvar(ncid,'absrho_T',(/vid/),(/absrhotide            /),start=(/vertind/))
   call ncinqputvar(ncid,'argrho_T',(/vid/),(/argrhotide*180/3.14159/),start=(/vertind/))
   call ncerr(NF90_close(ncid))


end subroutine correlation
end module m_correlation
