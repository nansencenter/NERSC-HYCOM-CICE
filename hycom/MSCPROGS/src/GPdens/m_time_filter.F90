module m_time_filter
contains
subroutine time_filter(cm2,gp2,cmfilt,gpfilt,cmtide,gptide,nrm,hours,chfilt,deep)
   use mod_data_new
   use m_mdian1
   use m_shfact
   use m_shfilt
   implicit none
   integer, intent(in) ::  nrm
   type(data_new), intent(in) ::  cm2(nrm),gp2(nrm)
   real, intent(in) :: hours
   type(data_new), intent(out) :: cmfilt(nrm),gpfilt(nrm)
   type(data_new), intent(out) :: cmtide(nrm),gptide(nrm)
   character(len=3), intent(in) :: chfilt
   integer, intent(in) :: deep

   real tmp(nrm)

   real, dimension(1000) :: gpu,gpv,cmu,cmv
   integer m,ih,ic,i
   integer ish
   real torad,todeg
   real sh(0:100)

   torad=pi/180.0
   todeg=1.0/torad

   ih=nint(hours/2.0)

   do m=1,nrm
      ic=0
      do i=max(1,m-ih),min(nrm,m+ih)
         if ( abs(cm2(i)%day-cm2(m)%day) <= hours ) then
            ic=ic+1
            cmu(ic)=cm2(i)%u
            cmv(ic)=cm2(i)%v
            gpu(ic)=gp2(i)%u
            gpv(ic)=gp2(i)%v
         endif
      enddo
      cmfilt(m)=cm2(m)
      gpfilt(m)=gp2(m)

      if (chfilt == 'MED') then
         call mdian1(cmu,ic,cmfilt(m)%u)
         call mdian1(cmv,ic,cmfilt(m)%v)
         call mdian1(gpu,ic,gpfilt(m)%u)
         call mdian1(gpv,ic,gpfilt(m)%v)
      elseif  (chfilt == 'AVE') then
         cmfilt(m)%u=sum(cmu(1:ic))/float(ic)
         cmfilt(m)%v=sum(cmv(1:ic))/float(ic)
         gpfilt(m)%u=sum(gpu(1:ic))/float(ic)
         gpfilt(m)%v=sum(gpv(1:ic))/float(ic)
      else
         STOP 'ERROR in m_time_filter: chfilt'
      endif
   enddo

   ish=1
   call shfact(ish,sh)
   do m=1,10
   tmp=cmfilt%u; call shfilt(ish,sh,nrm,tmp,1,cmfilt%u,1,100)
   tmp=cmfilt%v; call shfilt(ish,sh,nrm,tmp,1,cmfilt%v,1,100)
   tmp=gpfilt%u; call shfilt(ish,sh,nrm,tmp,1,gpfilt%u,1,100)
   tmp=gpfilt%v; call shfilt(ish,sh,nrm,tmp,1,gpfilt%v,1,100)
   enddo

   do m=1,nrm
      gpfilt(m)%speed=sqrt(gpfilt(m)%u**2+gpfilt(m)%v**2)
      cmfilt(m)%speed=sqrt(cmfilt(m)%u**2+cmfilt(m)%v**2)

      gpfilt(m)%dir=atan2(gpfilt(m)%v,gpfilt(m)%u)*todeg
      cmfilt(m)%dir=atan2(cmfilt(m)%v,cmfilt(m)%u)*todeg

      gptide(m)%u=gp2(m)%u-gpfilt(m)%u
      gptide(m)%v=gp2(m)%v-gpfilt(m)%v

      cmtide(m)%u=cm2(m)%u-cmfilt(m)%u
      cmtide(m)%v=cm2(m)%v-cmfilt(m)%v

      gptide(m)%speed=sqrt(gptide(m)%u**2+gptide(m)%v**2)
      cmtide(m)%speed=sqrt(cmtide(m)%u**2+cmtide(m)%v**2)

      gptide(m)%dir=atan2(gptide(m)%v,gptide(m)%u)*todeg
      cmtide(m)%dir=atan2(cmtide(m)%v,cmtide(m)%u)*todeg
   enddo


!   open(10,file='filtered.out',form='formatted',status='unknown',access='sequential')
!      do m=1,nrm
!         write(10,'(2(f9.4,2f9.3,4f8.2,tr2))')gpfilt(m),cmfilt(m)
!      enddo
!   close(10)
!
!
!   open(10,file='tides.out',form='formatted',status='unknown',access='sequential')
!      do m=1,nrm
!         write(10,'(f9.4,4f8.2))')cm2(m)%day,gptide(m)%u,gptide(m)%v,cmtide(m)%u,cmtide(m)%v
!      enddo
!   close(10)

end subroutine time_filter
end module m_time_filter
