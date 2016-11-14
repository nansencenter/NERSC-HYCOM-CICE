module m_merge
contains
subroutine merge(gp,cm,gp2,cm2,ngp,ncm,maxnr,nrm)
   use mod_modtype
   implicit none
   integer, intent(in) ::  ngp,ncm,maxnr
   type(data_new), intent(in)  ::  cm(maxnr),gp(maxnr)
   type(data_new), intent(out) ::  cm2(maxnr),gp2(maxnr)
   integer, intent(out) ::  nrm

   integer i,k,m,igp,icm,j,i1
   real diffa,diffb,adiffa,adiffb
   real torad,todeg
   torad=pi/180.0
   todeg=1.0/torad

   if (gp(1)%day < cm(1)%day) then
      do i=1,ngp
         if (gp(i)%day > cm(1)%day) exit
      enddo
      igp=i
      icm=1
   else
      do i=1,ncm
         if (cm(i)%day > gp(1)%day) exit
      enddo
      icm=i-1
      igp=1
   endif

!   print *,'icm=',icm
!   print *,'igp=',igp

   j=0
   i1=icm
   do i=igp,ngp
      if (gp(i)%day > max(gp(ngp)%day,cm(ncm)%day)) exit
      do k=i1,ncm-1

         diffa=cm(k)%day-gp(i)%day
         diffb=cm(k+1)%day-gp(i)%day
         adiffa=abs(diffa)
         adiffb=abs(diffb)

         if ((diffa <= 0.0).and.(diffb > 0.0).and.(min(adiffa,adiffb) <= 1.0/24.0)) then
            j=j+1
            cm2(j)%u=cm(k)%u+ ((cm(k+1)%u-cm(k)%u)/(cm(k+1)%day-cm(k)%day))*(gp(i)%day-cm(k)%day)
            cm2(j)%v=cm(k)%v+ ((cm(k+1)%v-cm(k)%v)/(cm(k+1)%day-cm(k)%day))*(gp(i)%day-cm(k)%day)
            cm2(j)%speed=sqrt(cm2(j)%u**2+cm2(j)%v**2)
            cm2(j)%dir=atan2(cm2(j)%v,cm2(j)%u)*todeg
            cm2(j)%day=gp(i)%day
            cm2(j)%s=cm(k)%s+ ((cm(k+1)%s-cm(k)%s)/(cm(k+1)%day-cm(k)%day))*(gp(i)%day-cm(k)%day)
            cm2(j)%t=cm(k)%t+ ((cm(k+1)%t-cm(k)%t)/(cm(k+1)%day-cm(k)%day))*(gp(i)%day-cm(k)%day)
            gp2(j)=gp(i)
!            if (abs(diffa)  < abs(diffb)) then
!               cm2(j)=cm(k)
!            else
!               cm2(j)=cm(k+1)
!            endif
            i1=k
         endif
      enddo
   enddo
   nrm=j
   print '(a,i7)','Merging done: Number of points is ',nrm
   if (nrm < 300) then
   !if (nrm < 1000) then
      call system('touch no_merged_points')
      stop 'number of merged points is less than 1000'
   endif

   open(10,file='merged.out',form='formatted',status='unknown',access='sequential')   
      do m=1,nrm
         write(10,'(2(f9.4,2f9.3,4f8.2,tr2))')gp2(m),cm2(m)
      enddo
   close(10)

end subroutine merge
end module m_merge
