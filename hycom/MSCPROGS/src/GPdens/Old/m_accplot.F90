module m_accplot
contains
subroutine accplot(cm2,gp2,nrm,deep)
   use mod_modtype
   implicit none
   integer, parameter :: nn=50

   integer, intent(in) ::  nrm
   type(data_new), intent(in) ::  cm2(nrm),gp2(nrm)
   integer, intent(in) :: deep
   real cmsort(nrm),gpsort(nrm)
   integer m,index(nn)
   real dx

! Sort speeds and compute exceedence plots
   open(10,file='cm.usort')
      do m=1,nrm
         write(10,'(f10.4)')cm2(m)%speed
      enddo
   close(10)
   call system('cat cm.usort | sort > cm.sort; rm cm.usort')

   open(10,file='gp.usort')
      do m=1,nrm
         write(10,'(f10.4)')gp2(m)%speed
      enddo
   close(10)
   call system('cat gp.usort | sort > gp.sort; rm gp.usort')

   open(10,file='cm.sort')
      do m=1,nrm
         read(10,'(f10.4)')cmsort(m)
      enddo
   close(10)

   open(10,file='gp.sort')
      do m=1,nrm
         read(10,'(f10.4)')gpsort(m)
      enddo
   close(10)

   dx=100.0/float(nn)
   do m=1,nn
      index(m)=nint(float(nrm)*(float(m)*dx/100.0))
   enddo

   open(10,file='accplot.dat')
      do m=1,nn
         write(10,'(f8.2,2f10.4)')float(m)*dx/100.0,gpsort(index(m)),cmsort(index(m))
      enddo
   close(10)

end subroutine accplot
end module m_accplot
