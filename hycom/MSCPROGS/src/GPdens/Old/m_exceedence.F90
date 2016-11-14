module m_exceedence
contains
subroutine exceedence(cm2,gp2,nrm,deep)
   use mod_modtype
   implicit none
   integer, intent(in) ::  nrm
   type(data_new), intent(in) ::  cm2(nrm),gp2(nrm)
   integer, intent(in) :: deep
   real cmsort(nrm),gpsort(nrm)
   integer m,index(5)

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

   index(1)=nint(float(nrm)*01.0/100.0)
   index(2)=nint(float(nrm)*10.0/100.0)
   index(3)=nint(float(nrm)*50.0/100.0)
   index(4)=nint(float(nrm)*90.0/100.0)
   index(5)=nint(float(nrm)*99.0/100.0)

   open(10,file='gpexceedence.dat')
      write(10,'(i5,tr1)',advance='no')deep
      do m=1,5
         write(10,'(f10.4,tr1)',advance='no')gpsort(index(m))
      enddo
      write(10,'(f10.4)')gpsort(index(5))
   close(10)

   open(10,file='cmexceedence.dat')
      write(10,'(i5,tr1)',advance='no')deep
      do m=1,5
         write(10,'(f10.4,tr1)',advance='no')cmsort(index(m))
      enddo
      write(10,'(f10.4)')cmsort(index(5))
   close(10)

   open(10,file='exceed_ratio.dat')
      write(10,'(i5,tr1)',advance='no')deep
      do m=1,5
         write(10,'(f10.4,tr1)',advance='no')gpsort(index(m))/cmsort(index(m))
      enddo
      write(10,'(f10.4)')gpsort(index(5))/cmsort(index(5))
   close(10)


end subroutine exceedence
end module m_exceedence
