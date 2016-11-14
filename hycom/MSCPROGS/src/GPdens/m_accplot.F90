module m_accplot
contains


subroutine accplot(cm2,gp2,nrm,deep)
   use mod_data_new
   use mod_netcdf_ops
   !use m_sort
   implicit none
   integer, parameter :: nn=50

   integer, intent(in) ::  nrm
   type(data_new), intent(in) ::  cm2(nrm),gp2(nrm)
   integer, intent(in) :: deep
   real cmsort(nrm),gpsort(nrm)
   integer sortindx(nrm)
   integer m,index(nn),tmp(nn)
   real dx
   integer :: ncid, vid, accid, vertind

   !KAL -- The above usually works but  is unsafe - Sort is lexical - not  numeric
   !KAL -- (besides its ugly)
   gpsort=gp2%speed
   cmsort=cm2%speed
   call sort(nrm,gpsort,sortindx)
   call sort(nrm,cmsort,sortindx)
!   open(10,file='gp.sort')
!   do m=1,nrm
!      write(10,'(f10.4)')gpsort(m)
!   enddo
!   close(10)
!   open(10,file='cm.sort')
!   do m=1,nrm
!      write(10,'(f10.4)')cmsort(m)
!   enddo
!   close(10)

   dx=100.0/float(nn)
   do m=1,nn
      index(m)=nint(float(nrm)*(float(m)*dx/100.0))
   enddo

   do m=1,nn
      tmp(m) = real(m)
   enddo



! KAL - added netcdf option - keeps all info in one file
   call ncopencreate(ncfile,ncid)
   call ncinqdefdim (ncid,'acc_levels',nn,accid)
   call ncinqdefdim (ncid,'vlevel'  ,NF90_UNLIMITED,vid)

   ! Inq and set vertical levels variable (id = vid). Also get vertial index 
   call ncinqdefvertvar(ncid,vid,deep,vertind)

   ! Inq&put variable for exceedence levels 
   call ncinqputvar(ncid,'acc_levels',(/accid/),real(tmp)*dx,(/1/))

   ! Inq&put variable for exceedence - gp
   call ncinqputvarslice(ncid,'acc_GP',(/accid,vid/),gpsort(index),(/1,vertind/))

   ! Inq&put variable for exceedence - cm
   call ncinqputvarslice(ncid,'acc_CM',(/accid,vid/),cmsort(index),(/1,vertind/))


   call ncputatt(ncid,'acc_levels','comment','exceedence levels ')
   call ncputatt(ncid,'acc_CM','comment','exceedence values for CM time series')
   call ncputatt(ncid,'acc_GP','comment','exceedence values for GP time series')

   call ncerr(NF90_close(ncid))

end subroutine accplot
end module m_accplot
