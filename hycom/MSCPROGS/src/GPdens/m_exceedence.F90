module m_exceedence
contains

subroutine exceedence(cm2,gp2,nrm,deep)
   use mod_data_new
   !use m_sort
   use netcdf
   use mod_netcdf_ops
   implicit none
   integer       ,  intent(in) ::  nrm
   type(data_new),  intent(in) ::  cm2(nrm),gp2(nrm)
   integer       ,  intent(in) :: deep

   real cmsort(nrm),gpsort(nrm)
   integer sortindx(nrm)
   integer m,index(5)
   real :: tmpr(5), percent(5)
   logical :: ex
   integer :: excid, vid, indx, rlen, varid,ncid,i, ierr
   integer, allocatable :: tmplev(:)

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

   percent(1)=1.0
   percent(2)=10.0
   percent(3)=50.0
   percent(4)=90.0
   percent(5)=99.0
   index(1)=nint(float(nrm)*percent(1)/100.0)
   index(2)=nint(float(nrm)*percent(2)/100.0)
   index(3)=nint(float(nrm)*percent(3)/100.0)
   index(4)=nint(float(nrm)*percent(4)/100.0)
   index(5)=nint(float(nrm)*percent(5)/100.0)


! KAL - added netcdf option - keeps all info in one file
   call ncopencreate(ncfile,ncid)
   call ncinqdefdim (ncid,'e_levels',5,excid)
   call ncinqdefdim (ncid,'vlevel'  ,NF90_UNLIMITED,vid)

   ! Inq and set vertical levels variable (id = vid). Also get vertial index 
   call ncinqdefvertvar(ncid,vid,deep,indx)

   ! Inq&put variable for exceedence levels 
   call ncinqputvar(ncid,'e_levels',(/excid/),percent,(/1/))

   ! Inq&put variable for exceedence - gp
   tmpr=gpsort(index)
   call ncinqputvarslice(ncid,'exceedence_GP',(/excid,vid/),tmpr,(/1,indx/))

   ! Get&put variable for exceedence - cm
   call ncinqputvarslice(ncid,'exceedence_CM',(/excid,vid/),cmsort(index),(/1,indx/))

   ! Get&put variable for exceedence - gp
   call ncinqputvarslice(ncid,'exceedence_ratio',(/excid,vid/),gpsort(index)/cmsort(indeX),(/1,indx/))


   print *
   call ncputatt(ncid,'e_levels','comment','exceedence levels')
   call ncputatt(ncid,'exceedence_CM','comment','exceedence values for CM time series')
   call ncputatt(ncid,'exceedence_GP','comment','exceedence values for GP time series')
   call ncputatt(ncid,'exceedence_ratio','comment','exceedence ratio  GP rel CM')
   
   call ncerr(NF90_close(ncid))


end subroutine exceedence
end module m_exceedence
