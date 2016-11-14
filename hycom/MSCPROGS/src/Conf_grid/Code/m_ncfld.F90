module m_ncfld
contains
subroutine ncfld(fname,nx,ny,lon,lat,deep,ires,weight)
 use netcdf
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: nx,ny,ires
   real, intent(in) :: lon(nx,ny)
   real, intent(in) :: lat(nx,ny)
   real, intent(in) :: deep(nx,ny)
   real, intent(in), optional :: weight(nx,ny)
   !real,slon((nx/ires)+1,(ny/ires)+1),slat((nx/ires)+1,(ny/ires)+1)
   !real,sdeep((nx/ires)+1,(ny/ires)+1),sweight((nx/ires)+1,(ny/ires)+1)
   real :: slon((nx/ires)+1,(ny/ires)+1),slat((nx/ires)+1,(ny/ires)+1)
   real :: sdeep((nx/ires)+1,(ny/ires)+1),sweight((nx/ires)+1,(ny/ires)+1)
   integer i,j,ii,jj,nsx,nsy
   integer :: nxid, nyid, var_id, ncid, ierr
   character(len=80) :: ncfile
   nsx=(nx/ires)+1
   nsy=(ny/ires)+1
   do j=1,ny,ires
      jj=(j/ires)+1
      do i=1,nx,ires
         ii=(i/ires)+1
         slon(ii,jj)=lon(i,j)
         slat(ii,jj)=lat(i,j)
         sdeep(ii,jj)=deep(i,j)
         if (present(weight)) then
          sweight(ii,jj)=weight(i,j)
         endif
      enddo
   enddo
   ncfile=fname//'.nc'
   if (NF90_CREATE(trim(ncfile),NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
        stop '(ncfld)'
   end if
   ierr=NF90_DEF_DIM(ncid,'nx',nsx,nxid)
   ierr=NF90_DEF_DIM(ncid,'ny',nsy,nyid)
   ierr=NF90_DEF_VAR(ncid,'lon',NF90_Float,(/nxid,nyid/),var_id)
   ierr=NF90_ENDDEF(ncid)
   ierr=NF90_PUT_VAR(ncid,var_id,slon)

   ierr=NF90_REDEF(ncid)
   ierr=NF90_DEF_VAR(ncid,'lat',NF90_Float,(/nxid,nyid/),var_id)
   ierr=NF90_ENDDEF(ncid)
   ierr=NF90_PUT_VAR(ncid,var_id,slat)
   
   ierr=NF90_REDEF(ncid)
   ierr=NF90_DEF_VAR(ncid,'weight',NF90_Float,(/nxid,nyid/),var_id)
   ierr=NF90_ENDDEF(ncid)
   ierr=NF90_PUT_VAR(ncid,var_id,sweight)

   ierr=NF90_REDEF(ncid)
   ierr=NF90_DEF_VAR(ncid,'depths',NF90_Float,(/nxid,nyid/),var_id)
   ierr=NF90_PUT_ATT(ncid,var_id,'missing_value',0e4)
   ierr=NF90_ENDDEF(ncid)
   ierr=NF90_PUT_VAR(ncid,var_id,real(sdeep,kind=4))
   ierr=NF90_CLOSE(ncid)
end subroutine ncfld

end module m_ncfld
