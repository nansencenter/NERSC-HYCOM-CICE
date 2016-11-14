module m_tecplot_dump
contains
subroutine tecplot_dump(depths,nx,ny,nz,fld,nfld,normal,n_ncrec,ltecplot,hfile)
   use mod_types
   use netcdf
   use m_handle_err
   use mod_hycomfile_io
   implicit none
   integer, intent(in)            :: nx,ny,nz
   integer, intent(in)            :: nfld
   type(fields),       intent(in) :: fld(nfld)
   logical, intent(in)            :: normal
   integer, intent(in)            :: n_ncrec
   logical, intent(in) :: ltecplot
   real, intent(in) :: depths(nx,ny)
   type(hycomfile), intent(in) :: hfile

   integer k,i,n,l,j, ix,jx
   real, dimension(nx,ny) :: twod1, dplayer
   logical :: levelmask(0:nz)
   logical found
   character(len=20) :: vartitle
   integer :: dimms3d(3)
   integer :: rdimid, jdimid, idimid,varid,ncid
   character(len=2) :: clay
   character(len=1) :: zoneinfo
   real*4, parameter :: undefr4=real(undef,kind=4)

   open(10,file='.zoneinfo')
   read(10,'(t5,a1)')zoneinfo
   close(10)

   if (NF90_OPEN('tmp1.nc',NF90_WRITE,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(restart2netcdf)'
   end if
   call handle_err (nf90_inq_dimid(ncid, 'idim', idimid))
   call handle_err (nf90_inq_dimid(ncid, 'jdim', jdimid))
   call handle_err (nf90_inq_dimid(ncid, 'rdim', rdimid))
   dimms3D=(/idimid,jdimid,rdimid/)

   ! Put time variable (floating-point year, and more?)
   if (zoneinfo=='=') then
      call handle_err(NF90_REDEF(ncid))
      call handle_err(NF90_DEF_VAR(ncid,'fyear',NF90_Float,rdimid,varid))
      call handle_err(NF90_ENDDEF(ncid))
   else
      call handle_err(Nf90_inq_varid(ncid,'fyear', varid))
   end if
   call handle_err(NF90_PUT_VAR(ncid,varid,hfile%fyear,start=(/n_ncrec/)))


   ! Find depth levels  to process
   levelmask=.false.
   do i=1,nfld
      levelmask(fld(i)%laya:fld(i)%layb) = .true.
   end do

   ! Extends layers up to surface
   do i=nfld-1,1,-1
      levelmask(i) = levelmask(i) .or. levelmask(i+1) 
   end do

   do k=0,nz
   if (levelmask(k)) then
      if (k>0) then
         ! Although not necessary for all files - we read layer interfaces
         ! to get layer mask. This way empty layers get undefed in output file
         call HFReadDPField_m(hfile,dplayer,nx,ny,k,1) ! _m always return meters
      end if

                                    
      do i=1,nfld
      if (fld(i)%option .and. (fld(i)%laya <= k) .and. (k <= fld(i)%layb) &
             .and. (.not.fld(i)%vecflag) ) then

         found=.false.
         call HFReadField(hfile,twod1,nx,ny,fld(i)%fextract,k,1)

         if (ltecplot) then
            where (abs(dplayer)<1) twod1=0.
            write(22,'(10(1x,e12.5))') ((twod1(min(l,nx),min(j,ny)),l=1,nx+1),j=1,ny+1)
         end if
         found=.true.

         ! Dump to netcdf file
         write(clay,'(i2.2)') k 
         vartitle=trim(adjustl(fld(i)%fextract))//clay
         !print *,vartitle
         if (zoneinfo=='=') then
            call handle_err(NF90_REDEF(ncid))
            call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,dimms3D,varid))
            call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))
            call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
            call handle_err(NF90_ENDDEF(ncid))
         else
            call handle_err(Nf90_inq_varid(ncid, trim(vartitle), varid))
         end if
         if (is3DVar(hfile,fld(i)%fextract,1)) then
            where (abs(dplayer)<1) twod1=undef
         end if
         where (depths<.1 .or. depths > 1e20) twod1=undef
         call handle_err(NF90_PUT_VAR(ncid,varid,twod1,start=(/1,1,n_ncrec/)))


      endif
      enddo
   end if
   enddo

   call handle_err(NF90_CLOSE(ncid))
end subroutine tecplot_dump
end module m_tecplot_dump
