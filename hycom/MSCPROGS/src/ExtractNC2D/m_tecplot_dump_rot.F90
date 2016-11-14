module m_tecplot_dump_rot
contains
subroutine tecplot_dump_rot(depthsin,nx,ny,nz,fld,nfld,plat,plon,normal,rot,sphere,n_ncrec,  &
                            ltecplot,hfile)
   use mod_types
   use netcdf
   use m_handle_err
   use mod_hycomfile_io
   implicit none
   integer, intent(in)            :: nx,ny,nz
   integer, intent(in)            :: nfld
   type(fields),       intent(in) :: fld(nfld)
   real,    intent(in)            :: depthsin(nx,ny)
   real,    intent(in)            :: plat(0:nx+1,0:ny+1)
   real,    intent(in)            :: plon(0:nx+1,0:ny+1)
   logical, intent(in)            :: normal
   logical, intent(in)            :: rot
   logical, intent(in)            :: sphere
   integer, intent(in) :: n_ncrec
   logical, intent(in) :: ltecplot
   type(hycomfile), intent(in) :: hfile
   logical :: levelmask(0:nz)

   real  :: unorm(0:nx+1,0:ny+1)
   real  :: vnorm(0:nx+1,0:ny+1)
   real  :: urot(0:nx+1,0:ny+1)
   real  :: vrot(0:nx+1,0:ny+1)
   real  :: depths(0:nx+1,0:ny+1)
   real  :: intf(0:nx+1,0:ny+1)
   real  :: tmp(0:nx+1,0:ny+1)
   real  :: oldintf(0:nx+1,0:ny+1)
   real  :: twod(nx,ny)

   integer k,ifld,n,i,j,ix,jx
   character(len=8) char8,vname,uname
   character(len=9) pre
   character(len=20) tst
   real rad
   integer :: varid,idimid,jdimid,rdimid,ncid
   integer :: dimms(3)
   character(len=2) :: clay
   character(len=1) :: zoneinfo
   real*4,parameter :: undefr4=real(undef,kind=4)
   logical :: vecflag

   depths=0.
   depths(1:nx,1:ny)=depthsin
   depths(1:nx,   0)=depths(1:nx, 1)
   depths(1:nx,ny+1)=depths(1:nx,ny)
   depths(   0,1:ny)=depths(1 ,1:ny)
   depths(nx+1,1:ny)=depths(nx,1:ny)

   !print '(a)','#############################################'
   if (normal) print '(a)','Processing index velocities'
   if (sphere) print '(a)','Processing sphere velocities'
   if (rot) print '(a)','Processing east/north velocities'

   rad=4.*atan(1.)/180.

   if (NF90_OPEN('tmp1.nc',NF90_WRITE,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(tecplot_dump_rot)'
   end if
   call handle_err (nf90_inq_dimid(ncid, 'idim', idimid))
   call handle_err (nf90_inq_dimid(ncid, 'jdim', jdimid))
   call handle_err (nf90_inq_dimid(ncid, 'rdim', rdimid))
   dimms=(/idimid,jdimid,rdimid/)

   open(10,file='.zoneinfo')
   read(10,'(t5,a1)')zoneinfo
   close(10)

   ! Find depth levels  to process
   levelmask=.false.
   do i=1,nfld
      levelmask(fld(i)%laya:fld(i)%layb) = .true.
   end do
   do i=nfld-1,1,-1
      levelmask(i) = levelmask(i) .or. levelmask(i+1)
   end do


   oldintf=0.
   intf=0.
   do k=0,nz
   if(levelmask(k)) then

      if (k>0) then
         !call HFReadDPField(hfile,twod,nx,ny,k,1)
         call HFReadDPField_m(hfile,twod,nx,ny,k,1) ! _m always returns meters
         oldintf=intf
         intf=intf+twod
         intf(1:nx,   0)=intf(1:nx, 1)
         intf(1:nx,ny+1)=intf(1:nx,ny)
         intf(   0,1:ny)=intf(1 ,1:ny)
         intf(nx+1,1:ny)=intf(nx,1:ny)
      end if


      do ifld=1,nfld-1
      if (fld(ifld)%vecflag .and. fld(ifld)%option) then
         uname=fld(ifld  )%fextract
         vname=fld(ifld+1)%fextract ; 
         pre=fld(ifld)%vecpost ; 
         !print *,len_trim(pre),'pre:',pre
         !pre=uname

         if ((fld(ifld)%laya <= k) .and. (k <= fld(ifld)%layb)) then
            urot=0.
            vrot=0.
            call HFReadField(hfile,urot(1:nx,1:ny),nx,ny,uname,k,1)
            call HFReadField(hfile,vrot(1:nx,1:ny),nx,ny,vname,k,1)

            ! Unrotated velocities
            unorm=urot
            vnorm=vrot

            ! Rotated velocities
            urot(1:nx,   0)=urot(1:nx, 1)
            urot(1:nx,ny+1)=urot(1:nx,ny)
            vrot(1:nx,   0)=vrot(1:nx, 1)
            vrot(1:nx,ny+1)=vrot(1:nx,ny)
            urot(   0,1:ny)=urot(1 ,1:ny)
            urot(nx+1,1:ny)=urot(nx,1:ny)
            vrot(   0,1:ny)=vrot(1 ,1:ny)
            vrot(nx+1,1:ny)=vrot(nx,1:ny)
            call rotate(urot,vrot,plat,plon,nx+2,ny+2,'m2l')
            do ix=1,nx
               urot(ix,   0)=urot(ix,ny)
               urot(ix,ny+1)=urot(ix,1)
               vrot(ix,   0)=vrot(ix,ny)
               vrot(ix,ny+1)=vrot(ix,1)
            end do
            do jx=1,ny
               urot(0   ,jx)=urot(nx,jx)
               urot(nx+1,jx)=urot( 1,jx)
               vrot(0   ,jx)=vrot(nx,jx)
               vrot(nx+1,jx)=vrot( 1,jx)
            end do

            write(clay,'(i2.2)') k
            !print *,clay

            if (zoneinfo=='=') then
               call handle_err(NF90_REDEF(ncid))
               if (normal) then
                  call handle_err(NF90_DEF_VAR(ncid,trim(uname)//clay,NF90_Float,dimms,varid))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))

                  call handle_err(NF90_DEF_VAR(ncid,trim(vname)//clay,NF90_Float,dimms,varid))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))
               end if

               if (rot) then
                  call handle_err(NF90_DEF_VAR(ncid,'UROT_'//trim(pre)//clay,NF90_Float,dimms,varid))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))

                  call handle_err(NF90_DEF_VAR(ncid,'VROT_'//trim(pre)//clay,NF90_Float,dimms,varid))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))
               end if
               if (sphere) then
                  call handle_err(NF90_DEF_VAR(ncid,'USPH_'//trim(pre)//clay,NF90_Float,dimms,varid))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))

                  call handle_err(NF90_DEF_VAR(ncid,'VSPH_'//trim(pre)//clay,NF90_Float,dimms,varid))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))

                  call handle_err(NF90_DEF_VAR(ncid,'WSPH_'//trim(pre)//clay,NF90_Float,dimms,varid))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
                  call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))
               end if
               call handle_err(NF90_ENDDEF(ncid))
            end if

            if (normal) then
               !write(6,'(a)',advance='yes') '|--->rotated velocities .. '
               if (ltecplot) then
                  write(22,'(10(1x,e12.5))') ((unorm(i,j),i=1,nx+1),j=1,ny+1)
                  write(22,'(10(1x,e12.5))') ((vnorm(i,j),i=1,nx+1),j=1,ny+1)
               end if

               where(abs(depths-oldintf)<1.) unorm=undef
               where(abs(depths-oldintf)<1.) vnorm=undef
               where(depths>1e20) unorm=undef
               where(depths>1e20) vnorm=undef

               call handle_err(nf90_inq_varid(ncid,trim(uname)//clay, varid))
               call handle_err(NF90_PUT_VAR(ncid,varid,unorm(1:nx,1:ny),start=(/1,1,n_ncrec/)))

               call handle_err(nf90_inq_varid(ncid,trim(vname)//clay, varid))
               call handle_err(NF90_PUT_VAR(ncid,varid,vnorm(1:nx,1:ny),start=(/1,1,n_ncrec/)))
            endif

            if (rot) then
               !write(6,'(a)',advance='yes') '|--->rotated velocities .. '
               if (ltecplot) then
                  write(22,'(10(1x,e12.5))') ((urot(i,j),i=1,nx+1),j=1,ny+1)
                  write(22,'(10(1x,e12.5))') ((vrot(i,j),i=1,nx+1),j=1,ny+1)
               end if

               where(abs(depths-oldintf)<1.) urot=undef
               where(abs(depths-oldintf)<1.) vrot=undef
               where(depths>1e20) urot=undef
               where(depths>1e20) vrot=undef

               !print *,'UROT_'//trim(pre)//clay
               call handle_err(nf90_inq_varid(ncid,'UROT_'//trim(pre)//clay, varid))
               call handle_err(NF90_PUT_VAR(ncid,varid,urot(1:nx,1:ny),start=(/1,1,n_ncrec/)))

               !print *,'VROT_'//trim(pre)//clay
               call handle_err(nf90_inq_varid(ncid,'VROT_'//trim(pre)//clay, varid))
               call handle_err(NF90_PUT_VAR(ncid,varid,vrot(1:nx,1:ny),start=(/1,1,n_ncrec/)))
            endif

            if (sphere) then
               !write(6,'(a)',advance='yes') '|--->3D sphere velocities .. '
               if (ltecplot) then
                  write(22,'(10(1x,e12.5))')((-urot(i,j)*sin(plon(i,j)*rad)&
                                 -vrot(i,j)*sin(plat(i,j)*rad)*cos(plon(i,j)*rad)&
                                ,i=1,nx+1),j=1,ny+1)
                  write(22,'(10(1x,e12.5))')(( urot(i,j)*cos(plon(i,j)*rad)&
                                 -vrot(i,j)*sin(plat(i,j)*rad)*sin(plon(i,j)*rad)&
                                ,i=1,nx+1),j=1,ny+1)
                  write(22,'(10(1x,e12.5))')(( vrot(i,j)*cos(plat(i,j)*rad),i=0,nx+1),j=0,ny+1)
               end if

               call handle_err(nf90_inq_varid(ncid,'USPH_'//trim(pre)//clay, varid))
               do jx=0,ny+1
               do ix=0,nx+1
                  tmp(ix,jx)=-urot(ix,jx)*sin(plon(ix,jx)*rad) &
                      -vrot(ix,jx)*sin(plat(ix,jx)*rad)*cos(plon(ix,jx)*rad)
               end do
               end do
               where(abs(depths-oldintf)<1.) tmp=undef
               where(depths>1e20) tmp=undef
               call handle_err(NF90_PUT_VAR(ncid,varid,tmp(1:nx,1:ny),start=(/1,1,n_ncrec/)))

               call handle_err(nf90_inq_varid(ncid,'VSPH_'//trim(pre)//clay, varid))
               do jx=0,ny+1
               do ix=0,nx+1
                  tmp(ix,jx)= urot(ix,jx)*cos(plon(ix,jx)*rad) &
                      -vrot(ix,jx)*sin(plat(ix,jx)*rad)*sin(plon(ix,jx)*rad)
               end do
               end do
               where(abs(depths-oldintf)<1.) tmp=undef
               where(depths>1e20) tmp=undef
               call handle_err(NF90_PUT_VAR(ncid,varid,tmp(1:nx,1:ny),start=(/1,1,n_ncrec/)))

               call handle_err(nf90_inq_varid(ncid,'WSPH_'//trim(pre)//clay, varid))
               do jx=0,ny+1
               do ix=0,nx+1
                  tmp(ix,jx)= vrot(ix,jx)*cos(plat(ix,jx)*rad)
               end do
               end do
               where(abs(depths-oldintf)<1.) tmp=undef
               where(depths>1e20) tmp=undef
               call handle_err(NF90_PUT_VAR(ncid,varid,tmp(1:nx,1:ny),start=(/1,1,n_ncrec/)))
            endif ! sphere


         endif ! layer test

      end if
      enddo

   end if
   enddo
   !print *,'DONE...'
   !print *
   call handle_err(NF90_CLOSE(ncid))
end subroutine tecplot_dump_rot
end module m_tecplot_dump_rot
