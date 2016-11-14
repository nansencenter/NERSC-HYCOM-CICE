module m_tecplot_zoneinfo
contains
subroutine tecplot_zoneinfo(nx,ny,nz,plot_time,qlon,qlat,dx,dy,depths,sphere,ltecplot)
   use netcdf
   use m_handle_err
   implicit none
   integer, intent(in) :: nx,ny,nz
   character(len=24), intent(in) ::  plot_time
   real, intent(in) :: qlat(nx+1,ny+1)
   real, intent(in) :: qlon(nx+1,ny+1)
   real, intent(in) :: dx(nx,ny)
   real, intent(in) :: dy(nx,ny)
   real, intent(in) :: depths(nx,ny)
   logical, intent(in) :: sphere
   logical, intent(in) :: ltecplot

   character(len=1) zoneinfo
   character(len=8)   :: dateinfo
   character(len=20) :: vartitle
   integer i,j
   real rad
   real*4, parameter :: undef = -1e14
   real :: tmp(nx,ny)


   integer :: lon_id, lat_id, depth_id, dx_id, dy_id
   integer :: xsphere_id, ysphere_id, zsphere_id,ncid
   integer :: idim_id, jdim_id, rdim_id, dimms3D(3), dimms2D(2)

   open(10,file='.zoneinfo')
      read(10,'(t5,a1)')zoneinfo
   close(10)

   rad=4.*atan(1.)/180.
   if (ltecplot) then
      open(22,file='outfile',status='unknown')
      write(22,103)plot_time,nx+1,ny+1,1 
   end if


   if (zoneinfo == '=') then


      ! Tecplot file
      if (ltecplot) then
         write(22,105)((i,i=1,nx+1),j=1,ny+1)
         write(22,105)((j,i=1,nx+1),j=1,ny+1)
         write(22,104)((qlon  (i,j),i=1,nx+1),j=1,ny+1)
         write(22,104)((qlat  (i,j),i=1,nx+1),j=1,ny+1)
         write(22,104)((-depths(min(i,nx),min(j,ny)),i=1,nx+1),j=1,ny+1)
      end if


      ! netcdf file
      ! This also creates the netcdf file
      if (NF90_CREATE('tmp1.nc',NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(restart2netcdf)'
      end if

      call handle_err(NF90_DEF_DIM(ncid,'idim',nx,idim_id))
      call handle_err(NF90_DEF_DIM(ncid,'jdim',ny,jdim_id))
      call handle_err(NF90_DEF_DIM(ncid,'rdim',nf90_unlimited,rdim_id))
      dimms2D=(/idim_id,jdim_id/)
      dimms3D=(/idim_id,jdim_id,rdim_id/)
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title','tecconv generated fields'))
      call date_and_time(date=dateinfo)
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Dataset_generation_date', &
         dateinfo(1:8)))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'history', &
         'Created '//dateinfo(1:8)//' by program tecconv'))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions','-'))

      ! Create Longitude
      vartitle='longitude'
      !print *,vartitle
      call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_Float,dimms2D,lon_id))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'units','degrees_east'))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'valid_range',(/-180.0,180.0/)))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'standard_name','longitude'))

      ! Latitude
      vartitle= 'latitude'
      call handle_err(NF90_DEF_VAR(ncid,'latitude',NF90_Float,dimms2D,lat_id))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'units','degrees_north'))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'valid_range',(/-90.0, 90.0/)))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'standard_name','latitude'))

      ! Latitude
      vartitle= 'dx'
      call handle_err(NF90_DEF_VAR(ncid,'dx',NF90_Float,dimms2D,dx_id))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'long_name','p-cell delta x'))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'units','m'))

      ! Latitude
      vartitle= 'dy'
      call handle_err(NF90_DEF_VAR(ncid,'dy',NF90_Float,dimms2D,dy_id))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'long_name','p-cell delta y'))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'units','m'))

      ! Latitude
      vartitle= 'depth'
      !print *,vartitle
      call handle_err(NF90_DEF_VAR(ncid,'depth',NF90_Float,dimms2D,depth_id))
      call handle_err(NF90_PUT_ATT(ncid,depth_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,depth_id,'units','meters'))
      call handle_err(NF90_PUT_ATT(ncid,depth_id,'valid_range',(/0.,14000./)))
      call handle_err(NF90_PUT_ATT(ncid,depth_id,'missing_value',undef))
      call handle_err(NF90_PUT_ATT(ncid,depth_id,'_FillValue',undef))
      call handle_err(NF90_PUT_ATT(ncid,depth_id,'standard_name','latitude'))
      !call handle_err(NF90_PUT_ATT(ncid,depth_id,'axis','Y'))


      ! End define mode
      call handle_err(NF90_ENDDEF(ncid))

      ! Put vars into netcdf files
      call handle_err(NF90_PUT_VAR(ncid,lon_id,qlon(1:nx,1:ny)))
      call handle_err(NF90_PUT_VAR(ncid,lat_id,qlat(1:nx,1:ny)))
      call handle_err(NF90_PUT_VAR(ncid,dx_id,dx))
      call handle_err(NF90_PUT_VAR(ncid,dy_id,dy))
      tmp=depths
      where(depths<1.) tmp=undef
      where(depths>1e20) tmp=undef
      call handle_err(NF90_PUT_VAR(ncid,depth_id,tmp))
      call handle_err(NF90_CLOSE(ncid))


   else
      if (ltecplot) then
         write(22,'(a)',advance='no')'D=(1,2,3,4,5'
      endif
   endif

   if (sphere) then
      if (zoneinfo == '=') then
         if (ltecplot) then
            write(22,104)((cos(qlat(i,j)*rad)*cos(qlon(i,j)*rad),i=1,nx+1),j=1,ny+1)
            write(22,104)((cos(qlat(i,j)*rad)*sin(qlon(i,j)*rad),i=1,nx+1),j=1,ny+1)
            write(22,104)((sin(qlat(i,j)*rad),i=1,nx+1),j=1,ny+1)
         end if

         if (NF90_OPEN('tmp1.nc',NF90_WRITE,ncid) /= NF90_NOERR) then
            print *,'An error occured when opening the netcdf file'
            stop '(tecplot_zoneinfo)'
         end if

         ! Redefine (add new data sets)
         call handle_err(NF90_REDEF(ncid))

         vartitle='Xsphere'
         !print *,vartitle
         call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,dimms2D,xsphere_id))
         call handle_err(NF90_PUT_ATT(ncid,xsphere_id,'long_name',trim(vartitle)))

         vartitle='Ysphere'
         !print *,vartitle
         call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,dimms2D,ysphere_id))
         call handle_err(NF90_PUT_ATT(ncid,ysphere_id,'long_name',trim(vartitle)))

         vartitle='Zsphere'
         !print *,vartitle
         call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,dimms2D,zsphere_id))
         call handle_err(NF90_PUT_ATT(ncid,zsphere_id,'long_name',trim(vartitle)))

         call handle_err(NF90_ENDDEF(ncid))
         call handle_err(NF90_PUT_VAR(ncid,xsphere_id,(cos(qlat(1:nx,1:ny)*rad)*cos(qlon(1:nx,1:ny)*rad))))
         call handle_err(NF90_PUT_VAR(ncid,ysphere_id,(cos(qlat(1:nx,1:ny)*rad)*sin(qlon(1:nx,1:ny)*rad))))
         call handle_err(NF90_PUT_VAR(ncid,zsphere_id,(sin(qlat(1:nx,1:ny)*rad))))
         call handle_err(NF90_CLOSE(ncid))

      else
         if (ltecplot) then
            write(22,'(a)',advance='no')',6,7,8)'
         end if
      endif
   else
      if (zoneinfo /= '=') then
         if (ltecplot) then
            write(22,'(a)')')'
         end if
      endif
   endif        

  103 format('ZONE T="',a24,'", F=BLOCK,',' I=',I3,',J=',I3,',K=',I3)
  104 format(10(1x,e12.5))
  105 format(30I4)
end subroutine tecplot_zoneinfo
end module m_tecplot_zoneinfo
