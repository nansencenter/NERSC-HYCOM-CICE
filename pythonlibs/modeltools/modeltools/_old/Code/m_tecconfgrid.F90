module m_tecconfgrid
contains
subroutine tecconfgrid(depths,dx,dy,modlon,modlat,dbodc,detopo5,dibcao,dgebco,dconman,nx,ny)
! --- write in tecplot format
   use netcdf
   implicit none
   integer, intent(in) :: nx,ny
   real,    intent(in) :: depths(nx,ny)        ! final model depths
   real,    intent(in) :: dx(nx,ny)            ! grid spacing
   real,    intent(in) :: dy(nx,ny)            ! grid spacing
   real,    intent(in) :: modlon(nx,ny)        ! model longitudes (p-point)
   real,    intent(in) :: modlat(nx,ny)        ! model latitudes
   real,    intent(in) :: dbodc(nx,ny)         ! model depths based on BODC data
   real,    intent(in) :: detopo5(nx,ny)       ! model depths based on ETOPO5 data
   real,    intent(in) :: dibcao(nx,ny)        ! model depths based on IBCAO data
   real,    intent(in) :: dgebco(nx,ny)        ! model depths based on GEBCO data
   real,    intent(in) :: dconman(nx,ny)        ! model depths based on GEBCO data

   real rad
   integer i,j
   character(len=80) :: ncfile
   integer :: idmid, jdmid, var_id, ncid, ierr

      rad=4.*atan(1.)/180.
      open(10,file='confgrid.dat')
      WRITE(10,*)'TITLE="Conformal grid"'
      write(10,'(A)')  &
      'VARIABLES = "i-index""j-index""Longitude""Latitude""depths""X""Y""Z"&
                   &"etopo""bodc""ibcao""gebco""conman""dx""dy"'

      WRITE(10,*)'ZONE I=',nx,',J=',ny,',F=BLOCK'

      WRITE(10,101)((i,i=1,nx),j=1,ny)
      WRITE(10,101)((j,i=1,nx),j=1,ny)

      write(10,402) ((modlon(i,j),i=1,nx),j=1,ny)
      write(10,402) ((modlat(i,j),i=1,nx),j=1,ny)

      write(10,402) ((depths(i,j),i=1,nx),j=1,ny)

      write(10,402) ((cos(modlat(i,j)*rad)*cos(modlon(i,j)*rad),i=1,nx),j=1,ny)
      write(10,402) ((cos(modlat(i,j)*rad)*sin(modlon(i,j)*rad),i=1,nx),j=1,ny)
      write(10,402) ((sin(modlat(i,j)*rad),i=1,nx),j=1,ny)

      write(10,402) ((detopo5(i,j),i=1,nx),j=1,ny)
      write(10,402) ((dbodc(i,j),i=1,nx),j=1,ny)
      write(10,402) ((dibcao(i,j),i=1,nx),j=1,ny)
      write(10,402) ((dgebco(i,j),i=1,nx),j=1,ny)
      write(10,402) ((dconman(i,j),i=1,nx),j=1,ny)
      write(10,402) ((dx(i,j),i=1,nx),j=1,ny)
      write(10,402) ((dy(i,j),i=1,nx),j=1,ny)
      close(10)


  101 FORMAT(30i5)
  402 format(15g13.5)
  403 format('ZONE T=".", F=BLOCK, I=',I3,',J=',I3)


  ! Also dump netcdf file
      ncfile='confgrid.nc'
      print *,'Dumping to netcdf file ',trim(ncfile)
      if (NF90_CREATE(trim(ncfile),NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(obsstats)'
      end if

      print *,'max/mean depth',maxval(depths),sum(depths)/(nx*ny)


      ierr=NF90_DEF_DIM(ncid,'idm',nx,idmid)
      ierr=NF90_DEF_DIM(ncid,'jdm',ny,jdmid)

      ierr=NF90_DEF_VAR(ncid,'modlon',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,modlon)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'modlat',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,modlat)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'depths',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_PUT_ATT(ncid,var_id,'missing_value',0e4)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,real(depths,kind=4))

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'detopo5',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_PUT_ATT(ncid,var_id,'missing_value',0e4)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,detopo5)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'dbodc',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,dbodc)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'dibcao',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,dibcao)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'dgebco',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,dgebco)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'dconman',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,dconman)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'dx',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,dx)

      ierr=NF90_REDEF(ncid)
      ierr=NF90_DEF_VAR(ncid,'dy',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,dy)

      ierr=NF90_CLOSE(ncid)

end subroutine tecconfgrid
end module m_tecconfgrid

