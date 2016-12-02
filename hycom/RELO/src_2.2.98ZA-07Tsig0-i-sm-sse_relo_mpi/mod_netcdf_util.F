      module  mod_netcdf_util
      implicit none


      integer function open_nc(
     &   filename,l_arch, c_arch, dim_arch,nfields,kkout,ncid)
      implicit none
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: nfields, kkout
      logical,          intent(in)  :: l_arch  (nfields)
      character(len=*), intent(in)  :: c_arch  (nfields)
      integer,          intent(in)  :: dim_arch(nfields)

      if (mnproc==1) then 
         mystatus = open_nc_master(
     &   filename,l_arch, c_arch, dim_arch,nfields,kkout,ncid)
      end if
      Communicate mystatus to mpi tasks

C --- Put initial data in file

C --- TODO: Possible to use parallel netcdf support here, but must specify start and end.
C --- This avoids the xcaget call
      call xcaget(depths,fld)
      if (mnproc==1) then 
         mystatus=nf90_inq_varid(ncid, "depths", varid)
      end if
      Communicate mystatus to mpi tasks
      Handle mpi status
      if (mnproc==1) then 
         mystatus=nf90_put_var(  ncid, varid,fld)
      end if
      Communicate mystatus to mpi tasks
      Handle mpi status
      call xcaget(plon  ,fld)
      if (mnproc==1) then 
         mystatus=nf90_inq_varid(ncid, "plon", varid)
      end if
      Communicate mystatus to mpi tasks
      Handle mpi status
      if (mnproc==1) then 
         mystatus=nf90_put_var(  ncid, varid,fld)
      end if
      Communicate mystatus to mpi tasks
      Handle mpi status
      call xcaget(plat  ,fld)
      if (mnproc==1) then 
         mystatus=nf90_inq_varid(ncid, "plat", varid)
      end if
      Communicate mystatus to mpi tasks
      Handle mpi status
      if (mnproc==1) then 
         mystatus=nf90_put_var(  ncid, varid,fld)
      end if
      Communicate mystatus to mpi tasks
      Handle mpi status

      end function open_netcdf_file



      integer function open_nc_master(
     &   filename,l_arch, c_arch, dim_arch,nfields,kkout,ncid)
      implicit none
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: nfields, kkout
      logical,          intent(in)  :: l_arch  (nfields)
      character(len=*), intent(in)  :: c_arch  (nfields)
      integer,          intent(in)  :: dim_arch(nfields)
      integer,          intent(out) :: ncid
      integer :: vdim(4),i,j,k,vdim2D(2),vdim2D_plus_time(3),
     &           vdim3D_plus_time(4)
      real, dimension(itdm,jtdm) :: fld
C
C --- Create file, dont overwrite (matches archv .ab file usage)
C --- No chunk size setA
C --- TODO: Possible to use parallel netcdf support here, but must set flags
      open_nc_master=nf90_create(trim(filename),NF90_NOCLOBBER,ncid)
      if (open_nc_master/= nf90_noerr) return
C
C --- Populate file with some useful attributes. Sky is the limit
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"iversn",iversn)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"iexpt" ,iexpt)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"nhybrd",nhybrd)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"nsigma",nsigma)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"dp00"  ,dp00)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"dp00x" ,dp00x)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"dp00f" ,dp00f)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"ds00"  ,ds00)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"ds00x" ,ds00x)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"ds00f" ,ds00f)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"isotop" ,isotop)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"locsig" ,logsig)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"kapref" ,kapref)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"thflag",thflag)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"thbase",thbase)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"vsigma",vsigma)
      if (open_nc_master/= nf90_noerr) return

      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"iniflg",iniflg)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"jerlv0",jerlv0)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"yrflag",yrflag)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"sshflg",sshflg)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"bnstfq",bnstfq)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"nestfq",nestfq)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"baclin",baclin)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_put_att(ncid,NF90_GLOBAL,"batrop",batrop)
      if (open_nc_master/= nf90_noerr) return
C --- TODO. Stop here for now ....

C --- Create dimensions
      open_nc_master=f90_def_dim(ncid,    "i",idm,vdim(1)))
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_dim(ncid,    "j",jdm,vdim(2)))
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_dim(ncid,"sigma",kdm,vdim(3)))
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_dim(ncid,"time" ,NF90_UNLIMITED,vdim(4)))
      if (open_nc_master/= nf90_noerr) return
      vdim2D=(vdim(1:2))
      vdim2D_plus_time=(/vdim(1),vdim(2),vdim(4)/)
      vdim3D_plus_time=vdim

C --- Create variables.  TODO: Put relevant attributes. TODO: Put lon, lat, depth?
      open_nc_master=f90_def_var(ncid,"sigma" ,NF90_REAL,vdim(3),varid) !  Coordinate variable
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_var(ncid,"time"  ,NF90_REAL,vdim(4),varid) !  Coordinate variable
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_var(ncid,"dp0k"  ,NF90_REAL,vdim(3),varid)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_var(ncid,"ds0k"  ,NF90_REAL,vdim(3),varid)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_var(ncid,"depths",NF90_REAL,vdim2D ,varid)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_var(ncid,"plon"  ,NF90_REAL,vdim2D ,varid)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=f90_def_var(ncid,"plat"  ,NF90_REAL,vdim2D ,varid)
      if (open_nc_master/= nf90_noerr) return
C
C --- Create requested variables
      do i=1,nfields
C ------ 3D field
         if (l_arch(i) .and. kkout > 1 .and dim_arch(i) == 3) then
            open_nc_master=nf90_def_var(ncid,trim(c_arch(i)),NF90_REAL,
     &                            vdim3D_plus_time,varid)
            if (open_nc_master/= nf90_noerr) return
         end if
C
C ------ 2D field
         if (l_arch(i) .and. kkout == 1) then
            open_nc_master=nf90_def_var(ncid,trim(c_arch(i)),NF90_REAL,
     &         vdim2D_plus_time,varid)
            if (open_nc_master/= nf90_noerr) return
         end if
      end do
C
      open_nc_master=f90_enddef(ncid)
      if (open_nc_master/= nf90_noerr) return

      open_nc_master=nf90_inq_varid(ncid, "sigma", varid)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=nf90_put_var(  ncid, varid,sigma)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=nf90_inq_varid(ncid, "ds0k", varid)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=nf90_put_var(  ncid, varid,ds0k)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=nf90_inq_varid(ncid, "dp0k", varid)
      if (open_nc_master/= nf90_noerr) return
      open_nc_master=nf90_put_var(  ncid, varid,dp0k)
      if (open_nc_master/= nf90_noerr) return
C
      end function open_nc_master
C
      subroutine write_nc(tfld,ncid,varname,level,time_index)
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: varname
      integer,          intent(in)  :: level
      integer,          intent(in)  :: time_index
      real,             intent(in)  :: 
     &   tfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
C
      real  ::  fld(itdm,jtdm)
      integer :: dimids(NF90_MAX_DIMS)
C
C --- TODO: Possible to use parallel netcdf support here, but must specify start and end.
C --- This avoids the xcaget call
      call xcaget(tfld  ,gfld)
C
      if (mnproc==1) then 
         mystatus = write_nc_master(gfld,ncid,varname,level,time_index)
      end if
      communicate mystatus to all cpus
      if (mystatus /= nf90_noerr) then
         call handle_err(mystatus)
      end if
      end subroutine write_nc


C
      integer function write_nc_master(gfld,ncid,varname,level,time_index)
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: varname
      integer,          intent(in)  :: level
      integer,          intent(in)  :: time_index
      real,             intent(in)  :: gfld(itdm,jtdm)
C ------ Get variable id matching name. 
      write_nc_master = nf90_inq_varid(ncid,trim(varname), varid)
      if (write_netcdf_master/=nf90_noerr) return

C
C ------ Get dimensions of variable
      write_nc_master=nf90_inquire_variable(
     &   ncid, varid, dimids = dimids, ndims = ndims))
      if (write_netcdf_master/=nf90_noerr) return
C
C --- Assume this is a 2D variable with a time index. To be 100% sure we could check
C --- dim names ...
      if (ndims==3) then
         write_nc_master=f90_put_var(ncid, varid,fld,
     &      start=(/time_index,1,1/)))
         if (write_netcdf_master/=nf90_noerr) return
      elseif (ndims==4) then
         write_nc_master=f90_put_var(ncid, varid,fld,
     &      start=(/time_index,level,1,1/)))
         if (write_netcdf_master/=nf90_noerr) return
      end if
      end subroutine
C
      subroutine close_netcdf_file(ncid)
      implicit none
      integer,          intent(in)  :: ncid
      call handle_error(nf90_close(ncid))
      end subroutine
C
      subroutine handle_err(status)
      implicit none
      integer,          intent(in)  :: status
      if (status/=nf90_noerr) then
         write(msg,'(a)') "mod_netcdf_util: "//
     &      trim(nf90_strerror(status))
      end subroutine

      end module


