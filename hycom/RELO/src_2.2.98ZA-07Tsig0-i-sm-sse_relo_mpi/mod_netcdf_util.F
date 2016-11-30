      module  mod_netcdf_util
      implicit none


      subroutine open_netcdf_file(filename,l_arch, c_arch, dim_arch, 
     &                            nfields,kkout,ncid)
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

C --- Create file, dont overwrite (matches archv .ab file usage)
C --- No chunk size setA
C --- TODO: Possible to use parallel netcdf support here, but must set flags
      call handle_error(nf90_create(trim(filename),NF90_NOCLOBBER,ncid))
C
C --- Populate file with some useful attributes. Sky is the limit
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"iversn",iversn)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"iexpt" ,iexpt)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"nhybrd",nhybrd)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"nsigma",nsigma)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"dp00"  ,dp00)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"dp00x" ,dp00x)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"dp00f" ,dp00f)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"ds00"  ,ds00)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"ds00x" ,ds00x)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"ds00f" ,ds00f)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"isotop" ,isotop)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"locsig" ,logsig)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"kapref" ,kapref)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"thflag",thflag)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"thbase",thbase)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"vsigma",vsigma)

      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"iniflg",iniflg)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"jerlv0",jerlv0)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"yrflag",yrflag)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"sshflg",sshflg)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"bnstfq",bnstfq)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"nestfq",nestfq)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"baclin",baclin)
      call handle_error(nf90_put_att(ncid,NF90_GLOBAL,"batrop",batrop)
C --- TODO. Stop here for now ....

C --- Create dimensions
      call handle_error(nf90_def_dim(ncid,    "i",idm,vdim(1)))
      call handle_error(nf90_def_dim(ncid,    "j",jdm,vdim(2)))
      call handle_error(nf90_def_dim(ncid,"sigma",kdm,vdim(3)))
      call handle_error(nf90_def_dim(ncid,"time" ,
     &                  NF90_UNLIMITED,vdim(4)))
      vdim2D=(vdim(1:2))
      vdim2D_plus_time=(/vdim(1),vdim(2),vdim(4)/)
      vdim3D_plus_time=vdim

C --- Create variables.  TODO: Put relevant attributes. TODO: Put lon, lat, depth?
#if defined(REAL4) 
      call handle_error(nf90_def_var(ncid,"sigma" ,NF90_REAL,vdim(3))  ,
     &                                varid) !  Coordinate variable
      call handle_error(nf90_def_var(ncid,"time"  ,NF90_REAL,vdim(4))  ,
     &                                varid) !  Coordinate variable
      call handle_error(nf90_def_var(ncid,"dp0k"  ,NF90_REAL,vdim(3))  ,
     &                                varid)
      call handle_error(nf90_def_var(ncid,"ds0k"  ,NF90_REAL,vdim(3))  ,
     &                               varid)
      call handle_error(nf90_def_var(ncid,"depths",NF90_REAL,vdim2D    ,
     &                               varid)
      call handle_error(nf90_def_var(ncid,"plon"  ,NF90_REAL,vdim2D    ,
     &                               varid)
      call handle_error(nf90_def_var(ncid,"plat"  ,NF90_REAL,vdim2D    ,
     &                               varid)
#else
      call handle_error(nf90_def_var(ncid,"sigma" ,NF90_DOUBLE,vdim(3)),
     &                               varid) !  Coordinate variable
      call handle_error(nf90_def_var(ncid,"time"  ,NF90_DOUBLE,vdim(4)),
     &                               varid) !  Coordinate variable
      call handle_error(nf90_def_var(ncid,"dp0k"  ,NF90_DOUBLE,vdim(3)),
     &                               varid)
      call handle_error(nf90_def_var(ncid,"ds0k"  ,NF90_DOUBLE,vdim(3)),
     &                               varid)
      call handle_error(nf90_def_var(ncid,"depths",NF90_DOUBLE,vdim2D  ,
     &                               varid)
      call handle_error(nf90_def_var(ncid,"plon"  ,NF90_DOUBLE,vdim2D  ,
     &                               varid)
      call handle_error(nf90_def_var(ncid,"plat"  ,NF90_DOUBLE,vdim2D  ,
     &                               varid)
#endif

C --- Create requested variables
      do i=1,nfields
C ------ 3D field
         if (l_arch(i) .and. kkout > 1 .and dim_arch(i) == 3) then
#if defined(REAL4) 
            call handle_error(nf90_def_var(
     &         ncid,trim(c_arch(i)),NF90_REAL  ,vdim3D_plus_time,varid)
#else
            call handle_error(nf90_def_var(
     &         ncid,trim(c_arch(i)),NF90_DOUBLE,vdim3D_plus_time,varid)
#endif
         end if
C
C ------ 2D field
         if (l_arch(i) .and. kkout == 1) then
#if defined(REAL4) 
            call handle_error(nf90_def_var(
     &         ncid,trim(c_arch(i)),NF90_REAL  ,vdim2D_plus_time,varid)
#else
            call handle_error(nf90_def_var(
     &         ncid,trim(c_arch(i)),NF90_DOUBLE,vdim2D_plus_time,varid)
#endif
         end if
      end do
C
      call handle_error(nf90_enddef(ncid)
C

C --- Put initial data in file
      call handle_error(nf90_inq_varid(ncid, "sigma", varid)
      call handle_error(nf90_put_var(  ncid, varid,sigma)
      call handle_error(nf90_inq_varid(ncid, "ds0k", varid)
      call handle_error(nf90_put_var(  ncid, varid,ds0k)
      call handle_error(nf90_inq_varid(ncid, "dp0k", varid)
      call handle_error(nf90_put_var(  ncid, varid,dp0k)
C --- TODO: Possible to use parallel netcdf support here, but must specify start and end.
C --- This avoids the xcaget call
      call xcaget(depths,fld)
      call handle_error(nf90_inq_varid(ncid, "depths", varid)
      call handle_error(nf90_put_var(  ncid, varid,fld)
      call xcaget(plon  ,fld)
      call handle_error(nf90_inq_varid(ncid, "plon", varid)
      call handle_error(nf90_put_var(  ncid, varid,fld)
      call xcaget(plat  ,fld)
      call handle_error(nf90_inq_varid(ncid, "plat", varid)
      call handle_error(nf90_put_var(  ncid, varid,fld)
      end subroutine open_netcdf_file
C
C
      subroutine write_netcdf_file(tfld,ncid,varname,level,time_index)
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: varname
      integer,          intent(in)  :: level
      integer,          intent(in)  :: time_index
      real,             intent(in)  :: 
     &   tfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
C
      integer :: dimids(NF90_MAX_DIMS)
C
C --- TODO: Possible to use parallel netcdf support here, but must specify start and end.
C --- This avoids the xcaget call
      call xcaget(plon  ,fld)
C
C --- Get variable id matching name. 
      call handle_error(nf90_inq_varid(ncid,trim(varname), varid)
C
C --- Get dimensions of variable
      call handle_err(nf90_inquire_variable(
     &   ncid, varid, dimids = dimids, ndims = ndims))
C
C --- Assume this is a 2D variable with a time index. To be 100% sure we could check
C --- dim names ...
      if (ndims==3) then
         call handle_error(nf90_put_var(ncid, varid,fld,
     &      start=(/time_index,1,1/)))
      elseif (ndims==4) then
         call handle_error(nf90_put_var(ncid, varid,fld,
     &      start=(/time_index,level,1,1/)))
      end if
      end subroutine


      subroutine close_netcdf_file(ncid)
      implicit none
      integer,          intent(in)  :: ncid
      call handle_error(nf90_close(ncid))

