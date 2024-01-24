      module  mod_netcdf_util
      use netcdf
      use mod_xc
      use mod_cb_arrays
      implicit none
      private
      real*4, parameter :: undef=-1e20

#if defined(NERSC_HYCOM_CICE) && defined(NETCDF_ARCHV)
      public :: write_nc_archv, open_nc, write_nc, close_nc, &
                write_current_time_nc
      contains
!
!
! --- Retrieve number of dimensions for this variable
      integer function num_dim_var(cstr)
      implicit none
      character(len=*), intent(in) :: cstr
!
! --- Translate field name to num dims. Defaults to 2
      select case (trim(cstr))
      case ('u-vel.')
         num_dim_var=3
      case ('v-vel.')
         num_dim_var=3
      case ('thknss')
         num_dim_var=3
      case ('temp')
         num_dim_var=3
      case ('salin')
         num_dim_var=3
      case ('t-diff' , 's-diff', 'viscty')
         num_dim_var=3
      case ('ke')
         num_dim_var=3
      case ('density')
         num_dim_var=3
      case ('tracer')
         num_dim_var=4
      case default 
         num_dim_var=2
      end select
      return
      end function
!
!
! --- Handle error by printing and aborting
      subroutine handle_err(myproc,i) 
      implicit none
      integer, intent(in) :: i, myproc
      real :: tmp(1)
      integer :: mystatus
      character(len=512) :: msg
!
! --- uces xcastr to broadcast integer to all tasks. Hacky, but easy
! --- to implement and is sure to follow hycom XC setup
      if (mnproc==myproc) tmp(1) = i
      call xcastr(tmp,myproc)
      mystatus=int(tmp(1))
      if (mystatus/=nf90_noerr) then
         if (mnproc==myproc) then
            write(msg,'(a)') "mod_netcdf_util: "// &
              trim(nf90_strerror(mystatus))
         write(lp,'(a)') trim(msg)
         end if
         call xcstop('(mod_netcdf_util:handle_err)')
         stop '(mod_netcdf_util:handle_err)'
      end if
      end subroutine
!
! --- Routine for opening and defining netcdf  file
      subroutine open_nc( &
         filename,l_arch, c_arch, nfields,kkout,ncid, &
         filetype,averaging)
      implicit none
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: nfields, kkout
      integer,          intent(out)  :: ncid
      logical,          intent(in)  :: l_arch  (nfields)
      character(len=*), intent(in)  :: c_arch  (nfields)
      character(len=*), intent(in)  :: filetype
      real            , intent(in), optional  :: averaging
!
      real, dimension(itdm,jtdm) :: gfld,gip
      integer :: varid, mystatus
      real    :: myave
!
      if (present(averaging)) then
         myave=averaging
      else 
         myave=0.
      end if
!
! --- For parallel netcdf, all members probably need to call this routine
      if (mnproc==1) then 
         mystatus = open_nc_one( &
            filename,l_arch, c_arch, nfields,kkout,ncid,filetype, &
            averaging=myave)
      end if
      call handle_err(1,mystatus)
!
! --- Put static data in file
      call write_nc(depths,ip,ncid,"depths",1,1,1)
      call write_nc(plon  ,ip,ncid,"plon"  ,1,1,1)
      call write_nc(plat  ,ip,ncid,"plat"  ,1,1,1)
      end subroutine open_nc
!
! --- Routine for opening and defining netcdf  file. At the
! --- moment this routine is designed to run on one task
      integer function open_nc_one( &
         filename,l_arch, c_arch, nfields,kkout,ncid, &
         filetype,averaging)
      implicit none
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: nfields, kkout
      logical,          intent(in)  :: l_arch  (nfields)
      character(len=*), intent(in)  :: c_arch  (nfields)
      integer,          intent(out) :: ncid
      character(len=*), intent(in)  :: filetype
      real            , intent(in), optional  :: averaging
      integer :: i,j,k,vdim2D(2),vdim2D_plus_time(3), &
                 vdim3D_plus_time(4), varid,  &
                 vdim3D_plus_tracer_plus_time(5),dimid_kdm, &
                 dimid_itdm, dimid_jtdm, dimid_time, dimid_kkout, &
                 dimid_ntrcr
      real, dimension(itdm,jtdm) :: fld
      real :: myave
      include 'stmt_fns.h'
!
      if (present(averaging)) then
         myave=averaging
      else 
         myave=0.
      end if
!
! --- Create file, dont overwrite (matches archv .ab file usage)
! --- No chunk size setA
! --- TODO: Possible to use parallel netcdf support here, but must set flags
      open_nc_one=nf90_create(trim(filename),NF90_NOCLOBBER,ncid)
      if (open_nc_one/= nf90_noerr) return
!
! --- Populate file with some useful attributes. Sky is the limit
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"iversn",iversn)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"iexpt" ,iexpt)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"nhybrd",nhybrd)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"nsigma",nsigma)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"dp00"  ,dp00)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"dp00x" ,dp00x)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"dp00f" ,dp00f)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"ds00"  ,ds00)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"ds00x" ,ds00x)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"ds00f" ,ds00f)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"isotop" ,isotop)
      if (open_nc_one/= nf90_noerr) return
      if (locsig) then
         open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"locsig" ,1)
      else
         open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"locsig" ,0)
      end if
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"kapref" ,kapref)
      if (open_nc_one/= nf90_noerr) return
      !open_nc_onen=f90_put_att(ncid,NF90_GLOBAL,"thflag",thflag)
      !if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"sigver",sigver)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"thbase",thbase)
      if (open_nc_one/= nf90_noerr) return
      if (vsigma) then
         open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"vsigma" ,1)
      else
         open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"vsigma" ,0)
      end if
      if (open_nc_one/= nf90_noerr) return

      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"iniflg",iniflg)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"jerlv0",jerlv0)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"yrflag",yrflag)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"sshflg",sshflg)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"bnstfq",bnstfq)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"nestfq",nestfq)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"baclin",baclin)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"batrop",batrop)
      if (open_nc_one/= nf90_noerr) return
! --- TODO. Stop here for now, but all blkdat params should be set in
! --- netcdf file
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL,"filetype",filetype)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,NF90_GLOBAL, &
         "averaging_period_in_days",myave)
      if (open_nc_one/= nf90_noerr) return
!
! --- Create dimensions
      open_nc_one=nf90_def_dim(ncid, "itdm",itdm,dimid_itdm)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_def_dim(ncid, "jtdm",jtdm,dimid_jtdm)
      if (open_nc_one/= nf90_noerr) return
!
! --- define sigma which will be used for 3D variables. if kkout > 1 
! --- define sigma_all, which will be used for 1D variables
      open_nc_one=nf90_def_dim(ncid,"sigma",kkout,dimid_kkout)
      dimid_kdm=dimid_kkout
      if (open_nc_one/= nf90_noerr) return
      if (kkout .ne. kdm) then
         open_nc_one=nf90_def_dim(ncid,"sigma_all",kdm,dimid_kdm)
         if (open_nc_one/= nf90_noerr) return
      end if
! 
      if (kkout.gt.1 .and. ntracr > 0 ) then
         open_nc_one=nf90_def_dim(ncid,"tracer" ,ntracr,dimid_ntrcr)
         if (open_nc_one/= nf90_noerr) return
      end if
!
! --- Time dimension
      open_nc_one=nf90_def_dim(ncid,"time" ,NF90_UNLIMITED,dimid_time)
      if (open_nc_one/= nf90_noerr) return
!
! --- Define convenient sets of dimensions
      vdim2D=(/dimid_itdm,dimid_jtdm/)
      vdim2D_plus_time=(/dimid_itdm,dimid_jtdm,dimid_time/)
      vdim3D_plus_time=(/dimid_itdm,dimid_jtdm,dimid_kkout,dimid_time/)
      if (kkout.gt.1 .and. ntracr > 0 ) then
         vdim3D_plus_tracer_plus_time =  &
           (/dimid_itdm,dimid_jtdm,dimid_kkout,dimid_ntrcr,dimid_time/)
      end if
!
! --- Create variables.  TODO: Put relevant attributes.
      open_nc_one=nf90_def_var(ncid,"sigma",NF90_REAL,dimid_kkout,varid)!  Coordinate variable
      if (open_nc_one/= nf90_noerr) return
      if (kkout.ne.kdm) then 
         open_nc_one=nf90_def_var(ncid,"sigma_all", &
              NF90_REAL,dimid_kdm,varid)
         if (open_nc_one/= nf90_noerr) return
      end if
!
      open_nc_one=nf90_def_var(ncid,"time",NF90_DOUBLE,dimid_time,varid)!  Coordinate variable
      if (open_nc_one/= nf90_noerr) return
      if (yrflag.eq.3) then
         open_nc_one=nf90_put_att(ncid,varid,"unit", &
            "seconds since 01-01-1900 00:00:00 UTC")
         if (open_nc_one/= nf90_noerr) return
      end if
!
      open_nc_one=nf90_def_var(ncid,"dp0k"  ,NF90_REAL,dimid_kdm,varid)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_def_var(ncid,"ds0k"  ,NF90_REAL,dimid_kdm,varid)
      if (open_nc_one/= nf90_noerr) return
!
! --- Static 2D Vars 
      open_nc_one=nf90_def_var(ncid,"depths",NF90_REAL,vdim2D ,varid)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,varid,"_FillValue",undef)
      if (open_nc_one/= nf90_noerr) return
!
      open_nc_one=nf90_def_var(ncid,"plon"  ,NF90_REAL,vdim2D ,varid)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,varid,"_FillValue",undef)
      if (open_nc_one/= nf90_noerr) return
!
      open_nc_one=nf90_def_var(ncid,"plat"  ,NF90_REAL,vdim2D ,varid)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_att(ncid,varid,"_FillValue",undef)
      if (open_nc_one/= nf90_noerr) return
!
! --- Create requested variables from c_arch
      do i=1,nfields
! ------ 4D field (tracer) if kkout > 1, ntracr > 0, dim=4
         if (l_arch(i) .and. ntracr > 0 .and. &
            num_dim_var(trim(c_arch(i))) == 4) then
            open_nc_one=nf90_def_var(ncid,trim(c_arch(i)),NF90_REAL, &
                                     vdim3D_plus_tracer_plus_time,varid)
            if (open_nc_one/= nf90_noerr) return
! --------- Set fillvalue
            open_nc_one=nf90_put_att(ncid,varid,"_FillValue",undef)
            if (open_nc_one/= nf90_noerr) return
!
! ------ 3D field if kkout > 1, dim=3
         elseif (l_arch(i) .and. num_dim_var(trim(c_arch(i))) == 3) then
            open_nc_one=nf90_def_var(ncid,trim(c_arch(i)),NF90_REAL, &
                                  vdim3D_plus_time,varid)
            if (open_nc_one/= nf90_noerr) return
! --------- Set fillvalue
            open_nc_one=nf90_put_att(ncid,varid,"_FillValue",undef)
            if (open_nc_one/= nf90_noerr) return
!
! ------ 2D field for the rest ...
         elseif (l_arch(i)) then
            open_nc_one=nf90_def_var(ncid,trim(c_arch(i)),NF90_REAL, &
               vdim2D_plus_time,varid)
            if (open_nc_one/= nf90_noerr) return
! --------- Set fillvalue
            open_nc_one=nf90_put_att(ncid,varid,"_FillValue",undef)
            if (open_nc_one/= nf90_noerr) return
         end if
      end do
!
! --- End of definitions
      open_nc_one=nf90_enddef(ncid)
      if (open_nc_one/= nf90_noerr) return
!
! --- Write some 1D fields. Sigma
      open_nc_one=nf90_inq_varid(ncid, "sigma", varid)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_var(  ncid, varid,sigma(1:kkout))
      if (open_nc_one/= nf90_noerr) return
!
! --- sigma for all layers
      if (kkout.ne.kdm) then 
         open_nc_one=nf90_inq_varid(ncid, "sigma_all", varid)
         if (open_nc_one/= nf90_noerr) return
         open_nc_one=nf90_put_var(  ncid, varid,sigma)
         if (open_nc_one/= nf90_noerr) return
      end if
!
! --- Shallow layer z-thickness
      open_nc_one=nf90_inq_varid(ncid, "ds0k", varid)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_var(  ncid, varid,ds0k)
      if (open_nc_one/= nf90_noerr) return
!
! --- Deep layer z-thickness
      open_nc_one=nf90_inq_varid(ncid, "dp0k", varid)
      if (open_nc_one/= nf90_noerr) return
      open_nc_one=nf90_put_var(  ncid, varid,dp0k)
      if (open_nc_one/= nf90_noerr) return
!
      end function open_nc_one
!
! --- Gets current time from hycom variables and writes it to a netcdf
! --  file
      subroutine write_current_time_nc(ncid,record) 
      implicit none
      integer,          intent(in)  :: ncid,record
      real*8  :: mydtime
      integer :: mystatus,varid
      integer*8 :: nts_day
!
      nts_day = nint(86400.0d0/baclin)
      mydtime=(nstep/nts_day)+mod(nstep,nts_day)*(baclin/86400.0d0)
      mydtime=mydtime*86400.0d0
!
! --- Get variable id matching name. 
      if (mnproc==1) mystatus = nf90_inq_varid(ncid,"time", varid)
      call handle_err(1,mystatus)
!
! --- Write time
      if (mnproc==1) mystatus=nf90_put_var(ncid, varid,mydtime)
      call handle_err(1,mystatus)
      end subroutine
!
!
! --- Writes variables to netcdf file
      subroutine write_nc(tfld,ip,ncid,varname,level,time_index,ktr)
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: varname
      integer,          intent(in)  :: level
      integer,          intent(in)  :: time_index
      integer,          intent(in)  :: ktr
      real,             intent(in)  ::  &
         tfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer,          intent(in)  ::  &
         ip  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
      real  ::  gfld(itdm,jtdm),gip(itdm,jtdm)
      integer :: dimids(NF90_MAX_VAR_DIMS), mystatus
!
! --- TODO: Possible to use parallel netcdf support here, but 
! --- all members must call, and must specify subarea in write call
! --- xcaget not needed in this case
      call xcaget(gfld  ,tfld,1)
! 
! --- Apply mask to land points
      call xcaget(gip   ,real(ip),1)
      where (gip<0.5)  gfld=undef
!
! --- Do write on master
      if (mnproc==1) then 
         mystatus = write_nc_one(gfld,ncid,varname,level,time_index,ktr)
      end if
      call handle_err(1,mystatus)
      end subroutine write_nc
!
!
! --- Routine for writing field to netcdf  file. At the
! --- moment this routine is designed to run on one MPI task
      integer function write_nc_one(gfld,ncid,varname,level, &
         time_index,ktr)
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: varname
      integer,          intent(in)  :: level
      integer,          intent(in)  :: time_index
      integer,          intent(in)  :: ktr
      real,             intent(in)  :: gfld(itdm,jtdm)
!
      integer :: varid, ndims
      integer :: dimids(NF90_MAX_VAR_DIMS)
!
! ------ Get variable id matching name. 
      write_nc_one = nf90_inq_varid(ncid,trim(varname), varid)
      if (write_nc_one/=nf90_noerr) return
!
! ------ Get dimensions of variable
      write_nc_one=nf90_inquire_variable( &
         ncid, varid, dimids = dimids, ndims = ndims)
      if (write_nc_one/=nf90_noerr) return
!
! --- Tracer fields with time dimension
      if (ndims==5) then   ! tracers
         write_nc_one=nf90_put_var(ncid, varid,gfld &
            ,start = (/1,1,level,ktr,time_index/)  )
!
! --- 3D field with time dimension
      elseif (ndims==4) then
         write_nc_one=nf90_put_var(ncid, varid,gfld &
            ,start = (/1,1,level,time_index/)  )
         if (write_nc_one/=nf90_noerr) return
!
! --- 2D field with time dimension
      elseif (ndims==3) then
         write_nc_one=nf90_put_var(ncid, varid,gfld &
            ,start = (/1,1,time_index/)  )
         if (write_nc_one/=nf90_noerr) return
! --- 2D field without time dimension
      elseif (ndims==2) then
         write_nc_one=nf90_put_var(ncid, varid,gfld &
            ,start = (/1,1/)  )
         if (write_nc_one/=nf90_noerr) return
      end if
      end function write_nc_one
!
      subroutine close_nc(ncid)
      implicit none
      integer,          intent(in)  :: ncid
      integer :: mystatus
      if (mnproc==1) then 
         mystatus=close_nc_one(ncid)
      end if
      call handle_err(1,mystatus)
      end subroutine close_nc
!
      integer function close_nc_one(ncid)
      implicit none
      integer,          intent(in)  :: ncid
      integer :: mystatus
      close_nc_one=nf90_close(ncid)
      return
      end function close_nc_one
!
!
! --- Write archv-like netcdf file. Works for archs and archv.
! --- Writes one time step only. NB: Does not use archs.input
      subroutine write_nc_archv(filename,kkout,n)
      implicit none
      character(len=*), intent(in) :: filename
      integer,          intent(in) :: kkout,n
!
      integer, parameter          :: nfields_nc=25
      character*6 :: c_arch_nc(nfields_nc)
      logical     :: l_arch_nc(nfields_nc)
      integer :: k,ktr, ncid
!
      c_arch_nc( 1) = 'montg1'
      c_arch_nc( 2) = 'srfhgt'
      c_arch_nc( 3) = 'steric'
      c_arch_nc( 4) = 'surflx'
      c_arch_nc( 5) = 'salflx'
      c_arch_nc( 6) = 'bldpth'
      c_arch_nc( 7) = 'mldpth'
      c_arch_nc( 8) = 'covice'
      c_arch_nc( 9) = 'thkice'
      c_arch_nc(10) = 'temice'
      c_arch_nc(11) = 'ubtrop'
      c_arch_nc(12) = 'vbtrop'
      c_arch_nc(13) = 'u-vel.'
      c_arch_nc(14) = 'v-vel.'
      c_arch_nc(15) = 'thknss'
      c_arch_nc(16) = 'temp  '
      c_arch_nc(17) = 'salin '
      c_arch_nc(18) = 'surtx '  !output after salflx
      c_arch_nc(19) = 'surty '
      c_arch_nc(20) = 'si_u  '
      c_arch_nc(21) = 'si_v  '
      c_arch_nc(22) = "t-diff"
      c_arch_nc(23) = "s-diff"
      c_arch_nc(24) = "viscty"
      c_arch_nc(25) = "tracer"
!
      l_arch_nc=.true.
      l_arch_nc(22) = difout
      l_arch_nc(23) = difout
      l_arch_nc(24) = difout
      l_arch_nc(25) = ntracr>0
!
      do k=1,nfields_nc
         if (k ==  3) then
            l_arch_nc(k)=sshflg.ne.0
         elseif (k==8 .or. k==9 .or. k==20 .or. k==21 ) then
            l_arch_nc(k)=iceflg.ne.0
         elseif (k >= 22) then
            l_arch_nc(k)=kkout.gt.1 .and. l_arch_nc(k)
         else 
            l_arch_nc(k)=kkout.gt.1 .or. l_arch_nc(k)
         end if
      end do
!
      write(lp,'(a)') "Writing to "//trim(filename)
      call open_nc(trim(filename), l_arch_nc,  &
            c_arch_nc, nfields_nc,kkout,ncid, &
            "HYCOM 3D Archive file", &
            averaging=0.)
      call write_current_time_nc(ncid,1) 
! 
! ------ 2D variables
      if (l_arch_nc(1)) call write_nc(montg1,ip,ncid, &
            trim(c_arch_nc(1)),1,1,1)
!
      if (l_arch_nc(2)) call write_nc(srfhgt,ip,ncid, &
            trim(c_arch_nc(2)),1,1,1)
!
      if (l_arch_nc(3)) call write_nc(steric,ip,ncid, &
            trim(c_arch_nc(3)),1,1,1)
!
      if (l_arch_nc(4)) call write_nc(surflx,ip,ncid, &
            trim(c_arch_nc(4)),1,1,1)
!
      if (l_arch_nc(5)) call write_nc(salflx,ip,ncid, &
            trim(c_arch_nc(5)),1,1,1)
!
      if (l_arch_nc(18)) call write_nc(surtx ,ip,ncid, &
            trim(c_arch_nc(18)),1,1,1)
!
      if (l_arch_nc(19)) call write_nc(surty ,ip,ncid, &
            trim(c_arch_nc(19)),1,1,1)
!
      if (l_arch_nc(7)) call write_nc(dpmixl(1-nbdy,1-nbdy,n), &
            ip,ncid,trim(c_arch_nc(7)),1,1,1)
!
      if (l_arch_nc(8)) call write_nc(covice,ip,ncid, &
            trim(c_arch_nc(8)),1,1,1)
!
      if (l_arch_nc(9)) call write_nc(thkice,ip,ncid, &
            trim(c_arch_nc(9)),1,1,1)
!
      if (l_arch_nc(20)) call write_nc(si_u  ,ip,ncid, &
            trim(c_arch_nc(20)),1,1,1)
!
      if (l_arch_nc(21)) call write_nc(si_v  ,ip,ncid, &
            trim(c_arch_nc(21)),1,1,1)
!
      if (l_arch_nc(11)) call write_nc(ubavg(1-nbdy,1-nbdy,n), &
            iu,ncid,trim(c_arch_nc(11)),1,1,1)
!
      if (l_arch_nc(12)) call write_nc(vbavg(1-nbdy,1-nbdy,n), &
            iv,ncid,trim(c_arch_nc(12)),1,1,1)
!
!        call write_nc(displd_mn, ip,ncid,"disp_ld",1,1,1)
!        call write_nc(dispqd_mn, ip,ncid,"disp_qd",1,1,1)
!        call write_nc(tidepg_mn,ip,ncid,"tidepg_mn",1,1,1)
!
! ------ 3D variables
      do k=1,kkout
      if (l_arch_nc(13)) call write_nc(u(1-nbdy,1-nbdy,k,n), &
            iu,ncid,trim(c_arch_nc(13)),k,1,1)
!
      if (l_arch_nc(14)) call write_nc(v(1-nbdy,1-nbdy,k,n), &
            iv,ncid,trim(c_arch_nc(14)),k,1,1)
!
      if (l_arch_nc(15)) call write_nc(dp(1-nbdy,1-nbdy,k,n), &
            ip,ncid,trim(c_arch_nc(15)),k,1,1)
!
      if (l_arch_nc(16)) call write_nc(temp(1-nbdy,1-nbdy,k,n), &
            ip,ncid,trim(c_arch_nc(16)),k,1,1)
!
      if (l_arch_nc(17)) call write_nc(saln(1-nbdy,1-nbdy,k,n), &
            ip,ncid,trim(c_arch_nc(17)),k,1,1)
!
! ------ Tracer fields
      do ktr= 1,ntracr
       if (l_arch_nc(25)) call write_nc(tracer(1-nbdy,1-nbdy,k,n,ktr), &
            ip,ncid,trim(c_arch_nc(25)),k,1,ktr)
      end do
!
      if (l_arch_nc(24)) call write_nc(vcty(1-nbdy,1-nbdy,k+1), &
            ip,ncid,trim(c_arch_nc(24)),k,1,1)
!
      if (l_arch_nc(22))  call write_nc(dift(1-nbdy,1-nbdy,k+1), &
            ip,ncid,trim(c_arch_nc(22)),k,1,1)
!
      if (l_arch_nc(23)) call write_nc(difs(1-nbdy,1-nbdy,k+1), &
            ip,ncid,trim(c_arch_nc(23)),k,1,1)
      end do
      call close_nc(ncid)
      end subroutine write_nc_archv
!
#endif /* defined(NERSC_HYCOM_CICE) && defined(NETCDF_ARCHV)*/
      end module


