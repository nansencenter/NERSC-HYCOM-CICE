module mod_netcdf_ops
use netcdf

   interface ncinqputvar
     module procedure      &
          ncinqputvar_r1d, &
          ncinqputvar_r2d
   end interface

   interface ncinqputvarslice
     module procedure            &
          ncinqputvarslice_r2d , &
          ncinqputvarslice_r3d
   end interface

   interface ncputatt
     module procedure            &
          ncputatt_r , &
          ncputatt_c
   end interface

   character(len=100) :: ncfile
! Some operations which are common among the gpdens subroutines are deefined here
contains



   ! Inquire if file exists - if it does open, otherwise create it
   subroutine ncopencreate(ncfile,ncid)
   implicit none
   character(len=*), intent(in) :: ncfile
   integer, intent(out) :: ncid
   logical :: ex

   !print *,'Diag dumped to netcdf file '//trim(ncfile)
   inquire(exist=ex,file=ncfile)
   if (.not. ex) then
      call ncerr(NF90_create(trim(ncfile),NF90_CLOBBER,ncid))
   else
      call ncerr(NF90_open(trim(ncfile),NF90_WRITE,ncid))
      call ncerr(NF90_redef(ncid))
   end if
   end subroutine


   ! First inquire for dimension - if unsuccessful define dimension
   subroutine ncinqdefdim(ncid,dname,dsize,dimid)
   implicit none
   character(len=*), intent(in)  :: dname
   integer         , intent(in)  :: ncid,dsize
   integer         , intent(out) :: dimid
   integer :: ierr
   ierr = nf90_inq_dimid(ncid, trim(dname), dimid)
   if (ierr/=NF90_NOERR) call ncerr(NF90_def_dim(ncid,dname,dsize , dimid))
   end subroutine


   ! First inquire for then retreve vertical  variable id. 
   ! If it does not exist - define it
   !
   ! On exit irec is vertical index to put data in , vid is vert var id
   ! netcdf file must be in define mode on entry. On exit netcdf file is
   ! in define mode
   subroutine ncinqdefvertvar(ncid,vid,depth,irec)
   implicit none
   integer         , intent(in)  :: ncid,depth,vid
   integer         , intent(out) :: irec

   integer :: ierr,rlen, indx, i,varid
   integer, allocatable :: tmplev(:)

   ! Get length of (vertical) dimension
   call ncerr(nf90_inquire_dimension(ncid, vid, len=rlen))

   ! Get vertical levels
   ierr = nf90_inq_varid(ncid, 'vlevel', varid)
   if (ierr==NF90_NOERR .and. rlen>0) then
      allocate(tmplev(rlen))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_get_var(ncid,varid,tmplev))
      call ncerr(NF90_redef(ncid))
      indx=-1
      do i=1,rlen
         if (tmplev(i)==depth) indx=i
      end do
   else if (ierr==NF90_NOERR ) then
      indx=-1
   else
      indx=-1
      rlen=0
      call ncerr(NF90_def_var(ncid,'vlevel',NF90_INT,(/vid/),varid)) 
   end if
   if (indx==-1) then
      irec=rlen+1
   else
      irec=indx
   end if
   call ncerr(NF90_enddef(ncid))
   call ncerr(NF90_put_var(ncid,varid,depth,start=(/irec/)))
   call ncerr(NF90_redef(ncid))
   end subroutine


 

   ! Put an entire or partial real 1-d data array into variable vname
   subroutine ncinqputvar_r1d(ncid,vname,dimid,var,start)
   implicit none
   character(len=*), intent(in)  ::  vname
   integer, intent(in)           ::  ncid, dimid(1)
   real   , intent(in)           :: var(:)
   integer, intent(in)           :: start(1)
   integer :: ierr, varid
   ierr = nf90_inq_varid(ncid, vname, varid)
   if (ierr/=NF90_NOERR) call ncerr(NF90_def_var(ncid,vname,NF90_FLOAT,dimid,varid)) 
   call ncerr(NF90_enddef(ncid))
   call ncerr(NF90_put_var(ncid,varid,var,start=start))
   call ncerr(NF90_redef(ncid))
   end subroutine


   ! Put an entire or partial real 1-d data array into variable vname
   subroutine ncinqputvar_r2d(ncid,vname,dimid,var,start)
   implicit none
   character(len=*), intent(in)  ::  vname
   integer, intent(in)           ::  ncid, dimid(2)
   real   , intent(in)           :: var(:,:)
   integer, intent(in)           :: start(2)
   integer :: ierr, varid
   ierr = nf90_inq_varid(ncid, vname, varid)
   if (ierr/=NF90_NOERR) call ncerr(NF90_def_var(ncid,vname,NF90_FLOAT,dimid,varid)) 
   call ncerr(NF90_enddef(ncid))
   call ncerr(NF90_put_var(ncid,varid,var,start=start))
   call ncerr(NF90_redef(ncid))
   end subroutine


   ! Put a 1-D slice into  real 2-d data array into variable vname
   subroutine ncinqputvarslice_r2d(ncid,vname,dimid,var,start)
   implicit none
   character(len=*), intent(in)  ::  vname
   integer, intent(in)           ::  ncid, dimid(2)
   real   , intent(in)           :: var(:)
   integer, intent(in)           :: start(2)
   integer :: ierr, varid
   ierr = nf90_inq_varid(ncid, vname, varid)
   if (ierr/=NF90_NOERR) call ncerr(NF90_def_var(ncid,vname,NF90_FLOAT,dimid,varid)) 
   call ncerr(NF90_enddef(ncid))
   call ncerr(NF90_put_var(ncid,varid,var,start=start))
   call ncerr(NF90_redef(ncid))
   end subroutine

   ! Put a 2-D slice into  real 3-d data array into variable vname
   subroutine ncinqputvarslice_r3d(ncid,vname,dimid,var,start)
   implicit none
   character(len=*), intent(in)  ::  vname
   integer, intent(in)           ::  ncid, dimid(3)
   real   , intent(in)           :: var(:,:)
   integer, intent(in)           :: start(3)
   integer :: ierr, varid
   ierr = nf90_inq_varid(ncid, vname, varid)
   if (ierr/=NF90_NOERR) call ncerr(NF90_def_var(ncid,vname,NF90_FLOAT,dimid,varid)) 
   call ncerr(NF90_enddef(ncid))
   call ncerr(NF90_put_var(ncid,varid,var,start=start))
   call ncerr(NF90_redef(ncid))
   end subroutine


   subroutine ncputatt_c(ncid,varname,attname,att)
   implicit none
   integer,          intent(in) :: ncid
   character(len=*), intent(in) :: varname, attname,att
   integer :: varid
   if (trim(varname)=='NF90_GLOBAL') then
      varid=NF90_GLOBAL
   else
      call ncerr(nf90_inq_varid(ncid, varname, varid))
   end if
   call ncerr(nf90_put_att(ncid, varid, attname, att))
   end subroutine


   subroutine ncputatt_r(ncid,varname,attname,att)
   implicit none
   integer,          intent(in) :: ncid
   character(len=*), intent(in) :: varname, attname
   real            , intent(in) :: att
   integer :: varid
   if (trim(varname)=='NF90_GLOBAL') then
      varid=NF90_GLOBAL
   else
      call ncerr(nf90_inq_varid(ncid, varname, varid))
   end if
   call ncerr(nf90_put_att(ncid, varid, attname, att))
   end subroutine


   subroutine ncerr(errcode)
   use netcdf
   implicit none
   integer, intent(in) :: errcode
   if (errcode/=NF90_NOERR) then
     write(6,'(a)') NF90_STRERROR(errcode)
      stop '(ncerr)'
   end if
   end subroutine


end module

