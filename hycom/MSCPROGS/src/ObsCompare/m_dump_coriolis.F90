module m_dump_coriolis
private :: errcond, getdimms,getvarstat
contains

   !Reads Coriolis data from a NetCDF file

   subroutine dump_coriolis(filename,varprefix,varname,var,ndepths,nprof,undef)
      use netcdf
      implicit none

      ! Dummy variables
      integer, intent(in) :: nprof, ndepths
      character(len=*),             intent(in)  :: filename,varname,varprefix
      real, dimension(ndepths,nprof) :: var
      real, intent(in) :: undef
      integer :: natts_orig,var_id_orig

      ! Local variables
      character(len=1), allocatable, dimension(:,:) :: &
         qcvar, qcpres,qcdeph
      character(len=1), allocatable, dimension(:) :: &
         qcpos
      integer, dimension(nf90_max_dims) ::                 &
         lonvdims,latvdims,vardims,qcvardims, qcpresvdims, &
         presvdims,juldvdims,qcposvdims,dephvdims,qcdephvdims, &
         pltfvdims
      integer ::    var_id,&
         lonvndim, latvndim, varndim, qcvarndim,       &
         qcpresvndim, presvndim, juldvndim,qcposvndim, &
         dephvndim, qcdephvndim,pltfvndim
      real, allocatable :: deph(:,:)

      integer :: profdimid,profdimlen,zlevdimid,zlevdimlen,s8dimid,s8dimlen
      integer :: &
         lonvarid,latvarid, presvarid, qcvarid, pltfvarid, &
         qcpresvarid, varid,juldvarid,qcposvarid,dephvarid,qcdephvarid
      integer :: sstvarid,sstvndim
      integer :: ncid
      integer :: i,stat,j
      logical :: ex
      real :: updiff,lwdiff,diff
      real :: add_offset,scale_factor, fillvalue
      logical :: lpres, ldeph, lfillvalue
      character(len=20) :: attname

      ! Check if variable type is correct
      !if (varname/='PSAL' .and. varname/='TEMP') then
      !   print *,'Variable name '//varname//'is unknown'
      !   print *,'valid names are PSAL or TEMP'
      !   stop '(read_coriolis)'
      !end if



      ! Check for existence of netcdf file
      inquire(file=trim(filename),exist=ex)
      if (.not.ex) then
         print *,'Can not find file '//trim(filename)
      end if

      print *,'dumping '//varprefix//varname//'  fields to '//trim(filename)
      stat = nf90_open(trim(filename),nf90_write,ncid) 
      if (stat /= nf90_noerr) then
         call ncerr( stat)
      end if

      if (.not. probe_var(ncid,varname)) then
         print *,'netcdf file does not have this fiels - Skipped !'
         return
         call ncerr(NF90_CLOSE(ncid))
      end if


      ! Dimension IDs and lengths (assuming the names are right)
      call getdimms(ncid,'N_PROF',profdimid,profdimlen)
      call getdimms(ncid,'N_LEVELS',zlevdimid,zlevdimlen)

      ! Naive security cehck
      if (nprof/=profdimlen .or. zlevdimlen/=ndepths) then
         print *,'Dimension mismatch in dump_coriolis'
         stop
      end if

      ! Define variable     
      call ncerr(nf90_redef(ncid))
      call ncerr(NF90_DEF_VAR(ncid,varprefix//varname,NF90_FLOAT,(/zlevdimid,profdimid/),var_id))

      ! Get # Attributes of original variable name
      call ncerr(nf90_inq_varid       (ncid, varname, var_id_orig))
      call ncerr(nf90_inquire_variable(ncid, var_id_orig, nAtts=natts_orig))

      ! Copy original attributes to new variable
      lfillvalue=.false.
      do i=1,natts_orig
         call ncerr(nf90_inq_attname(ncid, var_id_orig, i, attname))
         call ncerr(nf90_copy_att   (ncid, var_id_orig, attname, ncid, var_id))

         ! Retrieve Fillvalue
         if (trim(attname)=='_FillValue') then
            lfillvalue=.true.
            call ncerr(nf90_get_att(ncid, var_id_orig, trim(attname), fillvalue))
            where (var==undef) var=fillvalue
         end if
      end do
      !call ncerr(NF90_PUT_ATT(ncid,var_id,'_FillValue'   ,undef))
      !call ncerr(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
      call ncerr(nf90_enddef(ncid))
      call ncerr(NF90_PUT_VAR(ncid,var_id,var))
      call ncerr(NF90_CLOSE(ncid))

      ! Revert to undef values  for ascii dump
      if (lfillvalue)  where (var==fillvalue) var=undef



   end subroutine dump_coriolis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The following are auxillary routines 



   ! Give netcdf error message
   subroutine ncerr(error)
      use netcdf
      implicit none
      integer, intent(in) :: error
      if (error/=nf90_noerr) then
         print *,'read_coriolis: '//nf90_strerror(error)
         stop
      end if
   end subroutine ncerr




   subroutine getdimms(ncid,dname,dmid,dmlen)
      use netcdf
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: dname
      integer,          intent(out) :: dmid,dmlen
      call ncerr( nf90_inq_dimid(ncid,dname,dmid) )
      call ncerr( nf90_inquire_dimension(ncid,dmid,len=dmlen) )
   end subroutine


   logical function probe_var(ncid,varname)
      use netcdf
      implicit none
      character(len=*), intent(in) :: varname
      integer         , intent(in) :: ncid
      integer :: test, varid
      test=nf90_inq_varid(ncid,varname,varid) 
      probe_var = test==nf90_noerr
   end function probe_var

end module m_dump_coriolis
