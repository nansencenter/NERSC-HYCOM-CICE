module m_read_coriolis
private :: errcond, getdimms,getvarstat
contains

   !Reads Coriolis data from a NetCDF file

   subroutine read_coriolis(filename,var,pres,longitude,latitude,juld,floatid, &
                            errorflag,l_invalid,varname,undef)
      use netcdf
      implicit none

      ! Dummy variables
      character(len=*),                intent(in)  :: filename,varname
      real, pointer,   dimension(:)   :: longitude,latitude,juld
      real, pointer,   dimension(:,:) :: var,pres
      logical, pointer,   dimension(:,:) :: errorflag
      character(len=8),dimension(:) , pointer :: floatid
      logical         ,dimension(:) , pointer :: l_invalid
      real, intent(in) :: undef
      character(len=16) :: tmp
      character(len=80), dimension(:), allocatable :: depthname,varname2

      ! Local variables
      character(len=1), allocatable, dimension(:,:,:) ::  &
         sparam
      character(len=1), allocatable, dimension(:,:) :: &
         qcvar, qcpres,qcdeph,dstate,dmode
      character(len=1), allocatable, dimension(:) :: &
         qcpos,qcjuld
      logical         , allocatable, dimension(:) :: &
         l_deph, l_adj

      integer, dimension(nf90_max_dims) :: vdims ! Variable dimensions - not kept
      integer                           :: vndim ! number of variable dims - not kept
      integer ::   dstatvarid, dmodevarid,sparamvarid

      ! Dimension id - not kept
      integer :: dimid

      ! Dimension lengths
      integer :: profdimlen,zlevdimlen,s8dimlen,s4dimlen, nparamdimlen, s16dimlen 

      ! Variable ids - Kept !
      integer :: &
         lonvarid,latvarid, presvarid, qcvarid, pltfvarid, &
         qcpresvarid, varid,juldvarid,qcposvarid,qcjuldvarid,dephvarid,qcdephvarid
      integer :: ncid
      integer :: i,stat,j,k,k2
      logical :: ex
      real :: updiff,lwdiff,diff
      real :: add_offset,scale_factor, fillvalue
      logical :: lpres, ldeph,l_present
      integer :: ncresult
      real :: fillvalue_pres,fillvalue_var

      ! Check if variable type is correct
      if (varname/='PSAL' .and. varname/='TEMP') then
         print *,'Variable name '//varname//'is unknown'
         print *,'valid names are PSAL or TEMP'
         stop '(read_coriolis)'
      end if



      ! Check for existence of netcdf file
      inquire(file=trim(filename),exist=ex)
      if (.not.ex) then
         print *,'Can not find file '//trim(filename)
         stop '(read_coriolis)'
      end if

      print *,'reading Coriolis fields from '//trim(filename)
      stat = nf90_open(trim(filename),nf90_nowrite,ncid) 
      if (stat /= nf90_noerr) then
         call ncerr( stat)
      end if

      ! Check if sought variable is present
      l_present=probe_var(ncid,varname)

      !print *,'reading Coriolis fields from '//trim(filename)


      ! Dimension IDs and lengths (assuming the names are right)
      call getdimms(ncid,'N_PROF'  ,dimid,profdimlen  )
      call getdimms(ncid,'N_LEVELS',dimid,zlevdimlen  )
      call getdimms(ncid,'N_PARAM' ,dimid,nparamdimlen)
      call getdimms(ncid,'STRING16',dimid,s16dimlen   )
      call getdimms(ncid,'STRING8' ,dimid,s8dimlen    )
      call getdimms(ncid,'STRING4' ,dimid,s4dimlen    )
      print *,'# Profiles , # levels ',profdimlen, zlevdimlen
      !print *,profdimid, zlevdimid

      ! Better safe than sorry ...
      if (profdimlen<1 .or. zlevdimlen<1) then
         print *,'Dimensions are zero for this file ... '
         print *,'Exiting..'
         stop '(read_coriolis)'
      end if


      ! Get data state indicator, mode and station parameters
      call getvarstat(ncid,'DATA_STATE_INDICATOR',dstatvarid ,vndim  ,vdims)
      call getvarstat(ncid,'DATA_MODE'           ,dmodevarid ,vndim  ,vdims)
      call getvarstat(ncid,'STATION_PARAMETERS'  ,sparamvarid,vndim  ,vdims)
      allocate(sparam   (s16dimlen,nparamdimlen,profdimlen))
      allocate(dstate   (s4dimlen,profdimlen))
      allocate(dmode    (1       ,profdimlen))
      allocate(l_deph   (profdimlen))
      allocate(l_adj    (profdimlen))
      allocate(l_invalid(profdimlen))
      l_deph=.false.
      l_adj =.false.
      l_invalid =.false.
      do j=1,profdimlen
         
         ! Data state indicator
         call ncerr( nf90_get_var(ncid, dstatvarid, dstate(:,j), &
             start=(/1,j/),count=(/s4dimlen,1/)))
         !print *,dstate(:,j)

         ! Data mode (R)ealtime, (D)elayed or (A)djusted
         call ncerr( nf90_get_var(ncid, dmodevarid, dmode (:,j), &
             start=(/1,j/),count=(/1,1/)))
         !print *,dmode(1,j)

         ! Station parameters
         call ncerr( nf90_get_var(ncid, sparamvarid, sparam (:,:,j), &
             start=(/1,1,j/),count=(/s16dimlen,nparamdimlen,1/)))
         !print *,sparam(:,:,j)


         do k=1,nparamdimlen


            ! Check what variable to use for depth


            tmp=''
            do k2=1,16 
              tmp(k2:k2)=sparam(k2,k,j)
            end do
            !print *,sparam(:,k,j)
            if (trim(tmp)=='DEPH') then
               l_deph(j)=.true.
            elseif (trim(tmp)=='PRES') then
               l_deph(j)=.false.
            end if

            !! Check that variable we want is present
            !tmp=''
            !do k2=1,16 
            !  tmp(k2:k2)=sparam(k2,k,j)
            !end do
            !if (trim(tmp)==trim(varname)) then
            !   l_present(j)=.true.
            !end if
         end do


         ! Check if adjusted variables should be used
         tmp=''
         do k2=1,4 
           tmp(k2:k2)=dstate(k2,j)
         end do
         if (dmode(1,j)(1:1)=='A') then
            l_adj(j)=.true.
         elseif (dmode(1,j)(1:1)=='R') then
            if ( trim(tmp) == "2C" .or. trim(tmp)=="2B+") then
               l_adj(j)=.true.
            end if
         end if

         ! Check if profile is invalid
         if ( trim(tmp) == "0A") then
            l_invalid(j)=.true.
         end if

         !print *,tmp,l_deph(j),l_adj(j),l_invalid(j)
      end do


      ! platform, lon, lat and julian day
      call getvarstat(ncid,'PLATFORM_NUMBER',pltfvarid,vndim  ,vdims)
      call getvarstat(ncid,'LONGITUDE'      ,lonvarid ,vndim  ,vdims)
      call getvarstat(ncid,'LATITUDE'       ,latvarid ,vndim  ,vdims)
      call getvarstat(ncid,'JULD'           ,juldvarid,vndim  ,vdims)
      allocate(floatid  (profdimlen))
      allocate(longitude(profdimlen))
      allocate(latitude (profdimlen))
      allocate(juld     (profdimlen))
      call ncerr( nf90_get_var(ncid, juldvarid, juld ) )
      call ncerr( nf90_get_var(ncid, lonvarid, longitude) )
      call ncerr( nf90_get_var(ncid, latvarid, latitude ) )
      do j=1,profdimlen
         call ncerr( nf90_get_var(ncid, pltfvarid    ,floatid(j)(:),   &
                start=(/1,j/),count=(/s8dimlen,1/)) )
         !print *,floatid(j)
      end do

      ! QC flags for juld and position
      call getvarstat(ncid,'POSITION_QC',qcposvarid ,vndim  ,vdims)
      call getvarstat(ncid,'JULD_QC'    ,qcjuldvarid,vndim  ,vdims)
      allocate(qcpos    (profdimlen))
      allocate(qcjuld   (profdimlen))
      call ncerr( nf90_get_var(ncid, qcposvarid    , qcpos,   &
         start=(/1/),count=(/profdimlen/)) )
      call ncerr( nf90_get_var(ncid, qcjuldvarid    , qcjuld,   &
         start=(/1/),count=(/profdimlen/)) )

      !print *,'qcpos :',qcpos
      !print *,'qcjuld:',qcjuld

      ! Adjust invalid flag based on position and juld quality flags
      do j=1,profdimlen

         if(qcpos (j) /= '1' .and. qcpos(j) /='0') then
            l_invalid(j)=.true.
         end if

         if(qcjuld(j) /= '1' .and. qcjuld(j) /='0') then
            l_invalid(j)=.true.
         end if
      end do


      ! Set appropriate depth variable name
      allocate(depthname(profdimlen))
      allocate(varname2 (profdimlen))
      do j=1,profdimlen
         if (l_deph(j) ) then
            depthname(j)='DEPH'
         else
            depthname(j)='PRES'
         end if

         if (l_adj(j)) then
            depthname(j)=depthname(j)//'_ADJUSTED'
            varname2 (j)=trim(varname)//'_ADJUSTED'
         else
            varname2 (j)=trim(varname)
         end if

         !print *,j,trim(depthname(j))

         ! Security check - we assume that either DEPTH or PRES is present....
         ! If this becomes a problem we must gracefully handle it....
         !if (j>1) then
         !   if (depthname(j)(1:4)/=depthname(1)(1:4)) then
         !      print *,'Error - both PRES and DEPH are used in netcdf file!'
         !      stop
         !   end if
         !end if
         ! KAL - the routine should be able to handle this, depthname is now a
         ! profile vector
      end do

      

      ! Read profile depths and values
      allocate(var   (zlevdimlen,profdimlen))
      allocate(pres  (zlevdimlen,profdimlen))
      allocate(qcvar (zlevdimlen,profdimlen))
      allocate(qcpres(zlevdimlen,profdimlen))
      var=undef
      do j=1,profdimlen

         !print *,j,trim(varname2(j)),' ',trim(depthname(j))

         ! hmm.. overkill for var stats but...
         call getvarstat(ncid,trim(depthname(j))       ,presvarid  ,vndim,vdims )
         call getvarstat(ncid,trim(depthname(j))//'_QC',qcpresvarid    ,vndim,vdims )
         if (l_present) then
            call getvarstat(ncid,trim(varname2 (j))       ,varid      ,vndim,vdims )
            call getvarstat(ncid,trim(varname2 (j))//'_QC',qcvarid,vndim,vdims )
         end if


         
         call ncerr( nf90_get_var(ncid,presvarid    , pres(:,j), &
                                  start=(/1,j/),count=(/zlevdimlen,1/)))
         if (l_present) then
            call ncerr( nf90_get_var(ncid, varid    , var(:,j), &
                                     start=(/1,j/),count=(/zlevdimlen,1/)))

            ! get fillvalue attribute of pressure and variable name
            ncresult=nf90_get_att(ncid, varid    , '_FillValue', fillvalue_var)
            where (var (:,j)==fillvalue_var ) var (:,j)=undef
         end if

         ! set to undef (input) where variable and pressure equals fillvalues
         ncresult=nf90_get_att(ncid, presvarid, '_FillValue', fillvalue_pres)
         where (pres(:,j)==fillvalue_pres) pres(:,j)=undef
         

         !print *,var(:,j)
         call ncerr( nf90_get_var(ncid, qcpresvarid, qcpres(:,j), &
             start=(/1,j/),count=(/zlevdimlen,1/)))
         if (l_present) then
            call ncerr( nf90_get_var(ncid, qcvarid    , qcvar(:,j), &
                   start=(/1,j/),count=(/zlevdimlen,1/)))
         end if
         !print *,qcvar(:,j)

      end do
      !print *,'hei'

         

      ! Go through data and discard poor observations (errorflag=true)
      allocate(errorflag(zlevdimlen,profdimlen))
      errorflag=.false.
      do j=1,profdimlen
      do i=1,zlevdimlen
         if (l_present) then
            if(qcvar (i,j) /= '1' .and. qcvar(i,j)  /='0') errorflag(i,j)=.true.
         end if
         if(qcpres(i,j) /= '1' .and. qcpres(i,j) /='0') errorflag(i,j)=.true.
      end do
      end do


      do j=1,profdimlen
         if (l_invalid(j)) then
            errorflag(:,j)=.true.
            var(:,j)=undef
         end if
      end do



   end subroutine read_coriolis


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



   subroutine getvarstat(ncid,varname,varid,var_ndim,var_dimids)
      use netcdf
      implicit none
      integer,                          intent(in)  :: ncid
      character(len=*),                 intent(in)  :: varname
      integer,                          intent(out) :: varid,var_ndim
      integer, dimension(nf90_max_dims),intent(out) :: var_dimids
      call ncerr( nf90_inq_varid(ncid,varname,varid) )
      call ncerr( nf90_inquire_variable(ncid,varid,ndims=var_ndim,dimids=var_dimids) )
   end subroutine getvarstat


end module m_read_coriolis
