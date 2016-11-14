      program check_nest
      use mod_nest_io
      use m_parse_blkdat
      use netcdf
      implicit none
      character(len=80)           :: fndepths,nestfile,matfil, ncfile,
     &   tmparg
      character(len=20)           :: ctitle
      real, allocatable, dimension(:,:) :: tmpfld
      real, parameter  :: onem=9806.
      integer :: i1,i2,j1,j2,indbnd,imax,iday,iyear,i,j,k,l,k2, 
     &   kdm,ihour
      integer :: ncid, varid, idmid, jdmid, rdimid, ierr, kdmid
      real    :: rvar
      real, parameter :: undef=-1e14
      integer*4, external :: iargc
c
      if (iargc()/=3) then
         print *,'Program diagnoses nesting files for v 2.2'
         print *,'of NERSC-HYCOM. Produces a netcdf file containing '
         print *,'the nesting fields. Run this routine in the directory'
         print *,'containing the nesting files (in ./nest/)'
         print *
         print *,'Usage:  checknest <year> <julian day> <hour>'
         call exit(1)
      else
c ---    Get date to process
         call getarg(1,tmparg) ; read(tmparg,*) iyear
         call getarg(2,tmparg) ; read(tmparg,*) iday
         call getarg(3,tmparg) ; read(tmparg,*) ihour
      end if


c
c --- Construct file names
      write(ncfile,'(a,i4.4,a,i3.3,a,i2.2,a)')
     &  'nest_',iyear,'_',iday,'_',ihour,'.nc'
      write(nestfile,'(a,i4.4,a,i3.3,a,i2.2)')
     &  'nest/nest_',iyear,'_',iday,'_',ihour
      write(6,'(a)') 'Input file is '//trim(nestfile)
      if (NF90_CREATE(trim(ncfile),NF90_CLOBBER,ncid) /= NF90_NOERR) 
     &   then
         print *,'An error occured when opening the netcdf file'
         stop '(obsstats)'
      end if
c
c --- Read the depth matrix of the local grid 
      call parse_blkdat('kdm   ','integer',rvar,kdm)
      call readinnergrid('nest/')
      allocate(tmpfld(inner_idm,inner_jdm))
c
c --- Netcdf file def & put follows
      ierr=NF90_DEF_DIM(ncid,'idm',inner_idm,idmid)
      ierr=NF90_DEF_DIM(ncid,'jdm',inner_jdm,jdmid)
      ierr=NF90_DEF_DIM(ncid,'kdm',kdm,kdmid)
c --- Longitude
      ierr=NF90_DEF_VAR(ncid,'longitude',NF90_Float,(/idmid,jdmid/),
     &   varid)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,varid,inner_lon)
      ierr=NF90_REDEF(ncid)
c --- Latitude
      ierr=NF90_DEF_VAR(ncid,'latitude',NF90_Float,(/idmid,jdmid/),
     &   varid)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,varid,inner_lat)
      ierr=NF90_REDEF(ncid)
c --- Depth
      ierr=NF90_DEF_VAR(ncid,'depth',NF90_Float,(/idmid,jdmid/),varid)
      tmpfld=inner_depths
      where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,varid,tmpfld)
      ierr=NF90_REDEF(ncid)
c --- Temperature       
      print *,'temperature ...'
      ierr=NF90_DEF_VAR(ncid,'temp',NF90_Float,
     &   (/idmid,jdmid,kdmid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpfld=undef
      do k=1,kdm
         call readnest(trim(nestfile),'temp    ',k,tmpfld,
     &      inner_idm, inner_jdm,inner_inest)
         where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k/))
      end do
      ierr=NF90_REDEF(ncid)
c --- Salinity       
      print *,'salinity ...'
      ierr=NF90_DEF_VAR(ncid,'saln',NF90_Float,
     &   (/idmid,jdmid,kdmid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpfld=undef
      do k=1,kdm
         call readnest(trim(nestfile),'saln    ',k,tmpfld,
     &      inner_idm, inner_jdm,inner_inest)
         where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k/))
      end do
      ierr=NF90_REDEF(ncid)
c --- Interfaces       
      print *,'interface ...'
      ierr=NF90_DEF_VAR(ncid,'intf',NF90_Float,
     &   (/idmid,jdmid,kdmid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpfld=undef
      do k=1,kdm
         call readnest(trim(nestfile),'intf    ',k,tmpfld,
     &      inner_idm, inner_jdm,inner_inest)
         where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k/))
      end do
      ierr=NF90_REDEF(ncid)
c --- Baroclinic velocity       
      print *,'u ...'
      ierr=NF90_DEF_VAR(ncid,'u',NF90_Float,
     &   (/idmid,jdmid,kdmid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpfld=undef
      do k=1,kdm
         call readnest(trim(nestfile),'u       ',k,tmpfld,
     &      inner_idm, inner_jdm,inner_inest)
         where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k/))
      end do
      ierr=NF90_REDEF(ncid)
c --- Baroclinic velocity       
      print *,'v ...'
      ierr=NF90_DEF_VAR(ncid,'v',NF90_Float,
     &   (/idmid,jdmid,kdmid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpfld=undef
      do k=1,kdm
         call readnest(trim(nestfile),'v       ',k,tmpfld,
     &      inner_idm, inner_jdm,inner_inest)
         where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
         ierr=NF90_PUT_VAR(ncid,varid,tmpfld,start=(/1,1,k/))
      end do
      ierr=NF90_REDEF(ncid)
c --- Barotropic velocity       
      print *,'ubavg ...'
      ierr=NF90_DEF_VAR(ncid,'ub',NF90_Float, (/idmid,jdmid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpfld=undef
      call readnest(trim(nestfile),'ubavg   ',0,tmpfld,
     &      inner_idm, inner_jdm,inner_inest)
      where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
      ierr=NF90_PUT_VAR(ncid,varid,tmpfld)
      ierr=NF90_REDEF(ncid)
c --- Barotropic velocity       
      print *,'vbavg ...'
      ierr=NF90_DEF_VAR(ncid,'vb',NF90_Float, (/idmid,jdmid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpfld=undef
      call readnest(trim(nestfile),'vbavg   ',0,tmpfld,
     &      inner_idm, inner_jdm,inner_inest)
      where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
      ierr=NF90_PUT_VAR(ncid,varid,tmpfld)
      ierr=NF90_REDEF(ncid)
c --- Barotropic pressure component       
      print *,'pbavg ...'
      ierr=NF90_DEF_VAR(ncid,'pb',NF90_Float, (/idmid,jdmid/),varid)
      ierr=NF90_PUT_ATT(ncid,varid,'_FillValue',real(undef,kind=4))
      ierr=NF90_ENDDEF(ncid)
      tmpfld=undef
      call readnest(trim(nestfile),'pbavg   ',0,tmpfld,
     &      inner_idm, inner_jdm,inner_inest)
      where (inner_depths<.1 .or. inner_depths>1e20) tmpfld=undef
      ierr=NF90_PUT_VAR(ncid,varid,tmpfld)
      ierr=NF90_close(ncid)
C --- Done..
      print *,'Output in '//trim(ncfile)
      end program check_nest
