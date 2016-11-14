program p_postprocess_mersea
   use netcdf
   implicit none
#if defined (IARGC)
   integer*4, external :: iargc
#endif
   character(len=80) :: corfile
   character(len=nf90_max_name) :: varname
   integer :: ncid,varid,lparamdimid,cparamdimid,paramvarid
   integer :: lparamdimlen,cparamdimlen
   integer :: k,numnew,nvars,ncresult
   integer :: indx,indx1
   
   if (iargc()/=1) then
      print *,'Routine for converting a argo profile processed'
      print *,'with argocmp -ncdump to the correct format used'
      print *,'by the MERSEA project'
      print *,'Usage:'
      print *,'  postprocess_mersea modified-argo-file'
      stop
   end if
   call getarg(1,corfile)

   print *,'reading fields from '//trim(corfile)
   ncresult = nf90_open(trim(corfile),nf90_write,ncid) 
   if (ncresult /= nf90_noerr) then
      call ncerr( ncresult)
   end if

   ! Get variables in netcdf file
   call ncerr(nf90_inquire(ncid, nVariables=nvars))

   ! Get new variables
   numnew=0
   do k=1,nvars
      ! Get variable name
      call ncerr(nf90_inquire_variable(ncid, k, name=varname))

      if (mersea_check_new(varname)) then
         numnew=numnew+1
      end if
   end do

   if (numnew>0) then

      ! Put in define mode
      call ncerr( nf90_redef(ncid))

      ! Define new dimension "I_PARAMETERINT" and "C_PARAMETERINT"
      lparamdimlen=numnew
      cparamdimlen=12
      !print *,'def dim'
      call ncerr(nf90_def_dim(ncid, 'L_PARAMETERINT', lparamdimlen, lparamdimid))
      call ncerr(nf90_def_dim(ncid, 'C_PARAMETERINT', cparamdimlen, cparamdimid))

      ! Define new variable PARAMETERINT
      !print *,'def var'
      call ncerr(nf90_def_var(ncid, 'PARAMETERINT', nf90_char, (/cparamdimid,lparamdimid/), paramvarid))

      ! Put in datamode
      call ncerr(nf90_enddef(ncid))
      ! Get new variables
      numnew=1
      do k=1,nvars

         ! Get variable name
         call ncerr(nf90_inquire_variable(ncid, k, name=varname))

         if (mersea_check_new(varname)) then ! yes
            print *,'New variable :',trim(varname)
            call ncerr(nf90_put_var(ncid, paramvarid,varname(1:12),start=(/1,numnew/)))
            numnew=numnew+1
         end if
      end do
      numnew=numnew-1

      print '(a,i4,a)','File '//trim(corfile)//' modified; filled PARAMETERINT with ', &
              numnew,' fields'
   else
      print *,'No new variables found'
   end if

   ! Close file
   call ncerr(nf90_close(ncid))

   contains

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

   logical function mersea_check_new(varname)
   implicit none
   character(len=*), intent(in) :: varname

   ! List of new names
   integer, parameter :: nfields=18
   character(len=*), parameter, dimension(nfields) :: newfields =  &
      (/'TMLEV       ', &
        'SMLEV       ', &
        'WOA2001_PSAL', &
        'WOA2001_TEMP', &
        'HDCST_PSAL  ', &
        'HDCST_TEMP  ', &
        'ANA_PSAL    ', &
        'ANA_TEMP    ', &
        'FRCST6D_TEMP', &
        'FRCST6D_PSAL', &
        'PERS6D_TEMP ', &
        'PERS6D_PSAL ', &
        'FRCST3D_TEMP', &
        'FRCST3D_PSAL', &
        'PERS3D_TEMP ', &
        'PERS3D_PSAL ', &
        'FSTGS_TEMP  ', &
        'FSTGS_PSAL  '/)

   integer :: k

   mersea_check_new=.false.
   do k=1,nfields
      if (trim(varname)==trim(newfields(k))) mersea_check_new=.true.
   end do
   end function

end program p_postprocess_mersea
