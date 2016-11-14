      program ensstat_field
      use mod_xc
      use mod_za
      use mod_grid
      use mod_hycomfile_io
      use m_parse_blkdat
      use netcdf
      implicit none 
#if defined (IARGC)
      integer*4, external :: iargc
#endif
      character(len=2) ck1,ck2
      character(len=8) cfld1,cfld2
      integer :: k1, k2, kdm, xpoint, ypoint, iarg
      real, allocatable, dimension(:,:) :: fld2, fld1
      real*8, allocatable, dimension(:,:) :: var1, var2, &
         ave1, ave2, cov, corr 
      real*4, dimension(:,:), allocatable :: iofld4
      real*8 :: rscale,rscale1
      integer :: numrec
      real :: rdummy

      character(len=80) :: ncfile, fname, ftype, tmparg
      integer :: idmid, jdmid, var_id,ncid,ierr, pointid
      type(hycomfile) :: hfile


      ! Initialize IO
      call xcspmd()
      call zaiost()
      call get_grid()

      ! Get vertical dim from blkdat
      call parse_blkdat('kdm   ','integer',rdummy,kdm)

      
      allocate(fld1(idm,jdm))
      allocate(fld2(idm,jdm))
      allocate(ave1(idm,jdm))
      allocate(ave2(idm,jdm))
      allocate(var1(idm,jdm))
      allocate(var2(idm,jdm))
      allocate(cov (idm,jdm))
      allocate(corr(idm,jdm))
      allocate(iofld4(idm,jdm))

      ! Input arguments
      if (iargc()>4) then
         call getarg(1,cfld1)
         call getarg(2,tmparg) ; read(tmparg,*) k1
         call getarg(3,cfld2)
         call getarg(4,tmparg) ; read(tmparg,*) k2
         if ( k1<0 .or.  k1 > kdm .or. k2<0 .or. k2>kdm) then
            print *,'vertical coordinate out of range'
            print *,'k1=',k1
            print *,'k2=',k2
            print *,'kdm=',kdm
            call exit(1)
         end if
      else 
         print *
         print *
         print *
         print '(a)','**************** ensstat_field ********************'
         print '(a)','*This routine will read several hycom files, and  *'
         print '(a)','*calculate average, variance and correlation      *'
         print '(a)','*between fields specified    as input arguments   *'
         print '(a)','*To specify fields, specify its name (ex temp)    *'
         print '(a)','*and its vertical layer index (ex 1 for surface)  *'
         print '(a)','*Two such field/vertical index indications must be*'
         print '(a)','*given.                                           *'
         print '(a)','*                                                 *'
         print '(a)','*The routine will calculate the covariance/corr   *'
         print '(a)','*between field1 and field2, where the correlation *'
         print '(a)','*is calculated at the same horizontal position    *'
         print '(a)','**************** ensstat_field ********************'
         print '(a)','Usage:'
         print '(a)','  ensstat_field fld1 k1 fld2 k2 files'
         print '(a)','Example:'
         print '(a)','  ensstat_field saln 1 temp 1 FORDAILY_2008_210*.a '
         print '(a)','Will calculate correlation between salinity and '
         print '(a)','temperature, both in vertical layer number 1 (surf)'
         stop '(ensstat_field)'
      end if


      numrec=0
      ave1=0
      ave2=0
      var1=0
      var2=0
      cov=0
      corr=0
      do  iarg=5,iargc()
         numrec=numrec+1
         call getarg(iarg,fname)
         print '(a,i4,x,a)','File ',numrec,trim(fname)

         ! Init & read file
         ftype=getfiletype(trim(fname))
         call initHF(hfile,fname,trim(ftype))
         call HFReadField(hfile,fld1,idm,jdm,cfld1,k1,1)
         call HFReadField(hfile,fld2,idm,jdm,cfld2,k2,1)


         ! Mean - both fields
         ave1=ave1+fld1
         var1=var1+fld1**2
         ave2=ave2+fld2
         var2=var2+fld2**2
         cov = cov + fld1*fld2
      enddo


      ! Final fields
      if (numrec<=1) then
         print *,'Insufficient records from HYCOM restart files'
         print *,'read_state_records returned ',numrec
         stop '(p_hcorr)'
      else
         rscale=1.0/(numrec)
         rscale1=1.0/(numrec-1.0)
      end if
      print *
      var1   = max(0.,rscale1 * (var1-rscale*ave1*ave1))
      var2   = max(0.,rscale1 * (var2-rscale*ave2*ave2))
      cov    = rscale1 * (cov -rscale*ave1*ave2)
      ave1   = rscale  * ave1
      ave2   = rscale  * ave2
      where (sqrt(var1*var2) <= 1d-6 * abs(cov))
         corr   = 0.
      elsewhere
         corr   = cov/(sqrt(var1*var2))
      endwhere

      print '(a,2g10.2)','Max/min average '//trim(cfld1)//':', &
                          minval(ave1),maxval(ave1)
      print '(a,2g10.2)','Max/min std '//trim(cfld1)//':', &
                          maxval(sqrt(var1)),minval(sqrt(var1))
      print *
      print '(a,2g10.2)','Max/Min average '//trim(cfld2)//':', &
                          minval(ave2),maxval(ave2)
      print '(a,2g10.2)','Max/min std '//trim(cfld2)//':', &
                          maxval(sqrt(var2)),minval(sqrt(var2))
      print *
      print '(a,2g10.2)','Min/Max covariance ('//trim(cfld1)//','//trim(cfld2)//'):',&
                         minval(cov),maxval(cov)
      print '(a,2g10.2)','Min/Max correlation ('//trim(cfld1)//','//trim(cfld2)//'):', &
                         minval(corr),maxval(corr)
      print *


      ! Netcdf file creation
      ncfile='ensstat_field.nc'
      write(ck1,'(i2.2)') k1
      write(ck2,'(i2.2)') k2
      print *,'Dumping to netcdf file ',trim(ncfile)
      if (NF90_CREATE(trim(ncfile),NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(obsstats)'
      end if
      ierr=NF90_put_att(ncid,NF90_GLOBAL,'number_of_records',numrec)
      ierr=NF90_DEF_DIM(ncid,'idm',idm,idmid)
      ierr=NF90_DEF_DIM(ncid,'jdm',jdm,jdmid)

      iofld4=plon
      ierr=NF90_DEF_VAR(ncid,'lon',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4)
 
      ! Redef/enddef is inefficient with classic netcdf, but it means we dont have to define
      ! one id for each var. I'm lazy, I know, but the penalty is not great.
      iofld4=plat
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'lat',NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=plon
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'lon',NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=depths;  where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'depth',NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=ave1 ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'ave_a'//trim(cfld1)//ck1,NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=ave2 ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'ave_b'//trim(cfld2)//ck2,NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,ave2)

      iofld4=var1 ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'var_a'//trim(cfld1)//ck1,NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=var2 ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'var_b'//trim(cfld2)//ck2,NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=cov ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'cov_a'//trim(cfld1)//ck1//'_b'//trim(cfld2)//ck2, &
                        NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=corr ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'corr_a'//trim(cfld1)//ck1//'_b'//trim(cfld2)//ck2,&
                        NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 
      ierr=NF90_CLOSE(ncid)

      end program ensstat_field


