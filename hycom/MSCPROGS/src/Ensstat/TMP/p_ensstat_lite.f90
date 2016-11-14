      program ensstat_lite
      use mod_xc
      use mod_za
      use mod_read_rstab
      use mod_ice_io
      use m_parse_args
      use m_parse_blkdat
      use netcdf
      implicit none 




      integer*4, external :: iargc

      character(len=80) infile,outfile,filebase,cline, &
         rstfile,icefile,filebase2,icedriftfile,auxfile
      character(len=6) cvarin
      character(len=7) tag7
      character(len=3) cmem
      character(len=2) ck1,ck2
      logical periodic, ex, lok

      real rv,rh,alp,bet, sd_d
      integer iens,i,j, n1,n2, layer, level
      integer :: kdm
      integer :: ios,fab

      real :: hmina,hmaxa,hminb,hmaxb,rvar

      character(len=3) rungen
      character(len=50) :: tmparg
      integer k,l,iostat
      real, allocatable, dimension(:,:) :: fld2, fld1
      real*8, allocatable, dimension(:,:) :: var1, var2, &
         ave1, ave2, cov, corr 
      real, allocatable, dimension(:,:)  ::modlon, modlat, depths
      real*8, dimension(:,:), allocatable :: iofld
      real*4, dimension(:,:), allocatable :: iofld4
      integer, dimension(:,:), allocatable :: imask
      real*8 :: rscale,rscale1
      integer :: numrec
      logical :: exa, exb
      real, parameter :: undef=-1e14
      real :: amax,amin, rdummy

      integer :: iret





      integer :: status
      character(len=80) :: ncfile
      integer :: idmid, jdmid, var_id,ncid,ierr, pointid
      integer :: ifile

! -------------- parse blkdat ------------------------------
      call parse_blkdat('idm   ','integer',rdummy,idm)
      call parse_blkdat('jdm   ','integer',rdummy,jdm)
      call parse_blkdat('kdm   ','integer',rdummy,kdm)
! -------------- parse blkdat ------------------------------

      ! Initialize IO
      call xcspmd()
      call zaiost()
      
      allocate(modlon(idm,jdm)) ! Nice to have
      allocate(modlat(idm,jdm)) ! Nice to have
      allocate(depths(idm,jdm)) ! Nice to have
      allocate(imask(idm,jdm)) ! Nice to have
      allocate(fld1(idm,jdm))
      allocate(fld2(idm,jdm))
      allocate(ave1(idm,jdm))
      allocate(ave2(idm,jdm))
      allocate(var1(idm,jdm))
      allocate(var2(idm,jdm))
      allocate(cov (idm,jdm))
      allocate(corr(idm,jdm))
      allocate(iofld(idm,jdm))
      allocate(iofld4(idm,jdm))


!KAL - Shouldnt be used anymore
!      open(10,file='newpos.uf',form='unformatted',status='old')
!      read(10) modlat,modlon 
!      close(10)

      call zaiopf('regional.depth.a','old',11)
      call zaiord(depths,imask,.false.,amin,amax,11)
      call zaiocl(11)
      where(depths>1e20) depths=0.


      call zaiopf('regional.grid.a','old',11)
      call zaiord(modlon,imask,.false.,amin,amax,11)
      call zaiord(modlat,imask,.false.,amin,amax,11)
      call zaiocl(11)

      call parse_args(idm,jdm,kdm)



      ios=0
      numrec=0
      ave1=0
      ave2=0
      var1=0
      var2=0
      cov=0
      corr=0
      do  ifile=1,nfiles
         numrec=numrec+1

         ! Get base file name
         rstfile=filelist(ifile)(:)
         fab=max(index(rstfile,'.a'),index(rstfile,'.b'))
         filebase=rstfile(1:fab-1)
         icefile=trim(filebase)//'ICE.uf'
         icedriftfile=auxfile
         !print *,'File base is      :',trim(filebase)
         !print *,'ICE file is       :',trim(icefile)
         !print *,'Ice drift file is :',trim(icedriftfile)

         ! Check for .b and .a file
         inquire(exist=exa,file=trim(filebase)//'.a')
         inquire(exist=exb,file=trim(filebase)//'.b')


         if (exa .and. exb) then
            print *,'Member ',numrec
            print *,trim(filebase)//'.[ab]'

            ! Simple sub located at bottom of this file - shares scope with main
            ! program
            call readonefield(cfld1,k1,fld1) 
            call readonefield(cfld2,k2,fld2) 

            ! Mean - both fields
            ave1=ave1+fld1
            var1=var1+fld1**2

            if (trim(mode)=="--field") then
               ave2=ave2+fld2
               var2=var2+fld2**2
               cov = cov + fld1*fld2
            elseif (trim(mode)=="--point") then
               ave2=ave2+fld2(xpoint,ypoint)
               var2=var2+fld2(xpoint,ypoint)**2
               cov = cov + fld1*fld2(xpoint,ypoint)
            end if
            print *
         else 
            print *,'(could not find one of '//trim(filebase2)//'.[ab])'
            stop
         end if
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

      print '(a,2g10.2)','Max/min average '//trim(cfld1)//':', minval(ave1),maxval(ave1)
      print '(a,2g10.2)','Max/min std '//trim(cfld1)//':',maxval(sqrt(var1)),minval(sqrt(var1))
      print *
      print '(a,2g10.2)','Max/Min average '//trim(cfld2)//':', minval(ave2),maxval(ave2)
      print '(a,2g10.2)','Max/min std '//trim(cfld2)//':',maxval(sqrt(var2)),minval(sqrt(var2))
      print *
      print '(a,2g10.2)','Min/Max covariance ('//trim(cfld1)//','//trim(cfld2)//'):',&
         minval(cov),maxval(cov)
      print '(a,2g10.2)','Min/Max correlation ('//trim(cfld1)//','//trim(cfld2)//'):', &
         minval(corr),maxval(corr)
      print *


      ! Netcdf file creation
      ncfile='ensstat_'//trim(mode(3:))//'mode.nc'
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
      if (trim(mode)=="--point") then
         ierr=NF90_DEF_DIM(ncid,'point',1,pointid)
      end if

      iofld4=modlon
      ierr=NF90_DEF_VAR(ncid,'lon',NF90_Float,(/idmid,jdmid/),var_id)
      ierr=NF90_ENDDEF(ncid)
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4)
 
      ! Redef/enddef is inefficient with classic netcdf, but it means we dont have to define
      ! one id for each var. I'm lazy, I know, but the penalty is not great.
      iofld4=modlat
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'lat',NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=modlon
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

      if (trim(mode)=="--point") then
         iofld4=modlon
         ierr=NF90_REDEF(ncid) ;  
         ierr=NF90_DEF_VAR(ncid,'lon_point',NF90_Float,(/pointid/),var_id) ; 
         ierr=NF90_ENDDEF(ncid) ; 
         ierr=NF90_PUT_VAR(ncid,var_id,iofld4(xpoint,ypoint)) ; 

         iofld4=modlat
         ierr=NF90_REDEF(ncid) ;  
         ierr=NF90_DEF_VAR(ncid,'lat_point',NF90_Float,(/pointid/),var_id) ; 
         ierr=NF90_ENDDEF(ncid) ; 
         ierr=NF90_PUT_VAR(ncid,var_id,iofld4(xpoint,ypoint)) ; 

         ierr=NF90_REDEF(ncid) ;  
         ierr=NF90_DEF_VAR(ncid,'i_point',NF90_INT,(/pointid/),var_id) ; 
         ierr=NF90_ENDDEF(ncid) ; 
         ierr=NF90_PUT_VAR(ncid,var_id,xpoint) ; 

         ierr=NF90_REDEF(ncid) ;  
         ierr=NF90_DEF_VAR(ncid,'j_point',NF90_INT,(/pointid/),var_id) ; 
         ierr=NF90_ENDDEF(ncid) ; 
         ierr=NF90_PUT_VAR(ncid,var_id,ypoint) ; 
      end if

      iofld4=ave1 ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'ave_'//trim(cfld1)//ck1,NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      if (trim(mode)=="--field") then
         iofld4=ave2 ; where(depths<.1) iofld4=undef
         ierr=NF90_REDEF(ncid) ;  
         ierr=NF90_DEF_VAR(ncid,'ave_'//trim(cfld2)//ck2,NF90_Float,(/idmid,jdmid/),var_id) ; 
         ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
         ierr=NF90_ENDDEF(ncid) ; 
         ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 
      else if (trim(mode)=="--field") then
         ierr=NF90_REDEF(ncid) ;  
         ierr=NF90_DEF_VAR(ncid,'ave_'//trim(cfld2)//ck2,NF90_Float,(/pointid/),var_id) ; 
         ierr=NF90_ENDDEF(ncid) ; 
         ierr=NF90_PUT_VAR(ncid,var_id,ave2(1,1))
      end if

      iofld4=var1 ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'var_'//trim(cfld1)//ck1,NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      if (trim(mode)=="--field") then
         iofld4=var2 ; where(depths<.1) iofld4=undef
         ierr=NF90_REDEF(ncid) ;  
         ierr=NF90_DEF_VAR(ncid,'var_'//trim(cfld2)//ck2,NF90_Float,(/pointid/),var_id) ; 
         ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
         ierr=NF90_ENDDEF(ncid) ; 
         ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 
      else if (trim(mode)=="--field") then
         ierr=NF90_REDEF(ncid) ;  
         ierr=NF90_DEF_VAR(ncid,'ave_'//trim(cfld2)//ck2,NF90_Float,(/pointid/),var_id) ; 
         ierr=NF90_ENDDEF(ncid) ; 
         ierr=NF90_PUT_VAR(ncid,var_id,var2(1,1) )
      end if

      iofld4=cov ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'cov_'//trim(cfld1)//ck1//'_'//trim(cfld2)//ck2,NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      iofld4=corr ; where(depths<.1) iofld4=undef
      ierr=NF90_REDEF(ncid) ;  
      ierr=NF90_DEF_VAR(ncid,'corr_'//trim(cfld1)//ck1//'_'//trim(cfld2)//ck2,NF90_Float,(/idmid,jdmid/),var_id) ; 
      ierr=NF90_put_att(ncid,var_id,'_FillValue',real(undef,4))
      ierr=NF90_ENDDEF(ncid) ; 
      ierr=NF90_PUT_VAR(ncid,var_id,iofld4) ; 

      ierr=NF90_CLOSE(ncid)



contains

      subroutine readonefield(cfld,k,fld) 
      implicit none
      integer, intent(in) :: k
      character(len=*), intent(in)  ::  cfld
      real, dimension(idm,jdm), intent(out) :: fld

            ! Read field 2
            if (trim(cfld2)=='drfticex' .or. trim(cfld2)=='drfticey') then
!               if (k/=0) then
!                  print *,'Ice variables must have level=0'
!                  stop '(p_ensstat_lite)'
!               end if
!               call read_icedrift(trim(icedriftfile),fld,cfld,idm,jdm,numrec,undef)
stop '(ice read too dagerous at the moment)'
!            else if (trim(cfld2)=='hicem' .or. trim(cfld2) == 'ficem' .or. &
!                     trim(cfld2)=='hsnwm' ) then
!               if (k2/=0) then
!                  print *,'Ice variables must have level=0'
!                  stop '(p_ensstat_lite)'
!               end if
!               call read_restart_ice(trim(icefile),fld,cfld,idm,jdm,numrec)
!               write(lp,*) 'Got field '//trim(cfld2)
!stop '(ice read too dagerous at the moment)'
            else
               call read_rstfield2d(filebase,cfld,fld,idm,jdm,k,undef)
            endif 
            where(depths<.1) fld=undef
            !print *,minval(fld,fld/=undef),maxval(fld,fld/=undef)
      end subroutine 

      end program ensstat_lite


