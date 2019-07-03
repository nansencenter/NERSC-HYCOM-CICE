program tecconv
! A simpler version of the m2nc/m2t scripts and associated programs.
! Variables to plot are read from extract - files, and then put into
! 4D/3D netcdf variables with horizontal, vertical and time dimensions
! Everything is presented on the model grid.
! TODO: adhere to rotation/sphere etc
!
! Knut Liseter
   use mod_xc
   use mod_za
   use mod_grid
   use mod_year_info
   use netcdf
   use m_handle_err
   use m_fields_to_plot
   use mod_hycomfile_io
   implicit none
   real*4, parameter :: undefr4=real(undef,kind=4)
   type(fields)       fld(1000)
   integer i,j,k
   integer nfld
   logical sphere,rotatell,normal,ex
   integer :: rdimid,jdimid,idimid,ncid,nrec
#if defined(IARGC)
   integer*4, external :: iargc
#endif
   character(len=80) :: filename, ftype
   integer :: kdm
   real    :: rtime
   real, allocatable :: twod1(:,:), twod2(:,:), dp(:,:,:), lind(:), &
      urot(:,:), vrot(:,:)
   integer :: dimms2D(2),dimms3D(3),dimms4D(4)
   integer :: idim_id, jdim_id, kdim_id, rdim_id,maxk, varid, varid2, &
      varid3, varid4, varidtime
   integer :: ifld, ifile,irec
   logical :: first, vecflag
   character(len=8):: uname, vname, pre

   ! KAL - not all these are implemented yet
   logical :: ltecplot     ,   &   ! Input flag -- produce tecplot output
              forceweekly  ,   &   ! Input flag -- force weekly  type files
              forcedaily   ,   &   ! Input flag -- force daily   type files
              forcerestart ,   &   ! Input flag -- force restart type files
              forcepak             ! Input flag -- force pak     type files
   type(hycomfile) :: hfile

   
   ltecplot=.false.
   forceweekly=.false.
   forcedaily=.false.
   forcerestart=.false.
   forcepak=.false.

   first=.true.
   irec=1

   if (iargc()<1) then
      print *,'h2nc will extract data from hycom .[ab] files'
      print *,'and put these into netcdf file tmp1.nc.'
      print *,'The fields to extract are specified in extract'
      print *,'files which corresponds to the file type you '
      print *,'want to extract from (for example extract.daily, '
      print *,'extract.restart ...). Sample extract files can '
      print *,'be found in Input directory under MSCProgs directory'
      print *,' '
      print *,'Several files can be specified which results in '
      print *,'several time records in the netcdf file'
      print *,''
      print *,'Vectors can be rotated, depending on the flags in'
      print *,'extract-files (TODO!)'
      print *,''
      print *,'usage:  h2nc  file(s)'
      print *
      print *,'The main difference between this routine and m2nc'
      print *,'is that relevant fields are dumped into 3D netcdf'
      print *,'variables. m2nc will dump everything into 2D variables.'
      print *
      print *,'file name must be supplied '
      stop '(tecconv)'
   end if
   call xcspmd()
   call zaiost()

   ! Get model grid
   call get_grid


   do ifile=1,iargc()

      call getarg(ifile,filename) 

      ! Check file type
      ftype=getfiletype(filename)


      ! Read file headers to get dimensions etc
      call initHF(hfile,trim(filename),trim(ftype))
      rtime=hfile%fyear
      kdm=vDim(hfile) 


      ! determine number of and which fields to be plotted (stored in fld), 
      ! the number of fields (nfld) and options for grid projections.
      if (first) then

         call fields_to_plot(sphere,rotatell,normal,fld,nfld,hfile,kdm)
         allocate(lind(kdm))
         ! Get max k dim
         maxk=0
         do k=1,nfld
            maxk=max(maxk,max(fld(k)%layb,fld(k)%laya))
         end do

         ! Layer index
         do k=1,maxk
            lind(k)=k
         end do
      end if
      write(6,'(a)') 'Processing '//trim(filename)


      ! Allocate fields
      allocate(dp     (idm,jdm,kdm))
      allocate(twod1  (idm,jdm))
      allocate(twod2  (idm,jdm))
      allocate(urot   (idm,jdm))
      allocate(vrot   (idm,jdm))



      ! First pass we have to create the file - this sets up bathy, 
      ! lon, lat etc etc ...
      if ( first) then

         ! Start dumping to netcdf file
         if (NF90_create('tmp1.nc',NF90_NETCDF4,ncid) /= NF90_NOERR) then
            print *,'An error occured when opening the netcdf file'
            stop '(restart2netcdf)'
         end if

         ! Define dimensions
         call handle_err(NF90_DEF_DIM(ncid,'idim',idm            ,idim_id))
         call handle_err(NF90_DEF_DIM(ncid,'jdim',jdm            ,jdim_id))
         if (maxk>0) then
            call handle_err(NF90_DEF_DIM(ncid,'kdim',maxk          ,kdim_id))
         end if
         call handle_err(NF90_DEF_DIM(ncid,'rdim',nf90_unlimited,rdim_id))
         dimms2D=(/idim_id,jdim_id/)
         dimms3D=(/idim_id,jdim_id,rdim_id/)
         dimms4D=(/idim_id,jdim_id,kdim_id,rdim_id/)

         ! Time info ("floating-point year")
         call handle_err(NF90_DEF_VAR(ncid,'fyear',NF90_Float,rdim_id,varidtime))

         ! Define some basic vars
         call handle_err(NF90_DEF_VAR(ncid,'depth',NF90_Float,dimms2D,varid))
         call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',real(undefr4,kind=4)))
         call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
         call handle_err(NF90_ENDDEF(ncid))
         twod1=depths
         where(depths<.1) twod1=undefr4
         call handle_err(NF90_PUT_VAR(ncid,varid,twod1))
         call handle_err(NF90_REDEF(ncid))
          
         call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_Float,dimms2D,varid))
         call handle_err(NF90_ENDDEF(ncid))
         call handle_err(NF90_PUT_VAR(ncid,varid,plon(1:idm,1:jdm)))
         call handle_err(NF90_REDEF(ncid))
          
         call handle_err(NF90_DEF_VAR(ncid,'latitude',NF90_Float,dimms2D,varid))
         call handle_err(NF90_ENDDEF(ncid))
         call handle_err(NF90_PUT_VAR(ncid,varid,plat(1:idm,1:jdm)))
         call handle_err(NF90_REDEF(ncid))

         if (maxk>0) then
            call handle_err(NF90_DEF_VAR(ncid,'layer_index',NF90_INT,kdim_id,varid))
            call handle_err(NF90_ENDDEF(ncid))
            call handle_err(NF90_PUT_VAR(ncid,varid,lind(1:maxk)))
         else
            call handle_err(NF90_ENDDEF(ncid))
         end if

      end if

      ! Put time
      call handle_err(NF90_PUT_VAR(ncid,varidtime,rtime,start=(/irec/)))


      ! Read layer thickness
      do k=1,maxk
         !call HFReadDPField(hfile,dp(:,:,k),idm,jdm,k,1)
         call HFReadDPField_m(hfile,dp(:,:,k),idm,jdm,k,1) ! Returns dp in meters
         where(depths   <.1) dp(:,:,k)=0.
      end do


      ! Define variables, and put them in the netcdf file
      ! NB: This loop does not test on vectors (see below)
      do ifld=1,nfld
      if (fld(ifld)%option .and. .not.fld(ifld)%vecflag) then

         ! define/inquire the variables -- See subroutine at the end of this file
         call deforinqvar(hfile,ncid,varid,fld(ifld)%fextract,fld(ifld)%fextract,dimms3D,dimms4D,first)


         ! Cycle layers to be extracted
         do k=fld(ifld)%laya,fld(ifld)%layb
            call HFReadField(hfile,twod1,idm,jdm,fld(ifld)%fextract,k,1)

            ! Mask data
            where(depths   <.1) twod1=undefr4
            if (is3Dvar(hfile,fld(ifld)%fextract,1)) then
               if (k>0) then
                  where(dp(:,:,k)<0.1) twod1=undefr4
               end if
            end if

            ! Dump data into netcdf file
            if (is3Dvar(hfile,fld(ifld)%fextract,1)) then
               call handle_err(NF90_PUT_VAR(ncid,varid,twod1,start=(/1,1,k,irec/)))
            else
               call handle_err(NF90_PUT_VAR(ncid,varid,twod1,start=(/1,1,irec/)))
            end if

         end do ! k
      end if
      end do ! ifld



      ! Define variables, and put them in the netcdf file
      ! NB: This loop tests on vectors only
      do ifld=1,nfld-1
      if (fld(ifld)%option .and.  fld(ifld)%vecflag) then
         uname=fld(ifld  )%fextract
         vname=fld(ifld+1)%fextract
         pre=fld(ifld+1)%vecpost

         ! define/inquire the variables -- See subroutine at the end of this file
         if (rotatell) then
            call deforinqvar(hfile,ncid,varid ,uname,'UROT_'//trim(pre),dimms3D,dimms4D,first)
            call deforinqvar(hfile,ncid,varid2,vname,'VROT_'//trim(pre),dimms3D,dimms4D,first)
         end if

         ! define/inquire the variables -- See subroutine at the end of this file
         if (normal) then
            call deforinqvar(hfile,ncid,varid3,uname,uname,dimms3D,dimms4D,first)
            call deforinqvar(hfile,ncid,varid4,vname,vname,dimms3D,dimms4D,first)
         end if

         ! Cycle layers to be extracted
         do k=fld(ifld)%laya,fld(ifld)%layb
            call HFReadField(hfile,twod1,idm,jdm,uname,k,1)
            call HFReadField(hfile,twod2,idm,jdm,vname,k,1)

            urot=twod1
            vrot=twod2

            if (rotatell) then
               call rotate(urot,vrot,plat,plon,idm,jdm,'m2l')
            end if

            ! Mask data
            where(depths   <.1) 
               twod1=undefr4
               twod2=undefr4
               urot =undefr4
               vrot =undefr4
            end where
            if (k>0) then
               where(dp(:,:,k)<0.1)
                  twod1=undefr4
                  twod2=undefr4
                  urot =undefr4
                  vrot =undefr4
               end where
            end if

            ! Dump unrotated data into netcdf file
            if (normal) then
               if (is3Dvar(hfile,fld(ifld)%fextract,1)) then
                  call handle_err(NF90_PUT_VAR(ncid,varid3,twod1,start=(/1,1,k,irec/)))
                  call handle_err(NF90_PUT_VAR(ncid,varid4,twod2,start=(/1,1,k,irec/)))
               else
                  call handle_err(NF90_PUT_VAR(ncid,varid3,twod1,start=(/1,1,irec/)))
                  call handle_err(NF90_PUT_VAR(ncid,varid4,twod2,start=(/1,1,irec/)))
               end if
            end if

            ! Dump rotated data into netcdf file
            if (rotatell) then
               if (is3Dvar(hfile,fld(ifld)%fextract,1)) then
                  call handle_err(NF90_PUT_VAR(ncid,varid ,urot,start=(/1,1,k,irec/)))
                  call handle_err(NF90_PUT_VAR(ncid,varid2,vrot,start=(/1,1,k,irec/)))
               else
                  call handle_err(NF90_PUT_VAR(ncid,varid ,urot,start=(/1,1,irec/)))
                  call handle_err(NF90_PUT_VAR(ncid,varid2,vrot,start=(/1,1,irec/)))
               end if
            end if

         end do ! k
      end if ! vector flag vecflag
      end do ! ifld

      
      deallocate(dp   )
      deallocate(twod1)
      deallocate(twod2)
      deallocate(urot)
      deallocate(vrot)
      first=.false.
      irec=irec+1

   end do ! ifile
   call handle_err(NF90_close(ncid))
   print *,'Data dumped in tmp1.nc'



   contains 
      
      

      ! Subroutine defines variable on first pass, retrieves variable id on 
      ! subsequent passes
      subroutine deforinqvar(hfile,ncid,varid,vnameh,vnamenc,dimms3D,dimms4D,first)
      implicit none
      type(hycomfile),  intent(in)    :: hfile
      integer        ,  intent(in)    :: ncid, dimms3D(3), dimms4D(4)
      integer        ,  intent(out)   :: varid
      character(len=*), intent(in)    :: vnameh  ! name in hycom file
      character(len=*), intent(in)    :: vnamenc ! name in netcdf file
      logical         , intent(in)    :: first


      ! First file we define the variables
      if (first) then
         call handle_err(NF90_REDEF(ncid))
         if (is3Dvar(hfile,vnameh,1)) then
            call handle_err(NF90_DEF_VAR(ncid,trim(vnamenc),NF90_Float,dimms4D,varid))
         else
            call handle_err(NF90_DEF_VAR(ncid,trim(vnamenc),NF90_Float,dimms3D,varid))
         end if
         call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',real(undefr4,kind=4)))
         call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
         call handle_err(NF90_ENDDEF(ncid))

      ! On subsequent files we retrieve netcdf variable id
      else
         call handle_err(nf90_inq_varid(ncid, trim(vnamenc),varid))
      end if

      end subroutine








end program tecconv
