!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine creates NetCDF files from (various) hycom files. The netcdf
! files hold the data on a projection grid at specified depthlevels. Routine 
! works on several HYCOM files; restart, archive, nersc weekly average and nersc
! daily average.
!
! -------------------------------------------------------------------------
! Preparation for data processing:
!
! 1) Specify a projection grid in longitude latitude coordinates in the file named 
! "proj.in".  TODO: Examples ?
!
! 2) Specify depth levels to which you want to interpolate in the file
! "depthlevels.in". Here you specify the number of depth levels on the first
! line, then the actual depthlevels on the following lines. Depth levels must be
! increasing (down is positive direction).
!
! 3) Specify fields to extract from the HYCOM files in 'extract' - files. TODO:
! example?
!
! 4) Several other datafiles will be needed as well. These are files used when
! HYCOM runs, so if you run "hyc2proj" in the same catalogue as you run
! hycom, everything should be ok.
!
! 5) Run "hyc2proj" to generate netcdf files. Run without arguments for info on
! its usage.  For each file given on the command line, a netcdf file will be
! created.  ! To browse the data file you can use the program "ncview" installed on tre.
!
!
!NB: Sample infiles are given in the catalogue SAMPLE_INFILES. TODO: They are
!not
!
! KAL 20080515 - Added time info inside netcdf file as well. Also worked on path0
!                command line input.
! KAL 20090113 - Re-write of code. Should be easier to maintain/update.
! LB  20090209 - Added MLD and BSF variables
!---------------------------------------------------------------------------------


program p_hyc2proj
   use mod_parameters
   use mod_xc, only: idm,jdm, xcspmd
   use mod_grid, only : get_grid, depths, plon, plat, qlat, qlon
   use mod_year_info
   use mod_levitus
   use mod_toproj
   use mod_spline_calc
   use mod_hycomfile_io
   use netcdf
   use m_fields_to_plot
   use mod_rotate
   use m_strmf_eval
   use m_mixlayer_depths
   use m_bio_conversions
   use m_mersea_prepare
   use mod_za , only: zaiost, zaioempty
   use mod_netcdf_file
   use m_process_arguments
   use m_read_mean_ssh
   use m_construct_filename
   implicit none

   real, parameter :: thref=1e-3
   include 'stmt_funcs.H'
   ! Version info of this program
   character(len=*), parameter :: cver = 'V0.3'

   type(year_info) :: rt
   integer            :: k,ifile, i,j, i2, i3, i4, kdm
   real, allocatable, dimension(:,:,:) :: hy3d, regu3d, regusp3d, depthint, &
      hy3d2,pres, levint, temp, sal, dens
!AS06092011 - adding biological variables for MyOcean
   real, allocatable, dimension(:,:,:) :: fla, dia, nit, pho, oxy, pp, biovar,s1000
   real, allocatable, dimension(:,:,:) :: micro, meso, sil
   real, allocatable, dimension(:,:)   :: biovar2d
   real, allocatable, dimension(:,:)   :: hy2d, hy2d2, regu2d, strmf, &
      mld1, mld2, dplayer, meanssh, sla, ub, vb, mqlon, mqlat
   real, allocatable, dimension(:) :: tmpx, tmpy
   character(len=80)  :: ncfil, fname, ftype,vintmethod      ! Name of new netcdf file
   type(hycomfile) :: hfile

   type (file_state) :: ncstate
   logical :: usespline=.true., listmode
   logical :: lerr, msherr
   type(fields) :: fld(1000)
   integer :: nfld, ifld, filestart
   logical :: ex
   logical :: lev
#if defined (IARGC)
   integer*4, external :: iargc
#endif

   !Turn on/off the output of levitus climatology
   lev=.false.


   call process_arguments(vintmethod,listmode,filestart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff. 3 infiles are required:
! depthlevels.in  -- Depth levels used 
! proj.in         -- Specification of regular grid

   ! Init file io
   call xcspmd()
   call zaiost()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get bathymetry - lon/lat fields and sets up bigrid
   call get_grid ()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call proj_ini() ! Retrieves prjection from from proj.in and inits

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call spline_calc_ini_fromfile(trim(vintmethod))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ready to step through files to process

   print *
   do ifile=filestart,iargc()
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read raw data


      call getarg(ifile,fname)
      inquire(exist=ex,file=trim(fname))
      if (.not.ex) then
         print *,'File '//trim(fname)//' does not exist. Skipping this file'
         cycle
      end if

      ! What file is this? (daily? weekly? restart? archv?)
      ftype=trim(getfiletype(fname))

      ! Inits hycom file type
      call initHF(hfile,trim(fname),trim(ftype))
      !write (6,'(a)',advance='yes') 'Reading from '//trim(getHycomfileName(hfile))

      ! Vertical dim
      kdm=vDim(hfile)

      call fields_to_plot(fld,nfld,hfile,kdm)

      ! Start allocating fileds now..
      allocate(meanssh  (idm,jdm))        ! Holds meanssh on original grid
      allocate(sla      (idm,jdm))        ! Holds sla on original grid
      allocate(dplayer  (idm,jdm))        ! Holds layer thickness on original grid
      allocate(pres     (idm,jdm,kdm+1))  ! Holds pressure  on original grid
      allocate(hy3d     (idm,jdm,kdm))    ! Holds 3D vars on original grid
      allocate(hy3d2    (idm,jdm,kdm))    ! Holds pressure 
      allocate(regu3d   (nxp,nyp,kdm))    ! 3D vars on projection grid, original vertical grid
      allocate(regusp3d (nxp,nyp,ndeep))  ! 3D vars on projection grid, z-level vertical grid
      allocate(strmf    (idm,jdm))        ! 2D stream func 
      allocate(mld1     (idm,jdm))        ! MLD1
      allocate(mld2     (idm,jdm))        ! MLD2
      allocate(hy2d     (idm,jdm))        ! Holds 2D vars on original grid
      allocate(hy2d2    (idm,jdm))        ! Holds 2D vector vars on original grid
      allocate(regu2d   (nxp,nyp))        ! 2D vars on projection grid
      allocate(depthint (nxp,nyp,kdm))    ! Keeps bottom interface on regu grid


      call construct_filename(hfile,ncfil)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Start NC definition section 
      call OpenNCFile(ncstate,trim(ncfil))
 
      ! Define dimensions -- x,y,z  and possibly time
      call DefNCDim(ncstate,ndeep)

      ! Define header -- put some vital attributes
      call DefNCHeader(ncstate,hfile, cver)
 
      ! Set up projection vars
      call defNCProjectionVar(ncstate)
 
      ! Set up time variable
      call NCTimeVar(ncstate,hfile)
 
      ! Define space coordinates
      where (lons<-180.) lons=lons+360.
      where (lons> 180.) lons=lons-360.
      call NCSpaceVar(ncstate,deeps,ndeep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dump model bathymetry
      call to_proj(depths,regu2d)
      call putNCVar(ncstate,regu2d,nxp,nyp,1,'model_depth',2,.false.)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get layer interfaces
      pres(:,:,1)=0.

      do k=1,kdm
         !call HFReadDPField(hfile,dplayer,idm,jdm,k,1)

         ! Ensure dp is in pressure coordinates
         call HFReadDPField_p(hfile,dplayer,idm,jdm,k,1)

         pres(:,:,k+1)=pres(:,:,k)+dplayer
      end do

      ! At this stage, pressures (pres) should be in pressure coords, 
      call to_proj(pres(:,:,2:kdm+1)/onem,depthint,kdm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do ifld=1,nfld



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Special case (and rules) for vectors -
         ! 1) No scalar variables should begin with u,v,taux or tauy
         ! 2) v-component must come immediately after u-component
         ! See also checks in  m_fields_to_plot
         if (fld(ifld)%option .and. fld(ifld)%vecflag) then

            if (is3DVar(hfile,fld(ifld)%fextract,1)) then
               call HFReadField3D(hfile,hy3d ,idm,jdm,kdm,fld(ifld  )%fextract,1)
               call HFReadField3D(hfile,hy3d2,idm,jdm,kdm,fld(ifld+1)%fextract,1)

               print '(a)','Processing 3D Vector pair '// fld(ifld  )%fextract//' '//fld(ifld+1)%fextract


               if (trim(cprojection)/='native') then
                  if (gridrotate) then
                     do k=1,kdm
                        call rotate_general(hy3d(:,:,k),hy3d2(:,:,k),   yproj,xproj,idm,jdm,'m2l')
                     end do
                  else
                     do k=1,kdm
                        call rotate(hy3d(:,:,k),hy3d2(:,:,k),   plat,plon,idm,jdm,'m2l')
                     end do
                  end if
               end if

               ! NB - this prevents horizontal interpolation to thin layers
               ! TODO: This may mess up calculations near the sea floor
               do k=1,kdm
                  where (pres(:,:,k+1)-pres(:,:,k)<.1*onem)
                     hy3d(:,:,k)=undef
                     hy3d2(:,:,k)=undef
                  end where
               end do
!KAL 20151111  Fill up thin layers with data from above. Better, but not implemented yet (may influence mask)
!               do k=2,kdm
!                  where (pres(:,:,k+1)-pres(:,:,k)<.1*onem)
!                     hy3d(:,:,k)=hy3d(:,:,k-1)
!                     hy3d2(:,:,k)=hy3d2(:,:,k-1)
!               end do

               call to_proj(hy3d (:,:,1:kdm),regu3d,kdm)
               call spline_calc(regu3d,depthint,nxp,nyp,ongrid,regusp3d,ndeep,kdm)
               call putNCVar(ncstate,regusp3d,nxp,nyp,ndeep,fld(ifld)%fextract,4,gridrotate)
               call to_proj(hy3d2(:,:,1:kdm),regu3d,kdm)
               call spline_calc(regu3d,depthint,nxp,nyp,ongrid,regusp3d,ndeep,kdm)
               call putNCVar(ncstate,regusp3d,nxp,nyp,ndeep,fld(ifld+1)%fextract,4,gridrotate)

            ! 2D vector case
            else 
               print '(a)','Processing 2D Vector pair '// fld(ifld  )%fextract//' '//fld(ifld  )%fextract
               call HFReadField(hfile,hy2d ,idm,jdm,fld(ifld  )%fextract,0,1)
               call HFReadField(hfile,hy2d2,idm,jdm,fld(ifld+1)%fextract,0,1)
               if (trim(cprojection)/='native') then
                  if (gridrotate) then 
                     call rotate_general(hy2d,hy2d2,yproj,xproj,idm,jdm,'m2l')
                  else 
                     call rotate(hy2d,hy2d2,plat,plon,idm,jdm,'m2l')
                  end if
               end if
               call to_proj(hy2d,regu2d)
               call putNCVar(ncstate,regu2d,nxp,nyp,1,fld(ifld  )%fextract,3,gridrotate)
               call to_proj(hy2d2,regu2d)
               call putNCVar(ncstate,regu2d,nxp,nyp,1,fld(ifld+1)%fextract,3,gridrotate)
            end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! 3D Scalar  case. 
         else if (fld(ifld)%option .and. (.not. fld(ifld)%vecflag)) then
            if (is3DVar(hfile,fld(ifld)%fextract,1)) then
               print '(a)','Processing 3D Scalar  '// fld(ifld  )%fextract
!AS06092011 - adding biological variables for MyOcean
               if (trim(fld(ifld)%fextract)=='chla') then 
                  ! Compute chlorophyll a (kg m-3)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'DIA     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'FLA     ',1)
                  call chlorophyll(fla,dia,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)
               else if (trim(fld(ifld)%fextract)=='chla_eco') then    
                  ! Compute chlorophyll a (kg m-3)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'CHLDIA     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'CHLFLA     ',1)
                  call chlorophyll_eco(fla,dia,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)

               else if (trim(fld(ifld)%fextract)=='attcoef') then 
                  ! Compute attenuation coefficient (m-1)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'DIA     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'FLA     ',1)
                  call attenuation(dia,fla,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)
               else if (trim(fld(ifld)%fextract)=='attc_eco') then
                  ! Compute attenuation coefficient (m-1)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'CHLDIA     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'CHLFLA     ',1)
                  call attenuation_eco(dia,fla,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)

               else if (trim(fld(ifld)%fextract)=='nitrate') then 
                  ! Compute nitrate (mole m-3)
                  allocate(nit(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,nit,idm,jdm,kdm,'NIT     ',1)
                  call nitrate_conv(nit,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(nit,biovar)
               else if (trim(fld(ifld)%fextract)=='nit_eco') then
                  ! Compute nitrate (mole m-3)
                  allocate(nit(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,nit,idm,jdm,kdm,'NIT_eco     ',1)
                  call nitrate_conv_eco(nit,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(nit,biovar)

               else if (trim(fld(ifld)%fextract)=='phosphat') then 
                  ! Compute phosphate (mole m-3)
                  allocate(pho(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,pho,idm,jdm,kdm,'PHO     ',1)
                  call phosphate_conv(pho,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(pho,biovar)
               else if (trim(fld(ifld)%fextract)=='pho_eco') then
                  ! Compute phosphate (mole m-3)
                  allocate(pho(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,pho,idm,jdm,kdm,'PHO_eco     ',1)
                  call phosphate_conv_eco(pho,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(pho,biovar)

               else if (trim(fld(ifld)%fextract)=='sil_eco') then
                  ! Compute silicate (mole m-3)
                  allocate(sil(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,sil,idm,jdm,kdm,'SIL_eco     ',1)
                  call silicate_conv_eco(sil,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(sil,biovar)

               else if (trim(fld(ifld)%fextract)=='pbiomass') then 
                  ! Compute phytoplankton biomass (mole N m-3)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'DIA     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'FLA     ',1)
                  call pbiomass(dia,fla,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)
               else if (trim(fld(ifld)%fextract)=='pbio_eco') then
                  ! Compute phytoplankton biomass (mole C m-3)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'DIA_eco     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'FLA_eco     ',1)
                  call pbiomass_eco(dia,fla,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)

                else if (trim(fld(ifld)%fextract)=='zbiomass') then 
                  ! Compute zooplankton biomass (mole N m-3)
                  allocate(micro(idm,jdm,kdm))
                  allocate(meso(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,micro,idm,jdm,kdm,'MICRO   ',1)
                  call HFReadField3D(hfile,meso,idm,jdm,kdm,'MESO    ',1)
                  call pbiomass(micro,meso,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(micro,meso,biovar)
                else if (trim(fld(ifld)%fextract)=='zbio_eco') then
                  ! Compute zooplankton biomass (mole C m-3)
                  allocate(micro(idm,jdm,kdm))
                  allocate(meso(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,micro,idm,jdm,kdm,'MICR_eco   ',1)
                  call HFReadField3D(hfile,meso,idm,jdm,kdm,'MESO_eco    ',1)
                  call pbiomass_eco(micro,meso,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(micro,meso,biovar)

               else if (trim(fld(ifld)%fextract)=='oxygen') then 
                  ! Compute dissolved oxygen (mole m-3)
                  allocate(oxy(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,oxy,idm,jdm,kdm,'OXY     ',1)
                  call oxygen_conv(oxy,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(oxy,biovar)
               else if (trim(fld(ifld)%fextract)=='oxy_eco') then
                  ! Compute dissolved oxygen (mole m-3)
                  allocate(oxy(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,oxy,idm,jdm,kdm,'OXY_eco     ',1)
                  call oxygen_conv_eco(oxy,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(oxy,biovar)
               else if (trim(fld(ifld)%fextract)=='salt1000') then 
                  ! Compute salinity as parts psu/1000 (for the oil drift sumlations
                  allocate(sal(idm,jdm,kdm))
                  allocate(s1000(idm,jdm,kdm))
                  call HFReadField3D(hfile,sal,idm,jdm,kdm,'saln    ',1)
                  s1000=sal/1000.0
                  hy3d=s1000
                  deallocate(sal,s1000)
! _FABM__caglar_
               else if (trim(fld(ifld)%fextract)=='chl_fabm') then
                  ! Compute chlorophyll a (kg m-3)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'ECO_diac     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'ECO_flac     ',1)
                  call chlorophyll_fabm(fla,dia,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)
               else if (trim(fld(ifld)%fextract)=='nit_fabm') then
                  ! Compute nitrate (mole m-3)
                  allocate(nit(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,nit,idm,jdm,kdm,'ECO_no3     ',1)
                  call nitrate_conv_fabm(nit,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(nit,biovar)
               else if (trim(fld(ifld)%fextract)=='sil_fabm') then
                  ! Compute silicate (mole m-3)
                  allocate(sil(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,sil,idm,jdm,kdm,'ECO_sil     ',1)
                  call silicate_conv_fabm(sil,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(sil,biovar)
               else if (trim(fld(ifld)%fextract)=='pho_fabm') then
                  ! Compute phosphate (mole m-3)
                  allocate(pho(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,pho,idm,jdm,kdm,'ECO_pho     ',1)
                  call phosphate_conv_fabm(pho,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(pho,biovar)
               else if (trim(fld(ifld)%fextract)=='pbiofabm') then
                  ! Compute phytoplankton biomass (mmoleC m-3)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'ECO_dia     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'ECO_fla     ',1)
                  call pbiomass_fabm(dia,fla,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)
               else if (trim(fld(ifld)%fextract)=='zbiofabm') then
                  ! Compute zooplankton biomass (mmole C m-3)
                  allocate(micro(idm,jdm,kdm))
                  allocate(meso(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,micro,idm,jdm,kdm,'ECO_micr   ',1)
                  call HFReadField3D(hfile,meso,idm,jdm,kdm,'ECO_meso    ',1)
                  call pbiomass_fabm(micro,meso,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(micro,meso,biovar)
               else if (trim(fld(ifld)%fextract)=='oxy_fabm') then
                  ! Compute dissolved oxygen (mmole m-3)
                  allocate(oxy(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,oxy,idm,jdm,kdm,'ECO_oxy     ',1)
                  call oxygen_conv_fabm(oxy,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(oxy,biovar)
               else if (trim(fld(ifld)%fextract)=='prmpfabm') then
                  ! Compute gross primary production (mg m-3 d-1)
                  allocate(pp(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,pp,idm,jdm,kdm,'ECO_prim',1)
                  call pp_conv_fabm(pp,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(pp,biovar)
               else if (trim(fld(ifld)%fextract)=='attcfabm') then
                  ! Compute attenuation coefficient (m-1)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar(idm,jdm,kdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'ECO_diac     ',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'ECO_flac     ',1)
                  call attenuation_fabm(dia,fla,biovar,idm,jdm,kdm)
                  hy3d=biovar
                  deallocate(dia,fla,biovar)
! _FABM__caglar_
               else  ! LB normal case 
                  call HFReadField3D(hfile,hy3d,idm,jdm,kdm,fld(ifld)%fextract,1)
               end if

               ! NB - this prevents horizontal interpolation to thin layers
               ! TODO: This may mess up calculations near the sea floor
               do k=1,kdm
                  where (pres(:,:,k+1)-pres(:,:,k)<.1*onem) hy3d(:,:,k)=undef
               end do
!KAL 20151111  Fill up thin layers with data from above. Better, but not implemented yet (may influence mask)
!               do k=2,kdm
!                  where (pres(:,:,k+1)-pres(:,:,k)<.1*onem) hy3d(:,:,k)=hy3d(:,:,k-1)
!               end do

               call to_proj(hy3d(:,:,1:kdm),regu3d,kdm)
               call spline_calc(regu3d,depthint,nxp,nyp,ongrid,regusp3d,ndeep,kdm)
               call putNCVar(ncstate,regusp3d,nxp,nyp,ndeep,fld(ifld)%fextract,4,.false.)

            ! 2D Scalar case
            else 
               print '(a)','Processing 2D Scalar  '// fld(ifld  )%fextract
               if (trim(fld(ifld)%fextract)=='mld1' .or. trim(fld(ifld)%fextract)=='mld2') then 
                  ! Compute mix layer depths
                  allocate(temp(idm,jdm,kdm))
                  allocate(sal(idm,jdm,kdm))
                  call HFReadField3D(hfile,temp,idm,jdm,kdm,'temp    ',1)
                 !call HFReadField3D(hfile,sal ,idm,jdm,kdm,'saln    ',1)
                  call HFReadField3D(hfile,sal ,idm,jdm,kdm,'salin    ',1)
                  ! LB, pressure is expected in meters
                  call mixlayer_depths(temp,sal,pres/onem,mld1,mld2,idm,jdm,kdm)
                  deallocate(temp,sal)
                  if (trim(fld(ifld)%fextract)=='mld1') hy2d=mld1
                  if (trim(fld(ifld)%fextract)=='mld2') hy2d=mld2
                  
               elseif (trim(fld(ifld)%fextract)=='GS_MLD') then 
                  ! Compute the mixed layer depth interpolated
                  allocate(temp(idm,jdm,kdm))
                  allocate(sal(idm,jdm,kdm))
                  allocate(dens(idm,jdm,kdm))
                  allocate(biovar2D(idm,jdm))
                  call HFReadField3D(hfile,temp,idm,jdm,kdm,'temp',1)
                  call HFReadField3D(hfile,sal,idm,jdm,kdm,'salin',1)
                  !call HFReadField3D(hfile,sal,idm,jdm,kdm,'saln',1)
                  do i=1,idm
                    do j=1,jdm
                      do k=1,kdm
                        dens(i,j,k)=sig(temp(i,j,k),sal(i,j,k))
                      end do
                    end do
                  end do
                    print '(a)','calling gs_mld'
                  call gs_mld(-dens,pres/onem,biovar2d,idm,jdm,kdm,0.03)
                  hy2d=biovar2d
                  deallocate(temp,sal,dens,biovar2d)
               elseif (trim(fld(ifld)%fextract)=='bsf') then 
                  ! Compute barotropic streamfunction 
                  allocate(mqlon (0:idm+1,0:jdm+1))
                  allocate(mqlat (0:idm+1,0:jdm+1))
                  mqlon=qlon
                  mqlat=qlat
                  allocate(ub(idm,jdm))
                  allocate(vb(idm,jdm))
                  call HFReadField(hfile,ub,idm,jdm,'ubavg   ',0,1)
                  call HFReadField(hfile,vb,idm,jdm,'vbavg   ',0,1)
                  !call strmf_eval(idm,jdm,hy2d,ub,vb,depths,mqlat,mqlon)
                  call strmf_eval(idm,jdm,hy2d,ub,vb)
                  deallocate(mqlon,mqlat,ub,vb)
!AS06092011 - adding biological variables for MyOcean

               else if (trim(fld(ifld)%fextract)=='pp_depth') then 
                  ! Compute gross primary production (kg m-2 s-1)
                  allocate(pp(idm,jdm,kdm))
                  allocate(biovar2D(idm,jdm))
                  call HFReadField3D(hfile,pp,idm,jdm,kdm,'primprod',1)
                  call primary_production(pp,pres,biovar2d,onem,idm,jdm,kdm)
                  hy2d=biovar2d
                  deallocate(pp,biovar2d)
               else if (trim(fld(ifld)%fextract)=='pp_d_eco') then
                  ! Compute gross primary production (kg m-2 s-1)
                  allocate(pp(idm,jdm,kdm))
                  allocate(biovar2D(idm,jdm))
                  call HFReadField3D(hfile,pp,idm,jdm,kdm,'npp_eco',1)
                  call primary_production_eco(pp,pres,biovar2d,onem,idm,jdm,kdm)
                  hy2d=biovar2d
                  deallocate(pp,biovar2d)

              else if (trim(fld(ifld)%fextract)=='npp_euph') then 
                  ! Compute net primary production integrated over the euphotic depth (g m-2 day-1)
                  allocate(pp(idm,jdm,kdm))
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar2D(idm,jdm))
                  call HFReadField3D(hfile,pp,idm,jdm,kdm,'netpp',1)
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'DIA',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'FLA',1)
                  call net_primary_production(pp,dia,fla,pres,biovar2d,onem,idm,jdm,kdm)
                  hy2d=biovar2d
                  deallocate(pp,dia,fla,biovar2d)

              else if (trim(fld(ifld)%fextract)=='chl_opti') then 
                  ! Compute chlorophyll integrated over the euphotic depth(mg m-2)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar2D(idm,jdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'DIA',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'FLA',1)
                  call integrated_chlorophyll(dia,fla,pres,biovar2d,onem,idm,jdm,kdm)
                  hy2d=biovar2d
                  deallocate(dia,fla,biovar2d)
              else if (trim(fld(ifld)%fextract)=='chlo_eco') then
                  ! Compute chlorophyll integrated over the euphotic depth(mg m-2)
                  allocate(dia(idm,jdm,kdm))
                  allocate(fla(idm,jdm,kdm))
                  allocate(biovar2D(idm,jdm))
                  call HFReadField3D(hfile,dia,idm,jdm,kdm,'CHLDIA',1)
                  call HFReadField3D(hfile,fla,idm,jdm,kdm,'CHLFLA',1)
                  call integrated_chlorophyll_eco(dia,fla,pres,biovar2d,onem,idm,jdm,kdm)
                  hy2d=biovar2d
                  deallocate(dia,fla,biovar2d)

!KAL20151109 - Adding bottom temperature as a 2D field. Vertical interpolation to 10 meter above seabed
              else if (trim(fld(ifld)%fextract)=='btemp') then 
                  call HFReadField3D(hfile,hy3d,idm,jdm,kdm,'temp    ',1)
                  call spline_calc(hy3d,pres(:,:,2:kdm+1)/onem,idm,jdm,pres(:,:,kdm+1)/onem>10., hy2d,1,kdm,deepsin=(/-10./))
              else if (trim(fld(ifld)%fextract)=='bsaln') then 
                  call HFReadField3D(hfile,hy3d,idm,jdm,kdm,'saln    ',1)
                  call spline_calc(hy3d,pres(:,:,2:kdm+1)/onem,idm,jdm,pres(:,:,kdm+1)/onem>10., hy2d,1,kdm,deepsin=(/-10./))
!KAL20151109 - End adding bottom temperature
               else  ! LB normal case 
                  call HFReadField(hfile,hy2d,idm,jdm,fld(ifld)%fextract,0,1)
               endif

               ! Convention/unit issue for evaporation minus precipitation
               if (trim(fld(ifld)%fextract)=='emnp') hy2d=-hy2d*1000

               call to_proj(hy2d,regu2d)
               call putNCVar(ncstate,regu2d,nxp,nyp,1,fld(ifld)%fextract,3,.false.)
            end if
         end if


      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NEW - Levitus fields
     if (lev) then
         call levitus_setup('WOA2005')
         call forecastdate(hfile,rt)
         print '(a)','Levitus Salinity'
         allocate(levint(idm,jdm,ndeep))
         call levitus_interp3d('salinity',real(rt%idd),plon,plat,deeps, &
                                  levint, idm,jdm,ndeep)
         do k=1,ndeep
             call to_proj(levint(:,:,k),regu2d(:,:))
             call putNCVar(ncstate,regu2d,nxp,nyp,k,'levsaln',3,.false.)
         end do
         print '(a)','Levitus temperature'
         call levitus_interp3d('temperature',real(rt%idd),plon,plat,deeps, &
                                levint, idm,jdm,ndeep)
         do k=1,ndeep
            call to_proj(levint(:,:,k),regu2d(:,:))
            call putNCVar(ncstate,regu2d,nxp,nyp,k,'levtemp',3,.false.)
         end do
      end if

     call closeNCFile(ncstate)
     print *
     print '(a)','Generated '//trim(ncfil)
     print *
     print *


     ! Start deallocating fileds now..
     deallocate(sla     )
     deallocate(meanssh )
     deallocate(dplayer )
     deallocate(hy3d    )
     deallocate(hy3d2   )
     deallocate(pres   )
     deallocate(regu3d  )
     deallocate(regusp3d)
     deallocate(strmf   )
     deallocate(mld1    )
     deallocate(mld2    )
     deallocate(hy2d    )
     deallocate(hy2d2   )
     deallocate(regu2d  )
     deallocate(depthint)

     !call zaioempty ! deallocated za - fields - ready for next file

  enddo ! ifile-loop


end program p_hyc2proj
