!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine creates NetCDF files from (various) hycom files. The netcdf
! files hold the data on specified spatial locations. Routine 
! works on several HYCOM files; restart, archive, nersc weekly average and nersc
! daily average.
!
!
! 1) Specify positions to interpolate the model to. The positions are grouped 
! together by specifying a group name, which makes it possible to build points
! which are logically connected (e.g a section). This is done in the file
! "stations.in"
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
! hycom, everything should be ok
!
! 5) Run "hyc2stations" to generate netcdf files. At the command line you can
! specify files to process.  For each file,  netcdf file will be created. 
!
!
!NB: Sample infiles are given in the catalogue SAMPLE_INFILES : TODO
!
! 20090113 - KAL: Re-wrote the code. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program p_hyc2stations
   use mod_xc
   use mod_grid, only:  depths, plon, plat, get_grid
   use mod_year_info
   use mod_spline_calc
   use mod_netcdf_file
   use mod_station
   use mod_levitus
   use mod_rotate
   use mod_za , only: zaiost, zaioempty
   use mod_hycomfile_io
   use m_get_hyc_infiles
   use m_fields_to_plot
   use m_construct_filename
   implicit none

   ! Version info of this program
   character(len=*), parameter :: cver = 'V0.3'

   type(year_info)    :: rt
   integer            :: k,ifile,i,find,igroup, istat,nstat, i2, kdm

   real, allocatable, dimension(:,:,:) :: hy3d, hy3d2, pres
   real, allocatable, dimension(:,:) :: hy2d, hy2d2, dplayer, ubavg, vbavg
   real, allocatable, dimension(:,:,:) :: levint_3d
   real, allocatable, dimension(:,:,:,:) :: stat_3D, stat_3D_z
   real, allocatable, dimension(:,:,:)   :: stat_2D, stat_pres
   real, allocatable, dimension(:,:)     :: stat_lon,stat_lat,stat_depth

   character(len=80)  :: ncfil      ! Name of new netcdf file
   character(len=80)  :: ncfil2     ! Name of new netcdf file
   character(len=80)  :: fname     ! Name of new netcdf file
   character(len=80)  :: ftype     ! Name of new netcdf file

   real, allocatable, dimension(:,:)   :: tmpvar
   real, allocatable, dimension(:,:,:) :: stat_saln_lev, stat_temp_lev
   real, allocatable, dimension(:,:)   :: cspres
   real, allocatable, dimension(:)     :: csjuld
   type (file_state) :: ncstate
   type (hycomfile) :: hfile

   ! Input file with days
   type(fields) :: fld(1000)
   integer :: nfld, ifld
   logical :: ex
#if defined (IARGC)
   integer*4, external :: iargc
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Argument processing.  Exit if no arguments are specified
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff. 2 infiles are required:
! depthlevels.in -- Depth levels used 
! stations.in    -- Stations to use
   call spline_calc_ini_fromfile('spline') ! Retrieves depth levels too



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get bathymetry - lon/lat fields and sets up bigrid
   call xcspmd()
   call zaiost()
   call get_grid  ; where (depths > 1e10 ) depths=0.
   call get_stations()


                  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! This file will contain a list of data files produced by this routine
   open(96,file='hyc2stations.filelist',action='write',status='replace')

   ! Were ready...
   print *
   do ifile=1,iargc()
      
      call getarg(ifile,fname)
      inquire(exist=ex,file=trim(fname))
      if (.not.ex) then
         print *,'File '//trim(fname)//' does not exist. Skipping this file'
         cycle
      end if

      ! What file is this? (daily? weekly? restart? pak?)
      ftype=trim(getfiletype(fname))

      ! Inits hycom file type
      call initHF(hfile,trim(fname),trim(ftype))

      ! Check that grid sizes are in sync with those set by get_grid
      kdm=VDim(hfile)

      ! Init file io
      call zaiost()

      ! 
      call fields_to_plot(fld,nfld,hfile,kdm)

      !print *,ngroup,maxstat,kdm
      allocate(stat_pres (maxstat,kdm+1,ngroup))
      allocate(stat_lon  (maxstat,ngroup))
      allocate(stat_lat  (maxstat,ngroup))
      allocate(stat_depth(maxstat,ngroup))
      allocate(stat_3D  (maxstat,kdm,nfld,ngroup))
      allocate(stat_3D_z(maxstat,ndeep,nfld,ngroup))
      allocate(stat_2D (maxstat,nfld,ngroup))
      allocate(stat_saln_lev(maxstat,ndeep,ngroup))
      allocate(stat_temp_lev(maxstat,ndeep,ngroup))


      allocate(pres     (idm,jdm,kdm+1))      ! Holds 3D vars on original grid
      allocate(hy3d     (idm,jdm,kdm))        ! Holds 3D vars on original grid
      allocate(hy3d2    (idm,jdm,kdm))        ! Holds 3D vars on original grid
      allocate(hy2d     (idm,jdm))            ! Holds 2D vars on original grid
      allocate(hy2d2    (idm,jdm))            ! Holds 2D vars on original grid
      allocate(dplayer  (idm,jdm))            ! Holds 2D vars on original grid
      allocate(ubavg    (idm,jdm))            ! Holds 2D vars on original grid
      allocate(vbavg    (idm,jdm))            ! Holds 2D vars on original grid





      ! Days start from "0". Add 1 to get real ! dates
      call construct_filename(hfile,ncfil)

      ! Get layer interfaces
      pres(:,:,1)=0.
      do k=1,kdm
         call HFReadDPField(hfile,dplayer,idm,jdm,k,1)
         pres(:,:,k+1)=pres(:,:,k)+dplayer
      end do

      ! Interpolate pressures to station locations
      do igroup=1,ngroup
      do k=1,kdm+1
      do istat=1,stations_per_group(igroup)
         stat_pres(istat,k,igroup)=fld2d_to_station(pres(1,1,k),stations(igroup,istat))
      end do
      end do
      end do

      ! Interpolate positions to station locations
      do igroup=1,ngroup
      do istat=1,stations_per_group(igroup)
         stat_lon  (istat,igroup)=fld2d_to_station(plon  ,stations(igroup,istat))
         stat_lat  (istat,igroup)=fld2d_to_station(plat  ,stations(igroup,istat))
         stat_depth(istat,igroup)=fld2d_to_station(depths,stations(igroup,istat))
      end do
      end do


      ! Go through fields to extract and interpolate them to station locations
      ! TODO: rotate velocity
      do ifld=1,nfld

         ! Special case (and rules) for vectors -
         ! 1) No scalar variables should begin with u,v,taux or tauy
         ! 2) v-component must come immediately after u-component
         ! See also checks in  m_fields_to_plot
         if (fld(ifld)%fextract(1:1)=='u' .or. fld(ifld)%fextract(1:4)=='taux') then

            ! 3D vector to be rotated 
            if (is3Dvar(hfile,fld(ifld)%fextract,1)) then
               print *,fld(ifld)%fextract,'3D vector'
               call HFReadField3D(hfile,hy3d ,idm,jdm,kdm,fld(ifld  )%fextract,1)
               call HFReadField3D(hfile,hy3d2,idm,jdm,kdm,fld(ifld+1)%fextract,1)
               do k=1,kdm
                  print *,'Rotate ',fld(ifld)%fextract,fld(ifld+1)%fextract,k
                  call rotate(hy3d(:,:,k),hy3d2(:,:,k),   plat,plon,idm,jdm,'m2l')
               end do

               ! Interpolate onto model vertical grid
               do igroup=1,ngroup
               do k=1,kdm
               do istat=1,stations_per_group(igroup)
                  stat_3D(istat,k,ifld  ,igroup)=fld2d_to_station(hy3d (1,1,k),stations(igroup,istat))
                  stat_3D(istat,k,ifld+1,igroup)=fld2d_to_station(hy3d2(1,1,k),stations(igroup,istat))
               end do
               end do
               end do

               ! Spline to z-levels
               do igroup=1,ngroup
               do istat=1,stations_per_group(igroup)
               if (stat_depth(istat,igroup)>1.) then
                  call spline_calc_1d(stat_3D  (istat,:,ifld  ,igroup),stat_pres(istat,2:kdm+1,igroup)/onem, &
                                      stat_3D_z(istat,:,ifld  ,igroup),ndeep,kdm)
                  call spline_calc_1d(stat_3D  (istat,:,ifld+1,igroup),stat_pres(istat,2:kdm+1,igroup)/onem, &
                                      stat_3D_z(istat,:,ifld+1,igroup),ndeep,kdm)
               else
                  stat_3D_z(istat,:,ifld  ,igroup)=undef
                  stat_3D_z(istat,:,ifld+1,igroup)=undef
               end if
               end do
               end do

            ! 2D vector to be rotated 
            else
               print *,fld(ifld)%fextract,'2D vector'
               call HFReadField(hfile,hy2d ,idm,jdm,fld(ifld  )%fextract,0,1)
               call HFReadField(hfile,hy2d2,idm,jdm,fld(ifld+1)%fextract,0,1)
               call rotate(hy2d,hy2d2,plat,plon,idm,jdm,'m2l')
               do igroup=1,ngroup
               do istat=1,stations_per_group(igroup)
                  stat_2D(istat,ifld  ,igroup)=fld2d_to_station(hy2d (1,1),stations(igroup,istat))
                  stat_2D(istat,ifld+1,igroup)=fld2d_to_station(hy2d2(1,1),stations(igroup,istat))
               end do
               end do
            end if

         ! Scalar  case. Note that second vector component has been processed in vector case
         else if (fld(ifld)%fextract(1:1)/='v' .and. fld(ifld)%fextract(1:4)/='tauy') then


            ! 3D scalar to be interpolated
            if (is3Dvar(hfile,fld(ifld)%fextract,1)) then
               print *,fld(ifld)%fextract,'3D'
               call HFReadField3D(hfile,hy3d ,idm,jdm,kdm,fld(ifld  )%fextract,1)

               ! Interpolate from grid to station points
               do igroup=1,ngroup
               do k=1,kdm
               do istat=1,stations_per_group(igroup)
                  stat_3D(istat,k,ifld,igroup)=fld2d_to_station(hy3d(1,1,k),stations(igroup,istat))
               end do
               end do
               end do

               ! Spline to z-levels
               do igroup=1,ngroup
               do istat=1,stations_per_group(igroup)
               if (stat_depth(istat,igroup)>1.) then
                  call spline_calc_1d(stat_3D  (istat,:,ifld,igroup),stat_pres(istat,2:kdm+1,igroup)/onem, &
                                      stat_3D_z(istat,:,ifld,igroup),ndeep,kdm)
               else
                  stat_3D_z(istat,:,ifld,igroup)=undef
               end if
               end do
               end do

            !  2D scalar to be interpolated
            else 
               print *,fld(ifld)%fextract,'2D'
               call HFReadField(hfile,hy2d ,idm,jdm,fld(ifld  )%fextract,0,1)
               do igroup=1,ngroup
               do istat=1,stations_per_group(igroup)
                  stat_2D(istat,ifld,igroup)=fld2d_to_station(hy2d(1,1),stations(igroup,istat))
               end do
               end do
            end if
         end if
      end do




      ! Vertical interpolation of levitus data onto stations
      if (.true.) then
      call levitus_setup('WOA2005')
      stat_saln_lev=undef
      stat_temp_lev=undef
      do igroup=1,ngroup
         nstat=stations_per_group(igroup)
         allocate(csjuld(nstat)) ; 
         allocate(cspres(ndeep,nstat)) ; 
         allocate(tmpvar(ndeep,nstat)) ; 
         call forecastdate(hfile,rt)
         csjuld=rt%idd
         do k=1,ndeep
            cspres(k,:)=deeps(k)
         end do

         call stationsInterpLevitus('salinity', &
                             csjuld(1:nstat), &
                             stations(igroup,1:nstat)%lon, &
                             stations(igroup,1:nstat)%lat, &
                             cspres, &
                             tmpvar, &
                             nstat,ndeep)
         do k=1,ndeep
            stat_saln_lev(1:nstat,k,igroup)=tmpvar(k,:)
         end do
         !stop '(test)'

         call stationsInterpLevitus('temperature', &
                             csjuld(1:nstat), &
                             stations(igroup,1:nstat)%lon, &
                             stations(igroup,1:nstat)%lat, &
                             cspres, &
                             tmpvar, &
                             nstat,ndeep)
         do k=1,ndeep
            stat_temp_lev(1:nstat,k,igroup)=tmpvar(k,:)
         end do
         deallocate(csjuld)
         deallocate(cspres)
         deallocate(tmpvar)
      end do
      end if

!
         
      ! Fields are interpolated - produce files - one for each grouping
      do igroup = 1,ngroup

         find=index(ncfil,'.nc')
         ncfil2=ncfil(1:find-1)//'_group'//trim(groupname(igroup))//'.nc'
         nstat=stations_per_group(igroup)

         call openNCFile  (ncstate,trim(ncfil2))
         call DefNCHeader(ncstate,hfile, cver)
         call DefNCDimStat(ncstate,nstat,ndeep)
         call NCTimeVar(ncstate,hfile)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Start defining vars - put them as we move along
         call NCSpaceVarStat(ncstate,nstat,deeps,ndeep)


         ! longitude
         call putNCVar(ncstate, stations(igroup,1:nstat)%lon,nstat,1,1,'longitude',1,.false.)

         ! latitude
         call putNCVar(ncstate, stations(igroup,1:nstat)%lon,nstat,1,1,'latitude',1,.false.)

         ! model_depth
         call putNCVar(ncstate, stat_depth(1:nstat,igroup),nstat,1,1,'model_depth',1,.false.)
         
         ! Levitus salinity and temperature
         call putNCVar(ncstate, stat_temp_lev(1:nstat,:,igroup),nstat,ndeep,1,'levtemp',3,.false.)
         call putNCVar(ncstate, stat_saln_lev(1:nstat,:,igroup),nstat,ndeep,1,'levsaln',3,.false.)

         ! Loop through specified 2D and 3D variables
         do ifld=1,nfld
            if (is3Dvar(hfile,fld(ifld)%fextract,1)) then
               !print *,fld(ifld)%fextract,'3D'
               call putNCVar(ncstate, stat_3D_z(1:nstat,:,ifld,igroup),nstat,ndeep,1, &
                             fld(ifld)%fextract,3,.false.)
            else
               !print *,fld(ifld)%fextract,'2D'
               call putNCVar(ncstate, stat_2D(1:nstat,ifld,igroup),nstat,1,1, &
                             fld(ifld)%fextract,2,.false.)
            end if
         end do

         call closeNCFile(ncstate)


         ! Append to file list processed by this program
         write(96,'(a)') trim(ncfil2)


      end do ! igroup loop


     ! Start deallocating fileds now..
     deallocate(hy3d    )
     deallocate(hy3d2   )
     deallocate(hy2d    )
     deallocate(hy2d2   )
     deallocate(dplayer )
     deallocate(ubavg   )
     deallocate(vbavg   )
     deallocate(pres     )
     deallocate(stat_2D )
     deallocate(stat_3D )
     deallocate(stat_3D_z )
     deallocate(stat_lon)
     deallocate(stat_lat)
     deallocate(stat_depth)
     deallocate(stat_pres)

     call zaioempty ! deallocated za - fields - ready for next file

  enddo ! ifile-loop
  close(96)


end program p_hyc2stations
