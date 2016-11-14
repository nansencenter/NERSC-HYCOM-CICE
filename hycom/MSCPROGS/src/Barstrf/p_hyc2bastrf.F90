program p_hyc2bastrf
   !----------------------------------------------
   !Routine calculates barotropic stream function 
   !from the barotropic velocity components. If 
   !you supply more than one file, the average 
   !barotropic streamfunction is calculated 
   !
   !Created by knut Liseter, 2004
   !----------------------------------------------
   use mod_xc
   use mod_za
   use mod_grid
   use mod_hycomfile_io
   use mod_year_info
   use m_strmf_eval
   use m_handle_err
   use netcdf
   implicit none
#if defined(IARGC) 
   integer*4, external :: iargc
#endif
   character(len=80) :: fname, ncfil, ftype
   character(len=8)   :: dateinfo
   integer :: j,i, iarg, ios
   real    :: rtime
   real, allocatable, dimension(:,:)  :: ubavg, vbavg
   real, dimension(:,:),allocatable :: strf,strf_ave

   ! File, dimension and variable ids in netcdf file
   integer :: lat_id,  idm_dimension_id, ncid, &
              vars_2d(2), jdm_dimension_id, lon_id, stf_id
   type(hycomfile) :: hfile


   if (iargc() < 1 ) then
      print *,'Routine calculates barotropic stream function from '
      print *,'the barotropic velocity components. If you supply '
      print *,'More than one file, the average barotropic '
      print *,'streamfunction is calculated '
      print *, 'If several files are given, the average is calculated'
      print *, 'Usage: hyc2bastrf [file(s)]'
      stop
   end if

   ! use first file as base info
   call getarg(1,fname)
   ftype=getfiletype(trim(fname))
   call initHF(hfile,trim(fname),trim(ftype))

   ! Init IO, and Get model grid & Depths
   call xcspmd()
   call zaiost()
   call get_grid()

   ! Allocate arrays
   allocate(ubavg(idm,jdm))
   allocate(vbavg(idm,jdm))
   allocate(strf    (idm,jdm))
   allocate(strf_ave(idm,jdm))

   iarg=1
   strf_ave=0.
   do while (iarg<= iargc())

      call getarg(iarg,fname)
      ftype=getfiletype(trim(fname))
      call initHF(hfile,trim(fname),trim(ftype))
      rtime=hfile%fyear
      print *,trim(fname)

      ! Get barotropic velocities
      call HFreaduvbaro(hfile,ubavg,vbavg,idm,jdm,1)
      call strmf_eval(idm,jdm,strf,ubavg,vbavg)
      strf_ave= strf_ave+strf
      iarg = iarg + 1
    end do
    strf_ave=strf_ave/(iarg-1)
    strf_ave=strf_ave*1e-6  ! -> Sverdrup
    where(depths<.1) strf_ave=undef


   !Create netcdf file and fill it with data
   ncfil='BASTRF.nc'
   print '(a)','  Generating '//trim(ncfil)
   if (NF90_CREATE(trim(ncfil),NF90_CLOBBER,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(daily2mosf)'
   end if

   ! 
   call handle_err(NF90_DEF_DIM(ncid,'idm',idm ,idm_dimension_id   ))
   call handle_err(NF90_DEF_DIM(ncid,'jdm',jdm ,jdm_dimension_id     ))
   vars_2d       =(/idm_dimension_id,jdm_dimension_id/)
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title', 'BASTRF'))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution', 'NERSC, Thormoehlensgt 47, Bergen, Norway'))
   call date_and_time(date=dateinfo)
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Dataset_generation_date', dateinfo(1:8)))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'history',  'Created '//dateinfo(1:8)//' by program hyc2bastrf '))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'references', 'http://topaz.nersc.no'))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'field_type', trim(ftype)))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions', 'COARDS'))

   
   ! Latitude
   call handle_err(NF90_DEF_VAR(ncid,'latitude',NF90_Float,vars_2d,lat_id))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'long_name','latitude'))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'units','degrees_north'))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'valid_range',(/-90.0, 90.0/)))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'standard_name','latitude'))
   
   ! Longitude
   call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_Float,vars_2d,lon_id))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'long_name','longitude'))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'units','degrees_east'))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'valid_range',(/-90.0, 90.0/)))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'standard_name','latitude'))
   
   ! Longitude
   call handle_err(NF90_DEF_VAR(ncid,'bastrf',NF90_Float,vars_2d,stf_id))
   call handle_err(NF90_PUT_ATT(ncid,stf_id,'long_name','Barotropic SF'))
   call handle_err(NF90_PUT_ATT(ncid,stf_id,'units','Sverdrup'))
   call handle_err(NF90_PUT_ATT(ncid,stf_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,stf_id,'_FillValue',real(undef,kind=4)))
   call handle_err(NF90_ENDDEF(ncid))


   call handle_err(NF90_PUT_VAR(ncid,lon_id,plon    ))
   call handle_err(NF90_PUT_VAR(ncid,lat_id,plat    ))
   call handle_err(NF90_PUT_VAR(ncid,stf_id,strf_ave))
   call handle_err(NF90_CLOSE(ncid))


end program p_hyc2bastrf
