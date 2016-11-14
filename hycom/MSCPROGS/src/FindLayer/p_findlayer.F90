! Routine changed to low-memory-version for use on cluster machines

program p_mosf
   use mod_xc
   use mod_za
   use mod_grid
   use mod_hycomfile_io
   use m_handle_err
   use netcdf
   implicit none

   character(len=80) :: fname, ncfil,ftype
   CHARACTER(len=3)  :: runcode              
   character(len=8)   :: dateinfo
   integer :: j,i,ios,g_idm, g_jdm, g_kdm,k,ktop,kbot,kdm
   real    :: rtime

   real, allocatable, dimension(:,:)   ::  botintf, topintf,meanval, laythick, dp, &
      var3d,dpsum,tem,sal,ut,vt, v3dsum, temsum, salsum, utsum, vtsum
   integer, allocatable, dimension(:,:):: layerflag

   integer :: nlats, ndeeps

   integer :: z_id, lat_id, merht_id, strf_id, idm_dimension_id, ncid, &
              k_dimension_id, vars_2d(2), jdm_dimension_id, lon_id, top_id, &
              thk_id, val_id , thk2_id, dpt_id, ut_id, vt_id, tem_id, sal_id

   character(len=10) :: tmparg
   character(len= 8) :: variable
   character(len=40) :: varunit,cmaxval,cminval,standard_name
   real :: maxdpsum, cfac, minvalue, maxvalue
   type(hycomfile) :: hfile
#if defined(IARGC)
   integer*4, external :: iargc
#endif
   real :: thref
   include 'stmt_funcs.H'


   if (iargc() /= 4 ) then
      print '(a)', '*******************************************************************'
      print '(a)', '*Produces data for layers above and below the variable tresholds. *'
      print '(a)', '*Top interface of this layer and thickness of layer is dumped in  *'
      print '(a)', '*a netcdf file named LAYERTHK.nc.                                 *'
      print '(a)', '*Note that only the first layer from the surface which matches the*'
      print '(a)', '*criterion is dumped in the netcdf file.                          *'
      print '(a)', '*                                                                 *'
      print '(a)', '*Usage: findlayer VARIABLE minimumvalue maximumvalue file         *'
      print '(a)', '*  where "VARIABLE" is 3D var in  file (eg saln)                  *'
      print '(a)', '*  where minimumval is minimum of variable range                  *'
      print '(a)', '*  where maximumval is maximum of variable range                  *'
      print '(a)', '*  where VALUE is treshold value                                  *'
      print '(a)', '*NB: specifying th3d as variable will use density as variable     *'
      print '(a)', '*******************************************************************'
      print '(a)'
      print '(a)','Example: '
      print '(a)','   findlayer saln 35 35.2 FORDAILY_2008_302_2008_318.a'
      print '(a)'
      print '(a)','Outputs netcdf file with layer thickness and top surface of layer'
      stop
   else
      call getarg(1,variable)
      call getarg(2,tmparg) ; read(tmparg,*) minvalue
      call getarg(3,tmparg) ; read(tmparg,*) maxvalue
      call getarg(4,fname) 
   end if







   ! use first file as base info
   ftype=getfiletype(trim(fname))
   call initHF(hfile,trim(fname),trim(ftype))
   kdm=vDim(hfile)
   rtime=hfile%fyear
   print *,trim(fname)//' header read'

   ! Init hycom IO , init grid
   call xcspmd()
   call zaiost()
   call get_grid()
   print *,'dimensions are',idm,jdm,kdm


   ! Read UT, VT, Tem, sal, dp 
   allocate(layerflag(idm,jdm))
   allocate(topintf(idm,jdm))
   allocate(botintf(idm,jdm))
   allocate(laythick  (idm,jdm))
   allocate(var3d(idm,jdm))
   allocate(dp(idm,jdm))
   allocate(tem(idm,jdm))
   allocate(sal(idm,jdm))
   allocate(ut(idm,jdm))
   allocate(vt(idm,jdm))
   allocate(v3dsum(idm,jdm))
   allocate(dpsum(idm,jdm))
   allocate(temsum(idm,jdm))
   allocate(salsum(idm,jdm))
   allocate(utsum(idm,jdm))
   allocate(vtsum(idm,jdm))


   g_idm=idm
   g_jdm=jdm
   g_kdm=kdm


   print *,trim(fname)//' header read'

   ! Sanity check
   if (idm/=g_idm .or.  jdm/=g_jdm .or.  kdm/=g_kdm ) then
      print *,'Files have different dimensions !'
      stop
   end if

   topintf=1e14
   botintf=-1e14
   layerflag=0
   laythick=0
   salsum=0.
   utsum=0.
   vtsum=0.
   temsum=0.
   dpsum=0.
   v3dsum=0.
   do k=1,kdm
     

      !call HFReadDPField(hfile,dp,idm,jdm,k,1)
      !dp=dp/onem
      call HFReadDPField_m(hfile,dp,idm,jdm,k,1) ! Returns in meters
      call HFReadField(hfile,sal  ,idm,jdm,'saln    ',k,1) ! nb fails with archv
      call HFReadField(hfile,tem  ,idm,jdm,'temp    ',k,1) ! nb fails with archv
      call HFReadField(hfile,ut   ,idm,jdm,'utot    ',k,1) ! nb fails with archv/restart
      call HFReadField(hfile,vt   ,idm,jdm,'vtot    ',k,1) ! nb fails with archv/restart
      if (trim(variable)/='th3d') then
         call HFReadField(hfile,var3d,idm,jdm,variable  ,k,1)
      else
         do j=1,jdm
         do i=1,idm
            var3d(i,j)=sig(tem(i,j),sal(i,j))
         end do
         end do
      end if

      do j=1,jdm
      do i=1,idm
         ! Check that designated 3d var is within range
         if (var3d(i,j)>minvalue .and. var3d(i,j)<maxvalue .and. &
             layerflag(i,j)>=0) then
             if (layerflag(i,j)==0) layerflag(i,j)=1
             topintf(i,j)=min(topintf(i,j),dpsum(i,j))
             salsum(i,j)=salsum(i,j)+dp(i,j)*sal(i,j)
             temsum(i,j)=temsum(i,j)+dp(i,j)*tem(i,j)
             utsum (i,j)=utsum (i,j)+dp(i,j)*ut (i,j)
             vtsum (i,j)=vtsum (i,j)+dp(i,j)*vt (i,j)
             v3dsum(i,j)=v3dsum(i,j)+dp(i,j)*var3d(i,j)
             botintf(i,j)=max(botintf(i,j),dpsum(i,j)+dp(i,j))
             laythick(i,j)=botintf(i,j)-topintf(i,j)
             ! Skip more passes by setting layerflag < 0
          else if (layerflag(i,j)==1)  then
             layerflag(i,j)=-1
          end if
          dpsum (i,j)=dpsum (i,j)+dp(i,j)
       end do
       end do
   end do

   !correction factor if dp is in pressure coords
   cfac=maxval(dpsum)/maxval(depths)

   ! Average values now
   do j=1,jdm
   do i=1,idm
      if (laythick(i,j)/cfac>1) then
         salsum(i,j)= salsum(i,j)/laythick(i,j)
         temsum(i,j)= temsum(i,j)/laythick(i,j)
         utsum (i,j)= utsum (i,j)/laythick(i,j)
         vtsum (i,j)= vtsum (i,j)/laythick(i,j)
         v3dsum(i,j)= v3dsum(i,j)/laythick(i,j)
      else
         salsum(i,j)= undef
         temsum(i,j)= undef
         utsum (i,j)= undef
         vtsum (i,j)= undef
         v3dsum(i,j)= undef
      end if

      laythick(i,j)= laythick(i,j)/cfac
      topintf (i,j)= topintf (i,j)/cfac
      botintf (i,j)= botintf (i,j)/cfac

      if (topintf(i,j)>1e8) topintf(i,j)=0.
      botintf(i,j)=max(botintf(i,j),0.)
   end do
   end do



   ! guess variable units
   if     (trim(VARIABLE)=='saln') then
      varunit='1e-3'
      standard_name='sea_water_salinity'
   elseif (trim(VARIABLE)=='temp') then
      varunit='degrees_celsius'
      standard_name='sea_water_potential_temperature'
   elseif (trim(VARIABLE)=='utot') then
      varunit='m s-1'
      standard_name='NONE'
   elseif (trim(VARIABLE)=='th3d' ) then
      varunit='kg m-3'
      standard_name='sea_water_sigma_t'
   else
      varunit='UNKNOWN'
      standard_name='NONE'
   end if


         



   !Create netcdf file and fill it with data
   ncfil='LAYERTHK.nc'
   print '(a)','  Generating '//trim(ncfil)
   if (NF90_CREATE(trim(ncfil),NF90_CLOBBER,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(daily2mosf)'
   end if


   ! Define dimensions -- idm,jdm,kdm
   ! Following MERCATOR name definitions for netcdf var name
   !print *,'dims'
   call handle_err(NF90_DEF_DIM(ncid,'idm',idm ,idm_dimension_id   ))
   call handle_err(NF90_DEF_DIM(ncid,'jdm',jdm ,jdm_dimension_id     ))
   vars_2d       =(/idm_dimension_id,jdm_dimension_id/)

   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title', 'BASTRF'))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution',        &
      'NERSC, Edvard Griegs Vei 3A, Bergen, Norway'))
   call date_and_time(date=dateinfo)
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Dataset_generation_date', &
      dateinfo(1:8)))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'history', &
      'Created '//dateinfo(1:8)//' by program pak2bastrf '))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'references', &
      'http://topaz.nersc.no'))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'field_type', &
         trim(ftype)))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions', &
      'COARDS'))


   ! Latitude
   call handle_err(NF90_DEF_VAR(ncid,'latitude',NF90_REAL,vars_2d,lat_id))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'long_name','latitude'))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'units','degrees_north'))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'valid_range',(/-90.0, 90.0/)))
   
   ! Longitude
   call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_REAL,vars_2d,lon_id))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'long_name','longitude'))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'units','degrees_east'))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'valid_range',(/-180.0, 180.0/)))
   
   ! Depth
   call handle_err(NF90_DEF_VAR(ncid,'depths',NF90_REAL,vars_2d,dpt_id))
   call handle_err(NF90_PUT_ATT(ncid,dpt_id,'long_name','longitude'))
   call handle_err(NF90_PUT_ATT(ncid,dpt_id,'units','m'))

   !print *,'top'
   call handle_err(NF90_DEF_VAR(ncid,'top',NF90_REAL,vars_2d,top_id))
   call handle_err(NF90_PUT_ATT(ncid,top_id,'long_name','Layer top'))
   call handle_err(NF90_PUT_ATT(ncid,top_id,'units','m'))
   call handle_err(NF90_PUT_ATT(ncid,top_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,top_id,'_FillValue',real(undef,kind=4)))
   
   !print *,'thk'
   write(cminval,'(e12.4)')minvalue
   write(cmaxval,'(e12.4)')maxvalue
   call handle_err(NF90_DEF_VAR(ncid,'thk',NF90_REAL,vars_2d,thk_id))
   call handle_err(NF90_PUT_ATT(ncid,thk_id,'long_name', &
      'Thickness  of layers with '//trim(VARIABLE)//' in range '// &
      trim(cminval)//' '//trim(cmaxval)))
   call handle_err(NF90_PUT_ATT(ncid,thk_id,'units','m'))
   call handle_err(NF90_PUT_ATT(ncid,thk_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,thk_id,'_FillValue',real(undef,kind=4)))

   !print *,'val'
   call handle_err(NF90_DEF_VAR(ncid,'val',NF90_REAL,vars_2d,val_id))
   call handle_err(NF90_PUT_ATT(ncid,val_id,'long_name', &
      'Average value for layers with '//trim(VARIABLE)//' in range '// &
      trim(cminval)//' '//trim(cmaxval)))
   call handle_err(NF90_PUT_ATT(ncid,val_id,'units',trim(varunit)))
   call handle_err(NF90_PUT_ATT(ncid,val_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,val_id,'_FillValue',real(undef,kind=4)))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'standard_name', &
                                trim(standard_name)))

   !print *,'sal'
   call handle_err(NF90_DEF_VAR(ncid,'salinity',NF90_REAL,vars_2d,sal_id))
   call handle_err(NF90_PUT_ATT(ncid,sal_id,'long_name', &
      'Average value for layers with '//trim(VARIABLE)//' in range '// &
      trim(cminval)//' '//trim(cmaxval)))
   call handle_err(NF90_PUT_ATT(ncid,sal_id,'units',trim(varunit)))
   call handle_err(NF90_PUT_ATT(ncid,sal_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,sal_id,'_FillValue',real(undef,kind=4)))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'standard_name', &
                               'sea_water_salinity'))

   !print *,'temperature'
   call handle_err(NF90_DEF_VAR(ncid,'temperature',NF90_REAL,vars_2d,tem_id))
   call handle_err(NF90_PUT_ATT(ncid,tem_id,'long_name', &
      'Average value for layers with '//trim(VARIABLE)//' in range '// &
      trim(cminval)//' '//trim(cmaxval)))
   call handle_err(NF90_PUT_ATT(ncid,tem_id,'units',trim(varunit)))
   call handle_err(NF90_PUT_ATT(ncid,tem_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,tem_id,'_FillValue',real(undef,kind=4)))
   call handle_err(NF90_PUT_ATT(ncid,lon_id,'standard_name', &
                               'sea_water_potential_temperature'))

   !print *,'temperature'
   call handle_err(NF90_DEF_VAR(ncid,'ut',NF90_REAL,vars_2d,ut_id))
   call handle_err(NF90_PUT_ATT(ncid,ut_id,'long_name', &
      'Average value for layers with '//trim(VARIABLE)//' in range '// &
      trim(cminval)//' '//trim(cmaxval)))
   call handle_err(NF90_PUT_ATT(ncid,ut_id,'units',trim(varunit)))
   call handle_err(NF90_PUT_ATT(ncid,ut_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,ut_id,'_FillValue',real(undef,kind=4)))

   !print *,'temperature'
   call handle_err(NF90_DEF_VAR(ncid,'vt',NF90_REAL,vars_2d,vt_id))
   call handle_err(NF90_PUT_ATT(ncid,vt_id,'long_name', &
      'Average value for layers with '//trim(VARIABLE)//' in range '// &
      trim(cminval)//' '//trim(cmaxval)))
   call handle_err(NF90_PUT_ATT(ncid,vt_id,'units',trim(varunit)))
   call handle_err(NF90_PUT_ATT(ncid,vt_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,vt_id,'_FillValue',real(undef,kind=4)))

   !print *,'enddef...'
   call handle_err(NF90_ENDDEF(ncid))


   where (laythick<0.1) 
      topintf=undef
      laythick=undef
      v3dsum=undef
      salsum=undef
      temsum=undef
      utsum=undef
      vtsum=undef
   end where
   call handle_err(NF90_PUT_VAR(ncid,lon_id,plon  (1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,lat_id,plat  (1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,dpt_id,depths  (1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,top_id,topintf (1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,thk_id,laythick(1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,val_id,v3dsum (1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,sal_id,salsum (1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,tem_id,temsum (1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,ut_id ,utsum  (1:idm,1:jdm)))
   call handle_err(NF90_PUT_VAR(ncid,vt_id ,vtsum  (1:idm,1:jdm)))
   call handle_err(NF90_CLOSE(ncid))


end program p_mosf
