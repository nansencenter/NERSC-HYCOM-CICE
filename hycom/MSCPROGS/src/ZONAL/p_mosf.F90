program p_mosf
   use mod_xc
   use mod_za
   use mod_grid
   !use m_mosf
   !use m_mht
   use mod_overturning
   use mod_hycomfile_io
   use m_handle_err
   use netcdf
   implicit none

   character(len=80) :: pakfile, ncfil,fname,ftype
   CHARACTER(len=3)  :: runcode              
   CHARACTER(len=5)  :: char5              
   character(len=8)   :: dateinfo



   integer*4, external :: iargc
   integer :: n2d,n3d,ii1,jj1,j,i, kdm, l, k
   integer :: iday,iyear,ihour,imonth,iweek,layer,n
   real    :: rtime
   integer :: thisday,idayinmonth,endyear
   integer :: ios
   integer :: ilat,iarg

   real, allocatable, dimension(:,:,:) :: utot, vtot, &
      temp, saln, dp, th3d
   real, allocatable, dimension(:,:)   :: fice, hice, &
      uice, vice,twod, dpsum, ubavg, vbavg

   real, dimension(:),allocatable :: lats, merht, merht_ave
   real, dimension(:,:),allocatable :: mostrf,mostrf_ave, mostrf_sig0_ave

   integer :: z_id, lat_id, merht_id, strf_id, lat_dimension_id, ncid, &
              k_dimension_id, vars_2d(2), sig_id, sig_strf_id, sig_dimension_id

   integer :: g_idm, g_jdm, g_kdm
   type(hycomfile) :: hfile

   real, dimension(:), allocatable :: deep
   integer :: ndeep, nlat
   real, parameter  ::dstep=50.
   real, parameter :: dxfac=3
   real :: mindx, dlat, minlat, maxlat
   real, allocatable :: sigarray(:), mostrf_sig0(:,:)
   integer :: nsig
   real :: minsig, maxsig, dsig
   real, external :: spherdist
   real :: thref

   include 'stmt_funcs.H'


   if (iargc() < 1 ) then
      print *, 'Produces merid. overturning sf'
      print *, 'If several files are given, the average is calculated'
      print *, 'Usage: pak2mosf [file(s)]'
      stop
   end if

   ! use first file as base info
   call getarg(1,fname)
   ftype=getfiletype(trim(fname))
   call initHF(hfile,trim(fname),trim(ftype))
   kdm=vDim(hfile)

   ! Get model grid & Depths
   call xcspmd()
   call zaiost()
   call get_grid()


   ! Read UT, VT, Tem, sal, dp 
   allocate(utot(idm,jdm,kdm))
   allocate(vtot(idm,jdm,kdm))
   allocate(temp(idm,jdm,kdm))
   allocate(saln(idm,jdm,kdm))
   allocate(dp  (idm,jdm,kdm))
   allocate(th3d(idm,jdm,kdm))
   allocate(fice(idm,jdm))
   allocate(hice(idm,jdm))
   allocate(uice(idm,jdm))
   allocate(vice(idm,jdm))
   allocate(ubavg(idm,jdm))
   allocate(vbavg(idm,jdm))
   allocate(twod(idm,jdm))
   allocate(dpsum(idm,jdm))


   ! min grid size
   mindx=min(minval(scpy),minval(scpx))*dxfac ! factor subject to tuning

   ! Convert min dx to latitude increments
   dlat=mindx/spherdist(0.,0.,0.,1.)
   minlat=minval(plat)
   maxlat=maxval(plat)
   ! Security check 1
   if (mindx<1000.) then
      print *,'Too low min grid spacing :',mindx
      stop
   ! KAL - Added this test because some older regional files can have 1e30 for max lat
   else if (abs(minlat) > 90. .or. abs(maxlat) > 90) then
      print *,'Security check, |lat| > 90  ... Abandon ship'
      stop
   end if

   ! Allocate zonal bands
   nlat=(maxlat-minlat)/dlat
   allocate(lats(nlat))
   do l=1,nlat
      lats(l)=minlat + (l-1)*dlat
   end do
   print *,dlat,mindx,minlat,maxlat,nlat

   ! Allocate depth bands
   ndeep= 5000./dstep
   allocate(deep(ndeep))
   do k=1,ndeep
      deep(k)=k*dstep
   end do

   ! Allocate density bands
   minsig=22.
   maxsig=29.
   dsig=0.03
   nsig=(maxsig-minsig)/dsig 
   allocate(sigarray(nsig))
   do k=1,nsig
      sigarray(k)=minsig+(k-1)*dsig
   end do

   ! Allocate containers for heat transport and overturning stream functions
   allocate(mostrf     (nlat,ndeep))
   allocate(mostrf_sig0(nlat,nsig))
   allocate(mostrf_sig0_ave(nlat,nsig))
   allocate(merht(nlat))
   allocate(mostrf_ave(nlat,ndeep))
   allocate(merht_ave (nlat))


   g_idm=idm
   g_jdm=jdm
   g_kdm=kdm

   merht_ave=0.
   mostrf_ave=0.
   mostrf_sig0_ave=0.
   iarg=1

   do while (iarg<= iargc())

      call getarg(iarg,fname)
      ftype=getfiletype(trim(fname))
      call initHF(hfile,trim(fname),trim(ftype))
      rtime=hfile%fyear

      ! Sanity check
      if (idm/=g_idm .or.  jdm/=g_jdm .or.  kdm/=g_kdm ) then
         print *,'Files have different dimensions !'
         stop
      end if


      ! Get SALINIT, TEMPERATURE, VELOCITY and ICE variables
      call HFReaduvbaro(hfile,ubavg,vbavg,idm,jdm,1)
      do k=1,kdm
         !call HFReadDPField(hfile,dp(:,:,k),idm,jdm,k,1)
         call HFReadDPField_p(hfile,dp(:,:,k),idm,jdm,k,1) ! returns dp in pressure corrds
         call HFReaduvtot(hfile,utot(:,:,k),vtot(:,:,k),idm,jdm,k,1)
         where(dp  >0.5*huge) dp=0.
         where(utot>0.5*huge) utot=0.
         where(vtot>0.5*huge) vtot=0.
         utot(:,:,k)=utot(:,:,k)-ubavg
         vtot(:,:,k)=vtot(:,:,k)-vbavg

         call rotate(utot(:,:,k),vtot(:,:,k),plat,plon,idm,jdm,'m2l')

      end do
      dp=dp/onem
      call HFReadField3D(hfile,saln,idm,jdm,kdm,'saln    ',1)
      call HFReadField3D(hfile,temp,idm,jdm,kdm,'temp    ',1)
      call HFReadField(hfile,fice,idm,jdm,'ficem   ',0,1)
      call HFReadField(hfile,hice,idm,jdm,'hicem   ',0,1)
      call HFReadField(hfile,uice,idm,jdm,'uice    ',0,1)
      call HFReadField(hfile,vice,idm,jdm,'vice    ',0,1)

      ! Call mosf routine (meridional overturning streamfunction
      print *,'Calling mosf'
      call mosf(utot,vtot,dp,idm,jdm,kdm,mostrf,lats,nlat,deep,ndeep,undef) 
      print *,'Calling mht'
      call mht(utot,vtot,saln,temp,dp,qlon,qlat,idm,jdm,kdm,merht,lats,nlat)

      ! Densities
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
         th3d(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
      end do
      end do
      end do
      print *,'Calling mosf_sig0'
      print *,minval(th3d),maxval(th3d)
      call mosf_sig0(utot,vtot,th3d,dp,idm,jdm,kdm,mostrf_sig0,lats,nlat,sigarray,nsig,undef)
      print *,minval(mostrf_sig0),maxval(mostrf_sig0)

      mostrf_ave= mostrf_ave+mostrf
      mostrf_sig0_ave= mostrf_sig0_ave+mostrf_sig0
      merht_ave = merht_ave+merht

      iarg = iarg + 1
    end do

    mostrf_ave=mostrf_ave/(iarg-1)
    mostrf_sig0_ave=mostrf_sig0_ave/(iarg-1)
    merht_ave =merht_ave /(iarg-1)

   !Create netcdf file and fill it with data
   ncfil='MOSF.nc'
   print '(a)','  Generating '//trim(ncfil)
   if (NF90_CREATE(trim(ncfil),NF90_CLOBBER,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(daily2mosf)'
   end if


   ! Define dimensions -- idm,jdm,kdm
   ! Following MERCATOR name definitions for netcdf var name
   !print *,'dims'
   call handle_err(NF90_DEF_DIM(ncid,'latitude',nlat  ,lat_dimension_id   ))
   call handle_err(NF90_DEF_DIM(ncid,'depth'   ,ndeep ,k_dimension_id     ))
   call handle_err(NF90_DEF_DIM(ncid,'sig0'    ,nsig  ,sig_dimension_id   ))
   vars_2d       =(/lat_dimension_id,k_dimension_id/)

   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title', 'MOSF'))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution',        &
      'NERSC, Edvard Griegs Vei 3A, Bergen, Norway'))
   call date_and_time(date=dateinfo)
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Dataset_generation_date', &
      dateinfo(1:8)))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'history', &
      'Created '//dateinfo(1:8)//' by program pak2mosf '))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'references', &
      'http://topaz.nersc.no'))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'field_type', &
      trim(ftype)))
   call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions', &
      'COARDS'))
   !call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'field_date', &
   !   rtd%cyy//'-'//rtd%cmm//'-'//rtd%cdm))

   
   ! Latitude
   call handle_err(NF90_DEF_VAR(ncid,'latitude',NF90_REAL,lat_dimension_id,lat_id))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'long_name','latitude'))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'units','degrees_north'))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'valid_range',(/-90.0, 90.0/)))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'standard_name','latitude'))
   call handle_err(NF90_PUT_ATT(ncid,lat_id,'axis','X'))

   ! Depth levels
   call handle_err(NF90_DEF_VAR(ncid,'depth',NF90_REAL,k_dimension_id,z_id))
   call handle_err(NF90_PUT_ATT(ncid,z_id,'long_name','depth'))
   call handle_err(NF90_PUT_ATT(ncid,z_id,'units','m'))
   call handle_err(NF90_PUT_ATT(ncid,z_id,'valid_range',(/0., maxval(deep)+1/)))
   call handle_err(NF90_PUT_ATT(ncid,z_id,'standard_name','depth'))
   call handle_err(NF90_PUT_ATT(ncid,z_id,'positive','down'))
   call handle_err(NF90_PUT_ATT(ncid,z_id,'axis','Y'))


   ! Density levels
   call handle_err(NF90_DEF_VAR(ncid,'sig0',NF90_REAL,sig_dimension_id,sig_id))
   call handle_err(NF90_PUT_ATT(ncid,z_id,'axis','Y'))
   
   ! MO Streamfunction - Z
   call handle_err(NF90_DEF_VAR(ncid,'mosf',NF90_REAL,vars_2d,strf_id))
   call handle_err(NF90_PUT_ATT(ncid,strf_id,'long_name', &
      'Meridional overturning Streamfunction z levels'))
   call handle_err(NF90_PUT_ATT(ncid,strf_id,'units','million meter3 second-1'))
   call handle_err(NF90_PUT_ATT(ncid,strf_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,strf_id,'standard_name', &
                   'ocean_meridional_overturning_streamfunction'))
   call handle_err(NF90_PUT_ATT(ncid,strf_id,'_FillValue',real(undef,kind=4)))
   
   ! MO Streamfunction - sig0
   call handle_err(NF90_DEF_VAR(ncid,'mosf_sig0',NF90_REAL,(/lat_dimension_id,sig_dimension_id/),sig_strf_id))
   call handle_err(NF90_PUT_ATT(ncid,sig_strf_id,'long_name', &
      'Meridional overturning Streamfunction sig0'))
   call handle_err(NF90_PUT_ATT(ncid,sig_strf_id,'units','million meter3 second-1'))
   call handle_err(NF90_PUT_ATT(ncid,sig_strf_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,sig_strf_id,'standard_name', &
                   'ocean_meridional_overturning_streamfunction'))
   call handle_err(NF90_PUT_ATT(ncid,sig_strf_id,'_FillValue',real(undef,kind=4)))

   ! Meridional heat transport
   call handle_err(NF90_DEF_VAR(ncid,'merht',NF90_REAL,lat_dimension_id,merht_id))
   call handle_err(NF90_PUT_ATT(ncid,merht_id,'long_name', &
      'Meridional Heat Transport'))
   call handle_err(NF90_PUT_ATT(ncid,merht_id,'units','W'))
   call handle_err(NF90_PUT_ATT(ncid,merht_id,'missing_value',undef))
   call handle_err(NF90_PUT_ATT(ncid,merht_id,'standard_name', &
                   'northward_ocean_heat_transport'))
   !call handle_err(NF90_PUT_ATT(ncid,merht_id,'_FillValue',undef)) 
   call handle_err(NF90_ENDDEF(ncid))


   call handle_err(NF90_PUT_VAR(ncid,lat_id,lats))
   call handle_err(NF90_PUT_VAR(ncid,z_id,deep))
   call handle_err(NF90_PUT_VAR(ncid,sig_id,sigarray))
   call handle_err(NF90_PUT_VAR(ncid,merht_id,merht_ave))
   call handle_err(NF90_PUT_VAR(ncid,strf_id,mostrf_ave))
   call handle_err(NF90_PUT_VAR(ncid,sig_strf_id,mostrf_sig0_ave))
   call handle_err(NF90_CLOSE(ncid))


end program p_mosf
