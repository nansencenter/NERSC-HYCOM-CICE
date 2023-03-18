program p_icedrift_osisaf
! Calculates ice drift over a given number of days (2 days in previous)
!
! Uses data from EUMETSAT OSI SAF http://osisaf.met.no
!
! The point of this program is to calculate *model* ice drift in the
! *observation* space, so that the innovation needed for the analysis
! can be calculated.
!
! The program reads the observation nc file:
! ice_drift_nh_polstere-625_multi-oi_YYYYMMDD1200-YYYYMMDD1200.nc
! but only to find information about the observation grid and the from- and to-time.
! Otherwise the observations are processed in the EnKF-code
!
! From Topaz4/Hycom the program reads the model grid and the daily mean files. 
! From the daily mean (iceh) uice and vice is read. Using this the program calculates the model drift
! of the observation grid points in the model space. 
! Finally, these drift vectors are transformed to dX,dY displacements in the observation space.
!
! dX,dY are written to an output file with the same name as the input file, but with
! extension uf.
!
! These displacements are then compared to the observation displacements in the EnKF-code.
! Based on the previous version : 2012-10-26 Geir Arne Waagbø (met.no)
!
! Updated in February 2023 by Jiping Xie
!  1. Obseravtions replaced by 24h new product
!  2. Model daily u/vice reading from Netcdf files
!  3. The output from this programe supports two formats both binary and Netcdf
!  4. This program supports to prepare the ensemble drift with Netcdf file
!  5. To conviniently visialize the ice drift, the referred observation file wroten into a text file
!
  use netcdf, only: NF90_MAX_VAR_DIMS
  use mod_xc, only: xcspmd,idm,jdm
  use mod_za, only: zaiost
  use mod_grid, only: plon,plat,depths,scpx,scpy,get_grid
  use mod_confmap
  use mod_year_info
  use mod_hycomfile_io, only: undef
  use m_ncvar_dims
  use m_ncvar_read
  use m_rk2
  use m_handle_err
  use m_ncvar_dump
  implicit none

  !integer,parameter :: NDAYS_OSISAF_DRIFT=2    ! Two day interval for OSISAF product before
  !character(len=*),parameter :: OSISAF_FILE_NAME_PREFIX = 'ice_drift_nh_polstere-625_multi-oi_'
  integer,parameter :: NDAYS_OSISAF_DRIFT=1  ! Two day interval for OSISAF product before
  character(len=*),parameter :: OSISAF_FILE_NAME_PREFIX='ice_drift_nh_ease2-750_cdr-v1p0_24h-'
  character(len=*),parameter :: Def_Obsdir='/cluster/projects/nn2993k/TP4b0.12/idrft_osisaf/'

  ! Global variables derived from input arguments
  character(len=3)      :: rungen       ! Typically TP4/TP5
  character(len=80)    :: driftfile    ! Name of drift file from OSISAF
  character(len=80)    :: dailyfile    ! Name of daily nc-file from MODLE
  !character(len=4)      :: byear        ! Bulletin year (first date in daily mean file name)
  integer               :: enssize      ! Ensemble size

  call main()

  contains

  subroutine main()
      integer :: fyear,fjulday
      integer :: drnx,drny ! Dimensions of observation grid
      real,dimension(:,:),allocatable    :: drlon,drlat ! lon/lat-grid in observation space
      real,dimension(:,:),allocatable    :: x0,y0,x,y
      real,dimension(:,:),allocatable    :: dX,dY
      real,dimension(:,:),allocatable :: dflg  ! 30 -> vector was retrieved with nominal algorithm 

      integer :: imem

      rungen='Mod'
      call readArguments(fyear,fjulday)
      call getModelGrid()

      call getObservationGrid(drnx,drny,drlon,drlat,dflg)
      call getInitialPositions(drnx,drny,drlon,drlat,x0,y0)

      if (enssize>0) then
         do imem=1,enssize
            call calculateModelDrift(drnx,drny,fyear,fjulday,x0,y0,x,y,imem)
            call dXdYObservationSpace(drnx,drny,x0,y0,x,y,dX,dY)
            call writeOutputFile(drnx,drny,dX,dY,imem)
            call ncvar_dump2D(rungen,drnx,drny,drlon,drlat,dX,dY,dflg,imem)
         end do
      else
         call calculateModelDrift(drnx,drny,fyear,fjulday,x0,y0,x,y,0)
         call dXdYObservationSpace(drnx,drny,x0,y0,x,y,dX,dY)
         !call writeOutputFile(drnx,drny,dX,dY,imem)
         call ncvar_dump2D(rungen,drnx,drny,drlon,drlat,dX,dY,dflg,0)
      endif
      open(10,file='Obs_drift.txt',action='write')
         write(10,'(a,SP)') trim(driftfile)
      close(10)

  end subroutine

  subroutine readArguments(fyear,fjulday)
      integer,intent(out):: fyear,fjulday
      character(len=100) :: tmparg
      character(len=8) :: Sdate 
      integer :: pl,fmonth,fdm
      integer :: lyear,lmonth,ldm
      integer*4, external :: iargc
      logical             :: fexist

      if (iargc()>=2.and.iargc()<=3) then
         ! Rungen for model id (TP5)
         ! bulletin refyear from model
         ! call getarg(1,tmparg) ; read(tmparg,*) byear
         ! bulletin julian day from model
         !call getarg(2,tmparg) ; read(tmparg,*) bjuld

         ! Input ice drift file from model
         call getarg(1,tmparg) ; read(tmparg,*) dailyfile 

         ! Input ice drift file. Needed for initial grid positions,
         ! start and end times of model integrations
         call getarg(2,tmparg) ; !read(tmparg,*) driftfile
         
         pl=5
         if (len(dailyfile) .lt. pl .or. dailyfile(1:pl) .ne.'iceh.') then
            print *, 'icedriftnc_osisaf: ERROR'
            print *, 'Given filename does not follow with general prefix:'
            print *, 'iceh.????-??-??.nc'
            call exit(1)
         end if
         
         ! Get first 24h time to integrate over
         read(dailyfile(pl+1:pl+4),'(i4)') lyear
         read(dailyfile(pl+6:pl+7),'(i2)') lmonth
         read(dailyfile(pl+9:pl+10),'(i2)') ldm
         print *,'Last year, month, day'
         print *,lyear,lmonth,ldm

         ! Convert to julian days ref start year
         fjulday=datetojulian(lyear,lmonth,ldm,lyear,1,1)
         print *,'First year, month, day',fjulday
         if (fjulday<NDAYS_OSISAF_DRIFT) then
            fyear=lyear-1
            select case (NDAYS_OSISAF_DRIFT)
               case (2) 
                  fmonth=12
                  fdm=30
               case (3) 
                  fmonth=12
                  fdm=29
               case default
                  fmonth=12
                  fdm=31
            end select

         else
            call juliantodate(fjulday-NDAYS_OSISAF_DRIFT,fyear,fmonth,fdm,lyear,1,1)
         endif
         print *,fyear,fmonth,fdm,Sdate
         write(Sdate,'(i4.4,i2.2,i2.2)') fyear,fmonth,fdm
         write(dailyfile(pl+1:pl+4),'(i4)') fyear
         write(dailyfile(pl+6:pl+7),'(i2)') fmonth
         write(dailyfile(pl+9:pl+10),'(i2)') fdm
         fjulday=datetojulian(fyear,fmonth,fdm,fyear,1,1)

         driftfile=trim(OSISAF_FILE_NAME_PREFIX)//trim(Sdate)//'1200.nc'
         inquire(File=trim(tmparg)//trim(driftfile),EXIST=fexist)
         if (fexist) then
            !call system("ln -sf 'trim(tmparg)//trim(driftfile)//' .")
            call system("ln -sf "//trim(tmparg)//trim(driftfile)//" .")
            print *, 'Obs. Drift: ',trim(OSISAF_FILE_NAME_PREFIX)//trim(Sdate)//'1200.nc'
         else
            print *, 'Missing the drift Obs. Please check the input directory'
            print *,trim(tmparg)
            print *, 'Or try to use the default directory (Betzy): '
            print *, trim(Def_Obsdir)//' ?'
            call exit(1)
         endif

         enssize=0
         if (iargc()==3) then
            call getarg(3,tmparg) ; read(tmparg,*) enssize
            if (enssize<1.or.enssize>100) then
               print *, 'double check for the predefined index in ensemble!'
               call exit(1)
            endif
         endif 

       else
          print *,'icedriftnc_osisaf calculates ice drift based on a specified output file in HYCOM_CICE'
          print *
          print *,'Usage - 2/3 args:'
          print *,'icedriftnc_osisaf iceh_nc_file Obs_dir <size of ensemble>'
          print *
          call exit(1)
      end if
  end subroutine

  subroutine getObservationGrid(drnx,drny,drlon,drlat,drflg)
      integer,intent(out) :: drnx,drny
      real,dimension(:,:),allocatable,intent(inout) :: drlon,drlat
      real,dimension(:,:),allocatable,intent(inout) :: drflg
      integer, dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      integer :: recdim, numdim
      real,dimension(:,:),allocatable :: dt0,dt1

      ! Get dimensions of drift file
      print *,trim(driftfile)
      call ncvar_dims(driftfile,'dX',dimsizes,numdim,recdim)
      drnx=dimsizes(1)
      drny=dimsizes(2)

      ! get Quality flag: status_flag

      ! Read observation grid from drift file
      allocate(drlon(drnx,drny))
      allocate(drlat(drnx,drny))
      allocate(drflg(drnx,drny))
      call ncvar_read(driftfile,'lon',drlon,drnx,drny,1,1,1)
      call ncvar_read(driftfile,'lat',drlat,drnx,drny,1,1,1)
      call ncvar_read(driftfile,'status_flag',drflg,drnx,drny,1,1,1)
      !print *,minval(drlon),maxval(drlon)
      !print *,minval(drlat),maxval(drlat)

  end subroutine

!
! reading the iceh file as the daily output from HYCOM_CICE
  subroutine getIcedriftDaily(driu,driv)
      real,dimension(idm,jdm),intent(inout) :: driu,driv
      integer, dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      integer :: recdim, numdim
      integer :: drnx,drny
      ! skipping the grids outside of ice cover (aice_d<0.05)
      real,dimension(idm,jdm) :: varmask 

      ! Get dimensions of model file
      call ncvar_dims(dailyfile,'ULON',dimsizes,numdim,recdim)
      drnx=dimsizes(1)
      drny=dimsizes(2)

      call ncvar_read(dailyfile,'uvel_d',driu,drnx,drny,1,1,1)
      call ncvar_read(dailyfile,'vvel_d',driv,drnx,drny,1,1,1)
      call ncvar_read(dailyfile,'aice_d',varmask,drnx,drny,1,1,1)
      where(varmask<0.05.or.varmask>1.0) driv=0
      where(varmask<0.05.or.varmask>1.0) driu=0
      !print *,minval(driv),maxval(driv)
      !print *,minval(driu),maxval(driu)

  end subroutine



  subroutine getModelGrid()
      call xcspmd()
      call zaiost()
      call get_grid
      call initconfmap(idm,jdm)
  end subroutine

  subroutine getInitialPositions(drnx,drny,drlon,drlat,x0,y0)
      integer,intent(in) :: drnx,drny
      real,dimension(:,:),intent(in) :: drlon,drlat
      real,dimension(:,:),allocatable,intent(out) :: x0,y0
      integer :: i,j
      integer :: ip,jp ! Pivot points

      allocate(x0(drnx,drny))
      allocate(y0(drnx,drny))

      do j=1,drny
      do i=1,drnx
         call oldtonew(drlat(i,j),drlon(i,j),x0(i,j),y0(i,j))
         call pivotp(y0(i,j),x0(i,j),ip,jp)
         if (ip>0.and.ip<idm.and.jp>0 .and. jp<jdm-1) then
            call ll2gind(drlon(i,j),drlat(i,j),x0(i,j),y0(i,j))
         else
            x0(i,j)=undef
            y0(i,j)=undef
         end if
      end do
      end do

      !print *,'Initial Grid positions :'
      !print *,minval(x0,x0/=undef),maxval(x0,x0/=undef)
      !print *,minval(y0,y0/=undef),maxval(y0,y0/=undef)
  end subroutine

  subroutine calculateModelDrift(drnx,drny,fyear,fjulday,x0,y0,x,y,imem)
      integer,intent(in) :: drnx,drny,fyear,fjulday,imem
      real,dimension(:,:),intent(in) :: x0,y0
      real,dimension(:,:),allocatable,intent(out) :: x,y
      integer :: daystep,irecl
      logical :: first
      real :: time,delt
      type(year_info) :: forecastdate
      real,dimension(:,:,:),allocatable :: u,v

      allocate(u(idm,jdm,2))
      allocate(v(idm,jdm,2))

      allocate(x(drnx,drny))
      allocate(y(drnx,drny))

      ! We use 6 hours as the RK(2) time step
      delt=6*3600.

      time=fjulday
      first=.true.
      x=x0
      y=y0

      ! This cycles days
      do daystep=1,NDAYS_OSISAF_DRIFT
          !print *,'day step',daystep
          ! Initially ; Read new field into time slot 0 and 1
          if (first) then
            call year_day(time,fyear,forecastdate,'ecmwf')
            if (imem>0) then
               call readDailyMeanFile(forecastdate,u,v,1,imem)
            else
               call readDailyMeanFile(forecastdate,u,v,1)
            endif
            first=.false.
          else
            !print *,'switching'
            u(:,:,1)=u(:,:,2)
            v(:,:,1)=v(:,:,2)
          end if

          ! Note: year_day routine handles change from one year
          ! to the next. If time+daystep is 365 or larger it
          ! returns a day in January the following year
          call year_day(time+daystep,fyear,forecastdate,'ecmwf')
          if (imem>0) then
             call readDailyMeanFile(forecastdate,u,v,1,imem)
          else
             call readDailyMeanFile(forecastdate,u,v,2)
          endif
          ! One days worth of rk2 integration here
          call rk2(u,v,scpx,scpy,idm,jdm,x,y,drnx,drny,delt,undef)

      enddo

      !print *
      !print *,'Displacement max in grid coordinates:'
      !print *,'x:',maxval(abs(x-x0))
      !print *,'y:',maxval(abs(y-y0))
  end subroutine

  subroutine dXdYObservationSpace(drnx,drny,x0,y0,x,y,dX,dY)
      integer,intent(in) :: drnx,drny
      real,dimension(:,:),intent(in) :: x0,y0,x,y
      real,dimension(:,:),allocatable,intent(out) :: dX,dY
      integer :: i,j
      real :: lat0,lon0,lat1,lon1
      real :: xx0,yy0,xx1,yy1

      allocate(dX(drnx,drny))
      allocate(dY(drnx,drny))
      do j=1,drny
      do i=1,drnx
          ! From model space to lon/lat
          call gind2ll(x0(i,j),y0(i,j),lon0,lat0)
          call gind2ll(x (i,j),y (i,j),lon1,lat1)

          ! Put in same lon lat interval
          if (lon1 >  lon0+180.) lon1=lon1-360.
          if (lon1 <= lon0-180.) lon1=lon1+360.

          ! From lon/lat to observation space
          call calculateXY(lon0,lat0,xx0,yy0)
          call calculateXY(lon1,lat1,xx1,yy1)

          dX(i,j) = (xx1 - xx0)*1e-3 ! km
          dY(i,j) = (yy1 - yy0)*1e-3 ! km
      end do
      end do

      !print *
      print *,'Displacement max in km:'
      print *,maxval(sqrt( ((dX))**2  + ((dY))**2 ))
  end subroutine

  subroutine readDailyMeanFile(forecastdate,u,v,ind,indmem)
      type(year_info),intent(in) :: forecastdate
      real,dimension(:,:,:),intent(inout) :: u,v
      integer,intent(in) :: ind
      integer,optional,intent(in) :: indmem
      !real*4,dimension(:,:),allocatable :: iou,iov
      real,dimension(:,:),allocatable :: iou,iov
      character(len=80) :: hycfile
      integer :: irecl
      character(len=3)  :: cmem

      allocate(iou(idm,jdm))
      allocate(iov(idm,jdm))

      if (present(indmem)) then
         write(cmem,"(I3.3)") indmem
         dailyfile='iceh.'//forecastdate%cyy//'-'//forecastdate%cmm//  &
                '-'//forecastdate%cdm//"_"//trim(cmem)//'.nc' 
      else
         dailyfile='iceh.'//forecastdate%cyy// &
            '-'//forecastdate%cmm//'-'//forecastdate%cdm//'.nc' 
      endif

      print *,'reading '//trim(dailyfile),ind
      call getIcedriftDaily(iou,iov)
      u(:,:,ind) = iou
      v(:,:,ind) = iov
      deallocate(iou,iov)
  end subroutine

  subroutine writeOutputFile(drnx,drny,dX,dY,imem)
      integer,intent(in) :: drnx,drny
      real,dimension(:,:),intent(in) :: dX,dY
      integer,intent(in) :: imem
      integer :: irecl
      integer :: findnc
      character(len=80) :: drfilemodel
      real*4, dimension(:,:), allocatable :: dXio,dYio

      ! Name for resulting model drift file
      findnc=index(driftfile,'.nc')
      drfilemodel=driftfile(1:findnc-1)//'.uf'

      allocate(dXio(drnx,drny))
      allocate(dYio(drnx,drny))
      dXio=dX
      dYio=dY

      inquire(iolength=irecl)dXio,dYio
      open(10,file=trim(drfilemodel),access='direct',status='unknown',recl=irecl)
      write(10,rec=imem) dXio,dYio
      close(10)

   !   print *, 'Wrote: ',drfilemodel
  end subroutine

  subroutine calculateXY(lon,lat,x,y)
  ! Calculates x,y-coordinates in observation space
      real,intent(in)  :: lon,lat
      real,intent(out) :: x,y
      real,parameter :: TRUELAT1 = 70.
      real,parameter :: STDLON = -45.
      real,parameter :: E_DATUM  = .08181615300113177577
      real,parameter :: A_DATUM  = 6378273.
      real,parameter :: rad_per_deg = 4.0*atan(1.0)/180.

      ! Local variables
      real :: h, mc, tc, t, rho

      h = 1. ! h=-1 for South Polar Stereographic

      mc = cos(h*TRUELAT1*rad_per_deg)/sqrt(1.0-(E_DATUM*sin(h*TRUELAT1*rad_per_deg))**2.0)
      tc = sqrt(((1.0-sin(h*TRUELAT1*rad_per_deg))/(1.0+sin(h*TRUELAT1*rad_per_deg)))* &
                      (((1.0+E_DATUM*sin(h*TRUELAT1*rad_per_deg))/(1.0-E_DATUM*sin(h*TRUELAT1*rad_per_deg)))**E_DATUM ))

      t = sqrt(((1.0-sin(h*lat*rad_per_deg))/(1.0+sin(h*lat*rad_per_deg))) * &
                     (((1.0+E_DATUM*sin(h*lat*rad_per_deg))/(1.0-E_DATUM*sin(h*lat*rad_per_deg)))**E_DATUM))


      rho = A_DATUM * mc * t / tc
      x = h *  rho * sin((h*lon - h*STDLON)*rad_per_deg)
      y = h *(-rho)* cos((h*lon - h*STDLON)*rad_per_deg)
  end subroutine
end program
