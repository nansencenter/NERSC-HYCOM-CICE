program p_icedrift_osisaf
! Calculates ice drift over a given number of days (2 days currently)
!
! metno-version: This is a stripped down and cleaned up version of p_icedrift2.F90
!
! Uses data from EUMETSAT OSI SAF http://osisaf.met.no
! (instead of IFREMER data)
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
! From Topaz4/Hycom the program reads the model grid and the daily mean files. From the
! daily mean files uice and vice is read. Using this the program calculates the model drift
! of the observation grid points in the model space. Finally, these drift vectors are
! transformed to dX,dY displacements in the observation space.
!
! dX,dY are written to an output file with the same name as the input file, but with
! extension uf.
!
! These displacements are then compared to the observation displacements in the EnKF-code.
!
! 2012-10-26 Geir Arne Waagbø (met.no)
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
  implicit none


  integer,parameter :: NDAYS_OSISAF_DRIFT=2 ! Always two day interval for OSISAF product
  character(len=*),parameter :: OSISAF_FILE_NAME_PREFIX = 'ice_drift_nh_polstere-625_multi-oi_'


  ! Global variables derived from input arguments
  character(len=3)  :: rungen       ! Typically TP4
  character(len=80) :: driftfile    ! Name of drift file from OSISAF
  character(len=4)  :: byear        ! Bulletin year (first date in daily mean file name)
  character(len=3)  :: bjuld        ! Bulletin julian day
  integer           :: enssize      ! Ensemble size

  call main()

  contains

  subroutine main()
      integer :: fyear,fjulday
      integer :: drnx,drny ! Dimensions of observation grid
      real,dimension(:,:),allocatable :: drlon,drlat ! lon/lat-grid in observation space
      real,dimension(:,:),allocatable :: x0,y0,x,y
      real,dimension(:,:),allocatable :: dX,dY
      integer :: imem

      call readArguments(fyear,fjulday)
      call getModelGrid()
      call getObservationGrid(drnx,drny,drlon,drlat)
      call getInitialPositions(drnx,drny,drlon,drlat,x0,y0)
      do imem=1,enssize
          call calculateModelDrift(drnx,drny,fyear,fjulday,x0,y0,x,y,imem)
          call dXdYObservationSpace(drnx,drny,x0,y0,x,y,dX,dY)
          call writeOutputFile(drnx,drny,dX,dY,imem)
      end do
  end subroutine

  subroutine readArguments(fyear,fjulday)
      integer,intent(out):: fyear,fjulday
      character(len=80) :: tmparg
      integer :: pl,fmonth,fdm
#if defined(IARGC)
      integer*4, external :: iargc
#endif

      if (iargc()==5) then
          ! Rungen for model id (TP4)
          call getarg(1,tmparg) ; read(tmparg,*) rungen

          ! bulletin refyear from model
          call getarg(2,tmparg) ; read(tmparg,*) byear

          ! bulletin julian day from model
          call getarg(3,tmparg) ; read(tmparg,*) bjuld

          ! Input ice drift file. Needed for initial grid positions,
          ! start and end times of model integrations
          call getarg(4,tmparg) ; read(tmparg,*) driftfile

          pl = len(OSISAF_FILE_NAME_PREFIX)

          if (len(driftfile) .lt. pl .or. driftfile(1:pl) .ne. OSISAF_FILE_NAME_PREFIX) then
              print *, 'icedrift_osisaf: ERROR'
              print *, 'Given filename does not start with correct prefix:'
              print *, OSISAF_FILE_NAME_PREFIX
              call exit(1)
          end if

          ! Get first and last time to integrate over
          read(driftfile(pl+1:pl+4),'(i4)') fyear
          read(driftfile(pl+5:pl+6),'(i2)') fmonth
          read(driftfile(pl+7:pl+8),'(i2)') fdm

          !print *,'First year, month, day - last year, month, day'
          !print *,fyear,fmonth,fdm

          ! Convert to julian days ref start year
          fjulday=datetojulian(fyear,fmonth,fdm,fyear,1,1)

          call getarg(5,tmparg) ; read(tmparg,*) enssize

       else
          print *,'icedrift_osisaf calculates ice drift based on a specified input file'
          print *
          print *,'Usage - 5 args:'
          print *,'icedrift_osisaf rungen bull_year bull_day osisaf_file enssize'
          print *
          call exit(1)
      end if
  end subroutine

  subroutine getObservationGrid(drnx,drny,drlon,drlat)
      integer,intent(out) :: drnx,drny
      real,dimension(:,:),allocatable,intent(inout) :: drlon,drlat
      integer, dimension(NF90_MAX_VAR_DIMS) :: dimsizes
      integer :: recdim, numdim
      real,dimension(:,:),allocatable :: dt0,dt1

      ! Get dimensions of drift file
      call ncvar_dims(driftfile,'dX',dimsizes,numdim,recdim)
      drnx=dimsizes(1)
      drny=dimsizes(2)

      ! Read observation grid from drift file
      allocate(drlon(drnx,drny))
      allocate(drlat(drnx,drny))
      call ncvar_read(driftfile,'lon',drlon,drnx,drny,1,1,1)
      call ncvar_read(driftfile,'lat',drlat,drnx,drny,1,1,1)
      !print *,minval(drlon),maxval(drlon)
      !print *,minval(drlat),maxval(drlat)

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
            call readDailyMeanFile(forecastdate,u,v,1,imem)
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
          call readDailyMeanFile(forecastdate,u,v,2,imem)

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
      !print *,'Displacement max in km:'
      !print *,maxval(sqrt( ((dX))**2  + ((dY))**2 ))
  end subroutine

  subroutine readDailyMeanFile(forecastdate,u,v,ind,imem)
      type(year_info),intent(in) :: forecastdate
      real,dimension(:,:,:),intent(inout) :: u,v
      integer,intent(in) :: ind,imem
      real*4,dimension(:,:),allocatable :: iou,iov
      character(len=80) :: hycfile
      integer :: irecl

      allocate(iou(idm,jdm))
      allocate(iov(idm,jdm))

      hycfile=rungen//'DAILY_'//byear//'_'//bjuld// &
           '_'//forecastdate%cyy//'_'//forecastdate%cdd// &
           '_ICEDRIFT.uf'
      !print *,'reading '//trim(hycfile)

      inquire(iolength=irecl)iou,iov
      open(10,file=trim(hycfile),access='direct',status='old',recl=irecl)
      read(10,rec=imem) iou,iov
      close(10)
      u(:,:,ind) = iou
      v(:,:,ind) = iov
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

      !print *, 'Wrote: ',drfilemodel
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
