! --- -------------------------------------------------------------------
! --- River routine trip_riverflow
! --- -------------------------------------------------------------------
! --- Program to "connect the dots" in the trip05 river path database.
! --- For each TRIP point, data from a runoff database is assigned, and 
! --- a simple river transport model is used to calculate the flow from 
! --- the point towards the river outlets. All calculation is done on the 
! --- TRIP grid.
! ---
! --- For now this routine uses ERA5 data, but it can easily be changed 
! --- to other runoff products.
! ---
! --- Output from this routine is:
! --- netcdf file with runoff and river volume per grid cell.
! --- TODO: hycom river file
! --- -------------------------------------------------------------------
! --- Prerequisites:
! --- 1) river_weights must be called before running this 
! ---    routine - to get the mapping from runoff grid -> TRIP grid
! --- 2) ERA5 runoff data must be availabel in either the directory "./ERA5/",
! ---    or in the directory set in env variable ERA5_PATH
! --- 2) TRIP data base must be availabel in either the directory "./Data/",
! ---    or in the directory set in env variable TRIP_PATH
! --- -------------------------------------------------------------------

program trip_flow
   use mod_year_info
   use netcdf
   use m_handle_err
   use m_era40_fix
   use m_read_runoff_era40, only : nrolon_era40=>nlon, nrolat_era40=>nlat, &
                                   rolat_era40 => lat, rolon_era40 => lon, &
                                   init_runoff_era40, read_runoff_era40
   use m_read_runoff_erai, only : nrolon_erai=>nlon, nrolat_erai=>nlat, &
                                   rolat_erai => lat, rolon_erai => lon, &
                                   init_runoff_erai, read_runoff_erai
   use m_read_runoff_era5, only : nrolon_era5=>nlon, nrolat_era5=>nlat, &
                                   rolat_era5 => lat, rolon_era5 => lon, &
                                   init_runoff_era5, read_runoff_era5

   use mod_trip
   implicit none
!#if defined (TRIP05)
!   integer, parameter :: nx=720, ny=360 ! grid dimensions 0.5 X 0.5 grid
!   real,    parameter :: dx=0.5, dy=0.5 ! grid spacing
!   character(len=*), parameter :: tfile='trip05.txt'
!   character(len=720) txtline
!#elif defined (TRIP)
!   integer, parameter :: nx=360, ny=180 ! grid dimensions 1 X 1 grid
!   real,    parameter :: dx=1., dy=1. ! grid spacing
!   character(len=*), parameter :: tfile='trip.txt'
!   character(len=360) txtline
!#endif
!   character(len=200) :: pathtrip, cenv


   real, parameter :: rearth=6372.795477598 ! Quadratic mean radius (km)
   real, parameter :: radian=57.2957795
   real, parameter :: ueff  = 0.35

   character(len=4) :: cyy,cstartyy

   integer, dimension(nx,ny) :: destinationx, destinationy,lmask
!   real, dimension(nx) :: lon
!   real, dimension(ny) :: lat

   integer :: ios
   integer :: i,j,i2,j2,i3,j3, niter,k, itl,im1,ip1

   real    :: sumcatch(nx*ny)
   real    :: newvol(nx,ny),oldvol(nx,ny),riverflux(nx,ny), &
              riverflux_clim(nx,ny,12), &
              ro_clim       (nx,ny,12), &
              vol_clim      (nx,ny,12)
   integer :: numclim(12)
   real    :: triprunoff(nx,ny)
   real    :: tmplline(nx), tmp
   logical :: rmask(nx,ny)
   integer :: ocnmask(nx,ny)
   integer :: nriver

   integer              :: maxrocells,nxro,nyro
   real   , allocatable :: roweights(:,:,:)
   integer, allocatable :: nrocells(:,:), romapi(:,:,:), romapj(:,:,:)
   real   , allocatable :: ro(:,:)

   integer :: itrip,jtrip,destj,desti,nc_lastyear
   real    :: tmparea
   real    :: dt, dist, F

   integer :: ncid, varid, xdim, ydim, recdim, dims3D(3), irec, varidro, &
              varidriv, dims2D(2), varidrivcatch, rivdim, varidrec
   logical :: ex, asc_info=.false., excatch,firstclim=.true.,first_after_spinup=.true.
   real :: rtime, spinupdays,intdays
   integer :: startyear

    integer, parameter :: nriver_catch=30, maxtimes=5000
    real, allocatable, dimension(:,:) :: &
       river_outlets_flux
    real, allocatable, dimension(:) :: &
       river_outlets_lon,river_outlets_lat,river_outlets_catch
    integer, allocatable, dimension(:) :: &
       river_outlets_i,river_outlets_j

   integer :: iyear, imonth, idom,num_year

   real, external :: spherdist
   logical, parameter  :: ncdump=.true. ! put to netcdf file for illustrations
   integer  :: moddump

   integer  :: nrolon, nrolat
   real, dimension(:), allocatable :: rolon, rolat
   character(len=80) :: runoff_source,myfile
#if defined(IARGC)
   integer*4, external :: iargc
#endif
   if (iargc()>=1) then
      call getarg(1,runoff_source)
   else 
      runoff_source="era5" !default
   end if
   print *,"Runoff source: "//trim(runoff_source)


   ! Set up erai path and lon/lat
   if (trim(runoff_source) == "era40") then 
      call init_runoff_era40()
      nrolon = nrolon_era40
      nrolat = nrolat_era40
      allocate(rolon(nrolon))
      allocate(rolat(nrolat))
      rolon  = rolon_era40
      rolat  = rolat_era40
   elseif (trim(runoff_source) == "erai") then 
      call init_runoff_erai()
      nrolon = nrolon_erai
      nrolat = nrolat_erai
      allocate(rolon(nrolon))
      allocate(rolat(nrolat))
      rolon  = rolon_erai
      rolat  = rolat_erai
   elseif (trim(runoff_source) == "era5") then 
      call init_runoff_era5()
      nrolon = nrolon_era5
      nrolat = nrolat_era5
      allocate(rolon(nrolon))
      allocate(rolat(nrolat))
      rolon  = rolon_era5
      rolat  = rolat_era5
   else 

      print *,"Unknown runoff source "//trim(runoff_source)
      call exit(1)
   end if

   ! Initialize the trip dataset
   call init_trip()
      

!   ! "Flip" the data so that increasing j is northwards
!   do j=1,ny
!      tmp=lat(j) 
!      lat(j)=lat(ny-j+1)
!      lat(ny-j+1)=tmp
!
!      tmplline=direction(:,j)
!      direction(:,j)=direction(:,ny-j+1)
!      direction(:,ny-j+1)=tmplline
!   end do


   ! Read the cell info schtuff
   !open(10,file='rw_maxncells.asc',status='old')
   !read(10,*) maxrocells
   !close(10)

   inquire(exist=ex,file='rw_cellinfo.uf')
   if (.not. ex) then 
      print *,'Could not find rw_cellinfo.uf - run trip_weights first'
      stop
   end if


   open(10,file='rw_cellinfo.uf'   ,form='unformatted',status='old')
   read(10)  maxrocells,nxro,nyro
   close(10)
   print *,'Max number of cells, ronx, rony: ',maxrocells,nxro,nyro
   allocate(nrocells (nxro,nyro))
   allocate(romapi   (nxro,nyro,maxrocells))
   allocate(romapj   (nxro,nyro,maxrocells))
   allocate(roweights(nxro,nyro,maxrocells))

   ! Read the weights, ncells and mapping from a binary file
   open(10,file='rw_cellinfo.uf'   ,form='unformatted',status='old')
   read(10)  maxrocells,nxro,nyro       , &
             nrocells                   , &
             romapi   (:,:,1:maxrocells), &
             romapj   (:,:,1:maxrocells), &
             roweights(:,:,1:maxrocells)
    close(10)





!    ! Trip grid cell areas (approximate)
!    do j=1,ny
!    do i=1,nx
!       triparea(i,j)=sin(dx/radian)*sin(dy/radian)*sin(lat(j)/radian)*rearth**2*1e6
!    end do
!    end do
    allocate(ro   (nxro,nyro))



    ! Initialize potential river outlets
    nriver=0
    rmask=.false.
    do j=2,ny-1
    do i=1,nx
    if (direction(i,j)==9.or.direction(i,j)==0) then
       ip1=mod(i,nx)+1
       im1=mod(nx+i-2,nx)+1

       if ( (any(direction(im1,j-1:j+1)>0.and.direction(im1,j-1:j+1)<9)) .or.&
            (any(direction(i  ,j-1:j+1)>0.and.direction(i  ,j-1:j+1)<9)) .or.&
            (any(direction(ip1,j-1:j+1)>0.and.direction(ip1,j-1:j+1)<9)) ) then

          nriver=nriver+1
          rmask(i,j)=.true.
       end if
    end if
    end do
    end do
    print *,'Potential number of rivers:',nriver



    ! Read top 30 river outlets based on catchment data
    !Top rivers by catchment
    inquire(exist=excatch,file='catchment.asc') 
    if (excatch) then
       allocate(river_outlets_i    (nriver_catch))
       allocate(river_outlets_j    (nriver_catch))
       allocate(river_outlets_lon  (nriver_catch))
       allocate(river_outlets_lat  (nriver_catch))
       allocate(river_outlets_catch(nriver_catch))
       allocate(river_outlets_flux (maxtimes,nriver_catch))
       open(10,file='catchment.asc',status='old')
       do i=1,nriver_catch
          !read(10,'(2i5,2f14.5,e14.4)') basinx(i),basiny(i), &
          !   lon(basinx(i)),lat(basiny(i)), sumcatch(i)
          read(10,*) river_outlets_i(i),river_outlets_j(i), &
             river_outlets_lon(i),river_outlets_lat(i),     &
             river_outlets_catch(i)

          !print *,river_outlets_i(i),river_outlets_j(i),river_outlets_lat(i), &
          !        lat(river_outlets_j(i))

          ! Check for "flip" and river outlet error
          if (river_outlets_lat(i) /=  lat(river_outlets_j(i))) then
             print *,river_outlets_i(i),river_outlets_j(i),river_outlets_lat(i), &
                     lat(river_outlets_j(i))
             print *,'"Flip error" on lat variable'
             stop
          else if (.not.rmask(river_outlets_i(i),river_outlets_j(i)) ) then
             !This shouldnt happen either, so its a nice safety check
             print *,'river outlet apparently not in river mask?'
             stop
          end if
       end do
       close(10)
    end if


    !stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main loop starts here


    if (trim(runoff_source) == "era40") then 
       spinupdays=3*365  ! 1 years
       startyear=1958
       num_year=40
       intdays  =num_year*365  ! 40 years
       dt=6*3600                  ! Time step (6 hours)
    elseif (trim(runoff_source) == "erai") then 
       spinupdays=3*365  ! 1 years
       num_year=27
       intdays  =num_year*365  ! Up to and including 2015
       startyear=1989
       dt=6*3600                  ! Time step (6 hours)
    elseif (trim(runoff_source) == "era5") then 
       spinupdays=3*365           ! 3 years of spinup
       num_year=5
       intdays  =num_year*365
       startyear=1992
       dt=6*3600                  ! Time step (6 hours)
    else 
       print *,"Unknown runoff source "//trim(runoff_source)
       call exit(1)
    end if

    print *,'River integration starts '
    print '(a,f4.1)','Integration time in years: ',intdays/(365)
    print '(a,f4.1)','Spinup time before climatology starts to be estimated: ',spinupdays/(365)
    oldvol=0. 
    newvol=0. 
    riverflux=0.
    moddump = nint(86400./dt)  ! Dump every day 
    riverflux_clim=0.
    vol_clim=0.
    ro_clim =0.
    numclim=0
    nc_lastyear=-1
    write(cstartyy,'(i4.4)') startyear

    ! Time integration loop
    do itl = 1,floor(intdays*4)

       ! convert itl*dt to time (in days)
       rtime=itl*dt/86400.

       call juliantodate(floor(rtime),iyear,imonth, idom, startyear,1,1)
       if (mod(itl,100)==0) then
          print '(a,i5,i3,i3)','Date ',iyear, imonth, idom
       end if
       !print *,itl,rtime, iyear, imonth,floor(intdays)/365.


       ! Read one record of the runoff field
        if (trim(runoff_source) == "era40") then 
           call read_runoff_era40(startyear,rtime,ro,nxro,nyro)
           !Fanf: Introduce a fix to the precipitation bias
           call era40_fix('RO',ro,nxro,nyro,rolon(1),rolat(1),rolon(2)-rolon(1),rolat(2)-rolat(1))
        elseif (trim(runoff_source) == "erai") then 
           call read_runoff_erai(startyear,rtime,ro,nxro,nyro)
        elseif (trim(runoff_source) == "era5") then 
           call read_runoff_era5(startyear,rtime,ro,nxro,nyro)
        else 
           print *,"Unknown runoff source "//trim(runoff_source)
           call exit(1)
        end if
       !   print *,' value out from era',nxro,nyro,startyear,rtime,ro(242,29)
       !print *,"min runoff ",minval(ro)



       ! Assign runoff field to trip cells using mapping from trip_weights
       triprunoff=0.
       do j=1,nyro
       do i=1,nxro
          do k=1,nrocells(i,j)
             itrip   =romapi   (i,j,k)
             jtrip   =romapj   (i,j,k)
             tmparea =roweights(i,j,k)
             triprunoff(itrip,jtrip)= triprunoff(itrip,jtrip) +  &
             max(ro(i,j),0.0)*tmparea*1e6 ! triparea in km^2 units: m^3/s
                !if (tmparea<0 .or. ro(i,j)<0) then
                if (tmparea<0) then
                   print *,'This should not happen corect me',tmparea,ro(i,j),i,j
                   stop
                endif

          end do
       end do
       end do

       if (asc_info) then
          print *,'Max trip area           (m^2        ):',maxval(triparea)
          print *,'Max trip runoff         (m^3/s      ):',maxval(triprunoff)
          print *,'Max trip runoff per area(m^3/(s*m^2)):',maxval(triprunoff/triparea)
       end if


       ! Mass budget of each grid cell
       oldvol=newvol
       do j=2,ny-1
       do i=1,nx
       if (direction(i,j)/=9.and.direction(i,j)/=0) then ! land point with throughflow


          ! lat direction of flow
          destj=j
          if     (direction(i,j)==1.or.direction(i,j)==2.or.  direction(i,j)==8) then
             destj=j-1
          elseif (direction(i,j)==4.or.direction(i,j)==5.or.  direction(i,j)==6) then
             destj=j+1
          end if

          ! lon direction of flow
          desti=i
          if     (direction(i,j)==2.or.direction(i,j)==3.or.  direction(i,j)==4) then
             desti=mod(i,nx)+1
          elseif (direction(i,j)==6.or.direction(i,j)==7.or.  direction(i,j)==8) then
             desti=mod(nx+i-2,nx)+1
          end if

          ! Distance between grid cell centers
          !print *,i,j,desti,destj
          dist=spherdist(lon(i),lat(j),lon(desti),lat(destj))

          ! Flow from one cell to the downstream neighbour
          ! (Miller et al, J. Clim, 1994, "Continental-Scale River Flow in Climate models")
          F=oldvol(i,j)*ueff/dist

          ! Flow "limiter" - avoids negative volume but slows river flow. Only a
          ! problem in the Antarctic with the current value of "ueff". Could
          ! also reduce time step dt
          !F=min(oldvol(i,j)/dt,F)
          F=max(0.,min(newvol(i,j)/dt,F*9./10.))

          ! New mass in this cell
          newvol(i,j)=newvol(i,j)+(triprunoff(i,j)-F)*dt

          !if (triprunoff(i,j)>0.) print *,F,dt,triprunoff(i,j)
          ! OK then, this seems to happen from time to time even though the
          ! above should prevent it. Correct F again and reset newvol
          if (newvol(i,j)<0.) then
             print '(a,2i5,2f8.2,6e12.4)','Warn: Negative river volume in ', &
                i,j,lon(i),lat(j),F*dt,dist,oldvol(i,j),ueff*dt/dist,newvol(i,j),triprunoff(i,j)
             F= F + newvol(i,j)/dt
             newvol(i,j)=0.
          end if

          ! New mass in downstream cell (no runoff added here, see above)
          newvol(desti,destj)=newvol(desti,destj)+F*dt

       end if
       end do
       end do


       ! Calculate flux in river outlets
       do j=1,ny
       do i=1,nx
       if (direction(i,j)==9.or.direction(i,j)==0) then ! "end point"
          riverflux(i,j)=(newvol(i,j)-oldvol(i,j))/dt
          newvol(i,j)=0. ! NB 
       end if
       end do
       end do

       !! Error checks
       !do j=1,ny
       !do i=1,nx
       !   if (newvol(i,j)<0.) then
       !      print *,'Negative river volume in ',i,j,lon(i),lat(j)
       !   end if
       !end do
       !end do
       where (newvol<0.) newvol=0. ! Should not be necessary, but kept in just in case

       ! Climatology
       if (rtime>spinupdays) then
          if (firstclim) print *,'**Spinup time over - climatology calculation starts'
          do j=1,ny
          do i=1,nx
             if (direction(i,j)==9.or.direction(i,j)==0) then ! land point with throughflow
                riverflux_clim(i,j,imonth)=riverflux_clim(i,j,imonth) + &
                                            riverflux(i,j)
             end if
             vol_clim(i,j,imonth)= vol_clim(i,j,imonth)+ &
                newvol(i,j)
             ro_clim (i,j,imonth)= ro_clim (i,j,imonth)+ &
                triprunoff(i,j)
          end do
          end do
          numclim(imonth)=numclim(imonth)+1
          !print *,minval(newvol),maxval(newvol)
          firstclim=.false.
       end if


       if (asc_info) then
          print *,'Max vol diff            (m^3/s      ):',maxval(newvol-oldvol)/dt
          print *,'Max River volume        (m^3        ):',maxval(newvol)
          print *,'Max River flux          (m^3/s      ):',maxval(riverflux)
          print *,'--------------------- END OF LOOP ----------------'
       end if


       ! Put in a netcdf file at regular intervals - 2D maps
       if (rtime>spinupdays) then
       if (mod(itl,moddump)==0 .and. ncdump) then

          write(cyy,'(i4.4)') iyear
          if (trim(runoff_source) == "era40") then 
             myfile = 'trip_era40_'//cyy//'.nc'
          elseif (trim(runoff_source) == "erai") then 
             myfile = 'trip_erai_'//cyy//'.nc'
          elseif (trim(runoff_source) == "era5") then 
             myfile = 'trip_era5_'//cyy//'.nc'
          else 
             print *,"Unknown runoff source "//trim(runoff_source)
             call exit(1)
          end if


          !Replace file on first try
          ex=(itl==moddump) 

          ! Also replace on year changeover
          ex=ex .or. nc_lastyear /= iyear

          ! And replace if first time after spinup
          ex = ex .or. first_after_spinup
          first_after_spinup=.false.

          if (ex) then
             call handle_err(nf90_create(myfile,NF90_CLOBBER,ncid))
             call handle_err(nf90_def_dim(ncid,'lon'    ,nx,xdim))
             call handle_err(nf90_def_dim(ncid,'lat'    ,ny,ydim))
             call handle_err(nf90_def_dim(ncid,'time',nf90_unlimited,recdim))
             if (excatch) call handle_err(nf90_def_dim(ncid,'river_catchment',nriver_catch,rivdim))
             dims3d=(/xdim,ydim,recdim/)
             dims2d=(/rivdim,recdim/)

             call handle_err(NF90_DEF_VAR(ncid,'lon',NF90_Float,xdim,varid))
             call handle_err(NF90_ENDDEF(ncid))
             call handle_err(NF90_PUT_VAR(ncid,varid,lon))

             call handle_err(NF90_REDEF(ncid))
             call handle_err(NF90_DEF_VAR(ncid,'mask',NF90_INT,(/xdim,ydim/),varid))
             call handle_err(NF90_ENDDEF(ncid))
             call handle_err(NF90_PUT_VAR(ncid,varid,ocnmask))

             call handle_err(NF90_REDEF(ncid))
             call handle_err(NF90_DEF_VAR(ncid,'lat',NF90_Float,ydim,varid))
             call handle_err(NF90_ENDDEF(ncid))
             call handle_err(NF90_PUT_VAR(ncid,varid,lat))

             if (excatch) then
                call handle_err(NF90_REDEF(ncid))
                call handle_err(NF90_DEF_VAR(ncid,'riverbycatch_lon',NF90_Float,rivdim,varid))
                call handle_err(NF90_ENDDEF(ncid))
                call handle_err(NF90_PUT_VAR(ncid,varid,river_outlets_lon))

                call handle_err(NF90_REDEF(ncid))
                call handle_err(NF90_DEF_VAR(ncid,'riverbycatch_lat',NF90_Float,rivdim,varid))
                call handle_err(NF90_ENDDEF(ncid))
                call handle_err(NF90_PUT_VAR(ncid,varid,river_outlets_lat))
             end if

             call handle_err(NF90_REDEF(ncid))
             call handle_err(NF90_DEF_VAR(ncid,'runoff',NF90_Float,dims3D,varidro))
             call handle_err(NF90_PUT_ATT(ncid,varidro,'units','m3 s-1'))
             call handle_err(NF90_DEF_VAR(ncid,'volume',NF90_Float,dims3D,varid))
             call handle_err(NF90_PUT_ATT(ncid,varid,'units','m3'))
             call handle_err(NF90_DEF_VAR(ncid,'river',NF90_Float,dims3D,varidriv))
             call handle_err(NF90_PUT_ATT(ncid,varid,'units','m3 s-1'))
             if (excatch) then
                call handle_err(NF90_DEF_VAR(ncid,'riverbycatch',NF90_Float,dims2D,varidrivcatch))
                call handle_err(NF90_PUT_ATT(ncid,varidrivcatch,'units','m3 s-1'))
             end if
             call handle_err(NF90_DEF_VAR(ncid,'time',NF90_Float,recdim,varidrec))
             call handle_err(NF90_PUT_ATT(ncid,varidrec,'units', &
                'days since '//trim(cstartyy)//'-01-01 00:00:00'))
             call handle_err(NF90_ENDDEF(ncid))
             irec=1
          else

             call handle_err(nf90_open(myfile,NF90_WRITE,ncid))
             call handle_err(nf90_inq_dimid(ncid, 'time', recdim))
             if (excatch) then 
                call handle_err(Nf90_inq_varid(ncid, 'riverbycatch', varidrivcatch))
             end if
             call handle_err(nf90_Inquire_Dimension(ncid, recdim, len=irec))
             call handle_err(Nf90_inq_varid(ncid, 'runoff', varidro ))
             call handle_err(Nf90_inq_varid(ncid, 'volume', varid   ))
             call handle_err(Nf90_inq_varid(ncid, 'river' , varidriv))
             call handle_err(Nf90_inq_varid(ncid, 'time'  , varidrec))

             irec=irec+1
          end if
          call handle_err(NF90_PUT_VAR(ncid,varidrec,rtime, &
                                       start=(/irec/)))
          call handle_err(NF90_PUT_VAR(ncid,varid,newvol(1:nx,1:ny), &
                                       start=(/1,1,irec/)))
          call handle_err(NF90_PUT_VAR(ncid,varidro,triprunoff(1:nx,1:ny), &
                                       start=(/1,1,irec/)))
          call handle_err(NF90_PUT_VAR(ncid,varidriv,riverflux(1:nx,1:ny), &
                                       start=(/1,1,irec/)))
          if (excatch) then
          do i=1,nriver_catch
             river_outlets_flux (min(max(1,itl/moddump),maxtimes),i)= &
                riverflux(river_outlets_i(i),river_outlets_j(i))
             !if (i==1) then
             !   print *, river_outlets_flux (min(max(1,itl/moddump),maxtimes),i)
             !end if
             call handle_err(NF90_PUT_VAR(ncid,varidrivcatch,  &
             river_outlets_flux(min(max(1,itl/moddump),maxtimes),:), start=(/1,irec/)))
          end do
          end if


          call handle_err(nf90_close(ncid))
          nc_lastyear=iyear
          

       end if
       end if !if (rtime>spinupdays+ 366) then
    enddo ! main loop
    !call handle_err(nf90_close(ncid))






    ! Dump climatology if calculated
    print *
    print  '(a)','River integration finished. Sampled rivers in trip_era??_YYYY.nc'
    if (rtime>spinupdays+ 366) then
       do imonth=1,12
          riverflux_clim(:,:,imonth)=riverflux_clim(:,:,imonth)/numclim(imonth)
          ro_clim(:,:,imonth)=ro_clim(:,:,imonth)/numclim(imonth)
          vol_clim(:,:,imonth)=vol_clim(:,:,imonth)/numclim(imonth)
       enddo

       if (trim(runoff_source) == "era40") then 
          call handle_err(nf90_create('trip_era40_clim.nc',NF90_CLOBBER,ncid))
       elseif (trim(runoff_source) == "erai") then 
          call handle_err(nf90_create('trip_erai_clim.nc',NF90_CLOBBER,ncid))
       elseif (trim(runoff_source) == "era5") then
          call handle_err(nf90_create('trip_era5_clim.nc',NF90_CLOBBER,ncid))
       else 
          print *,"Unknown runoff source "//trim(runoff_source)
          call exit(1)
       end if
       call handle_err(nf90_def_dim(ncid,'lon'    ,nx,xdim))
       call handle_err(nf90_def_dim(ncid,'lat'    ,ny,ydim))
       call handle_err(nf90_def_dim(ncid,'record',12,recdim))
       dims3d=(/xdim,ydim,recdim/)
       call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,"climatology_first_year",startyear))
       call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,"climatology_last_year", startyear+num_year))

       call handle_err(NF90_DEF_VAR(ncid,'lon',NF90_Float,xdim,varid))
       call handle_err(NF90_ENDDEF(ncid))
       call handle_err(NF90_PUT_VAR(ncid,varid,lon))

       call handle_err(NF90_REDEF(ncid))
       call handle_err(NF90_DEF_VAR(ncid,'mask',NF90_INT,(/xdim,ydim/),varid))
       call handle_err(NF90_ENDDEF(ncid))
       call handle_err(NF90_PUT_VAR(ncid,varid,ocnmask))

       call handle_err(NF90_REDEF(ncid))
       call handle_err(NF90_DEF_VAR(ncid,'rtst',NF90_Float,(/xdim,ydim/),varid))
       call handle_err(NF90_ENDDEF(ncid))
       call handle_err(NF90_PUT_VAR(ncid,varid,ro_clim(:,:,1)))

       call handle_err(NF90_REDEF(ncid))
       call handle_err(NF90_DEF_VAR(ncid,'lat',NF90_Float,ydim,varid))
       call handle_err(NF90_ENDDEF(ncid))
       call handle_err(NF90_PUT_VAR(ncid,varid,lat))

       call handle_err(NF90_REDEF(ncid))
       call handle_err(NF90_DEF_VAR(ncid,'volume',NF90_Float,dims3D,varid))
       call handle_err(NF90_PUT_ATT(ncid,varid,'units','m3'))
       call handle_err(NF90_ENDDEF(ncid))
       call handle_err(NF90_PUT_VAR(ncid,varid ,vol_clim,start=(/1,1,1/)))

       call handle_err(NF90_REDEF(ncid))
       call handle_err(NF90_DEF_VAR(ncid,'runoff',NF90_Float,dims3D,varidro))
       call handle_err(NF90_PUT_ATT(ncid,varidro,'units','m3 s-1'))
       call handle_err(NF90_DEF_VAR(ncid,'river',NF90_Float,dims3D,varidriv))
       call handle_err(NF90_PUT_ATT(ncid,varidriv,'units','m3 s-1'))
       call handle_err(NF90_ENDDEF(ncid))

       do k=1,12 
          call handle_err(NF90_PUT_VAR(ncid,varidro,ro_clim(:,:,k),start=(/1,1,k/)))
          call handle_err(NF90_PUT_VAR(ncid,varidriv,riverflux_clim(:,:,k),start=(/1,1,k/)))
       end do
       call handle_err(NF90_CLOSE(ncid))
       print  '(a)','River climatology in trip_era??_clim.nc'
    end if





   end program
