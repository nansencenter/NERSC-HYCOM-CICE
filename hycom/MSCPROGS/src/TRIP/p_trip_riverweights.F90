! --- -------------------------------------------------------------------
! --- River routine trip_riverweights
! --- -------------------------------------------------------------------
! --- Program to map ERAI grid cells onto TRIP grid cells. Required
! --- by trip_riverflow.
! ---
! --- For now this routine uses ERAI data, but it can easily be changed 
! --- to other runoff products.
! ---
! --- Output from this routine is:
! --- unformatted file containing mapping from ERAI runoff grid -> TRIP grid
! --- -------------------------------------------------------------------
! --- Prerequisites:
! --- 1) ERAI landmask must be available in the path set in env variable ERAI_PATH
! --- 2) TRIP data base must be available in the path set in env variable TRIP_PATH
! --- -------------------------------------------------------------------

program trip_riverweights
   use netcdf
   use m_read_runoff_era40, only : nrolon_era40=>nlon, nrolat_era40=>nlat, &
                                   rolat_era40 => lat, rolon_era40 => lon, &
                                   init_runoff_era40
   use m_read_runoff_erai, only : nrolon_erai=>nlon, nrolat_erai=>nlat, &
                                   rolat_erai => lat, rolon_erai => lon, &
                                   init_runoff_erai
   use mod_trip
   use m_handle_err
   implicit none
   real, parameter :: rearth=6372.795477598 ! Quadratic mean radius (km)
   real, parameter :: radian=57.2957795
   integer :: ivar,i,j
   real, allocatable :: var(:)
   integer :: itrip,jtrip,iro,jro,im1,ip1
   real :: tripeast, tripwest, tripnorth, tripsouth
   real :: roeast, rowest, ronorth, rosouth
   real :: dlon_overlap,dlat_overlap,area_overlap,dro
   integer :: rominj,romaxj, ntmp, ierr, flip_lat
   integer, parameter :: itest=90,jtest=90
   integer, parameter :: maxrocells=20
   integer, allocatable, dimension(:,:,:) :: romapi,romapj ! trip cells map to these ro cells
   integer, allocatable, dimension(:,:) :: tmpi
   real   , allocatable, dimension(:,:,:) :: roweights     ! area of trip cells occ. ro cells
   real   , allocatable, dimension(:,:) :: tmp     ! area of trip cells occ. ro cells
   integer, allocatable, dimension(:,:)    :: nrocells      ! n ro cells per trip cell
   character(len=80) :: runoff_source

   integer  :: nrolon, nrolat
   real, dimension(:), allocatable :: rolon, rolat
#if defined(IARGC)
   integer*4, external :: iargc
#endif
   if (iargc()>=1) then
      call getarg(1,runoff_source)
   else 
      runoff_source="erai"
   end if

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
   else 
      print *,"Unknown runoff source "//trim(runoff_source)
      call exit(1)
   end if

   print *,"Runoff source: "//trim(runoff_source)



   ! Set up TRIP path and lon/lat
   call init_trip()


   ! Flip it when we are working with it here so that latitude is
   ! increasign with j
   call flip_trip()
   print *,'lat runof(1),tp(1),dx-ro,dx-tp',rolat(1),lat(1),rolat(2)-rolat(1),dx,rolat(nrolat),lat(ny)
   print *,'lon runof(1),tp(1),dx-ro,dx-tp',rolon(1),lon(1),rolon(2)-rolon(1),dx,rolon(nrolon),lon(nx)

   allocate(romapi   (nrolon,nrolat,maxrocells))
   allocate(romapj   (nrolon,nrolat,maxrocells))
   allocate(roweights(nrolon,nrolat,maxrocells))
   allocate(nrocells (nrolon,nrolat           ))
   flip_lat=0
   if (rolat(2)<rolat(1)) then
   !Fanf: quick fix of the function.
   !if lat is decreasing , flip them, and flip the weight at the end
   flip_lat=1
   dro=rolat(1)-rolat(2)
   do j= 1,nrolat
      rolat(j)=-90+(j-1)*dro
   enddo
   end if

   ! Go through each TRIP grid cell, and weight/assign it to the runoff grid
   ! cells
   romapi=-1
   romapj=-1
   nrocells=0
   roweights=0.
   print '(a)','Calculating how much each runoff cell covers each  TRIP cell'
   print '(a)','ru = runoff grid index, j1=low index, j2=high index '
   print '(a,3a6,a,3a12)','j',' jtrip',' ro_j1',' ro_j2','-->', &
                             'lat(jtrip)','rolat(ro_j1)','rolat(ro_j2)'
   print '(a)','-------------------------------------------------------------'
   do jtrip=2,ny-1

      ! Narrow the search range on runoff grid. NB - assumes runoff lat is 
      ! increasing
      rominj=1
      romaxj=nrolat
      !find the indices in runoff including [-2dx ; 2dx]
      do jro=1,nrolat
         if (rolat(jro)<lat(jtrip)-2*dx) rominj=jro
         if (rolat(nrolat-jro+1)>lat(jtrip)+2*dx) romaxj=nrolat-jro+1
      end do
      rominj=max(2,rominj)
      romaxj=min(nrolat-1,romaxj)

      if (mod(jtrip,20)==0) &
      print '(a,3i6,a,3f12.2)','j',jtrip,rominj,romaxj,'-->', &
                             lat(jtrip),rolat(rominj),rolat(romaxj)

      do itrip=1,nx
         ! west, east, south and north boundary of trip cell
         ip1=mod(itrip,nx)+1
         im1=mod(itrip+nx-2,nx)+1
         tripwest =(lon(itrip)+newloncenter(lon(itrip),lon(im1)))/2.
         tripeast =(newloncenter(lon(itrip),lon(ip1))+lon(itrip))/2.
         tripnorth=(lat(jtrip+1)+lat(jtrip))/2.
         tripsouth=(lat(jtrip)+lat(jtrip-1))/2.

         !print *,itrip,jtrip,im1,ip1,tripeast,tripwest,lon(itrip),lon(ip1),lon(im1)
         !print *,tripeast,tripwest,lon(itrip)



         ! Make sure tripeast is within -180 to 180 degrees from  tripwest
         tripeast=newloncenter(tripwest,tripeast) 


         ! Now go through the runoff cells and map them to this TRIP grid along with
         ! their areal coverage of this TRIP cell.
         do jro=rominj,romaxj
         do iro=1,nrolon

            ! west, east, south and north boundary of runoff cell
            ip1=mod(iro,nrolon)+1
            im1=mod(iro+nrolon-2,nrolon)+1
            roeast =(rolon(iro)+newloncenter(rolon(iro),rolon(ip1)))/2.
            rowest =(newloncenter(rolon(iro),rolon(im1))+rolon(iro))/2.
            ronorth=(rolat(jro+1)+rolat(jro))/2.
            rosouth=(rolat(jro)+rolat(jro-1))/2.

            ! Make sure roeast/rowest is within -180 to 180 degrees from  rowest
            rowest=newloncenter(tripwest,rowest) 
            roeast=newloncenter(tripwest,roeast) 

            ! Overlap of ro cells on trip grid
            dlon_overlap=                                &
               max( min(roeast,tripeast) , rowest ) -  &
               min( max(rowest,tripwest) , roeast ) 
            dlat_overlap=                                   &
               max( min(ronorth,tripnorth) , rosouth ) -  &
               min( max(rosouth,tripsouth) , ronorth ) 

!KAL     print *,'sn',tripsouth,tripnorth,rosouth, ronorth
!KAL     print *,'ew',tripwest,tripeast,rowest, roeast
!KAL     print *,max( min(rowest,tripwest) , tripeast ) &
!KAL            ,min( max(roeast,tripeast) , tripwest ) &
!KAL            ,max( min(ronorth,tripnorth) , tripsouth ) &
!KAL            ,min( max(rosouth,tripsouth) , tripnorth ) 

            ! Approximate area of overlap
            area_overlap=cos(rolat(jro)/radian)* &
                         sin(dlat_overlap/radian)* &
                         sin(dlon_overlap/radian)* &
                         rearth**2

            if (area_overlap>0.001.and.iro==itest.and.jro==jtest) then
               print '(a,2i5,3e14.2)','new:',itrip,jtrip,area_overlap, dlon_overlap,dlat_overlap
               print '(a,2i5)'   ,'--->',iro,jro
               print '(a,4f14.2)','--->sn',tripsouth,tripnorth,rosouth, ronorth
               print '(a,4f14.2)','--->we',tripwest,tripeast,rowest, roeast
               print *,'--------------------'
            end if

            if (area_overlap>0.001) then
               nrocells(iro,jro)=nrocells(iro,jro)+1
               ntmp=nrocells(iro,jro)

               if (ntmp>maxrocells) then
                  stop '(too many cells in mapping)'
               end if
               romapi     (iro,jro,ntmp)=itrip         ! mapping runoff -> TRIP
               ! NB - re-flip of TRIP cell
               romapj     (iro,jro,ntmp)=ny-jtrip+1    ! mapping runoff -> TRIP
               roweights  (iro,jro,ntmp)=area_overlap  ! Area runoff covers in TRIP cell
            end if

            !if (rolon(iro)>50.) stop

         end do
         end do

         !stop

      end do
   end do

   print *,'max number of TRIP cells pr runoff cell           :',maxval(nrocells)
   print *,'max area overlap (km^2) runoff cell/TRIP cell     :',minval(roweights),maxval(roweights)


   open(10,file='rw_maxncells.asc',status='replace')
   write(10,*) maxval(nrocells)
   close(10)
!     open(10,file='weight.uf'   ,form='unformatted',status='replace')
!    write(10) roweights   (:,:,1)
!    close(10)

   
   if (flip_lat==1) then
   ! Flip back the weight
     allocate(tmp(nrolon,maxval(nrocells)))
     allocate(tmpi(nrolon,maxval(nrocells)))
     do j=1,nrolat/2
        tmpi(:,:)=romapi   (:,j,1:maxval(nrocells))
        romapi(:,j,1:maxval(nrocells))=romapi(:,nrolat-j,1:maxval(nrocells))
        romapi(:,nrolat-j,1:maxval(nrocells))=tmpi(:,:)

        tmpi(:,:)=romapj   (:,j,1:maxval(nrocells))
        romapj(:,j,1:maxval(nrocells))=romapj(:,nrolat-j,1:maxval(nrocells))
        romapj(:,nrolat-j,1:maxval(nrocells))=tmpi(:,:)
      
        tmp(:,:)=roweights(:,j,1:maxval(nrocells))
        roweights(:,j,1:maxval(nrocells))=roweights(:,nrolat-j,1:maxval(nrocells))
        roweights(:,nrolat-j,1:maxval(nrocells))=tmp(:,:)
     end do
   endif
!     open(10,file='ro_i.uf'   ,form='unformatted',status='replace')
!    write(10) romapi   (:,:,1)
!    close(10)
!     open(10,file='ro_j.uf'   ,form='unformatted',status='replace')
!    write(10) romapj   (:,:,1)
!    close(10)
!    print *, 'check the boundary',romapi(1,1,1),romapj(1,1,1),roweights(1,1,1)
!     open(10,file='weight.uf'   ,form='unformatted',status='replace')
!    write(10) roweights   (:,:,1)
!    close(10)
    !stop





   
   ! Dump the weights, ncells and mapping into a binary file
   print '(a)','Weights and mappings dumped in rw_cellinfo.uf'
   open(10,file='rw_cellinfo.uf'   ,form='unformatted',status='replace')
   write(10) maxval(nrocells), nrolon, nrolat,  &
             nrocells,                          &
             romapi   (:,:,1:maxval(nrocells)), &
             romapj   (:,:,1:maxval(nrocells)), &
             roweights(:,:,1:maxval(nrocells))
    close(10)




   contains

      real function newloncenter(center,lon)
      implicit none
      real, intent(in) :: center, lon

      newloncenter=lon

      do while (abs(center-newloncenter)>180)
!         print *,center,newloncenter
         if (newloncenter>center) then
            newloncenter=newloncenter-360.
         else if (newloncenter<=center) then
            newloncenter=newloncenter+360.
         end if
      end do
      end function 
   end program


         



