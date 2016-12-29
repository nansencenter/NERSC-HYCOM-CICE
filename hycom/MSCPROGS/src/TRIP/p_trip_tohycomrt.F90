! --- -------------------------------------------------------------------
! --- River routine trip_tohycomrt
! --- -------------------------------------------------------------------
! --- Program to convert netcdf files produced by "trip_riverflow" to
! --- hycom forcing files
! ---
! --- Output from this routine is:
! ---   forcing.rivers.[ab]  
! ---   forcing.precip.[ab]  
! --- -------------------------------------------------------------------
! --- Prerequisites:
! --- You must have run trip_riverweight + trip_riverflow before this routine
! --- For now only the climatology is produced
! --- -------------------------------------------------------------------


module utilities
   implicit none
   contains 

      !tmpflag=1 if point within the circle of size radius and wet, else 0
      integer function myflag(ip,ncells,dist)
      implicit none
      real, intent(in)   :: dist
      integer,intent(in) :: ip,ncells
      myflag=nint(ip*0.5*(1.+sign(1.,ncells-dist)))
      end function myflag

      real function myweight(landdist,rdist,radius,landradius,myflag)
      implicit none
      real, intent(in)   :: landdist, landradius
      real, intent(in)   :: rdist, radius
      integer,intent(in) :: myflag
      myweight=exp(- (landdist/landradius)**2) * exp(- (rdist/radius)**2) * myflag
      end function myweight


      ! This was  done in two locations, so I collected it here to make it
      ! consistent..
      subroutine spread_rivers( &
         river_i,river_j,river_ipiv,river_jpiv,river_darea,river_dweight,nriver, &
         riv_flux,nxd,nyd, &
         myip,tmp3,landdist,idm,jdm, &
         radius, landradius)
      use mod_grid
      implicit none
      integer, intent(in)  :: nriver, nxd, nyd, idm, jdm
      integer, intent(in)  :: river_i   (nriver)
      integer, intent(in)  :: river_j   (nriver)
      integer, intent(in)  :: river_ipiv(nriver)
      integer, intent(in)  :: river_jpiv(nriver)
      real   , intent(in)  :: river_darea(nriver)
      real   , intent(in)  :: river_dweight(nriver)
      real   , intent(in)  :: riv_flux  (nxd,nyd)
      integer, intent(in)  :: myip        (idm,jdm)
      real   , intent(in)  :: landdist      (idm,jdm)
      real   , intent(out) :: tmp3      (idm,jdm)
      real, intent(in) :: landradius, radius
      integer :: iriver,i_river,j_river,i,j,ix,jy,ix2, ncells, tmpflag
      real    :: unif_flux, dist, sum_flux, tmpweight, darea,rdist,sum_weight
      real, external :: spherdist

      tmp3=0.

      do iriver=1,nriver
         ! Location of river in TRIP dataset
         i_river = river_i(iriver)
         j_river = river_j(iriver)

         ! Mapping to hycom
         i = river_ipiv(iriver)
         j = river_jpiv(iriver)

         ! flux from TRIP calc, spread over area
         unif_flux=riv_flux(i_river,j_river)/river_darea(iriver)
         sum_flux=0.
         sum_weight=0.

         ! Spread discharge over area
         ncells=ceiling(radius/scpx(i,j))
         do jy=max(j-ncells,1),min(j+ncells,jdm)
         do ix=i-ncells,i+ncells
            if (periodic) then
               ix2=mod(ix-1+idm,idm)+1
               dist=sqrt(real(ix-i)**2 + real(jy-j)**2)
            else
               ix2=min(idm,max(ix,1))
               dist=sqrt(real(ix2-i)**2 + real(jy-j)**2)
            end if
            if (myip(ix2,jy)==1) then 

               !tmpflag=1 if point within the circle of size radius and wet, else 0
               !tmpflag  =nint(myip(ix2,jy)*0.5*(1.+sign(1.,ncells-dist)))
               !tmpweight=exp(- (landdist(ix2,jy)/landradius)**2) * tmpflag / river_dweight(iriver)
               rdist=spherdist(plon(i,j),plat(i,j),plon(ix2,jy),plat(ix2,jy))

               ! Flag
               tmpflag=myflag(myip(ix2,jy),ncells,dist)

               ! Weights
               tmpweight=myweight(landdist(ix2,jy),rdist,radius,landradius,tmpflag)  / river_dweight(iriver)

               ! Total area
               darea =darea + scpx(ix2,jy)*scpy(ix2,jy)*tmpflag

               ! Cumulative across rivers
               tmp3(ix2,jy)=tmp3(ix2,jy)+tmpweight*unif_flux
               sum_flux=sum_flux + tmpweight*unif_flux
               sum_weight=sum_weight + tmpweight

            end if
         end do ! jy
         end do !ix
         !if (iriver==1) then
         !   print *,"flux check",iriver,sum_weight,river_dweight(iriver)
         !   print *,"wgth check",iriver,unif_flux,sum_flux
         !end if
      end do !iriver
   end subroutine


end module

program trip_tohycom
   use utilities
   use mod_xc
   use mod_za
   use mod_grid
   use mod_trip , nxd=>nx , nyd => ny
   use mod_confmap
   use netcdf
   use m_datetojulian
   use m_handle_err
   use m_nearestpoint
   implicit none
   real :: riv_flux(nxd,nyd,12)
   real,    allocatable ::  mod_riv_flux(:,:), &
      landdist(:,:), meanflux(:,:), tmp1(:,:), tmp2(:,:), tmp3(:,:),&
      fld1(:,:),fld2(:,:)
   integer, allocatable, dimension(:) :: river_i, river_j, river_ipiv, river_jpiv
   real,    allocatable, dimension(:) :: river_darea, river_dweight
   integer, allocatable, dimension(:) :: river_dpnts
   integer nriver,iriver
   integer :: tmpflag
   real    :: tmpweight

   integer :: i,j,ix,ix2,jy,k,ncid,varid,ncells, &
      ipiv,jpiv,ipiv0,jpiv0,ia,ja,ib,jb, sumpoints, nfailed, &
      topi(11), topj(11),nfailed2
   logical :: ass
   real    :: ri, rj, darea, dist, hmin, hmax, a1, a2, a3, a4, mindist, &
      unif_flux,topflx(11),buff,tmp, sum_flux
   real, external :: spherdist

   real, parameter :: d_searchradius=300000. ! search radius
   !We assume that the delta size is of 60 km 
   real, parameter :: d_radius      =60000. ! river radius in meters
   real, parameter :: d_landradius  =200000. ! river radius in meters (alongshore)
   character(len=80) :: tmparg
   real radius,landradius,searchradius
   character(len=80) :: runoff_source

   ! For synoptic calcs
   integer i_river,j_river
   integer max_ncells
   real :: dtime1, dtime2,ddtime, dtime
   
   integer :: day1, year1, hour1
   integer :: day2, year2, hour2
   integer :: numfiles
   integer, dimension(:), allocatable :: nchandle
   character(len=4) :: c4
   character(len=256) :: fname

   integer nrec,irec,recdim, rec1, ind, ref_year, ref_month, ref_day, jday_ref,i2, &
      varid_river, fld1_index, fld2_index
   integer :: xdim,ydim,dims2d(2),dims3d(3)
   real :: w0, w1, dw, rdist
   integer :: new_fld1_index, new_fld2_index
   integer, dimension(:), allocatable :: ncid_handles,time_index, varid_handles
   real   , dimension(:), allocatable :: all_times
   character(len=80) :: time_unit
   character(len=*), parameter :: unit_pattern="days since "
   logical, parameter :: diag_netcdf=.true.
   real*4, parameter :: fillval=-1e20


#if defined(IARGC)
   integer*4, external :: iargc
#endif

   print *,'Routine will calculate river discharge from TRIP+ERAi-derived'
   print *,'river fields. Input is across-shore radius and alongshore-shore radius.  Both in km'
   if (iargc()==6) then
      call getarg(1,runoff_source) 
      call getarg(2,tmparg) ; read(tmparg,*) radius
      call getarg(3,tmparg) ; read(tmparg,*) landradius
      call getarg(4,tmparg) ; read(tmparg,*) dtime1 ! start time in hycom dtime
      call getarg(5,tmparg) ; read(tmparg,*) dtime2 ! start time in hycom dtime
      call getarg(6,tmparg) ; read(tmparg,*) ddtime ! delta time in hycom dtime
      radius=radius*1000.
      landradius=landradius*1000.
   else 
      !runoff_source="erai"
      !radius=d_radius
      !landradius=d_landradius
      print *,"p_trip_tohycomrt.F90 Not correct number of arguments..."
      call exit(1)
   end if
   searchradius=d_searchradius
   print *,"Runoff source: "//trim(runoff_source)
   print *,"Alongshore   radius is ",landradius
   print *,"Across-shore radius is ",radius


   call init_trip()

   ! Read climatologies
   if (trim(runoff_source) == "era40") then 
      print '(a)',"Opening river climatology file trip_era40_clim.nc"
      call handle_err(nf90_open('trip_era40_clim.nc',NF90_NOWRITE,ncid))
   elseif (trim(runoff_source) == "erai") then 
      print '(a)',"Opening river climatology file trip_erai_clim.nc"
      call handle_err(nf90_open('trip_erai_clim.nc',NF90_NOWRITE,ncid))
   else 
      print *,"Unknown runoff source "//trim(runoff_source)
      call exit(1)
   end if

   call handle_err(NF90_INQ_VARID(ncid,'river',varid))
   call handle_err(NF90_GET_VAR(ncid,varid,riv_flux))
   call handle_err(NF90_CLOSE(ncid))
   !Fanf: here we match the average total such that it fits Perry et al. estimate
   !of P-E=1,268 Sv
   !Here with Erai  we obtain 1.3870Sv, we therefore multiply by 0.9142
   tmp=sum(sum(sum(riv_flux,dim=3),dim=2),dim=1)/12.
   print *, 'the world river run off is [m^3.s¯1]:', tmp
   riv_flux=riv_flux*0.9142
   tmp=sum(sum(sum(riv_flux,dim=3),dim=2),dim=1)/12.
   print *, 'the world river run off  after scaling is:', tmp

   ! Init HYCOM fields
   call xcspmd()
   call zaiost()
   call get_grid()
   call initconfmap(idm,jdm)

   ! First step Go through TRIP points, and place rivers on map (rivers placed
   ! on one mdoel pixel). Use climatology to investigate rivers
   allocate(river_i(nxd*nyd))
   allocate(river_j(nxd*nyd))
   allocate(river_ipiv(nxd*nyd))
   allocate(river_jpiv(nxd*nyd))
   iriver=0
   nfailed=0
   do j=1,nyd  !going through Trip grid
   do i=1,nxd
   if (any(riv_flux(i,j,:)>1e-4)) then

      call ll2gind(lon(i),lat(j),ri,rj)
      ipiv0=floor(ri)
      jpiv0=floor(rj)
      ! If off grid neglect altogether
      if (ipiv0>=1.and.ipiv0<idm.and.jpiv0>=1.and.jpiv0<jdm) then
         ! Approximate local search distance
         ncells=ceiling(searchradius/scpx(ipiv0,jpiv0))
         ! If on land search within searchradius for a ocean point
         if (ip(ipiv0,jpiv0)==0) then
            mindist=ncells+1
            ipiv=-1
            jpiv=-1
            do jy=max(jpiv0-ncells,1),min(jpiv0+ncells,jdm)
            do ix=ipiv0-ncells,ipiv0+ncells

               if (periodic) then
                  ix2=mod(ix-1+idm,idm)+1
                  dist=sqrt(real(ix-ipiv0)**2 + real(jy-jpiv0)**2)
               else
                  ix2=min(idm,max(ix,1))
                  dist=sqrt(real(ix2-ipiv0)**2 + real(jy-jpiv0)**2)
               end if

               if (dist<mindist .and. ip(ix2,jy)==1) then
                  ipiv=ix2
                  jpiv=jy
                  mindist=dist
               end if
            end do
            end do

            if (ipiv/=-1  .and. jpiv/=-1) then
               iriver=iriver+1
               river_i(iriver)=i
               river_j(iriver)=j
               river_ipiv(iriver)=ipiv
               river_jpiv(iriver)=jpiv
            else
               !print '(a,2f14.2,a)','Failed to place river at', lon(i),lat(j),' to model'
               nfailed=nfailed+1
            end if
         ! on ocean...
         else
            ! Make sure river is placed on a coastal cell
            nfailed2=0
            mindist=ncells+1
            ipiv=-1
            jpiv=-1
            do jy=max(jpiv0-ncells,2),min(jpiv0+ncells,jdm-1)
            do ix=ipiv0-ncells,ipiv0+ncells


               if (periodic) then
                  ix2=mod(ix-1+idm,idm)+1
                  ia=mod(ix2-2+idm,idm)+1 ! ipiv0-1
                  ib=mod(ix2+idm,idm)+1   ! ipiv0+1
                  dist=sqrt(real(ix-ipiv0)**2 + real(jy-jpiv0)**2)
               else
                  ix2=min(idm,max(ix,1))
                  ia=min(idm,max(ix2-1,1))
                  ib=min(idm,max(ix2+1,1))
                  dist=sqrt(real(ix2-ipiv0)**2 + real(jy-jpiv0)**2)
               end if
               ja=min(jdm,max(jy-1,1))
               jb=min(jdm,max(jy+1,1))

               ! Sum of neighbours and local point
               sumpoints=sum(ip(ia,ja:jb))+sum(ip(ix2,ja:jb))+sum(ip(ib,ja:jb))

               ! Accept point if it is ocean and has three or more land eighbours
               if (ip(ix2,jy)==1 .and. sumpoints <=6 .and. dist<mindist) then
                  mindist=dist
                  ipiv=ix2
                  jpiv=jy
               end if
            end do
            end do

            if (ipiv/=-1 .and. jpiv/=-1) then
               iriver=iriver+1
               river_i(iriver)=i
               river_j(iriver)=j
               river_ipiv(iriver)=ipiv
               river_jpiv(iriver)=jpiv
            else
               print '(a,2f14.2,a)','Failed to connect river at', lon(i),lat(j),' to model coastal point'
               nfailed2=nfailed2+1
            end if
         end if ! ocean/not ocean cell test
      end if ! on model grid test
   end if ! has river flux > 0 test
   end do
   end do
   nriver=iriver
   print '(a,i5,a)', 'Found ',nriver,' rivers inside domain'
   print '(a,i5,a)','Failed to place ',nfailed,' river points dry in model'
   print '(a,i5,a)','Failed to place ',nfailed2,' river points wet in model'
   !print *, 'the global river run off within the model grid is [m^3.s¯1]:', tmp


   !Fanf: ensure that the value entered are larger than the grid size
   print *, 'river delta set to (in meter):',radius
   print *, 'river plume diffuse with a Gaussian until (in meter):',landradius
   if (landradius>searchradius .or. radius>searchradius .or.&
   landradius<maxval(scpx) .or. radius< maxval(scpx)) then
     print *, 'radius < grid size or search radius too small'
     stop
   endif

   ! 2nd step - Find land distances
   allocate(landdist(idm,jdm))
   print '(a)','Calculating land distances '
   print '(a,i5)','Max=',jdm
   do j=1,jdm
   if (mod(j,jdm/8)==0.and.j/=jdm) then
      write (6,'(a,i5.5)',advance='no') '...',j
      call flush(6)
   end if
   do i=1,idm


      ! If on ocean search within searchradius for a land point
      if (ip(i,j)==1) then
         ! Approximate local search distance within 
         ncells=ceiling(searchradius/scpx(i,j))
         mindist=100*landradius ! initialize to a very large value
         do jy=max(j-ncells,1),min(j+ncells,jdm)
         do ix=i-ncells,i+ncells

            if (periodic) then
               ix2=mod(ix-1+idm,idm)+1
            else
               ix2=min(idm,max(ix,1))
            end if
            if (ip(ix2,jy)==0) then
               dist=spherdist(plon(i,j),plat(i,j),plon(ix2,jy),plat(ix2,jy))
               mindist=min(mindist,dist)
            end if
         end do
         end do
         landdist(i,j)=mindist
      else
         landdist(i,j)=0.
      end if

   end do
   end do
   write (6,'(a,i5.5)',advance='yes') '...',jdm
   !print *,'Calculated land distance', minval(landdist),maxval(landdist)
   print *
   call zaiopf('landdist.a', 'replace', 909) 
   call zaiowr(landdist,ip,.false.,hmin,hmax,909,.true.)
   call zaiocl(909)



   ! Calculate max cell distance possible
   max_ncells=maxval(radius/scpx)
   print '(a,i5)', 'Max cell radius for spreading=',max_ncells


   ! 2nd step - spread from river point to neighbouring cells
   print '(a)', 'Spreading discharge over ocean grid cells (land-distance approach)'
   print '(a,i5)', 'Max=',jdm
   allocate(river_dweight(nriver))
   allocate(river_darea(nriver))
   allocate(river_dpnts(nriver))
   river_darea=0.
   river_dpnts=0
   print '(a,i5)','Max=',nriver
   do iriver=1,nriver
      if (mod(iriver,nriver/8)==0.and.iriver/=nriver) then
         write (6,'(a,i5.5)',advance='no') '...',iriver
         call flush(6)
      end if

      ! Location of river in TRIP dataset
      i_river = river_i(iriver)
      j_river = river_j(iriver)

      ! Mapping to hycom
      i = river_ipiv(iriver)
      j = river_jpiv(iriver)


      ! This is a simple estimate based on the local grid distance
      ! spread over neighbour ocean cells in radius of 200 km (for now)
      ! uniform spread(for now)
      ncells=ceiling(radius/scpx(i,j))
      ! 1) tag cells in tmpflag - discharge area  in darea
      darea=0. ! discharge area

      do jy=max(j-ncells,1),min(j+ncells,jdm)
      do ix=i-ncells,i+ncells
         if (periodic) then
            ix2=mod(ix-1+idm,idm)+1
            dist=sqrt(real(ix-i)**2 + real(jy-j)**2)
         else
            ix2=min(idm,max(ix,1))
            dist=sqrt(real(ix2-i)**2 + real(jy-j)**2)
         end if
         if (ip(ix2,jy)==1) then 

            ! Flag
            !tmpflag=nint(ip(ix2,jy)*0.5*(1.+sign(1.,ncells-dist)))
            tmpflag=myflag(ip(ix2,jy),ncells,dist)
            rdist=spherdist(plon(i,j),plat(i,j),plon(ix2,jy),plat(ix2,jy))


            ! Weights
            !tmpweight=exp(- (landdist(ix2,jy)/landradius)**2) * tmpflag
            tmpweight=myweight(landdist(ix2,jy),rdist,radius,landradius,tmpflag)

            ! Total area
            darea =darea + scpx(ix2,jy)*scpy(ix2,jy)*tmpflag

            ! Sum of discharge area, weights and points for river
            river_dpnts  (iriver) = river_dpnts  (iriver) + tmpflag
            river_darea  (iriver) = river_darea  (iriver) + scpx(ix2,jy)*scpy(ix2,jy)*tmpflag
            river_dweight(iriver) = river_dweight(iriver) + tmpweight
         end if
      end do
      end do
   end do
   write (6,'(a,i5)',advance='yes') '...',nriver









   print '(a)', 'Spreading discharge over ocean grid cells for climatology'
   open (unit=909, file='forcing.rivers.b',  status='replace', action='write')
   write(909,'(a)') 'River mass fluxes from TRIP+ERAI  climatology'
   write(909,'(a)') ''
   write(909,'(a)') ''
   write(909,'(a)') ''
   write(909,'(a,2i5)') 'i/jdm = ',idm,jdm
   call zaiopf('forcing.rivers.a', 'replace', 909) 
   allocate(mod_riv_flux(idm,jdm))
   do k=1,12

      ! Its big, its ugly, but at least its all in the same place
      call spread_rivers( &
         river_i(1:nriver),river_j(1:nriver),river_ipiv(1:nriver),river_jpiv(1:nriver), &
         river_darea(1:nriver),river_dweight(1:nriver),nriver, &
         riv_flux(:,:,k),nxd,nyd, &
         ip,mod_riv_flux,landdist,idm,jdm, &
         radius, landradius)

      call zaiowr(mod_riv_flux(:,:),ip,.false.,hmin,hmax,909,.true.)
      write(909,'(" rivers:month,range = ",i2.2,2e16.8)')  k,hmin,hmax
   end do !k
   close(909)
   call zaiocl(909)


   !KAL - now process synoptic time series. dtime1 to dtime2 in steps of ddtime
   ! calculate start and end time in years.
   !dtime1 = 40000
   !dtime2 = 40400
   !ddtime = 10. 
   call forday(dtime1,3,year1,day1,hour1)
   call forday(dtime2,3,year2,day2,hour2)
   !print *,year1,day1,hour1
   !print *,year2,day2,hour2

   ! Open netcdf files in that time range. For this we do not copy stuff in here
   ! (files may be huge)
   numfiles=year2-year1+1
   !allocate(nchandle(numfiles))
   allocate(all_times   (numfiles*366)) ! Collect all times from netcdf files here
   allocate(ncid_handles(numfiles*366)) ! For each time, a handle to the netcdff file containing this time 
   allocate(varid_handles (numfiles*366)) ! For each time, a handle to the netcdff file containing this time 
   allocate(time_index  (numfiles*366)) ! For each time, a handle to the netcdff file containing this time 


   rec1=1
   do i=1,numfiles
      write(c4,'(i4)') year1+(i-1)
      fname=trim(trip_path0)//"/"//trim(runoff_source)//"/trip_"//trim(runoff_source)//"_"//c4//".nc"
      print '(a)',"Opening synoptic river forcing file "//trim(fname)
      call handle_err(nf90_open(trim(fname),NF90_NOWRITE,ncid))
      call handle_err(nf90_inq_dimid(ncid, 'time', recdim))
      call handle_err(Nf90_inq_varid(ncid, 'river', varid_river ))
      call handle_err(Nf90_inq_varid(ncid, 'time', varid ))
      call handle_err(nf90_Inquire_Dimension(ncid, recdim, len=nrec))
      call handle_err(nf90_get_var(ncid, varid,all_times(rec1:rec1+nrec-1)))
      call handle_err(nf90_get_att(ncid, varid,"units",time_unit))
      ncid_handles(rec1:rec1+nrec-1)=ncid
      varid_handles(rec1:rec1+nrec-1)=varid_river
      do i2=rec1,rec1+nrec-1
         time_index(i2) = i2-rec1+1
      end do

      ! This is the assumed form of the cf time string ... All of the below is
      ! very sensitive to how the unit is written...
      ind = index(time_unit,unit_pattern)
      if (ind ==0 ) then
         print '(a)',"Unable to parse time string from time unit "//trim(time_unit)
         call exit(1)
      end if

      !print *,"time unit = "//time_unit
      !print *,ind
      ind=len(unit_pattern)
      print '(a)',"Ref time:"//trim(time_unit(ind:ind+20))

      ! Year 
      !print *,time_unit(ind+1:ind+4)
      read(time_unit(ind+1:ind+4),*) ref_year

      ! Month 
      !print *,time_unit(ind+6:ind+7)
      read(time_unit(ind+6:ind+7),*) ref_month

      ! Day of month 
      !print *,time_unit(ind+9:ind+10)
      read(time_unit(ind+9:ind+10),*) ref_day



      ! Day of reference date rel to 1st of that year
      jday_ref = datetojulian(ref_year,ref_month,ref_day,ref_year,1,1)

      ! Convert times to dtime
      call dayfor(dtime,3,ref_year,1,0)
      do i2 =rec1,rec1+nrec-1
         all_times(i2)=dtime+jday_ref+all_times(i2)
      end do
      !print  *, all_times(rec1:rec1+nrec-1)

      rec1=rec1+nrec
   end do
   nrec=rec1-1

   ! Test increasing ts
   do i=2,nrec
      if (.not. all_times(i) > all_times(i-1)) then
         print *,"Time series not increasing !"
         print *,all_times(i),all_times(i-1)
         call exit(1)
      end if
   end do




   ! Create netcdf output for diagnostics
   if (diag_netcdf) then
      call handle_err(nf90_create("precip.nc",NF90_CLOBBER,ncid))
      call handle_err(nf90_def_dim(ncid,'idm'    ,idm,xdim))
      call handle_err(nf90_def_dim(ncid,'jdm'    ,jdm,ydim))
      call handle_err(nf90_def_dim(ncid,'time',nf90_unlimited,recdim))
      dims2d=(/xdim,ydim/)
      dims3d=(/xdim,ydim,recdim/)

      call handle_err(NF90_DEF_VAR(ncid,'lon',NF90_Float,dims2d,varid))
      call handle_err(NF90_ENDDEF(ncid))
      call handle_err(NF90_PUT_VAR(ncid,varid,plon))

      call handle_err(NF90_REDEF(ncid))
      call handle_err(NF90_DEF_VAR(ncid,'lat',NF90_Float,dims2d,varid))
      call handle_err(NF90_ENDDEF(ncid))
      call handle_err(NF90_PUT_VAR(ncid,varid,plat))

      call handle_err(NF90_REDEF(ncid))
      call handle_err(NF90_DEF_VAR(ncid,'depths',NF90_Float,dims2d,varid))
      call handle_err(NF90_ENDDEF(ncid))
      call handle_err(NF90_PUT_VAR(ncid,varid,depths))

      call handle_err(NF90_REDEF(ncid))
      call handle_err(NF90_DEF_VAR(ncid,'riverflux',NF90_Float,dims3d,varid))
      call handle_err(NF90_PUT_ATT(ncid,varid,'units','m s-1'))
      call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',fillval))
      call handle_err(NF90_ENDDEF(ncid))
   end if


   ! Now loop over dtimes, Find appropriate times from open netcdf files, and
   ! create forcing. Fake precipitation
   print '(a)', 'Spreading discharge over ocean grid cells for realtime'
   open (unit=909, file='forcing.precip.b',  status='replace', action='write')
   write(909,'(a)') 'River mass fluxes from TRIP+'//trim(runoff_source)//' implemented as precip'
   write(909,'(a)') 'precip (m s**-1)'
   write(909,'(a)') ''
   write(909,'(a)') ''
   write(909,'(a,2i5)') 'i/jdm = ',idm,jdm
   call zaiopf('forcing.precip.a', 'replace', 909)

   dtime=dtime1
   allocate(tmp1(idm,jdm))
   allocate(tmp2(idm,jdm))
   Allocate(tmp3(idm,jdm))
   allocate(fld1(idm,jdm))
   allocate(fld2(idm,jdm))
   fld1_index=-1
   fld2_index=-1
   irec=1
   do while (dtime <= dtime2)

      print '(a,f14.4)','Calculating "river precip" for time ',dtime

      !print *,"fld_indexes:",fld1_index,fld2_index

      ! Find times in data. times are assumed increasing ...
      i2=-1
      do i=1,nrec-1
         if (dtime < all_times(i+1) .and. dtime >= all_times(i)) then
            i2 = i
         end if
      end do
      if (i2 < 0 ) then
         print *,"Unable to find index in time range"
         print '(i5,3f20.4)',i2,dtime,all_times(i2),all_times(i2+1)
         call exit(1)
      !else 
      !   print '(i4,3f20.4)',i2,dtime,all_times(i2),all_times(i2+1)
      end if

      ! Linear weights
      dw = all_times(i2+1)-all_times(i2)
      w1 = (dtime - all_times(i2  ))/dw
      w0 = (all_times(i2+1)-dtime  )/dw

      ! Read fields ( if necessary)
      do k=1,2
         ind = i2+k-1

         ! Already read
         if     (ind == fld1_index) then 
            if (k==1) then
               tmp1           = fld1
               new_fld1_index = fld1_index
            else 
               tmp2           = fld1
               new_fld2_index = fld1_index
            end if
         ! Already read
         elseif (fld2_index == ind) then
            if (k==1) then
               tmp1           = fld2
               new_fld1_index = fld2_index
            else 
               tmp2           = fld2
               new_fld2_index = fld2_index
            end if
         ! Not read. Read it ...
         else 

            !print '(a,i8,a,i8)',"Reading k=",k,"index=",ind
            tmp3=0.
            call handle_err(NF90_GET_VAR(ncid_handles(ind),varid_handles(ind),riv_flux(:,:,k),start=(/1,1,time_index(ind)/)))
            !print *,"riv_flux(:,:,k):",minval(riv_flux(:,:,k)),maxval(riv_flux(:,:,k))


            ! Its big, its ugly, but at least its all in the same place
            call spread_rivers( &
               river_i(1:nriver),river_j(1:nriver),river_ipiv(1:nriver), &
               river_jpiv(1:nriver),river_darea(1:nriver),river_dweight(1:nriver),nriver, &
               riv_flux(:,:,k),nxd,nyd, &
               ip,tmp3,landdist,idm,jdm, &
               radius, landradius)


            if (k==1) then
               tmp1           = tmp3
               new_fld1_index = ind
            else 
               tmp2           = tmp3
               new_fld2_index = ind
            end if

         end if
      end do !k

      ! update handles
      fld1=tmp1
      fld2=tmp2
      !print *,"fld1:",minval(fld1),maxval(fld1)
      !print *,"fld2:",minval(fld2),maxval(fld2)
      fld1_index=new_fld1_index
      fld2_index=new_fld2_index

      ! Final river field
      mod_riv_flux = w0*fld1 + w1*fld2

      call zaiowr(mod_riv_flux(:,:),ip,.false.,hmin,hmax,909,.true.)
      write(909,'("precip:dtime1,range = ",3e16.8)')  dtime,hmin,hmax

      where(ip==0) mod_riv_flux=-1e20
      if (diag_netcdf) call handle_err(NF90_PUT_VAR(ncid,varid,mod_riv_flux,start=(/1,1,irec/)))

      dtime = dtime + ddtime
      irec=irec+1
   end do
   if (diag_netcdf)  call handle_err(NF90_CLOSE(ncid))
   close(909)
   call zaiocl(909)

   



         

!   ! Diagnostics: 10 top rivers - Locate them
!   topi  =-1
!   topj  =-1
!   topflx=0.
!
!   do j=1,jdm
!   do i=1,idm
!
!      ! "Stupid sort" - but we only keep the top 10
!      if (any(meanflux(i,j) > topflx)) then
!         do k=1,10
!            if (topflx(k)<=meanflux(i,j)) then
!               topflx(k+1:11) = topflx(k:10)
!               topi  (k+1:11) = topi  (k:10)
!               topj  (k+1:11) = topj  (k:10)
!               topflx(k) = meanflux(i,j)
!               topi  (k) = i
!               topj  (k) = j
!               exit
!            end if
!         end do
!      end if
!
!   end do
!   end do
!   print *
!   print '(a)','Top Rivers by month and discharge -- Location:'
!   print 100,'i     ',topi(1:10)
!   print 100,'j     ',topj(1:10)
!   do j=1,10
!      topflx(j)=plon(topi(j),topj(j))
!   end do
!   print 200,'lon   ',topflx(1:10)
!   do j=1,10
!      topflx(j)=plat(topi(j),topj(j))
!   end do
!   print 200,'lat   ',topflx(1:10)
!   print '(a)','Top Rivers by month and discharge -- discharge [m^3 s^-1]'
!   do k=1,12
!      do j=1,10
!         topflx(j)=mod_riv_flux(topi(j),topj(j),k)
!      end do
!      print 300,k,topflx(1:10)
!   end do
!   print '(a)','Top Rivers annual average discharge -- discharge [km^3 year^-1]'
!   do j=1,10
!      print 400,j,meanflux(topi(j),topj(j))*3600*24*364/1000000000
!   end do
      
100 format(a8,  10i14)
200 format(a8,  10f14.3)
300 format("Month ",i2,  10f14.3)
400 format("River ",i2,  10f14.3)






      
   





end program


   



