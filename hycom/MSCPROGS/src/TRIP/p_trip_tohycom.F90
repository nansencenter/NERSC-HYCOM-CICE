! --- -------------------------------------------------------------------
! --- River routine trip_tohycom
! --- -------------------------------------------------------------------
! --- Program to convert netcdf files produced by "trip_riverflow" to
! --- hycom forcing files
! ---
! --- Output from this routine is:
! ---   forcing.rivers.[ab]  
! --- -------------------------------------------------------------------
! --- Prerequisites:
! --- You must have run trip_riverweight + trip_riverflow before this routine
! --- For now only the climatology is produced
! --- -------------------------------------------------------------------

program trip_tohycom
   use mod_xc
   use mod_za
   use mod_grid
   use mod_trip , nxd=>nx , nyd => ny
   use mod_confmap
   use netcdf
   use m_handle_err
   use m_nearestpoint
   implicit none
   real :: riv_flux(nxd,nyd,12)
   real,    allocatable :: mod_riv_flux(:,:,:), mod_riv_flux_spread(:,:,:), &
      landdist(:,:), tmpweight(:,:), meanflux(:,:)
   integer, allocatable :: tmpflag(:,:)

   integer :: i,j,ix,ix2,jy,k,ncid,varid,ncells, &
      ipiv,jpiv,ipiv0,jpiv0,ia,ja,ib,jb, sumpoints, nfailed, &
      topi(11), topj(11),nfailed2
   logical :: ass
   real    :: ri, rj, darea, dist, hmin, hmax, a1, a2, a3, a4, mindist, &
      unif_flux,topflx(11),buff,tmp
   real, external :: spherdist

   real, parameter :: d_searchradius=300000. ! search radius
   !We assume that the delta size is of 60 km 
   real, parameter :: d_radius      =60000. ! river radius in meters
   real, parameter :: d_landradius  =200000. ! river radius in meters (alongshore)
   character(len=80) :: tmparg
   real radius,landradius,searchradius
   character(len=80) :: runoff_source


#if defined(IARGC)
   integer*4, external :: iargc
#endif

   print *,'Routine will calculate river discharge from TRIP+ERAi-derived'
   print *,'river fields. Input is across-shore radius and alongshore-shore radius.  Both in km'
   if (iargc()==3) then
      call getarg(1,runoff_source) 
      call getarg(2,tmparg) ; read(tmparg,*) radius
      call getarg(3,tmparg) ; read(tmparg,*) landradius
      radius=radius*1000.
      landradius=landradius*1000.
   else 
   !   runoff_source="erai"
   !   radius=d_radius
   !   landradius=d_landradius
      print *,"p_trip_tohycom.F90 Not correct number of arguments..."
      call exit(1)
   end if
   searchradius=d_searchradius
   print *,"Runoff source: "//trim(runoff_source)
   print *,"Alongshore   radius is ",landradius
   print *,"Across-shore radius is ",radius


   call init_trip()

   ! Read climatology (pointwise)
   if (trim(runoff_source) == "era40") then 
      call handle_err(nf90_open('trip_era40_clim.nc',NF90_CLOBBER,ncid))
   elseif (trim(runoff_source) == "erai") then 
      call handle_err(nf90_open('trip_erai_clim.nc',NF90_CLOBBER,ncid))
   else 
      print *,"Unknown runoff source "//trim(runoff_source)
      call exit(1)
   end if

   call handle_err(NF90_INQ_VARID(ncid,'river',varid))
   call handle_err(NF90_GET_VAR(ncid,varid,riv_flux))
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
   ! on one mdoel pixel)
   allocate(mod_riv_flux(idm,jdm,12))
   mod_riv_flux=0.
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
               mod_riv_flux(ipiv,jpiv,:)=mod_riv_flux(ipiv,jpiv,:)+riv_flux(i,j,:)
            else
               !print '(a,2f14.2,a)','Failed to place river at', lon(i),lat(j),' to model'
               !mod_riv_flux(ipiv0,jpiv0,:)=mod_riv_flux(ipiv0,jpiv0,:)+riv_flux(i,j,:)
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
               mod_riv_flux(ipiv,jpiv,:)=mod_riv_flux(ipiv,jpiv,:)+riv_flux(i,j,:)
            else
               print '(a,2f14.2,a)','Failed to connect river at', lon(i),lat(j),' to model coastal point'
               mod_riv_flux(ipiv0,jpiv0,:)=mod_riv_flux(ipiv0,jpiv0,:)+riv_flux(i,j,:)
               nfailed2=nfailed2+1
            end if
         end if ! ocean/not ocean cell test
      end if ! on model grid test
   end if ! has river flux > 0 test
   end do
   end do
   print '(a,i5,a)','Failed to place ',nfailed,' river points dry in model'
   print '(a,i5,a)','Failed to place ',nfailed2,' river points wet in model'
   tmp=sum(sum(sum(mod_riv_flux,dim=3),dim=2),dim=1)/12.
   print *, 'the global river run off within the model grid is [m^3.s¯1]:', tmp


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









   ! 2nd step - spread from river point to neighbouring cells
   print '(a)', 'Spreading discharge over ocean grid cells (land-distance approach)'
   print '(a,i5)', 'Max=',jdm
   allocate(mod_riv_flux_spread(idm,jdm,12))
   allocate(tmpflag(idm,jdm))
   allocate(tmpweight(idm,jdm))
   mod_riv_flux_spread=0.
   do j=1,jdm
   if (mod(j,jdm/8)==0.and.j/=jdm) then
      write (6,'(a,i5.5)',advance='no') '...',j
      call flush(6)
   end if
   do i=1,idm
   if (any(mod_riv_flux(i,j,:)>1e-4)) then
      ! This is a simple estimate based on the local grid distance
      ! spread over neighbour ocean cells in radius of 200 km (for now)
      ! uniform spread(for now)
      ncells=ceiling(radius/scpx(i,j))
      ! 1) tag cells in tmpflag - discharge area  in darea
      tmpflag=0
      tmpweight=0
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

         !tmpflag=1 if point within the circle of size radius and wet, else 0
         tmpflag(ix2,jy)=nint(ip(ix2,jy)*0.5*(1.+sign(1.,ncells-dist)))
         darea =darea + scpx(ix2,jy)*scpy(ix2,jy)*tmpflag(ix2,jy)
         ! Weights
         tmpweight(ix2,jy)=exp(- (landdist(ix2,jy)/landradius)**2) * tmpflag(ix2,jy)
      end do
      end do
      !if (darea<1e-4) then
      !   print *,'Zero discharge area for river '
      !   print *,mod_riv_flux(i,j,:)
      !   print *,darea*1e-6
      !end if


      if (darea>1e-4) then

         tmpweight=sum(tmpflag)*tmpweight/sum(tmpweight)
         ! Here we ensure that sum(tmpweight)=nb of weight points in the
         ! neighborood
         do k=1,12
            unif_flux=mod_riv_flux(i,j,k)/darea
            ! Use weights
            mod_riv_flux_spread(:,:,k)=mod_riv_flux_spread(:,:,k)+tmpweight*unif_flux
         end do
      end if

   end if
   end do
   end do
   write (6,'(a,i5)',advance='yes') '...',jdm
   print *,minval(mod_riv_flux_spread),maxval(mod_riv_flux_spread)
   ! 3rd step - Save to forcing files
   open (unit=909, file='forcing.rivers.b',  &
         status='replace', action='write')
   write(909,'(a)') 'River mass fluxes from TRIP+ERAI '
   write(909,'(a)') ''
   write(909,'(a)') ''
   write(909,'(a)') ''
   write(909,'(a,2i5)') 'i/jdm = ',idm,jdm
   call zaiopf('forcing.rivers.a', 'replace', 909) 
   do k=1,12
      call zaiowr(mod_riv_flux_spread(:,:,k),ip,.false.,hmin,hmax,909,.true.)
      write(909,'(" rivers:month,range = ",i2.2,2e16.8)')  k,hmin,hmax
   end do
   close(909)
   call zaiocl(909)


         

   ! Diagnostics: 10 top rivers - Locate them
   allocate(meanflux(idm,jdm))
   meanflux=sum(mod_riv_flux,dim=3)/12
   topi  =-1
   topj  =-1
   topflx=0.

   do j=1,jdm
   do i=1,idm

      ! "Stupid sort" - but we only keep the top 10
      if (any(meanflux(i,j) > topflx)) then
         do k=1,10
            if (topflx(k)<=meanflux(i,j)) then
               topflx(k+1:11) = topflx(k:10)
               topi  (k+1:11) = topi  (k:10)
               topj  (k+1:11) = topj  (k:10)
               topflx(k) = meanflux(i,j)
               topi  (k) = i
               topj  (k) = j
               exit
            end if
         end do
      end if

   end do
   end do
   print *
   print '(a)','Top Rivers by month and discharge -- Location:'
   print 100,'i     ',topi(1:10)
   print 100,'j     ',topj(1:10)
   do j=1,10
      topflx(j)=plon(topi(j),topj(j))
   end do
   print 200,'lon   ',topflx(1:10)
   do j=1,10
      topflx(j)=plat(topi(j),topj(j))
   end do
   print 200,'lat   ',topflx(1:10)
   print '(a)','Top Rivers by month and discharge -- discharge [m^3 s^-1]'
   do k=1,12
      do j=1,10
         topflx(j)=mod_riv_flux(topi(j),topj(j),k)
      end do
      print 300,k,topflx(1:10)
   end do
   print '(a)','Top Rivers annual average discharge -- discharge [km^3 year^-1]'
   do j=1,10
      print 400,j,meanflux(topi(j),topj(j))*3600*24*364/1000000000
   end do
      
100 format(a8,  10i14)
200 format(a8,  10f14.3)
300 format("Month ",i2,  10f14.3)
400 format("River ",i2,  10f14.3)






      
   





end program


   



