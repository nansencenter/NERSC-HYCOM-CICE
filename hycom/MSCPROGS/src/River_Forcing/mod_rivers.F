      module mod_rivers

      private
      real   :: rradius = 80000. ! Radius in meters
      real   :: alongshoreradius = 200000. ! Radius in meters

      type river_data
        logical banned ! 
        logical active
        character(len=60) string
        real annual
        real flux(12)
        real lat,lon
        integer ip,jp
        integer ip_land,jp_land
        real ipr,jpr
        real discharge_area
        integer ip_baro,jp_baro
      end type river_data
      type(river_data), save, allocatable :: river(:)

!#if defined(RIVER_PLUME)
!      integer, allocatable, save :: river_ip(:)
!      integer, allocatable, save :: river_jp(:)
!      real,    allocatable, save :: river_weight(:)
!      logical, allocatable, save :: river_direction(:,:)
!#endif

      integer, save :: riv_nr=0

      character(len=*) , parameter :: flriver = './Data/rivers.dat'

      public :: rivers_to_hycom, rradius, alongshoreradius

      contains




!#######################################################################
!##############       SUBROUTINE RIVERS_TO_HYCOM        ################
!############## --------------------------------------- ################
!############## Only routine visible from the outside   ################
!############## called from hycom to get riverine flux. ################
!#######################################################################
      subroutine rivers_to_hycom(modlon,modlat,scpx,scpy,depths)
      use mod_xc
      use mod_za
      use netcdf
      use m_handle_err
      implicit none
      real, intent(in), dimension(idm,jdm) :: modlon, modlat, depths,
     &   scpx,scpy

      integer :: lgth, mo
      real:: river_flux(1:idm,1:jdm,12)
      real :: hmin,hmax
      real, dimension(idm,jdm) :: tmp
      integer :: i,j,k,new_riv_nr
      integer, dimension(idm,jdm) :: ip
      logical :: ldiagriver=.true.

      integer :: varid, idmid, jdmid, rdimid,ncid, ierr, rivdim,
     &   N_2
      real    :: fillv

      ! Get numbers of rivers 
      riv_nr=rivernr()
      write(lp,*) 'Number of rivers=',riv_nr
      if (riv_nr<=0) then
         print *,'Error reading rivers or no rivers in rivers.dat'
         call exit(1)
      end if
      allocate(river(riv_nr))

      ! Read rivers
      write(lp,*) 'calling readrivers'
      call readrivers(modlon,modlat)

      ! readrivers sets riv_nr
      if (riv_nr<=0) then
         print *,'No rivers in domain - no forcing produced'
         call exit(1)
      end if


      ! Accumulate rivers
      write(lp,*) 'calling accrivers'
      river_flux=0.
      call accrivers(modlon,modlat,scpx,scpy,depths,river_flux)


      ! Dump rivers to hycom-style files
      !lgth = len_trim(flnmfor)
      print *,'Dumping hycom forcing fields'
      if (mnproc==1) then
         open (unit=909,
     &   file='forcing.rivers.b', 
     &         status='replace', action='write')
         write(909,'(a)') 'River mass fluxes '
         write(909,'(a)') ''
         write(909,'(a)') ''
         write(909,'(a)') ''
         write(909,'(a,2i5)') 'i/jdm = ',idm,jdm
      end if

      where (depths>.1) 
         ip=1
      elsewhere
         ip=0
      endwhere

      call zaiopf('forcing.rivers.a', 'replace', 909) 
      do mo=1,12
         call zaiowr(river_flux(:,:,mo),ip,.false.,hmin,hmax,909,.true.)
         if(mnproc==1) 
     &     write(909,'(" rivers:month,range = ",i2.2,2e16.8)') 
     &     mo,hmin,hmax
      end do
      if (mnproc==1) close(909)
      call zaiocl(909)

      if (riv_nr<=0) then
         print *,'No rivers in river.dat in model area - I will quit '
         stop 
      end if


      ! Diagnose rivers
      print *,'Dumping tecplot diagnostic file'
      if (ldiagriver) then
         if (mnproc==1) then
            OPEN(10,FILE='rivers.tec',STATUS='UNKNOWN')
            WRITE(10,'(''TITLE= "River Flux fields"'')')
            write(10,'(a)')
     &        'VARIABLES="i" "j" "lon" "lat" "depths" "flux[mm/day]" '
         end if
         do mo=1,12
            if (mo > 1.and. mnproc==1) then
               WRITE(10,'(''ZONE I='',I3,'', J='',I3,'', F=BLOCK'')')
     &            idm,jdm
               write(10,'(a)')'D=(1,2,3,4,5)'
            elseif (mnproc==1) then
               WRITE(10,'(''ZONE I='',I3,'', J='',I3,'', F=BLOCK'')')
     &            idm,jdm
               WRITE(10,99)((i,i=1,idm),j=1,jdm)
               WRITE(10,99)((j,i=1,idm),j=1,jdm)
               WRITE(10,100)((modlon(i,j),i=1,idm),j=1,jdm)
               WRITE(10,100)((modlat(i,j),i=1,idm),j=1,jdm)
               WRITE(10,100)((depths(i,j),i=1,idm),j=1,jdm)
            end if
            !call xcaget(tmp,river_flux(:,:,mo),0)
            tmp=river_flux(:,:,mo)
            tmp=tmp*86400.*1e3
            if (mnproc==1) WRITE(10,100)((tmp(i,j),i=1,idm),j=1,jdm)
         end do
         ! Place river pivot points markings on the grid
         if (mnproc==1) then
         do k=1,riv_nr
            if (river(k)%active) then
            write(10,400) river(k)%ipr,river(k)%jpr
            write(10,'(i4)') 1
            write(10,500) river(k)%ip_land,river(k)%jp_land
            write(10,'(i4)') 1
            write(10,300) river(k)%ip_land+2,river(k)%jp_land,
     &                    trim(river(k)%string(1:20))
            end if
         end do
         close(10)
         end if

         ! Diagnose rivers take 2 - netcdf file
         print *,'Dumping netcdf diagnostic file'
         if (NF90_CREATE("rivers.nc",NF90_CLOBBER,ncid)
     &       /= NF90_NOERR) then
            print *,'An error occured when opening the netcdf file'
            stop '(obsstats)'
         end if
         ierr=NF90_DEF_DIM(ncid,'idm',idm,idmid)
         ierr=NF90_DEF_DIM(ncid,'jdm',jdm,jdmid)
         ierr=NF90_DEF_DIM(ncid,'rdm',NF90_UNLIMITED,rdimid)
         ierr=NF90_DEF_DIM(ncid,'rivdim',riv_nr,rivdim)
         ierr=NF90_DEF_DIM(ncid,'N_2',2,N_2)

         ierr=NF90_DEF_VAR(ncid,'modlon',NF90_Float,
     &                     (/idmid,jdmid/),varid)
         call handle_err(NF90_ENDDEF(ncid))
         call handle_err(NF90_PUT_VAR(ncid,varid,modlon))
         call handle_err(NF90_REDEF(ncid))

         call handle_err(NF90_DEF_VAR(ncid,'modlat',NF90_Float,
     &                     (/idmid,jdmid/),varid))
         call handle_err(NF90_ENDDEF(ncid))
         call handle_err(NF90_PUT_VAR(ncid,varid,modlat))
         call handle_err(NF90_REDEF(ncid))

         tmp=depths
         where (depths<1 .or. depths > 1e26) tmp=0.
         call handle_err(NF90_DEF_VAR(ncid,'depths',NF90_Float,
     &                     (/idmid,jdmid/),varid))
         call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',
     &                   real(0.,kind=4)))
         call handle_err(NF90_ENDDEF(ncid))
         call handle_err(NF90_PUT_VAR(ncid,varid,tmp))
         call handle_err(NF90_REDEF(ncid))

         fillv=-1e14
         call handle_err(NF90_DEF_VAR(ncid,'riverpos_init',NF90_Float,
     &                     (/rivdim,N_2/),varid))
         call handle_err(NF90_PUT_ATT(ncid,varid,'comment',
     &                   'data file river placement (lon lat)'))
         call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',
     &                   real(fillv,kind=4)))
         call handle_err(NF90_ENDDEF(ncid))
         do k=1,riv_nr
            if (river(k)%active) then
               call handle_err(NF90_PUT_VAR(ncid,varid,
     &            river(k)%lon, start=(/k,1/)))
               call handle_err(NF90_PUT_VAR(ncid,varid,
     &            river(k)%lat, start=(/k,2/)))
            else
               call handle_err(NF90_PUT_VAR(ncid,varid,
     &            fillv,start=(/k,1/)))
               call handle_err(NF90_PUT_VAR(ncid,varid,
     &            fillv,start=(/k,2/)))
            end if
         end do
         call handle_err(NF90_REDEF(ncid))


         fillv=-1e14
         call handle_err(NF90_DEF_VAR(ncid,'riverpos_final',NF90_Float,
     &                     (/rivdim,N_2/),varid))
         call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',
     &                   real(fillv,kind=4)))
         call handle_err(NF90_PUT_ATT(ncid,varid,'comment',
     &                   'repositioned river placement (lon lat)'))
         call handle_err(NF90_ENDDEF(ncid))
         do k=1,riv_nr
            if (river(k)%active) then
               call handle_err(NF90_PUT_VAR(ncid,varid,
     &         modlon(river(k)%ip_land,river(k)%jp_land),start=(/k,1/)))
               call handle_err(NF90_PUT_VAR(ncid,varid,
     &         modlat(river(k)%ip_land,river(k)%jp_land),start=(/k,2/)))
            else
               call handle_err(NF90_PUT_VAR(ncid,varid,
     &            fillv,start=(/k,1/)))
               call handle_err(NF90_PUT_VAR(ncid,varid,
     &            fillv,start=(/k,2/)))
            end if
         end do
         call handle_err(NF90_REDEF(ncid))

         fillv=-1e14
         call handle_err(NF90_DEF_VAR(ncid,'riverflux',NF90_Float,
     &                     (/idmid,jdmid,rdimid/),varid))
         call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',
     &      real(fillv,kind=4)))
         call handle_err(NF90_PUT_ATT(ncid,varid,'units','mm day-3'))
         call handle_err(NF90_ENDDEF(ncid))
         do mo=1,12
            tmp=river_flux(:,:,mo)
            tmp=tmp*86400.*1e3
            where (depths<1 .or. depths > 1e26) tmp=fillv
            call handle_err(NF90_PUT_VAR(ncid,varid,tmp,
     &                                   start=(/1,1,mo/)))
         end do
      end if
      call handle_err(NF90_CLOSE(ncid))




      if (mnproc==1) then
         print *,'rivers calculated ok, see diagnostics in rivers.tec'
      end if

  99  FORMAT(30I4) 
 100  FORMAT(10(1x,e10.4)) 
 300  format('TEXT CS=GRID,   X=',
     &       i4,', Y=',i4,', T="',a,'"')
  400 format('GEOMETRY CS=GRID, X=',f10.2,', Y=',f10.2,
     &       ', T=CIRCLE, FC=BLACK')
  500 format('GEOMETRY CS=GRID, X=',i5,', Y=',i5,
     &       ', T=CIRCLE, C=BLACK, FC=RED')
      end subroutine rivers_to_hycom







!#######################################################################
!##############         FUNCTION RIVERNR()              ################
!############## --------------------------------------- ################
!############## Gives number of rivers Specified in the ################
!############## input file.                             ################
!#######################################################################
! Input:
!  
! Output:
!   rivernr (function) -- Number of rivers in input file
!

      function rivernr()
      use mod_xc
      implicit none
      integer rivernr
      integer, parameter :: rivermax=10000000
      logical active,ex
      integer k
      real dummy, dummy2
      character(len=60) string

      inquire(exist=ex,file=flriver)
      if (.not.ex) then
         if(mnproc==1) write(6,*) 'Can not find file '//flriver
         call xcstop('(rivernr)')
         stop '(rivernr)'
      end if
      open(18,file=flriver,STATUS='OLD')
      rivernr=0
      do k=1,rivermax
         read(18,'(l1,a60)',end=100,err=100)active,string
         read(18,'(f10.1,12f8.5)')dummy
         read(18,'(2f9.2)')dummy2
         if (active) rivernr=rivernr+1
      enddo
 100  close(18)

      if (rivernr == rivermax) then
         if(mnproc==1) write(*,*)'rivernr: rivernr=rivermax'
         if(mnproc==1) write(*,*)'rivernr: all rivers may not be used'
         call xcstop('(rivernr)')
         stop '(rivernr)'
      endif
      end function rivernr




!#######################################################################
!##############        SUBROUTINE READRIVERS            ################
!############## --------------------------------------- ################
!##############      Reads river data from file         ################
!##############                                         ################
!#######################################################################
! Input:
!   -riv_nr   -- Number of rivers in file 
!  
! Output:
!   -river    -- Contains flux ++ of each river
!
      subroutine readrivers(modlon,modlat)
      ! --- ------------------------------------------------------------------
      ! --- Include river data. 3 lines are specifying each river.
      ! --- Line 1.: Active, Name and other info (l1,t3,a60)
      ! --- Line 2.: Annual mean (km^3/year), monthly weights (f10.1,12f8.5)
      ! --- Line 3.: Lat, Lon position , lat,lon (2f9.1)
      ! --- ------------------------------------------------------------------
      use mod_xc
      use mod_confmap
      implicit none
      real, intent(in), dimension(idm,jdm) :: modlon, modlat

      logical riv_active,ex, riv_banned
      character(len=80) riv_string
      real riv_annual
      real riv_flux(12)
      real riv_lat,riv_lon
      real g_latriv,g_lonriv
      integer riv_ip,riv_jp
      real fl,annual
      integer i,j,k,m
      real latnew, lonnew          ! New transformed lat, lon
      real lat1,lat2,lon1,lon2, riv_ipr,riv_jpr
      type(river_data), allocatable :: tmpriver(:)

      if (mnproc==1) write(lp,*) 'read_rivers'


      inquire(file=flriver,exist=ex)
      if (.not.ex) then
         call xcstop('Data/rivers.dat file does not exist')
         stop 'Data/rivers.dat file does not exist'
      end if
      open(18,file=flriver,STATUS='OLD')

      k=1
      do m=1,10000
         read(18,'(l1,t3,a80)',end=100,err=100)riv_active,riv_string
        ! read(18,'(f10.1,12f8.5)')riv_annual,(riv_flux(i),i=1,12)
         read(18,'(f12.5,12f8.5)')riv_annual,(riv_flux(i),i=1,12)
         read(18,'(2f9.2)')riv_lat,riv_lon
         riv_banned=.false.


         if (riv_active) then
            if (mnproc==1) then
               write(lp,*)
               write(lp,'(l1,t3,a,a60,g13.2)')riv_active,riv_string,
     &            ' Annual discharge (m**3)= ',riv_annual
            end if

           ! Consistency check flux must add up to 1!
           if (abs(sum(riv_flux)-1.0) >0.03 ) then
              if (mnproc==1)  then
                 write(lp,*) 'fluxes for river does not add to 1 !!'
                 write(lp,*) 'sum river_flux:',sum(riv_flux)
              end if
              call xcstop('(readrivers)')
              stop '(readrivers)'
           end if

           ! Transform lat, lon to grid points  - pivot point
           call initconfmap(idm,jdm)
           call oldtonew(riv_lat,riv_lon,latnew,lonnew)
           call pivotp (lonnew,latnew,riv_ip,riv_jp)

           


           ! Convert from km^3/y to m^3/s
           riv_annual=riv_annual*1.e09/(365.*24.*3600.) 

           if (mnproc==1) then
              write(lp,*) 'Pivot point for river location       =',
     &             riv_ip,riv_jp
              write(lp,*) 'Latitude longitude for river location=',
     &             riv_lat, riv_lon
           end if

           if (riv_ip<1 .or. riv_ip>idm-1 .or.
     &         riv_jp<1 .or. riv_jp>jdm-1 ) then
              print *,'This river is outside the model domain - skipped'
              riv_active=.false.
              riv_banned=.true.
              cycle
           else

           ! real pivot point (useful for dumping river data to tecplot file)
              call oldtonew(modlat(riv_ip+1,riv_jp+1),
     &                      modlon(riv_ip+1,riv_jp+1),
     &                      lat1,lon1)
              call oldtonew(modlat(riv_ip,riv_jp),modlon(riv_ip,riv_jp),
     &                      lat2,lon2)
              riv_ipr=riv_ip+(lonnew-lon2)/(lon1-lon2)
              riv_jpr=riv_jp+(latnew-lat2)/(lat1-lat2)

           end if

    
          do i=1,12
             ! Monthly flux (m^3/s)
             riv_flux(i)=riv_annual*riv_flux(i)*12.0

             ! Corresponding level change (m/s)
             !riv_flux(i)=riv_flux(i)/darea
          enddo

          river(k)%active=riv_active
          river(k)%banned=riv_banned
          river(k)%string=riv_string
          river(k)%annual=riv_annual
          river(k)%flux(:)=riv_flux(:)
          river(k)%lat=riv_lat
          river(k)%lon=riv_lon
          river(k)%ip=riv_ip
          river(k)%jp=riv_jp
          river(k)%ipr=riv_ipr
          river(k)%jpr=riv_jpr

          k=k+1
        endif
123     continue
        riv_active=.false.
        riv_annual=0.0
        riv_flux(:)=0.0
        riv_lat=0.0
        riv_lon=0.0
        riv_ip=0
        riv_jp=0

      enddo

100   CONTINUE
      close(18)


      allocate(tmpriver(riv_nr)) 
      tmpriver=river
      deallocate(river)
      riv_nr=k-1
      print *,'new riv_nr=',riv_nr
      if (riv_nr>0) then
         allocate(river(riv_nr))
         river(1:riv_nr)=tmpriver(1:riv_nr)
      end if
      deallocate(tmpriver)


      do k=1,riv_nr
         if (mnproc==1) then
            write(lp,'(l1,t3,a60,g13.2)') 
     &         river(k)%active,river(k)%string,river(k)%annual
            write(lp,'(12g10.2,tr3,g10.2)') river(k)%flux(1:12),
     &         sum(river(k)%flux(1:12))
            write(lp,'(2f9.2,2i5)') river(k)%lat,river(k)%lon,
     &         river(k)%ip,river(k)%jp
         end if
      enddo
      end subroutine readrivers


!#######################################################################
!##############        SUBROUTINE ACCRIVERS             ################
!############## --------------------------------------- ################
!############## Calculates vertical speed of sea level  ################
!##############       due to river mass flux            ################
!#######################################################################
! Input:
!   -riv_nr     -- Number of rivers in file 
!   -river      -- River data
!  
! Output:
!   -river_flux -- Monthly river fluxes (m/s)
!
!
! KAL : New: let river flux be weighted by distance from origin

      subroutine accrivers(modlon,modlat,scpx,scpy,depths,river_flux)
      use mod_xc
      implicit none
      real, intent(in), dimension(idm,jdm) :: modlon, modlat, depths,
     &                                        scpx,scpy
      real,             intent(out) ::  river_flux(1:idm,1:jdm,12)

      real,dimension(idm,jdm) :: tmpriver, weight
      real, dimension(idm,jdm) :: landdist
      integer :: k,i,j,im1,ip1,jm1,jp1,mo
      logical :: isopen
      real    :: q, dist,tmpu,tmpv
      real*8  :: wsum

      real :: mindist, norm,maxdist
      integer :: i2,j2,irad,iedge,ind,imod,jmod
      integer, dimension(idm,jdm) :: ip
      logical, dimension(idm,jdm) :: connected
      logical :: cont_loop, anyban
      integer :: mini, minj, numrad
      real, external :: spherdist

      !include 'common_blocks.h'

      river_flux=0.


      where (depths>.1) 
         ip=1
      elsewhere
         ip=0
      endwhere


      ! New scheme
      do k=1,riv_nr

         print *
         print *
         print *,'Doing river ',river(k)%string

         ! Find ocean point bordering land which is closest to the
         ! specified pivot point
         mindist=1e8
         do j=2,jdm-1
         do i=2,idm-1
!KAL        if (ip(i,j)==1 .and. (ip(i-1,j)==0 .or. ip(i+1,j)==0  .or.
!KAL &                            ip(i,j-1)==0 .or. ip(i,j+1)==0)) then

            ! Accept point if it has three or more land neighbours
            if (ip(i,j)==1 .and. sum(ip(i-1:i+1,j-1:j+1))<=6) then

               dist= spherdist(river(k)%lon,river(k)%lat,
     &                         modlon(i,j),modlat(i,j))
               
               if (dist<mindist ) then
                  mindist=dist
                  mini=i
                  minj=j
               end if
            end if
         end do
         end do

         river(k)%ip_land=mini
         river(k)%jp_land=minj

         print *,'River lat lon   (from rivers.dat):',
     &      river(k)%lat, river(k)%lon
         print *,'old pivot point (from rivers.dat):',
     &      river(k)%ip,river(k)%jp
         print *,'new pivot point (land boundary  ):',
     &      mini,minj,depths(mini,minj)
         print *,'lon lat of new pivot point (land boundary  ):',
     &      modlat(mini,minj),modlon(mini,minj)
         river(k)%ip_land=mini
         river(k)%jp_land=minj
         !print *,ip(mini-1,minj)
         !print *,ip(mini+1,minj)
         !print *,ip(mini  ,minj-1)
         !print *,ip(mini  ,minj+1)

         ! Find where this point is "land connected" to the other 
         ! land points
         irad=1
         connected=.false.
         if (ip(mini-1,minj)==0) connected(mini-1,minj)=.true.
         if (ip(mini+1,minj)==0) connected(mini+1,minj)=.true.
         if (ip(mini,minj-1)==0) connected(mini,minj-1)=.true.
         if (ip(mini,minj+1)==0) connected(mini,minj+1)=.true.
         maxdist=0.
         cont_loop=.true.
         do while (cont_loop)

            do j=max(2,minj-irad),min(jdm-1,minj+irad)
            do i=max(2,mini-irad),min(idm-1,mini+irad)

               ! Set connectedness if this is land
               if (ip(i,j)==0) then
                  connected(i,j) = connected(i,j)  .or. connected(i-1,j)
                  connected(i,j) = connected(i,j)  .or. connected(i+1,j)
                  connected(i,j) = connected(i,j)  .or. connected(i,j-1)
                  connected(i,j) = connected(i,j)  .or. connected(i,j+1)
               endif

               ! In any case, calculate mindist for this irad iteration
               if (i==2 .or.j==2 .or. i==idm-1 .or. j==jdm-1) then

                  ! Border points do not count towards max dist calc
                  maxdist=maxdist ! doh

               else

                  !Maximal distance for other points
                  maxdist= max(maxdist,
     &                 spherdist(modlon(mini,minj),modlat(mini,minj),
     &                           modlon(   i,   j),modlat(   i,   j)))
               end if


            end do
            end do


            ! Conditions for proceeding while loop
            cont_loop=maxdist<max(rradius,alongshoreradius)
            irad=irad+1
!KAL        print *,'Radius from pivot point, max distance:',
!KAL &         irad,maxdist
         end do
         print *,'Connectedness calculated'
         print *,'Number of  land points connected to river:',
     &      count(connected)

         ! Keep irad!
         numrad=irad-1


         ! Find distance from all ocean points within radius
         ! to connected land points

         landdist=2*max(rradius,alongshoreradius)
         do j=max(1,minj-numrad),min(jdm,minj+numrad)
         do i=max(1,mini-numrad),min(idm,mini+numrad)

            if (depths(i,j)>.1) then

               irad=1
               maxdist=0.
               mindist=1e8
               do while (maxdist<max(rradius,alongshoreradius) .and. 
     &                   irad<numrad)

                  do j2=max(1,j-irad),min(jdm,j+irad)
                  do i2=max(1,i-irad),min(idm,i+irad)
   
                     if(connected(i2,j2)) then
                        maxdist= max(maxdist,
     &                           spherdist(modlon(i2,j2),modlat(i2,j2),
     &                                     modlon(i ,j ),modlat(i ,j )))
                        mindist= min(mindist,
     &                           spherdist(modlon(i2,j2),modlat(i2,j2),
     &                                     modlon(i ,j ),modlat(i ,j )))
                     end if
                  end do
                  end do
                  irad=irad+1
                  !print *,i2,j2,irad,mindist
               end do

            ! Atthis stage the routine failed to find a point or
            ! we have a min distance  < rradius and alongshore
            landdist(i,j)=mindist
            end if
            !print *,i,j,landdist(i,j)

         end do
         end do
         landdist(mini,minj)=0. ! Ahemm
         print *,'Ocean distances from connected land calculated'

         tmpriver=0.
         weight=0.
         mindist=1e9
         do j=1,jdm
         do i=1,idm
         if (depths(i,j)>.1) then

            ! Get distance from this point to river origin
            dist= spherdist(modlon(mini,minj),modlat(mini,minj),
     &                      modlon(   i,   j),modlat(   i,   j))
         
            ! weight according to  alongshore radius
            if (dist<alongshoreradius) then
               weight(i,j) = exp(-dist/alongshoreradius)
            else
               weight(i,j) = 0.
            end if

            ! the following puts more weight on points close to land
            ! (closer to land than river-radius)
            ! Only connected land points are considered (see above)
            if (landdist(i,j)<rradius.and.ip(i,j)==1.and.
     &          dist<alongshoreradius) then
               weight(i,j)=weight(i,j)*
     &            exp(-max(1.,2*(landdist(i,j)/rradius)))
            else
               weight(i,j) = 0.
            end if
               
 
            ! Set temporary river flux in this grid cell  - NB
            ! weights will generally not sum up to one (see  below)
            ! and actual flux not used yet
            !tmpriver(i,j) = ip(i,j)*weight(i,j)/(scpx(i,j)*scpy(i,j))
            tmpriver(i,j) = ip(i,j)*weight(i,j)

         end if
         enddo
         enddo
         print *,'Alongshore and radius weights calculated'

         ! Make sure tmpriver equals volume flux
         wsum=sum(tmpriver,ip==1)
         !print *,wsum

         ! wsum is the sum of the weights
         if (wsum<.1) then 
            if (mnproc==1) then
               write(lp,*) 'No discharge for river no', k,
     &                     trim(river(k)%string)
               write(lp,*) 'Land distance for river pivot point is :',
     &                     landdist(mini,minj)
               write(lp,*) 'Lat/lon        for river pivot point is :',
     &                     modlat(mini,minj),modlon(mini,minj)
               write(lp,*) 'Land distance treshold is ',
     &                     rradius
            end if
            call xcstop('(accrivers)')
            stop '(accrivers)'
         end if

         ! This correction mean that the weight used in loop above sum up to
         ! one
         weight=tmpriver/wsum
         wsum=sum(weight,ip==1)
         !print *,wsum

         do mo=1,12
!KAL -      q should be volume now, after changes in readrivers
            ! Interpolated sea level change due to river flux (m/s)
            q=river(k)%flux(mo)

            tmpriver=q*weight/(scpx*scpy)

            ! Add river contribution to final field
            river_flux(:,:,mo)=river_flux(:,:,mo)+tmpriver

            ! Double check - this should be equal to the volume flux from the river
            tmpriver=tmpriver*scpx*scpy
            !call xcsum(wsum,tmpriver,ip)
            wsum=sum(tmpriver,ip==1)

            if (abs(wsum-river(k)%flux(mo))>wsum*1e-4) then
               print *,'Calculated river flux does not match input'
               print *,wsum
               print *,river(k)%flux(mo)
               call xcstop('(accrivers)')
               stop '(accrivers)'
            end if

            ! Print river diagnostics
            if (mo==1) then
               print '(a,f10.2,a,f10.2,a)', 
     &            '---- Annual average flux ---- ',river(k)%annual,     
     &                 ' m^3 s^-1, or ', 
     &                 river(k)%annual*365*86400*1e-9, 
     &                 ' km^3 y^-1'
               print *,'----------------------------------------------'
            end if
            print '(a,i2.2,a,f10.2,a,f10.2,a)','Month ', mo,
     &         ' Monthly average flux ',wsum,     
     &         ' m^3 s^-1, or ',
     &         wsum*365*86400*1e-9, 
     &         ' km^3 y^-1'


         enddo ! month

         ! Diagnostic output
      enddo ! rivers

      print *
      print *

      anyban=.false.
      do k=1,riv_nr
      if (river(k)%banned) then
         if (.not. anyban) then
            Print '(a)','***Warning - The following rivers were banned 
     &      (outside of domain)'
            anyban=.true.
         end if
         print '(a,2f10.2,a)','In position ',river(k)%lon,river(k)%lat,
     &      river(k)%string
      end if
      end do

        

      print *
      print *

      end subroutine accrivers

 
      end module mod_rivers
 
