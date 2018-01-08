      module mod_biorivers

      private
      real   :: rradius = 80000. ! Radius in meters
      real   :: alongshoreradius = 200000. ! Radius in meters

      type bioriver_data
        logical banned ! 
        logical active
        character(len=60) string
        real annual
        real annual_nit
        real annual_pho
        real annual_sil
        real flux(12)
        real lat,lon
        integer ip,jp
        real nitflux(12)
        real phoflux(12)
        real silflux(12)
        integer ip_land,jp_land
        real ipr,jpr
        real discharge_area
        integer ip_baro,jp_baro
      end type bioriver_data
      type(bioriver_data), save, allocatable :: river(:)

!#if defined(RIVER_PLUME)
!      integer, allocatable, save :: river_ip(:)
!      integer, allocatable, save :: river_jp(:)
!      real,    allocatable, save :: river_weight(:)
!      logical, allocatable, save :: river_direction(:,:)
!#endif

      integer, save :: riv_nr=0

      character(len=*) , parameter :: flriver = './Data/biorivers.dat'

      public :: biorivers_to_hycom, rradius, alongshoreradius

      contains




!#######################################################################
!##############       SUBROUTINE BIORIVERS_TO_HYCOM     ################
!############## --------------------------------------- ################
!############## Only routine visible from the outside   ################
!############## called from hycom to get riverine flux. ################
!#######################################################################
      subroutine biorivers_to_hycom(modlon,modlat,scpx,scpy,depths)
      use mod_xc
      use mod_za
      use netcdf
      use m_handle_err
      implicit none
      real, intent(in), dimension(idm,jdm) :: modlon, modlat, depths,
     &   scpx,scpy

      integer :: lgth, mo
      real:: river_flux(1:idm,1:jdm,12)
      real:: nit_flux(1:idm,1:jdm,12)
      real:: pho_flux(1:idm,1:jdm,12)
      real:: sil_flux(1:idm,1:jdm,12)
      real :: hmin,hmax
      real, dimension(idm,jdm) :: tmp
      integer :: i,j,k,new_riv_nr
      integer, dimension(idm,jdm) :: ip
      logical :: ldiagriver=.true.

      integer :: varid, idmid, jdmid, rdimid,ncid, ierr, rivdim,
     &   N_2
      real    :: fillv

      ! Get numbers of rivers 
      riv_nr=biorivernr()
      write(lp,*) 'Number of biorivers=',riv_nr
      if (riv_nr<=0) then
         print *,'Error reading rivers or no rivers in rivers.dat'
         call exit(1)
      end if
      allocate(river(riv_nr))

      ! Read rivers
      write(lp,*) 'calling readbiorivers'
      call readbiorivers(modlon,modlat)

      ! readrivers sets riv_nr
      if (riv_nr<=0) then
         print *,'No rivers in domain - no forcing produced'
         call exit(1)
      end if


      ! Accumulate rivers
      write(lp,*) 'calling accbiorivers'
      river_flux=0.
      nit_flux=0.
      pho_flux=0.
      sil_flux=0.
      call accbiorivers(modlon,modlat,scpx,scpy,depths,river_flux,
     &                    nit_flux,pho_flux,sil_flux)
      ! Dump rivers to hycom-style files
      !lgth = len_trim(flnmfor)
      print *,'Dumping hycom forcing fields'
      if (mnproc==1) then
         open (unit=87,
     &   file='forcing.rivnitr.b', 
     &         status='replace', action='write')
         write(87,'(a)') 'River nitrate fluxes '
         write(87,'(a)') ''
         write(87,'(a)') ''
         write(87,'(a)') ''
         write(87,'(a,2i5)') 'i/jdm = ',idm,jdm
      end if

      ! Dump rivers to hycom-style files
      !lgth = len_trim(flnmfor)
      print *,'Dumping hycom forcing fields'
      if (mnproc==1) then
         open (unit=88,
     &   file='forcing.rivphos.b', 
     &         status='replace', action='write')
         write(88,'(a)') 'River phosphate fluxes '
         write(88,'(a)') ''
         write(88,'(a)') ''
         write(88,'(a)') ''
         write(88,'(a,2i5)') 'i/jdm = ',idm,jdm
      end if

      ! Dump rivers to hycom-style files
      !lgth = len_trim(flnmfor)
      print *,'Dumping hycom forcing fields'
      if (mnproc==1) then
         open (unit=89,
     &   file='forcing.rivsili.b', 
     &         status='replace', action='write')
         write(89,'(a)') 'River silicate fluxes '
         write(89,'(a)') ''
         write(89,'(a)') ''
         write(89,'(a)') ''
         write(89,'(a,2i5)') 'i/jdm = ',idm,jdm
      end if


      where (depths>.1) 
         ip=1
      elsewhere
         ip=0
      endwhere

      call zaiopf('forcing.rivnitr.a', 'replace', 87) 
      do mo=1,12
         call zaiowr(nit_flux(:,:,mo),ip,.false.,
     &               hmin,hmax,87,.true.)
         if(mnproc==1) 
     &     write(87,'(" river nitrate: month,range = ",i2.2,2e16.8)') 
     &     mo,hmin,hmax
      end do
      if (mnproc==1) close(87)
      call zaiocl(87)

      call zaiopf('forcing.rivphos.a', 'replace', 88) 
      do mo=1,12
         call zaiowr(pho_flux(:,:,mo),ip,.false.,
     &               hmin,hmax,88,.true.)
         if(mnproc==1) 
     &     write(88,'(" river phosphate: month,range = ",i2.2,2e16.8)') 
     &     mo,hmin,hmax
      end do
      if (mnproc==1) close(88)
      call zaiocl(88)

      call zaiopf('forcing.rivsili.a', 'replace', 89) 
      do mo=1,12
         call zaiowr(sil_flux(:,:,mo),ip,.false.,
     &               hmin,hmax,89,.true.)
         if(mnproc==1) 
     &     write(89,'(" river silicate: month,range = ",i2.2,2e16.8)') 
     &     mo,hmin,hmax
      end do
      if (mnproc==1) close(89)
      call zaiocl(89)

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
      end subroutine biorivers_to_hycom







!#######################################################################
!##############         FUNCTION BIORIVERNR()           ################
!############## --------------------------------------- ################
!############## Gives number of rivers Specified in the ################
!############## input file.                             ################
!#######################################################################
! Input:
!  
! Output:
!   biorivernr (function) -- Number of rivers (with nutrients) in input file
!

      function biorivernr()
      use mod_xc
      implicit none
      integer biorivernr
      integer, parameter :: biorivermax=10000000
      logical active,ex
      integer k
      real dummy, dummy2
      character(len=60) string

      inquire(exist=ex,file=flriver)
      if (.not.ex) then
         if(mnproc==1) write(6,*) 'Can not find file '//flriver
         call xcstop('(biorivernr)')
         stop '(biorivernr)'
      end if
      open(18,file=flriver,STATUS='OLD')
      biorivernr=0
      do k=1,biorivermax
         read(18,'(l1,a60)',end=100,err=100)active,string
         read(18,'(f10.1,12f8.5)')dummy
         read(18,'(f10.1,12f8.5)')dummy
         read(18,'(f10.1,12f8.5)')dummy
         read(18,'(f10.1,12f8.5)')dummy
         read(18,'(2f9.2)')dummy2
         if (active) biorivernr=biorivernr+1
      enddo
 100  close(18)

      if (biorivernr == biorivermax) then
         if(mnproc==1) write(*,*)'biorivernr: biorivernr=rivermax'
         if(mnproc==1) write(*,*)'biorivernr: all rivers may not 
     &                             be used'
         call xcstop('(biorivernr)')
         stop '(biorivernr)'
      endif
      end function biorivernr




!#######################################################################
!##############        SUBROUTINE READBIORIVERS         ################
!############## --------------------------------------- ################
!##############      Reads bioriver data from file      ################
!##############                                         ################
!#######################################################################
! Input:
!   -riv_nr   -- Number of rivers in file 
!  
! Output:
!   -river    -- Contains flux ++ of each river
!
      subroutine readbiorivers(modlon,modlat)
      ! --- ------------------------------------------------------------------
      ! --- Include river data. 6 lines are specifying each river.
      ! --- Line 1.: Active, Name and other info (l1,t3,a60)
      ! --- Line 2.: Annual mean (km^3/year), monthly weights (f10.1,12f8.5)
      ! --- Line 3.: Annual mean nitrate   (mg/m^3y), monthly weights (f10.1,12f8.5)
      ! --- Line 4.: Annual mean phosphate (mg/m^3y), monthly weights (f10.1,12f8.5)
      ! --- Line 5.: Annual mean silicate  (mg/m^3y), monthly weights (f10.1,12f8.5)
      ! --- Line 6.: Lat, Lon position , lat,lon (2f9.1)
      ! --- ------------------------------------------------------------------
      use mod_xc
      use mod_confmap
      implicit none
      real, intent(in), dimension(idm,jdm) :: modlon, modlat

      logical riv_active,ex, riv_banned
      character(len=60) riv_string
      real riv_annual
      real riv_flux(12)
      real riv_lat,riv_lon
      real g_latriv,g_lonriv
      real nit_annual
      real nit_flux(12)
      real pho_annual
      real pho_flux(12)
      real sil_annual
      real sil_flux(12)
      integer riv_ip,riv_jp
      real fl,annual
      integer i,j,k,m
      real latnew, lonnew          ! New transformed lat, lon
      real lat1,lat2,lon1,lon2, riv_ipr,riv_jpr
      type(bioriver_data), allocatable :: tmpriver(:)

      if (mnproc==1) write(lp,*) 'read_rivers'


      inquire(file=flriver,exist=ex)
      if (.not.ex) then
         call xcstop('Data/biorivers.dat file does not exist')
         stop 'Data/biorivers.dat file does not exist'
      end if
      open(18,file=flriver,STATUS='OLD')

      k=1
      do m=1,10000
         read(18,'(l1,t3,a60)',end=100,err=100)riv_active,riv_string
         read(18,'(f10.1,12f8.5)')riv_annual,(riv_flux(i),i=1,12)
         read(18,'(f10.1,12f8.5)')nit_annual,(nit_flux(i),i=1,12)
         read(18,'(f10.1,12f8.5)')pho_annual,(pho_flux(i),i=1,12)
         read(18,'(f10.1,12f8.5)')sil_annual,(sil_flux(i),i=1,12)
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
              call xcstop('(readbiorivers)')
              stop '(readbiorivers)'
           end if

          ! Consistency check flux must add up to 1!
           if (abs(sum(nit_flux)-1.0) >0.03 ) then
              if (mnproc==1)  then
                 write(lp,*) 'nitrate fluxes for river 
     &                        does not add to 1 !!'
                 write(lp,*) 'sum rivnit_flux:',sum(nit_flux)
              end if
              call xcstop('(readbiorivers)')
              stop '(readbiorivers)'
           end if

         ! Consistency check flux must add up to 1!
           if (abs(sum(pho_flux)-1.0) >0.03 ) then
              if (mnproc==1)  then
                 write(lp,*) 'phosphate fluxes for river 
     &                        does not add to 1 !!'
                 write(lp,*) 'sum rivpho_flux:',sum(pho_flux)
              end if
              call xcstop('(readbiorivers)')
              stop '(readbiorivers)'
           end if

         ! Consistency check flux must add up to 1!
           if (abs(sum(sil_flux)-1.0) >0.03 ) then
              if (mnproc==1)  then
                 write(lp,*) 'silicate fluxes for river does 
     &                        not add to 1 !!'
                 write(lp,*) 'sum rivsil_flux:',sum(sil_flux)
              end if
              call xcstop('(readbiorivers)')
              stop '(readbiorivers)'
           end if

           ! Transform lat, lon to grid points  - pivot point
           call initconfmap(idm,jdm)
           call oldtonew(riv_lat,riv_lon,latnew,lonnew)
           call pivotp (lonnew,latnew,riv_ip,riv_jp)

           


           ! Convert from km^3/y to m^3/s
           riv_annual=riv_annual*1.e09/(365.*24.*3600.) 
          ! Convert from mg/m3y to  mg/m3s
           nit_annual=nit_annual/(365.*24.*3600.)
           pho_annual=pho_annual/(365.*24.*3600.)
           sil_annual=sil_annual/(365.*24.*3600.)

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
             ! Nutrient flux (mg/s)
             nit_flux(i)=nit_annual*
     &                        nit_flux(i)*12.0
             pho_flux(i)=pho_annual*
     &                        pho_flux(i)*12.0
             sil_flux(i)=sil_annual*
     &                        sil_flux(i)*12.0

             ! Corresponding level change (m/s)
             !riv_flux(i)=riv_flux(i)/darea
          enddo

          river(k)%active=riv_active
          river(k)%banned=riv_banned
          river(k)%string=riv_string
          river(k)%annual=riv_annual
          river(k)%annual_nit=nit_annual
          river(k)%annual_pho=pho_annual
          river(k)%annual_sil=sil_annual
          river(k)%flux(:)=riv_flux(:)
          river(k)%lat=riv_lat
          river(k)%lon=riv_lon
          river(k)%ip=riv_ip
          river(k)%jp=riv_jp
          river(k)%ipr=riv_ipr
          river(k)%jpr=riv_jpr
          river(k)%nitflux(:)=nit_flux(:)
          river(k)%phoflux(:)=pho_flux(:)
          river(k)%silflux(:)=sil_flux(:)

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
        nit_annual=0.0
        nit_flux(:)=0.0
        pho_annual=0.0
        pho_flux(:)=0.0
        sil_annual=0.0
        sil_flux(:)=0.0

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
            write(lp,'(12g10.2,tr3,g10.2)') river(k)%nitflux(1:12),
     &         sum(river(k)%nitflux(1:12))
            write(lp,'(12g10.2,tr3,g10.2)') river(k)%phoflux(1:12),
     &         sum(river(k)%phoflux(1:12))
            write(lp,'(12g10.2,tr3,g10.2)') river(k)%silflux(1:12),
     &         sum(river(k)%silflux(1:12))
         end if
      enddo
      end subroutine readbiorivers


!#######################################################################
!##############        SUBROUTINE ACCBIORIVERS          ################
!############## --------------------------------------- ################
!############## Calculates vertical added nutrients     ################
!##############       due to river mass flux            ################
!#######################################################################
! Input:
!   -riv_nr     -- Number of rivers in file 
!   -river      -- River data
!  
! Output:
!   -river_flux -- Monthly river fluxes (m/s)
!   -nit_flux   -- Monthly river fluxes (mg/m^3s)
!   -pho_flux   -- Monthly river fluxes (mg/m^3s)
!   -sil_flux   -- Monthly river fluxes (mg/m^3s)
!
!
! KAL : New: let river flux be weighted by distance from origin

      subroutine accbiorivers(modlon,modlat,scpx,scpy,depths,river_flux, 
     &                    nit_flux,pho_flux,sil_flux)
      use mod_xc
      implicit none
      real, intent(in), dimension(idm,jdm) :: modlon, modlat, depths,
     &                                        scpx,scpy
      real,             intent(out) ::  river_flux(1:idm,1:jdm,12)
      real,             intent(out) ::  nit_flux(1:idm,1:jdm,12)
      real,             intent(out) ::  pho_flux(1:idm,1:jdm,12)
      real,             intent(out) ::  sil_flux(1:idm,1:jdm,12)

      real,dimension(idm,jdm) :: tmpriver, weight
      real,dimension(idm,jdm) :: tmpriver_nit
      real,dimension(idm,jdm) :: tmpriver_pho
      real,dimension(idm,jdm) :: tmpriver_sil
      real, dimension(idm,jdm) :: landdist
      integer :: k,i,j,im1,ip1,jm1,jp1,mo
      logical :: isopen
      real    :: q, dist,tmpu,tmpv
      real    :: q_nit,q_pho,q_sil  
      real*8  :: wsum, wsum_nit, wsum_pho, wsum_sil

      real :: mindist, norm,maxdist
      integer :: i2,j2,irad,iedge,ind,imod,jmod
      integer, dimension(idm,jdm) :: ip
      logical, dimension(idm,jdm) :: connected
      logical :: cont_loop, anyban
      integer :: mini, minj, numrad
      real, external :: spherdist

      !include 'common_blocks.h'

      river_flux=0.
      nit_flux=0.
      pho_flux=0.
      sil_flux=0.

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
         tmpriver_nit=0.; 
         tmpriver_pho=0.; 
         tmpriver_sil=0.; 
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
            call xcstop('(accbiorivers)')
            stop '(accbiorivers)'
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
            q_nit=river(k)%nitflux(mo)
            q_pho=river(k)%phoflux(mo)
            q_sil=river(k)%silflux(mo)

            tmpriver=q*weight/(scpx*scpy)
            tmpriver_nit=q_nit*weight!/(scpx*scpy)
            tmpriver_pho=q_pho*weight!/(scpx*scpy)
            tmpriver_sil=q_sil*weight!/(scpx*scpy)

            ! Add river contribution to final field
            river_flux(:,:,mo)=river_flux(:,:,mo)+tmpriver
            nit_flux(:,:,mo)=nit_flux(:,:,mo)+tmpriver_nit
            pho_flux(:,:,mo)=pho_flux(:,:,mo)+tmpriver_pho
            sil_flux(:,:,mo)=sil_flux(:,:,mo)+tmpriver_sil

            ! Double check - this should be equal to the volume flux from the river
            tmpriver=tmpriver*scpx*scpy
            tmpriver_nit=tmpriver_nit!*scpx*scpy
            tmpriver_pho=tmpriver_pho!*scpx*scpy
            tmpriver_sil=tmpriver_sil!*scpx*scpy

            !call xcsum(wsum,tmpriver,ip)
            wsum=sum(tmpriver,ip==1)
            wsum_nit=sum(tmpriver_nit,ip==1)
            wsum_pho=sum(tmpriver_pho,ip==1)
            wsum_sil=sum(tmpriver_sil,ip==1)

            if (abs(wsum-river(k)%flux(mo))>wsum*1e-4) then
               print *,'Calculated river flux does not match input'
               print *,wsum
               print *,river(k)%flux(mo)
               call xcstop('(accbiorivers)')
               stop '(accbiorivers)'
            end if

            if (abs(wsum_nit-river(k)%nitflux(mo))>wsum_nit*1e-4) then
               print *,'Calculated river flux does not match input'
               print *,wsum_nit
               print *,river(k)%nitflux(mo)
               call xcstop('(accbiorivers)')
               stop '(accbiorivers)'
            end if

            if (abs(wsum_pho-river(k)%phoflux(mo))>wsum_pho*1e-4) then
               print *,'Calculated river flux does not match input'
               print *,wsum_pho
               print *,river(k)%phoflux(mo)
               call xcstop('(accbiorivers)')
               stop '(accbiorivers)'
            end if

            if (abs(wsum_sil-river(k)%silflux(mo))>wsum_sil*1e-4) then
               print *,'Calculated river flux does not match input'
               print *,wsum_sil
               print *,river(k)%silflux(mo)
               call xcstop('(accbiorivers)')
               stop '(accbiorivers)'
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

            if (mo==1) then
               print '(a,f10.2,a,f10.2,a)', 
     &     '---- Annual average nitrate flux ---- ',river(k)%annual_nit,     
     &     ' mg/m^3 s^-1, or ', 
     &     river(k)%annual_nit*365*86400, 
     &                 ' mg/m^3 y^-1'
               print *,'----------------------------------------------'
            end if
            print '(a,i2.2,a,f10.2,a,f10.2,a)','Month ', mo,
     &         ' Monthly average flux ',wsum_nit,     
     &         ' mg/m^3 s^-1, or ',
     &         wsum_nit*365*86400, 
     &         ' mg/m^3 y^-1'

             if (mo==1) then
               print '(a,f10.2,a,f10.2,a)', 
     &     '---- Annual average phosph. flux ---- ',river(k)%annual_pho,     
     &     ' mg/m^3 s^-1, or ', 
     &     river(k)%annual_pho*365*86400, 
     &                 ' mg/m^3 y^-1'
               print *,'----------------------------------------------'
            end if
            print '(a,i2.2,a,f10.2,a,f10.2,a)','Month ', mo,
     &         ' Monthly average flux ',wsum_pho,     
     &         ' mg/m^3 s^-1, or ',
     &         wsum_pho*365*86400, 
     &         ' mg/m^3 y^-1'

             if (mo==1) then
               print '(a,f10.2,a,f10.2,a)', 
     &     '---- Annual average silica. flux ---- ',river(k)%annual_sil,     
     &     ' mg/m^3 s^-1, or ', 
     &     river(k)%annual_sil*365*86400, 
     &                 ' mg/m^3 y^-1'
               print *,'----------------------------------------------'
            end if
            print '(a,i2.2,a,f10.2,a,f10.2,a)','Month ', mo,
     &         ' Monthly average flux ',wsum_sil,     
     &         ' mg/m^3 s^-1, or ',
     &         wsum_sil*365*86400, 
     &         ' mg/m^3 y^-1'


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

      end subroutine accbiorivers

 
      end module mod_biorivers
 
