! KAL -- This routine creates "NERSC" nesting condition from hycom-files
! KAL -- For now it works with weekly average files, but in the future
! KAL -- it should support daily and restart files as well
! KAL
! KAL -- Program is designed to have a low memory footprint, this
! KAL -- means it will run slightly slower, but the generation of nesting
! KAL -- conditions by this method isn't done that often, so this should be ok.
! KAL -- Another sacrifice is in "readibility" of the program -- I will
! KAL -- try to document it better later on
! KAL 
! KAL -- Some things remain, but it should create reasonable nesting fields now.
! KAL -- Usable/Inital version: 30.10.2006 - Knut Lisæter

! KAL -- Cleaned up for new nesting setup in 2.2. This routine remains spaghetti-like...

      program nest_offline
      use mod_xc
      use mod_za
      use mod_grid
      use mod_year_info
      use mod_hycomfile_io
      use m_filenesting
      use mod_confmap
      use m_nearestpoint
      use m_rotate2
      use m_parse_blkdat
      use m_bavg_flds_from_ave
      implicit none

      character(len= *),parameter :: fnestloc='Nest/nestloc.uf'
      character(len=80)           :: fndepths

      integer :: nvar     ! Number of 2D flds for loop (velocity treated as one)
      integer :: nvar_tot ! Number of 2D flds for offset, velocity treated as two


      character(len=3 ) :: char3
      character(len=3 ) :: char3_2
      character(len=5 ) :: char5
      character(len=15) :: char15
      character(len=7 ) :: tag7
      character(len=80) :: cfile,nestingfile,tmparg

      ! Arrays for holding data on the "global" model grid
      real*8, dimension(:,:), allocatable ::  io1, io2
      real  , dimension(:,:), allocatable ::   &
         fld1, fld2, dp  , intf, intfold
      real  , dimension(:,:), allocatable ::   &
         ssh, ubavg, pbavg, vbavg
      logical, allocatable, dimension(:,:) ::  &
         gmsk

      ! Dimensions of inner grid
      integer :: inner_idm, inner_jdm, inner_inest,  &
         inner_iidm, inner_jjdm

      ! Arrays holding data on inner grid (where we need nesting
      ! conditions)
      real*8, allocatable, dimension(:,:) ::  &
         inner_ior8
      real*4, allocatable, dimension(:,:) ::  &
         inner_ior4
      real, allocatable, dimension(:,:) ::  &
         inner_lon, inner_lat, inner_tmp1, inner_tmp2, &
         inner_depths, inner_a1,inner_a2,inner_a3,inner_a4, &
         old_inner_tmp1, old_inner_tmp2
      integer, allocatable, dimension(:,:) ::  &
         inner_ipiv, inner_jpiv

      type(year_info)   :: rti,rtd

      real :: dpsum,basum
      real :: ba1,ba2,ba3,ba4,lon,lat,lon_n,lat_n, realvar
      integer :: i,j,k,ios,ipiv,jpiv,irec,klevel, irec_offset,ivar, &
         kdm, intvar,ipib,jpib
      integer :: igrace,jgrace
      logical :: ex, isvelocity,ass
      integer :: nrmem

      integer :: maxdpind(1)
      integer :: hycfile
      integer :: iyear, imonth, iweek, iday
      character(len=80) :: fbase,rstfile, ftype
      logical :: lperiodic,useold,nestice, samegrid

      integer :: indxrst

      integer :: ngrace=3
      real :: amin, amax
#if defined (IARGC)
      integer*4, external :: iargc
#endif
      type(hycomfile) :: hfile


      ! Program creates nesting files from .[ab] files (DAILY for now)
      !
      ! Needs:
      ! blkdat.input         - blkdat from global model
      ! grid.info            - grid info from global model
      ! regional.grid.[ab]   - regional grid   from global model
      ! regional.depth.[ab]  - regional depths from global model
      ! A bunch of daily average files (for now) from the global model
      ! Nest/nestloc.uf         - positions of local model
      ! Nest/ndepthsXXXXXXXX.uf - depths of local model
      ! nestoffline.in       - Contains files to nest, and their type
      ! botprop.in       - Contains name of file with "bottom properties".
      !
      ! Spits out:
      ! nesting files ;-) in "Nest/" directory

      ! Note that this routine dumps nesting conditions with the global
      ! model vertical discretization. If you use a different
      ! discretization in the local model you should run
      ! nest_layer_remap as well - this procedure is still  "experimental",
      ! however



      if (iargc()==0) then
         samegrid=.false.
      else if (iargc()==1) then
         call getarg(1,tmparg)
         samegrid = trim(tmparg) == 'samegrid'
         print *,'Same grid!'
      else
         print *,'Only one (optional) argument, "samegrid"'
         call exit(1)
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set up global  grids ...  There is redundant information here, blkdat.input
! contains idm,jdm,kdm. regional.depths.b contains idm jdm

      ! 1 - start with xcspmd, which accesses regional.grid.b
      ! Initialize Arrai IO
      call xcspmd()
      call zaiost()
      call get_grid()

      ! 2 - get kdm from blkdat.input - cross check idm/kdm with existing
      ! Search for idm/jdm/kdm
      print *,'parsing blkdat (consistency check  + kdm)'
      
      call parse_blkdat('idm   ','integer',realvar,intvar)
      if (intvar/=idm) then
         print *,'idm: blkdat.input and regional.depths.b mismatch'
         stop '(nest_offline)'
      end if
      call parse_blkdat('jdm   ','integer',realvar,intvar)
      if (intvar/=jdm) then
         print *,'jdm: blkdat.input and regional.depths.b mismatch'
         stop '(nest_offline)'
      end if
      call parse_blkdat('kdm   ','integer',realvar,intvar)
      kdm=intvar


      call initconfmap(idm,jdm)

      


      ! This does not incorporate ice (yet)
      nvar=4*kdm+2     ! Number of 2D flds for loop (velocity treated as one)
      nvar_tot=5*kdm+3 ! Number of 2D flds for offset, velocity treated as two
      nestice=.true.
      if (nestice) then
         nvar    = nvar+3
         nvar_tot=nvar_tot+4 ! Number of 2D flds for offset, velocity treated as two
      end if



      ! Allocate variables for global grid
      allocate(io1   (idm,jdm))
      allocate(io2   (idm,jdm))
      allocate(fld1  (idm,jdm))
      allocate(fld2  (idm,jdm))
      allocate(dp    (idm,jdm))
      allocate(intf  (idm,jdm))
      allocate(intfold  (idm,jdm))
      allocate(gmsk  (idm,jdm))

      allocate(ssh   (idm,jdm))
      allocate(ubavg (idm,jdm))
      allocate(pbavg (idm,jdm))
      allocate(vbavg (idm,jdm))



      where (depths>1e12) depths=0.
      gmsk=depths>.1

      ! Flag indicates if grid is periodic in 1st dimension
      if (any(depths(1,:)>1.) .or. any(depths(idm,:)>1.)) then 
         lperiodic=.true.
      else
         lperiodic=.false.
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read geo positions and depth from local grid
      inquire(file=fnestloc,exist=ex)
      if (.not.ex) then 
         print *,fnestloc//' does not exist!'
         stop '(nest_offline)'
      end if


      ! Read nesting positions for the internal grid
      open(10,file=fnestloc,form='unformatted')

      ! Read grid dimensions
      read(10)inner_idm,inner_jdm,inner_inest
      inner_iidm=inner_idm-inner_inest+1
      inner_jjdm=inner_jdm-inner_inest+1

      write(lp,'(a,5i5)')  &
         fnestloc//'  has dimensions: ', inner_idm,inner_jdm, &
         inner_iidm,inner_jjdm,inner_inest

      ! Allocate grid for inner model
      allocate(inner_lon(inner_idm,inner_jdm))
      allocate(inner_lat(inner_idm,inner_jdm))
      allocate(inner_ior8 (inner_idm,inner_jdm))
      allocate(inner_ior4(inner_idm,inner_jdm))

      ! Read grid
      read(10) inner_ior8 ; inner_lon=inner_ior8
      read(10) inner_ior8 ; inner_lat=inner_ior8
      close(10)

      write(lp,'(a)')   fnestloc//'  is read'

      ! Allocate temporary arrays for interpolating to inner grid
      ! + depth matrix
      allocate(inner_tmp1  (inner_idm,inner_jdm))
      allocate(inner_tmp2  (inner_idm,inner_jdm))
      allocate(old_inner_tmp1  (inner_idm,inner_jdm))
      allocate(old_inner_tmp2  (inner_idm,inner_jdm))
      allocate(inner_depths(inner_idm,inner_jdm))


      ! Read the depth matrix of the local grid 
! --- tag7 has 11 chars to accomodate for huge grids in future
      if (inner_idm>999 .or. inner_jdm > 999) then
         write(tag7,'(i5.5,a,i5.5)')inner_idm,'x',inner_jdm
      else
         write(tag7,'(i3.3,a,i3.3)')inner_idm,'x',inner_jdm
      end if

      fndepths='Nest/ndepths'//trim(tag7)//'.uf'

      inquire(file=trim(fndepths),  exist=ex)
      if (.not.ex) then
          write(lp,'(a)') 'nesting depths file for local grid does not exist:', &
                 'ndepths'//trim(tag7)//'.uf'
         stop '(nest_offline)'
      else
         open(10,file=trim(fndepths),form='unformatted')
         read(10) inner_ior8
         close(10)
         inner_depths=inner_ior8
      endif
         
      ! Extend local grid to boundary
      where(inner_depths>1e28) inner_depths=0.
      where (inner_depths(2          ,          :)>0.) &
         inner_depths(1,:)=inner_depths(2,:)
      where (inner_depths(inner_idm-1,          :)>0.) &
         inner_depths(inner_idm,:)=inner_depths(inner_idm-1,:)
      where (inner_depths(          :,          2)>0.) &
         inner_depths(:,1)=inner_depths(:,2)
      where (inner_depths(          :,inner_jdm-1)>0.) &
         inner_depths(:,inner_jdm)=inner_depths(:,inner_jdm-1)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setup pivot points for interpolation -- layer thickness not considered here

      ! These are the pivot points (on global grid) of points on the local grid
      ! Also bilinear interpolation coeffs
      allocate(inner_ipiv(inner_idm,inner_jdm))
      allocate(inner_jpiv(inner_idm,inner_jdm))
      allocate(inner_a1  (inner_idm,inner_jdm))
      allocate(inner_a2  (inner_idm,inner_jdm))
      allocate(inner_a3  (inner_idm,inner_jdm))
      allocate(inner_a4  (inner_idm,inner_jdm))

      print *,'Assigning pivot points'
      inner_ipiv=-1
      inner_jpiv=-1
      inner_a1=0.
      inner_a2=0.
      inner_a3=0.
      inner_a4=0.
      if (samegrid) then ! KAL should add security checks for lon, lat, depths here

         if (inner_idm/=idm .or. inner_jdm/=jdm) then
            print *,'This is not the same grid !'
            call exit(1)
         end if

         do j=1,inner_jdm
         do i=1,inner_idm
         if (inner_depths(i,j)>0.) then 
            inner_ipiv(i,j)=min(max(i,2),inner_idm-1)
            inner_jpiv(i,j)=min(max(j,2),inner_jdm-1)
            inner_a1(i,j)=1.
         end if
         end do
         end do

      else


         do j=1,inner_jdm
         do i=1,inner_idm
         if (inner_depths(i,j)>0.) then 

            ! Get pivot points of internal grid positions
            lon=inner_lon(i,j)
            lat=inner_lat(i,j)
            call oldtonew(lat,lon,lat_n,lon_n)
            call pivotp(lon_n,lat_n,ipiv,jpiv)

            ! For periodic grids -  make sure pivot point is in range
            if (lperiodic) ipiv=mod(idm+ipiv-1,idm)+1

            ! grace for i
            igrace=min(idm-ipiv,ipiv-1) ! negative when ipiv < 1 or ipiv > idm

            ! grace for j
            jgrace=min(jdm-jpiv,jpiv-1) ! negative when jpiv < 1 or jpiv > idm





            ! Check if this point is on global grid
            if (ipiv>=1.and.ipiv<idm.and.jpiv>=1.and.jpiv<jdm) then

               if (lperiodic) then
                  ipib=mod(ipiv,idm)+1
               else
                  ipib=ipiv+1
               end if

               ! check 
               if (ipiv+1 /= ipib) then
                  print *,'periodic grid needs testing !'
                  stop
               end if

               jpib=min(jdm,jpiv+1)

               ! Check if surrounded by sea cells on the global  grid
               if (gmsk(ipiv,jpiv).and.gmsk(ipib,jpiv).and. &
                   gmsk(ipiv,jpib).and.gmsk(ipib,jpib)) then

                  ! TODO -- bilincoeff should handle periodic grids
                  call bilincoeff(plon,plat,idm,jdm,lon,lat, &
                      ipiv,jpiv,ba1,ba2,ba3,ba4)

               ! Otherwise search for nearest point
               else
                  ! TODO -- nearestpoint should handle periodic grids
                  call nearestpoint(plon,plat,idm,jdm,lon,lat, &
                                    ipiv,jpiv,ba1,ba2,ba3,ba4, &
                                    gmsk,ass,lperiodic)

               endif

               inner_ipiv(i,j)=ipiv
               inner_jpiv(i,j)=jpiv
               inner_a1  (i,j)=ba1
               inner_a2  (i,j)=ba2
               inner_a3  (i,j)=ba3
               inner_a4  (i,j)=ba4
            else if (igrace>-ngrace .and. jgrace > -ngrace ) then
               print *,'Warning: Pivot point outside model domain - saved by grace'
               print *,'Warning: Pivot point outside model domain - saved by grace'
               print *,'Warning: Pivot point outside model domain - saved by grace'
               print *,'Warning: Pivot point outside model domain - saved by grace'
               print *,'Warning: Pivot point outside model domain - saved by grace'
               print *,'Warning: Pivot point outside model domain - saved by grace'
               print *,'Warning: Pivot point outside model domain - saved by grace'
               print *,ipiv,jpiv
               ipiv=max(1,min(ipiv,idm))
               jpiv=max(1,min(jpiv,jdm))
               ! TODO -- nearestpoint should handle periodic grids
               call nearestpoint(plon,plat,idm,jdm,lon,lat, &
                                 ipiv,jpiv,ba1,ba2,ba3,ba4, &
                                 gmsk,ass,lperiodic)
               inner_ipiv(i,j)=ipiv
               inner_jpiv(i,j)=jpiv
               inner_a1  (i,j)=ba1
               inner_a2  (i,j)=ba2
               inner_a3  (i,j)=ba3
               inner_a4  (i,j)=ba4
            else
               print *,'Warning: Pivot point outside model domain - not saved by grace'
               print *,'Warning: Pivot point outside model domain - not saved by grace'
               print *,'Warning: Pivot point outside model domain - not saved by grace'
               print *,'Warning: Pivot point outside model domain - not saved by grace'
               print *,'Warning: Pivot point outside model domain - not saved by grace'
               print *,'Warning: Pivot point outside model domain - not saved by grace'
               print *,'Warning: Pivot point outside model domain - not saved by grace'
               print *,ipiv,jpiv
               stop
            end if
         end if
         end do
         end do

      end if


      ! Store for later
      write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      print *,'Assigning points - done'
      print *, 'Storing to '//'Nest/pivots'//trim(tag7)//'.uf'
      open(10,file='Nest/pivots'//trim(tag7)//'.uf',  form='unformatted')
      write(10)inner_ipiv,inner_jpiv
      write(10)inner_a1  ,inner_a2  
      write(10)inner_a3  ,inner_a4  
      close(10)

      ! For tecplot diagnostics
      open(10,file='Nest/pivots'//trim(tag7)//'.tec', form='formatted',status='replace')
      do j=1,inner_jdm
      do i=1,inner_idm
      if (inner_depths(i,j)>0) then
         write(10,'(2i6,2f14.4,2i6,4f14.4)') i,j, &
             inner_lon(i,j),inner_lat(i,j), &
             inner_ipiv(i,j),inner_jpiv(i,j), &
             inner_a1(i,j), inner_a2(i,j), &
             inner_a3(i,j), inner_a4(i,j)
      end if
      end do
      end do
      close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We are set .. Now we must read the daily average files from the
! global model, interpolate the global daily average fields to the
! inner  model, then dump the data in suitable formats

      inquire(exist=ex,file='nestoffline.in') ! Contains daily average files, one per line

      if (.not. ex) then
         print *,'File "nestoffline.in" is not present'
         stop '(nest_offline)'
      else
         open(118,file='nestoffline.in')
      end if


      ! Get files to read from nestoffline.in
      read(118,*,iostat=ios) cfile
      do while( ios==0)

         ftype=getfiletype(trim(cfile))
         call initHF(hfile,trim(cfile),trim(ftype))
         print *,'File name is ',trim(cfile)
         print *,'File type is ',trim(ftype)

         iyear=hfile%iyear
         iday =hfile%iday

         ! For weekly/daily average - get ubavg/vbavg/pbavg 
         if (trim(ftype)=='nersc_daily' .or. trim(ftype)=='nersc_weekly') then
            ! We need a restart file to do this properly....
            ! Locate it in "botprop.in"
            inquire(exist=ex,file='botprop.in')
            if (ex) then
               open(10,file='botprop.in')
               read(10,*) rstfile
            else
               print *,'For hycom files of type weekly average, you'
               print *,'need to specify a file with the bottom properties'
               print *,'used in the model. This is any restart file from '
               print *,'the model run. Specify in botprop.in'
               stop
            end if
            ! Some work is needed .... put in this routine
            call bavg_flds_from_ave(ubavg,vbavg,pbavg,depths, trim(rstfile),trim(cfile),idm,jdm,kdm)
         else
            call HFReadField(hfile,pbavg,idm,jdm,'pbavg   ',0,1)
            call HFReadField(hfile,ubavg,idm,jdm,'ubavg   ',0,1)
            call HFReadField(hfile,vbavg,idm,jdm,'vbavg   ',0,1)
         end if


         nestingfile=filenesting(rtd,'Nest/')

         irec_offset = 1 ! One per day, so offset is one
         irec=irec_offset-1


         ! Open nesting files
         write(lp,'(a,a,a)') 'saving records ',  &
            ' to files=',trim(nestingfile)//'_xx'
         inquire(iolength=j)inner_ior4(1:inner_inest,1:inner_jdm)
         open(10,file=trim(nestingfile)//'_i1',form='unformatted', &
              access='direct',recl=j)
         inquire(iolength=j)inner_ior4(1:inner_inest,1:inner_jdm)
         open(11,file=trim(nestingfile)//'_ii',form='unformatted', &
              access='direct',recl=j)
         inquire(iolength=j)inner_ior4(1:inner_idm,1:inner_inest)
         open(12,file=trim(nestingfile)//'_j1',form='unformatted', &
              access='direct',recl=j)
         inquire(iolength=j)inner_ior4(1:inner_idm,1:inner_inest)
         open(13,file=trim(nestingfile)//'_jj',form='unformatted', &
              access='direct',recl=j)
         open(16,file=trim(nestingfile)//'.hdr',form='formatted', &
              access='direct',recl=100,status='unknown')
         open(15,file='testfld.tec', &
              form='formatted',status='replace')


         ! Go through variables in global model, read the necessary
         ! fields...
         intf=0.
         do ivar=1,nvar! Ivar needs to be rubbed out

            if     (ivar>=1 .and. ivar <=kdm) then ! Salinity
               klevel=ivar
               call HFReadField(hfile,fld1,idm,jdm,'saln    ',klevel,1)
               !call HFReadDPField(hfile,dp,idm,jdm,klevel,1) 
               call HFReadDPField_p(hfile,dp,idm,jdm,klevel,1) ! dp in pressure coords
               isvelocity=.false.
               write(char15,'(a11,i4)') 'saln level ',klevel
               write(char5,'(a3,i2.2)') 'SAL',klevel
               write(char3,'(a3)') 'SAL'
            elseif (ivar>=kdm+1 .and. ivar <=2*kdm) then ! Temperature
               klevel=ivar-kdm
               call HFReadField(hfile,fld1,idm,jdm,'temp    ',klevel,1)
               call HFReadDPField_p(hfile,dp,idm,jdm,klevel,1) ! dp in pressure coords
               isvelocity=.false.
               write(char15,'(a11,i4)') 'temp level ',klevel
               write(char5,'(a3,i2.2)') 'TEM',klevel
               write(char3,'(a3)') 'TEM'
            elseif (ivar>=2*kdm+1 .and. ivar <=3*kdm) then ! Interface values
               klevel=ivar-2*kdm
               call HFReadDPField_p(hfile,dp,idm,jdm,klevel,1) ! dp in pressure coords
               if (klevel==1) then
                  intf=dp
               else
                  intf=intf+dp
               end if

               fld1=intf
               isvelocity=.false.
               write(char15,'(a11,i4)') 'intf level ',klevel
               write(char5,'(a3,i2.2)') 'INT',klevel
               write(char3,'(a3)') 'INT'
               !print *,minval(dp,dp>1e-8),maxval(dp,dp>1e-8)
            elseif (ivar>=3*kdm+1 .and. ivar <=4*kdm) then ! Velocities
               klevel=ivar-3*kdm

               call HFReaduvtot(hfile,fld1,fld2,idm,jdm,k,1)
               call HFReaduvbaro(hfile,ubavg,vbavg,idm,jdm,1)
               !call HFReadDPField (hfile,dp,idm,jdm,klevel,1)
               call HFReadDPField_p(hfile,dp,idm,jdm,klevel,1) ! dp in pressure coords
               fld1=fld1-ubavg
               fld2=fld2-ubavg

               isvelocity=.true.
               call rotate2(fld1,fld2, plat,plon,idm,jdm,'m2l')
               write(char15,'(a11,i4)') 'Vel  level ',klevel
               write(char5,'(a3,i2.2)') 'VEL',klevel
               write(char3,'(a3)') 'UT '
               write(char3_2,'(a3)') 'VT '
            elseif (ivar==4*kdm+1) then ! Baro Pressure
               dp=onem ! Must be set prior to interpolation
               fld1=pbavg
               isvelocity=.false.
               char15='pbavg'
               write(char5,'(a3,i2.2)') 'PBA',1
               write(char3,'(a3)') 'PB '
            elseif (ivar==4*kdm+2) then ! Baro Velocities
               dp=onem ! Must be set prior to interpolation
               fld1=ubavg
               fld2=vbavg
               call rotate2(fld1,fld2, plat,plon,idm,jdm,'m2l')
               isvelocity=.true.
               char15='Bar vel'
               write(char5,'(a3,i2.2)') 'VBA',1
               write(char3,'(a3)') 'UB '
               write(char3_2,'(a3)') 'VB '
            elseif (ivar==4*kdm+3) then  ! Ice velocities
               dp=onem ! Must be set prior to interpolation
               klevel=0
               call HFReadfield(hfile,fld1,idm,jdm,'uice    ',0,1)
               call HFReadfield(hfile,fld2,idm,jdm,'vice    ',0,1)
               call rotate2(fld1,fld2, plat,plon,idm,jdm,'m2l')
               isvelocity=.true.
               char15='ice vel'
               write(char5,'(a3,i2.2)') 'VI ',1
               write(char3,'(a3)') 'UI '
               write(char3_2,'(a3)') 'VI '
            elseif (ivar==4*kdm+4) then  ! Ice Thickness
               dp=onem ! Must be set prior to interpolation
               klevel=0
               call HFReadfield(hfile,fld1,idm,jdm,'hicem   ',0,1)
               isvelocity=.false.
               char15='hicem'
               write(char5,'(a3,i2.2)') 'HI ',1
               write(char3,'(a3)') 'HI '
            elseif (ivar==4*kdm+5) then  ! Ice Concentration
               dp=onem ! Must be set prior to interpolation
               klevel=0
               call HFReadfield(hfile,fld1,idm,jdm,'ficem   ',0,1)
               isvelocity=.false.
               char15='ficem'
               write(char5,'(a3,i2.2)') 'FI ',1
               write(char3,'(a3)') 'FI '
            else

               if (mnproc==1) then
                  print *,'You shouldnt be here ...'
                  print *,'4*kdm+5  is ',4*kdm+5
                  print *,'ivar     is ',ivar
                  print *,'irec     is ',irec
                  print *,'nvar_tot is ',nvar_tot
                  print *,'nvar     is ',nvar
               end if
               call xcstop('nestcondsave')
               stop '(nestcondsave)'
            end if

            ! Use mask
            where (depths<.1) 
               fld1=0.
               fld2=0.
            end where

            if (isvelocity) then

               do j=1,jdm-1
               do i=1,idm-1

                  if (depths(i,j)<.1 .and. depths(i+1,j)<.1) then
                     fld1(i,j)=0.
                  end if

                  if (depths(i,j)<.1 .and. depths(i,j+1)<.1) then
                     fld2(i,j)=0.
                  end if
               end do
               end do
            end if

            if (count(fld1>1e25)>0) then
               print *,'Error in nestoffline'
               print *,'hugecount: ',count(fld1>1e25)
               stop
            end if


            ! Initialize inner model output arrays
            inner_tmp1=0.
            inner_tmp2=0.

! Do actual interpolation using what we calculated earlier 
!print *,'Warning - interpolation not ok for periodic grid!'
!print *,'Warning - interpolation not ok for periodic grid!'
!print *,'Warning - interpolation not ok for periodic grid!'
!print *,'Warning - interpolation not ok for periodic grid!'
!print *,'Warning - interpolation not ok for periodic grid!'
!print *,'Warning - interpolation not ok for periodic grid!'
!print *,'Warning - interpolation not ok for periodic grid!'
!print *,'Warning - interpolation not ok for periodic grid!'
inner_tmp1=0.
inner_tmp2=0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!$OMP PARALLEL DO PRIVATE(i,j,lon,lat,lon_n,lat_n,ipiv,jpiv,  &
!$OMP                     ba1,ba2,ba3,ba4,ass,k)
            do j=1,inner_jdm
            do i=1,inner_idm
            if (inner_depths(i,j)>.1) then
            !if (inner_tile_mask(i,j)) then 
               !ipiv=inner_ipiv(i,j)-i0 ! On local grid
               !jpiv=inner_jpiv(i,j)-j0 ! On local grid


               ipiv=inner_ipiv(i,j)
               jpiv=inner_jpiv(i,j)

               if (lperiodic) then
                  ipib=mod(ipiv,idm)+1
               else
                  ipib=ipiv+1
               end if

               ! check 
               if (ipiv+1 /= ipib) then
                  print *,'periodic grid needs testing !'
                  stop
               end if

               useold=.false.
               if (trim(char3)/='INT') then
                  ba1 =inner_a1  (i,j) * dp(ipiv  ,jpiv  )
                  ba2 =inner_a2  (i,j) * dp(ipib  ,jpiv  )
                  ba3 =inner_a3  (i,j) * dp(ipib  ,jpiv+1)
                  ba4 =inner_a4  (i,j) * dp(ipiv  ,jpiv+1)

                  dpsum=dp(ipiv,jpiv)+dp(ipib,jpiv+1)+dp(ipib,jpiv+1)+dp(ipiv,jpiv+1)
                  basum=ba1+ba2+ba3+ba4

                  !if (dpsum<.1 .and. inner_depths(i,j)>1.) then
                  !   print *,'dpsum error !!',i,j,inner_depths(i,j),klevel,dpsum
                  !   print *,depths(ipiv  ,jpiv  ),inner_a1(i,j)
                  !   print *,depths(ipib  ,jpiv  ),inner_a2(i,j)   
                  !   print *,depths(ipib  ,jpiv+1),inner_a3(i,j)
                  !   print *,depths(ipiv  ,jpiv+1),inner_a4(i,j)
                  !   stop
                  !end if
                  !print *,dpsum



                  if (basum>onem*1e-2) then
                     ba1=ba1/basum
                     ba2=ba2/basum
                     ba3=ba3/basum
                     ba4=ba4/basum
                  elseif (dpsum >onem) then
                     !use value at max dp point 
                     ba1=0.; ba2=0. ; ba3=0. ; ba4=0.
                     maxdpind=maxloc((/dp(ipiv,jpiv),dp(ipib,jpiv),dp(ipib,jpiv+1),dp(ipib,jpiv+1)/))
                     if (maxdpind(1)==1) then
                        ba1=1.
                     else if (maxdpind(1)==2) then
                        ba2=1.
                     else if (maxdpind(1)==3) then
                        ba3=1.
                     else if (maxdpind(1)==4) then
                        ba4=1.
                     else
                        print *,'Excuse me?'
                        stop
                     end if
                  elseif (klevel==1) then
                     print *,'threshold error - you may adjust this in nest_offline'
                     print *,'Sum of weights*dp is   ',basum
                     print *,'Sum of neighbour dp is ',dpsum
                     print *,'i,j =                  ',i,j
                     print *,minval(dp,dp>1e-8),maxval(dp,dp>1e-8)
                     print *,minval(inner_ipiv),maxval(inner_ipiv)
                     print *,minval(inner_jpiv),maxval(inner_jpiv)
                     print *,char3
                     call zaiopf('no_dptest.a','replace',99)
                     call zaiowr(dp,ip,.false.,amin,amax,99,.false.) !plon
                     dp=float(inner_ipiv)
                     call zaiowr(dp,ip,.false.,amin,amax,99,.false.) !plon
                     dp=float(inner_jpiv)
                     call zaiowr(dp,ip,.false.,amin,amax,99,.false.) !plon
                     call zaiowr(plon,ip,.false.,amin,amax,99,.false.) !plon
                     call zaiocl(99)



                     stop
                  else
                     ba1=0.
                     ba2=0.
                     ba3=0.
                     ba4=0.
                     useold=.true.
                  end if
               else
                  ba1 =inner_a1  (i,j)
                  ba2 =inner_a2  (i,j)
                  ba3 =inner_a3  (i,j)
                  ba4 =inner_a4  (i,j)
               end if

               ! Use old only if klevel > 1
               if (useold .and. klevel==1) then
                  print *,'Fatal error - can not use old if klevel=1'
                  print *,dpsum
                  print *,klevel
                  stop
               end if


               ! We should have a representation of this point now...
               if (useold) then
                  inner_tmp1(i,j)=old_inner_tmp1(i,j)
               else
                  inner_tmp1(i,j)=ba1*fld1(ipiv  ,jpiv  )      &
                                 +ba2*fld1(ipiv+1,jpiv  ) &
                                 +ba3*fld1(ipiv+1,jpiv+1)  &
                                 +ba4*fld1(ipiv  ,jpiv+1)
               end if

               ! If velocity, then tmp2 will contain second comp
               if (isvelocity)  then
                  if (useold) then
                     inner_tmp2(i,j)=old_inner_tmp2(i,j)
                  else
                     inner_tmp2(i,j)=ba1*fld2(ipiv  ,jpiv  )      &
                                    +ba2*fld2(ipiv+1,jpiv  ) &
                                    +ba3*fld2(ipiv+1,jpiv+1)  &
                                    +ba4*fld2(ipiv  ,jpiv+1)
                  end if
               end if
            end if
            end do
            end do

            ! Set old_inner_tmp1 old_inner_tmp2 - used when encountering
            ! empty layers _below_ layer 1
            old_inner_tmp1=inner_tmp1
            old_inner_tmp2=inner_tmp2

            ! Rotate if velociy
            if (isvelocity) then
               call rotate2(inner_tmp1,inner_tmp2,inner_lat,inner_lon, &
                           inner_idm,inner_jdm,'l2m')
            end if

            ! Dump it to nest files!
            irec=irec+1
            inner_ior4=inner_tmp1
            write(lp,'(a,i5,a,2e14.2)') 'Writing record ',irec, &
               ' '//char15,minval(inner_ior4),maxval(inner_ior4)
            write(10,rec=irec) inner_ior4(1:inner_inest,1:inner_jdm)
            write(11,rec=irec)  &
               inner_ior4(inner_iidm:inner_idm,1:inner_jdm)
            write(12,rec=irec) inner_ior4(1:inner_idm,1:inner_inest)
            write(13,rec=irec)  &
               inner_ior4(1:inner_idm,inner_jjdm:inner_jdm)
            write(16,103,rec=irec)char3,klevel,irec,irec_offset, &
                  rtd%iyy,rtd%idd,rtd%ihh,kdm,minval(inner_ior4), &
                  maxval(inner_ior4)


            ! Dump 2nd vector component if this is a velocity field!
            if (isvelocity) then
               irec=irec+1
               inner_ior4=inner_tmp2
               write(lp,'(a,i5,a,2e14.2)') 'Writing record ',irec, &
                  ' '//char15,minval(inner_ior4),maxval(inner_ior4)
                write(10,rec=irec) inner_ior4(1:inner_inest,1:inner_jdm)
                write(11,rec=irec)  &
                   inner_ior4(inner_iidm:inner_idm,1:inner_jdm)
                write(12,rec=irec) inner_ior4(1:inner_idm,1:inner_inest)
                write(13,rec=irec)  &
                   inner_ior4(1:inner_idm,inner_jjdm:inner_jdm)
                write(16,103,rec=irec)char3_2,klevel,irec,irec_offset, &
                     rtd%iyy,rtd%idd,rtd%ihh,kdm,minval(inner_ior4), &
                     maxval(inner_ior4)
            end if
         end do ! Variable loop
         close(10)
         close(11)
         close(12)
         close(13)
         close(16)
         close(15)

         print *,'Finished processing nesting conditions from '//trim(cfile)
         print *
         print *

         ! Read next file to process
         read(118,*,iostat=ios) cfile
      end do !  file loop

      print *,'Successful exit'

  101 format(30i5)      
  102 format(10e14.3)
  103 format(a3," level=",i4," record=",i5," offset=",i5," date=",i4, &
             " ",i3," ",i2," kdm=",i3," min/max:",2e12.2)


      end program nest_offline





         


