      module mod_nest_io
      integer, save :: inner_idm, inner_jdm
      integer, save :: inner_iidm, inner_jjdm
      integer, parameter :: inner_inest=20
      real, allocatable, save :: inner_lon(:,:)
      real, allocatable, save :: inner_lat(:,:)
      real, allocatable, save :: inner_depths(:,:)
      contains
c
      subroutine readnest(nestname, cvar,k,fld,nx,ny,inest)
      use mod_xc
      implicit none
      integer         , intent(in)  :: inest, nx, ny
      character(len=*), intent(in)  :: nestname
      character(len=8), intent(in)  :: cvar
      integer         , intent(in)  :: k
      real            , intent(out) :: fld(nx,ny)
c
      integer :: i1,i2,j1,j2, ibnd, j, tmp_kdm, tmp_inest, k2, i,
     &   ios,irec
      real*4 :: inner_ior4(nx, ny), hmin,hmax
      real*8 :: dtime
      character(len=2) :: cbnd
      character(len=8) :: cvarf
c
c --- Open header
      open(10,file=trim(nestname)//'.hdr',form='formatted', 
     &        access='sequential',status='old')
      read(10,104) dtime,tmp_kdm,tmp_inest
      if (inner_inest/=tmp_inest ) then
         write(6,'(a)') 'error reading nest file'
         call xcstop('(readnest)')
         stop '(readnest)'
      endif
c --- Get sought-for variable
      ios=0 ; irec=0 ; i=1
      do while (ios==0 .and. irec==0)
         read (10,103,iostat=ios) cvarf, k2, hmin,hmax
         !print *,cvarf,k2, hmin, hmax
         if (trim(cvar)==trim(cvarf) .and. k==k2) irec=i
         i=i+1
      end do
      close(10)
      if (ios/=0 .or. irec==0) then
         write(6,'(a)') 'Could not find variable '//cvar,k
         call xcstop('(readnest)')
         stop '(readnest)'
      endif
c
c --- Read 4 edges
      do ibnd=1,4
         if     (ibnd==1)  then
            cbnd='i1'
            i1=1    ; i2=inest
            j1=1    ; j2=ny
         elseif (ibnd==2)  then
            cbnd='j1'
            i1=1    ; i2=nx
            j1=1    ; j2=inest
         elseif (ibnd==3)  then
            cbnd='ii'
            i1=nx-inest+1   ; i2=nx
            j1=1    ; j2=ny
         elseif (ibnd==4)  then
            cbnd='jj'
            i1=1    ; i2=nx
            j1=ny-inest+1   ; j2=ny
         end if
c
         inquire(iolength=j)inner_ior4(i1:i2,j1:j2)
         open(11,file=trim(nestname)//'_'//cbnd,form='unformatted',
     &           access='direct',recl=j,status='old')
         read(11,rec=irec) inner_ior4(i1:i2,j1:j2)
         fld(i1:i2,j1:j2)=inner_ior4(i1:i2,j1:j2)
         close(11)
      end do
      !where(inner_depths<.1) fld=undef
  103 format(a8,": level=",i4," min/max:",2e14.6)
  104 format("dtime, kdm, inest :",f16.4,i3,i3) 
      end subroutine


c
      subroutine readinnergrid(nestname)
      use mod_xc
      implicit none 
      character(len=*), intent(in)  ::  nestname
      real*4, allocatable :: inner_io1(:)
      real*4, allocatable :: inner_io2(:,:)
      logical :: ex
      integer :: rl, n2drec
c
CKAL  Yes, I killed nestloc
c --- Read grid size from regional.grid.b
      inquire(file=trim(nestname)//'regional.grid.b', exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(*,'(a)')
     &   'nesting depths file for local grid does not exist:',
     &           trim(nestname)//'regional.grid.b'
         call xcstop('readinnergrid')
         stop '(readinnergrid)'
      else
         open(10,file=trim(nestname)//'regional.grid.b', status='old')
         read(10,*) inner_idm
         read(10,*) inner_jdm
         close(10)
      end if
c --- inest hardcoded for now
      inner_iidm=inner_idm-inner_inest+1
      inner_jjdm=inner_jdm-inner_inest+1
      if (mnproc==1) then
         write(*,'("inest     = ",i)') inner_inest
         write(*,'("inner_idm = ",i)') inner_idm
         write(*,'("inner_jdm = ",i)') inner_jdm
      end if
c --- Read cell center positions from regional.grid.a
      n2drec=((inner_idm*inner_jdm+4095)/4096)*4096
      allocate(inner_io1(n2drec))
      allocate(inner_io2   (inner_idm,inner_jdm))
      allocate(inner_lon   (inner_idm,inner_jdm))
      allocate(inner_lat   (inner_idm,inner_jdm))
      allocate(inner_depths(inner_idm,inner_jdm))
      inquire(file=trim(nestname)//'regional.grid.a', exist=ex)
      if (.not.ex) then
         if (mnproc==1) write(*,'(a)')
     &   'nesting depths file for local grid does not exist:',
     &           trim(nestname)//'regional.grid.a'
         call xcstop('readinnergrid')
         stop '(readinnergrid)'
      else
         inquire(iolength=rl) inner_io1
         open(10,file=trim(nestname)//'regional.grid.a',
     &      form='unformatted', status='old',recl=rl,access='direct')
         read(10,rec=1) inner_io2 ; inner_lon=inner_io2
         read(10,rec=2) inner_io2 ; inner_lat=inner_io2
         close(10)
      end if
c
c
c ---    Read the local grid - now uses regional.depth.a
CKAL  n2drec = ((inner_idm*inner_jdm+4095)/4096)*4096
      inquire(file=trim(nestname)//'regional.depth.a', exist=ex)
      if (.not.ex) then
         write(6,'(a)')
     &   'nesting depths file for local grid does not exist:',
     &   trim(nestname)//'regional.depth.a'
         call xcstop('readinnergrid')
         stop '(readinnergrid)'
      else
         inquire(iolength=rl) inner_io2
         open(710,file=trim(nestname)//'regional.depth.a',
     &      status='old',access='direct',recl=rl,action='read')
         read(710,rec=1) inner_io2
         close(710)
         inner_depths=inner_io2
      endif
c --- Correct huge values
      where (inner_depths>0.5*2.0**99) inner_depths=0.
c
c --- Extend local grid to boundary
      where (inner_depths(2,:)>0.) 
     &   inner_depths(1,:)=inner_depths(2,:)
      where (inner_depths(inner_idm-1,:)>0.) 
     &   inner_depths(inner_idm,:)=inner_depths(inner_idm-1,:)
      where (inner_depths(:,2)>0.) 
     &   inner_depths(:,1)=inner_depths(:,2)
      where (inner_depths(:,inner_jdm-1)>0.)
     &   inner_depths(:,inner_jdm)=inner_depths(:,inner_jdm-1)
      deallocate(inner_io1)
      deallocate(inner_io2)
      end subroutine 

c --- NB: ipiv,jpiv refers to points on the Global (outer) grid
c --- Assign points on inner grid to points on the global grid of this model (if possible)
      subroutine assignpivots(nestname, gmsk, modlon, modlat)
      use mod_xc
      use mod_confmap
      use m_nearestpoint
      use netcdf
      implicit none
      character(len=*), intent(in)    :: nestname
      real,    dimension(idm,jdm), intent(in) :: modlon, modlat
      logical, dimension(idm,jdm), intent(in) :: gmsk
c
      integer :: idmid, jdmid, varid, ncid
      real :: lon,lat, lon_n, lat_n, ba1, ba2, ba3, ba4
      integer :: i,j, ipiv, jpiv, iol, ioli, iol4
      logical :: ex, ass
      character(len=200) :: ncfil
      character(len=11) :: tag7g
      real, parameter :: undef=-1e14
      integer, allocatable, dimension(:,:) :: 
     &   inner_ipiv, inner_jpiv
      real,    allocatable, dimension(:,:) :: 
     &   inner_a1, inner_a2, inner_a3, inner_a4
      real*4   ,  allocatable, dimension(:,:) :: iofld4
      integer*4,  allocatable, dimension(:,:) :: iofldi

c
      write(tag7g,'(i5.5,a,i5.5)')idm,'x',jdm
c
      allocate(inner_ipiv(inner_idm,inner_jdm))
      allocate(inner_jpiv(inner_idm,inner_jdm))
      allocate(inner_a1  (inner_idm,inner_jdm))
      allocate(inner_a2  (inner_idm,inner_jdm))
      allocate(inner_a3  (inner_idm,inner_jdm))
      allocate(inner_a4  (inner_idm,inner_jdm))
      allocate(iofld4    (inner_idm,inner_jdm))
      allocate(iofldi    (inner_idm,inner_jdm))
c
      if (mnproc==1) 
     &   write(*,'(a)') 'Assigning pivot points '


c --- Initialize conformal mapping
      call initconfmap(idm,jdm)
      print *,'Start Loop'
c
      inner_ipiv=-1
      inner_jpiv=-1
      inner_a1=0.
      inner_a2=0.
      inner_a3=0.
      inner_a4=0.
      do j=1,inner_jdm
      do i=1,inner_idm
      if (inner_depths(i,j)>0.) then 
c
c ---    Get pivot points of internal grid positions
         lon=inner_lon(i,j)
         lat=inner_lat(i,j)
         call oldtonew(lat,lon,lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)
c
c ---    Check if this point is on global grid
         if (ipiv>=1.and.ipiv<idm.and.jpiv>=1.and.jpiv<jdm) then
c
c ---       Check if surrounded by sea cells on the global  grid
            if (gmsk(ipiv  ,jpiv  ).and.gmsk(ipiv+1,jpiv  ).and.
     &          gmsk(ipiv  ,jpiv+1).and.gmsk(ipiv+1,jpiv+1))  then
               call bilincoeff(modlon,modlat,idm,jdm,lon,lat,
     &            ipiv,jpiv,ba1,ba2,ba3,ba4)
c
c ---       Otherwise search for nearest point
            else
               call nearestpoint(modlon,modlat,idm,jdm,lon,lat,
     &                           ipiv,jpiv,ba1,ba2,ba3,ba4,
     &                           gmsk,ass,.false.)
               if (.not.ass) then
                  print *,'Beware Could not find pivot point!!!!!'
                  print *,'Verify the netcdf diagmostic file'
                  print *,inner_depths(i,j),lon,lat
               endif
            endif
            inner_ipiv(i,j)=ipiv
            inner_jpiv(i,j)=jpiv
            inner_a1  (i,j)=ba1
            inner_a2  (i,j)=ba2
            inner_a3  (i,j)=ba3
            inner_a4  (i,j)=ba4
         end if
      end if
      end do
      end do
c --- Store pivot points for later
c     print '(a)', 'Storing pivots in            '//trim(nestname)
c    &         //'pivots'//trim(tag7g)//'.uf'
      print '(a)', 'Storing pivots in            '//trim(nestname)
     &         //'pivots.uf'
      inquire(iolength=iol4) iofldi
      inquire(iolength=ioli) iofld4
      !print *,iol4,ioli
c     open(10,file=trim(nestname)//'pivots'//trim(tag7g)
c    &     //'.uf',  form='unformatted', access='direct',recl=iol4)
      open(10,file=trim(nestname)//'pivots.uf',  
     &   form='unformatted', access='direct',recl=iol4)
      iofldi=inner_ipiv ; write(10,rec=1) iofldi
      iofldi=inner_jpiv ; write(10,rec=2) iofldi
      iofld4=inner_a1   ; write(10,rec=3) iofld4
      iofld4=inner_a2   ; write(10,rec=4) iofld4
      iofld4=inner_a3   ; write(10,rec=5) iofld4
      iofld4=inner_a4   ; write(10,rec=6) iofld4
      close(10)
c
c --- Netcdf Diagnostics
      ncfil=trim(nestname)//'pivots.nc'
      print '(a)', 'Storing pivot diagnostics in '//trim(ncfil)
      iofld4=inner_depths; where(inner_depths<0.1) iofld4=undef
      call ncerr(NF90_create(ncfil,NF90_CLOBBER,ncid))
      call ncerr(NF90_def_dim(ncid,'idm',inner_idm,idmid))
      call ncerr(NF90_def_dim(ncid,'jdm',inner_jdm,jdmid))
      call ncerr(NF90_def_var(ncid,'depths',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,iofld4))
c
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'lon',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,inner_lon))
c
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'lat',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,inner_lat))
c
      iofld4=inner_ipiv; where(inner_depths<0.1) iofld4=undef
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'ipiv',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,real(iofld4)))
c
      iofld4=inner_jpiv; where(inner_depths<0.1) iofld4=undef
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'jpiv',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,real(iofld4)))
c
      iofld4=inner_a1; where(inner_depths<0.1) iofld4=undef
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'a1',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,real(iofld4)))
c
      iofld4=inner_a2; where(inner_depths<0.1) iofld4=undef
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'a2',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,real(iofld4)))
c
      iofld4=inner_a3; where(inner_depths<0.1) iofld4=undef
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'a3',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,real(iofld4)))
c
      iofld4=inner_a4; where(inner_depths<0.1) iofld4=undef
      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'a4',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'missing_value',undef))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,real(iofld4)))
      call ncerr(NF90_close(ncid))
c
      end subroutine

      subroutine ncerr(errcode)
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      if (errcode/=NF90_NOERR) then
         write(*,'(a)') NF90_STRERROR(errcode)
         stop '(ncerr)'
      end if
      end subroutine

c
      function filenesting2(dtime,yrflag,nestname)
      implicit none
      character(len=200) filenesting2
      real*8, intent(in) :: dtime
      character(len=*), intent(in) :: nestname
      integer, intent(in) :: yrflag
      character(len=4) :: cyy
      character(len=3) :: cdd
      character(len=2) :: chh
      integer l,trimlen, iyear, iday, ihour
c
      call forday(dtime, yrflag, iyear,iday,ihour)
      write(cyy,'(i4.4)') iyear
      write(cdd,'(i3.3)') iday-1
      write(chh,'(i2.2)') ihour
c
      filenesting2='nest_'//cyy//'_'//cdd//'_'//chh
      if (len_trim(nestname)/=0) then
         filenesting2=trim(nestname)//trim(filenesting2)
      end if
      end function filenesting2

      end module

