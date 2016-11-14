      PROGRAM csr2mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program reads information about the model grid (lon-lat-depth-dimensions)
! and amplitudes and phases for 8 constituents from a global grid, and interpolates
! these to the model boundary grid points and save the result to a file
! CSRobc_tide.dat which is read by HYCOM to compute the tidal boundary conditions.
!
! It reads the following files:
!    regional.grid.[ab]  --- contains model lon lat positions.
!    regional.depth.[ab] --- model bathymetry.
!    ortho_csr_4.0.ot    --- tidal amplitudes and phases for 8 constituents
!    input.usr           --- user specified constituents to include (DO NOT ALTER)
!
! Path to ortho_csr_4.0.ot is specified by environment variable CSR_PATH
!
! Grid dimensions are retrieved from regional.grid.b.
!
! Output is dumped to CSRobc_tide.dat
!
! Diagnostic output of phase and amplitude is dumped to csrtide.nc
!
! -Originally developed by Knud Simonsen, and later adapted and modified
!     by Geir Evensen,   Geir.Evensen@nersc.no
! -Modified for regional.grid/depth files by Knut Liseter. Grid sizes
!    are determined by values in regional.grid.b. Added netcdf diafgnostics
!    Generic makefile with arch-dependent includes from Make.Inc/make.$machine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use mod_xc , jj => idm,  kk=> jdm ! Uses mod_xc for dims
      use mod_za
      use m_constituents
      use m_harmonics
      use m_har_mdomain
      implicit none

      !integer jj,kk

      integer, parameter :: nconst1=8
      integer nconst

      real, allocatable :: mlon(:,:), mlat(:,:), amppha(:,:)
      integer, allocatable :: iutil(:,:)
      real, allocatable :: depths(:,:)
      character(len=7) tag7

      integer i,j,k,m,ii,i1,i2,itec
      real dlon,dlat
      real, parameter :: xmjd=0.0, xmjd0=0.0
      integer ifreq(nconst1)

      character*4 tname2(nconst1)
      character answ*1,dummy*100
      character(len=200) :: tecvar=' '
      real freq2(nconst1)

      logical ex
      logical pseudo,radiation,isdata
      integer mdyhms(6)
      real u(3,2), v(3,2)
      real, dimension(100) :: x,y,h1,h2,amp,pha 

      ! Nev vars introduced by Knut
      real, parameter :: huge=2.0**99
      real :: hmin, hmax
      integer, allocatable :: imsk(:,:)
      ! Change this to use old depths and newpos files
      logical, parameter :: use_old_files=.false. 
      character(len=200) :: csrpath, csrfile

      if (use_old_files) then
         write(*,'(a)',advance='no')
     &      ' Choose x and y dimensions for grid: '
         read(*,*)jj,kk
      else
         ! KAL - this uses regional.grid to retrieve grid dims
         inquire(file='regional.grid.b',exist=ex)
         if (.not. ex) then
            print *,'regional.grid.b not present - needed to get dims'
            call exit(1)
         end if
         call xcspmd()
         call zaiost()
      end if


      ! Not needed anymore
      allocate(mlon(jj,kk), mlat(jj,kk), amppha(2*(jj+kk),2*nconst1) )
      allocate(iutil(jj,kk) )
      allocate(depths(jj,kk))
      inquire(file='input.CSR',exist=ex)
      if (.not.ex) stop 'input.CSR file does not exist'
      OPEN(10,FILE='input.CSR',STATUS='OLD')
      READ(10,'(a)')dummy
      if (dummy(1:3) /= 'CSR') then
         print *,'error in input.CSR'
         stop
      endif
      READ(10,'(a)')dummy
      DO i=1,nconst1
       ifreq(i)=0
      ENDDO
      ii=0
      DO i=1,100
         READ(10,'(a1)',end=999)answ
         IF(answ.EQ.'y')THEN
           ii=ii+1
           if (ii > nconst1) stop 'nconst1 must be increased'
           ifreq(ii)=i
         ENDIF
      ENDDO
 999  CLOSE(10)

      IF (ii <= nconst1) nconst=ii
      print *,'nconst=',nconst

      CALL constituents(ifreq,nconst,freq2,tname2)
       
! ---
! --- ----------------- Load lon,lat of the model grid points ------------
! ---

      ! Deprecated for new versions  - set use old files in case of old model files
      if (use_old_files) then
         inquire(file='newpos.uf',exist=ex)
         if (.not.ex) stop 'newpos.uf file does not exist'
         OPEN (unit=10,file='newpos.uf',form='unformatted')
            read(10)mlat,mlon
         close(10)
!!GE     where (mlon < 0.0) mlon=mlon+360.0

         write(tag7,'(i3.3,a,i3.3)')jj,'x',kk
         inquire(file='depths'//tag7//'.uf',exist=ex)
         if (.not.ex) stop 'depths file does not exist'
         open(unit=10,file='depths'//tag7//'.uf',form='unformatted')
         read(10)depths
         close(10)
      else
         allocate(imsk(jj,kk)) ! dummy array - not used but must be allocated
         ! NB: No check on values here (.b - files)
         call zaiopf('regional.depth.a', 'old', 701)
         call zaiord(depths,imsk,.false.,hmin,hmax,701)
         call zaiocl(701)
         where (depths>0.5*huge) depths=0.
         call zaiopf('regional.grid.a', 'old', 701)
         call zaiord(mlon,imsk,.false.,hmin,hmax,701)
         call zaiord(mlat,imsk,.false.,hmin,hmax,701)
         call zaiocl(701)
      end if



! --- --------------------------------------------------------------------

      ! KAL - needs environment variable CSR_PATH as input
      call getenv('CSR_PATH',csrpath)
      if (trim(csrpath)=="") then
         print *,'environment variable CSR_PATH must be set to location'
     &    //' of CSR data'
         call exit(1)
      end if

      csrfile=trim(csrpath)//'ortho_csr_4.0.ot'
      inquire(file=trim(csrfile),exist=ex)
      if (.not.ex) then
         print *,trim(csrfile)//' not found...'
         call exit(1)
      end if
      OPEN (51,file=trim(csrfile))

C     inquire(file='/home/parallab/nersc/geirev/Progs/Tide_CSR/Data/orth
C    &o_csr_4.0.ot',exist=ex)
C     if (.not.ex) stop 'NO /home/parallab/nersc/geirev/Progs/Tide_CSR/D
C    &ata/ortho_csr_4.0.ot'
C     OPEN (51,file='/home/parallab/nersc/geirev/Progs/Tide_CSR/Data/ort
C    &ho_csr_4.0.ot') 

! --- ---------------------------------------- Initialize array ----------
      amppha=0.0

! --- ---------------------------------------- Estimate boundary values --
! --- Data are generated mlong for k=2, and k=kk-1 mlong the
! --- 'EW' boundaries and mlong j=2, j=jj-1 mlong the 'NS' boundaries.
! --- If this is changed, the calculations on index m below, and in 
! --- subroutine nest_barotp.F90  in the model code must be modeified.
! ---

      DO j=1,jj
         m=j                           !South  1 <m<jj
         k=2
 
         dlon=mlon(j,k)
         dlat=mlat(j,k)

         CALL harmonics(jj,kk,nconst,dlon,dlat,xmjd0,m,k,ifreq,amppha,
     .                  depths(j,k))

         m=2*jj+kk+1-j                    !North  jj+kk+1<m<2*jj+kk 
         k=kk-1
         dlon=mlon(j,k)
         dlat=mlat(j,k)

         CALL harmonics(jj,kk,nconst,dlon,dlat,xmjd0,m,k,ifreq,amppha,
     .                  depths(j,k))

        
      ENDDO

      DO k=1,kk
         m=jj+k                  !East     jj+1 <m<jj+kk
         j=jj-1
         dlon=mlon(j,k)
         dlat=mlat(j,k)
         CALL harmonics(jj,kk,nconst,dlon,dlat,xmjd0,m,k,ifreq,amppha,
     .                  depths(j,k))

         m=2*(jj+kk)+1-k
         j=2                    !West    2*jj+kk+1<m<2*(jj+kk)
         dlon=mlon(j,k)
         dlat=mlat(j,k)
         CALL harmonics(jj,kk,nconst,dlon,dlat,xmjd0,m,k,ifreq,amppha,
     .                  depths(j,k))
        
      ENDDO

      tecvar(1:13)='VARIABLES="m"'
      DO i=1,nconst
         i1=14+(i-1)*17
         tecvar(i1:i1+16)=',"axxxx","pxxxx"'
         WRITE(tecvar(i1+3:i1+6),'(a4)')tname2(i)
         WRITE(tecvar(i1+11:i1+14),'(a4)')tname2(i)
      ENDDO

 200  FORMAT('TEXT X=80, Y=',i2,', H=2,T="Model Dim.  ',2I4,'"')
 201  FORMAT('TEXT X=80, Y=',i2,', H=2,T="Constituents',I2,'"')
 202  FORMAT('TEXT X=80, Y=',i2,', H=2,T="Name Freq (cyc/d)"')
 203  FORMAT('TEXT X=80, Y=',i2,', H=2,T="',A4,1x,F11.9,'"')
 204  FORMAT('ZONE I=',I4,', F=POINT')

      OPEN(10,FILE='./CSRobc_elev.dat',STATUS='UNKNOWN')
      WRITE(10,'(200a)')trim(tecvar)
      itec=96
      WRITE(10,200)itec,jj,kk
      itec=itec-3
      WRITE(10,201)itec,nconst
      itec=itec-3
      WRITE(10,202)itec

      DO i=1,nconst
         itec=itec-3
         write(10,203)itec,tname2(i),freq2(i)
      ENDDO
      WRITE(10,204)2*(jj+kk)

      DO m=1,2*(jj+kk)
         write(10,'(i4,16f8.2)')m,(amppha(m,i),i=1,2*nconst)
      ENDDO
       
      WRITE(*,'(''SW  m=   1'')')
      WRITE(*,'(''SE  m='',I4)')jj
      WRITE(*,'(''NE  m='',I4)')jj+kk
      WRITE(*,'(''NW  m='',I4)')2*jj+kk

      WRITE(10,'(''GEOMETRY X=0, Y=0, T=LINE C=BLACK, CS=GRID'')')
      m=3
      WRITE(10,'(I2)')m             !Number of lines
      m=2
      i1=0                          !On axis value
      i2=1000                       !A max amp/phase value
      WRITE(10,'(I2)')m             !Points in line
      WRITE(10,'(2I5)')jj,i1
      WRITE(10,'(2I5)')jj,i2
      WRITE(10,'(I2)')m             !Points in line
      WRITE(10,'(2I5)')jj+kk,i1
      WRITE(10,'(2I5)')jj+kk,i2
      WRITE(10,'(I2)')m             !Points in line
      WRITE(10,'(2I5)')2*jj+kk,i1
      WRITE(10,'(2I5)')2*jj+kk,i2
      CLOSE(10)
      WRITE(*,*)
      WRITE(*,'(''OBS data are dumped to file CSRobc_elev.dat'')')

      CALL har_mdomain(jj,kk,nconst,mlon,mlat,depths,
     &   xmjd,ifreq,tname2)
       
      WRITE(*,*)
      WRITE(*,'(''Note,  *.dat files are in Tecplot format'')')

      CLOSE(51)
      STOP
      END  
