      module mod_gridp
c --- TODO: Complete module info
C --- Module contains tools for generating gridpoint data
C ---
C --- subroutine gridpinp - initialize gridpoint information arrays
C --- subroutine gridpdef - Allocate arrays used 
      use mod_xc
      implicit none
c --- Shield everything by default
      private

      type grid_points
         integer   station_index ! index to be used for file naming
         integer*4 i
         integer*4 j
         real*4 z(kdm)
         real*4 u(kdm)
         real*4 v(kdm)
         real*4 t(kdm)
         real*4 s(kdm)
         real*4 sealevel
         real*4 uwind
         real*4 vwind
#if defined (ICE)
         real*4 fice
         real*4 hice
         real*4 uice
         real*4 vice
#endif
      end type grid_points



c --- Array of grid point information. Actually only i/j points are kept
c --- in this variable (see type grid_points)
      type(grid_points), allocatable :: gp(:)      
c
c --- Array of stored variables. Keeps 2D and 3D fields at the
c --- gridpoints at different times
      type(grid_points), allocatable :: gpstore(:,:) ! Keeps info for diff times and points
c
c --- This variable holds information for each time stored in gpstore
      integer, allocatable :: itime(:)          
c
c --- Total nuber of grid points to be saved - set in gridpdef
      integer gpdim  
c
c --- Total times in temporary array. 31 is max number of days in a month
      integer, parameter ::  gptdim=31*24  
c
c --- Public flag for activating gridp - set in m_limits
      logical, save, public :: lgridp
c ---
      integer, save :: gpdim_global

c
c --- These are made available for hycom 
      public :: gridpinit, gridp_process

      contains

c --- -----------------------------------------------------------------
c --- Initialization of gridp routines
c --- -----------------------------------------------------------------
      subroutine gridpinit()
      use mod_diagnostics ! contains iday1, dday
      implicit none
      if (lgridp) then
         call gridpinp()
         call gridpalloc()
         call gridpclear()
      else
         if (mnproc==1)
     &      write(lp,*) 'NO GRIDP information stored:  lgridp=',lgridp
      endif
      end subroutine


c --- -----------------------------------------------------------------
c --- gridp processing logic, moved from hycom in here
c --- -----------------------------------------------------------------
      subroutine gridp_process(rt,dtime,m,n)
      use mod_diagnostics
      use mod_year_info, only: year_info, year_day,refyear
      implicit none
      type(year_info) , intent(in) :: rt
      integer         , intent(in) :: m,n
      real*8          , intent(in) :: dtime
      type(year_info) :: rttmp
      include 'common_blocks.h'
c
c --- Grid points
      if (lgridp.and.rt%iss < nint(baclin)) then
         call gridpstore(n,m,rt)
      endif
c
      if (nstep>nstep1+1) then
         call year_day(dtime+baclin/86400,refyear,rttmp,yrflag)
c ---    Gridp dump at diagnostic times
         if (lgridp.and.loutput) then
            if (mnproc==1) write(lp,*)'Gridp output on diagnostic day'
            call gridpdump(n,rt)
            call gridpclear()
c ---    Always dump at end of month
         else if (lgridp .and.rttmp%imm /= rt%imm) then
            if (mnproc==1) write(lp,*)'Gridp output on month changeover'
            call gridpdump(n,rt)
            call gridpclear()
         endif
      endif
      end subroutine


c --- -----------------------------------------------------------------
C --- This routine reads the input files, and sets a mask to determine
C --- which points are dumped later on. The array "gp" (see above) gets
C --- initialized with grid indices corresponding to lon/lat positions
c --- -----------------------------------------------------------------
      subroutine gridpinp()
      implicit none
      integer i,j,m,n,nrgps
      real mind,maxd,minlat,maxlat,minlon,maxlon
      integer igpskip, jgpskip, ibnd,k
      logical ex
      logical :: gpmask(itdm,jtdm)
      real, dimension(itdm,jtdm) :: gdepths,modlon,modlat
      integer, dimension(itdm,jtdm) :: stindex
      character*3 cthread
      include 'common_blocks.h'

      if (mnproc==1) write(lp,*) 'GRIDP activated'
      call xcaget(gdepths,depths,0)
      call xcaget(modlon,plon,0)
      call xcaget(modlat,plat,0)
C$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jtdm
      do i=1,itdm
         gpmask(i,j)=.false.
      end do
      end do
C$OMP END PARALLEL DO
c
c --- Read infile_gp.in 
      inquire(file='infile_gp.in',exist=ex)
      if (.not.ex)  then
         if (mnproc==1) 
     &   write(lp,*) 'infile_gp.in file does not exist'
         call xcstop('(gridpinp)')
         stop '(gridpinp)'
      end if
      open(11,file='infile_gp.in')
      read(11,*)mind    ! Minimum depth
      read(11,*)maxd    ! Max depth
      read(11,*)minlat  ! Minimum lat (lonlat box)
      read(11,*)maxlat  ! Maximum Lat (lonlat box)
      read(11,*)minlon  ! Minimum lon (lonlat box)
      read(11,*)maxlon  ! inximum lon (lonlat box)
      read(11,*)igpskip ! Index skip
      read(11,*)jgpskip ! Index skip
c
c --- counting number of grid points for grid
C$OMP PARALLEL DO PRIVATE(i,j)
      do j=2,jtdm-1,jgpskip
      do i=2,itdm-1,igpskip
         if     (mind  <= gdepths(i,j).and.gdepths(i,j)<=maxd 
     &     .and.minlon <= modlon(i,j) .and.modlon(i,j) <=maxlon
     &     .and.minlat <= modlat(i,j) .and.modlat(i,j) <=maxlat ) then
           gpmask(i,j)=.true.
         end if
      enddo
      enddo
C$OMP END PARALLEL DO
      if (mnproc==1)print *,'Number of regular points : ',count(gpmask)
c
c --- From now on we read specified stations. Input is now i,j points
c --- (NB: used to  be lon lat values)
      do m=1,1000
         read(11,*,end=100,err=100) i,j
         if (mnproc==1)  print '(a,2i5)','gridpinp: ',i,j
         if (1 < i .and. i < itdm .and. 1 < j .and. j < jtdm) then
            gpmask(i  ,j  )=gdepths(i  ,j  ) > 0.0
            gpmask(i+1,j  )=gdepths(i+1,j  ) > 0.0
            gpmask(i  ,j+1)=gdepths(i  ,j+1) > 0.0
            gpmask(i+1,j+1)=gdepths(i+1,j+1) > 0.0
         else
            if (mnproc==1) write(lp,*) 'Station NOT on grid !!'
            call xcstop ('(gridpinp)')
            stop '(gridpinp)'
         endif
      enddo
100   close(11)

c --- Stop if there are no points:
      if (all(.not.gpmask)) then
         if (mnproc==1) print *,'Error: There are no gp points!'
         call xcstop('(gridpinp)')
         stop '(gridpinp)'
      end if

      gpdim=0
C$OMP PARALLEL DO PRIVATE (i,j) REDUCTION(+:gpdim)
      do j=1,jtdm
      do i=1,itdm
         stindex(i,j)=0
         if (gpmask(i,j)) then 
            gpdim=gpdim+1
            stindex(i,j)=gpdim
         end if
      enddo
      enddo
      gpdim_global= gpdim
C$OMP END PARALLEL DO
      if (mnproc==1) then
      write(lp,*) 'gridpinp: Number of grid points is          : ',
     &   gpdim
      write(lp,*) 'gridpinp: Size of each file per year is (Mb): ',
     &   (kdm*5+3)*24*4*365*1.0E-06
      write(lp,*) 'gridpinp: Total archive size per year (Gb)) : ',
     &   (kdm*5+3)*24*4*365*1.0E-09*gpdim
      end if


c --- Now calculate for tiles:
      gpdim=0
C$OMP PARALLEL DO PRIVATE (i,j) REDUCTION(+:gpdim)
      do j=1,jj
      do i=1,ii
         if (gpmask(i0+i,j0+j)) then
            gpdim=gpdim+1
         end if
      end do
      end do
C$OMP END PARALLEL DO
      if (gpdim>0) allocate(gp(gpdim))

      call xcsync(flush_lp) 
      write(lp,'(a,2i5,a,i5)') 
     &' gridpinp: Number of grid points for tile ', mproc,nproc,
     &' is ', gpdim
      call xcsync(flush_lp) 

c --- define local indexes for gpdata
      if (gpdim>0) then
         m=0
         do j=1,jj
         do i=1,ii
            if (gpmask(i0+i,j0+j)) then
               m=m+1
               gp(m)%i=i
               gp(m)%j=j
               gp(m)%station_index=stindex(i+i0,j+j0)
            endif
         enddo
         enddo
      end if
c
c --- Save to tecplot file - one for each tile
c --- write(cthread,'(i3.3)') mnproc
c --- open(10,file='gridploc'//cthread//'.tec')
c --- Save to tecplot file - one file
      do k=1,ijpr
c ---    one MPI task at a time
         if (mnproc==k) then
            if (k==1) then
               open(10,file='gridploc.tec', status='replace')
            else
               open(10,file='gridploc.tec',position='append')
            end if
            do m=1,gpdim
               i=gp(m)%i
               j=gp(m)%j
               write(10,'(2f10.5,f14.2,3i4,i8)')plon(i,j),plat(i,j),
     &            depths(i,j),i+i0,j+j0,k,gp(m)%station_index
            enddo
            close(10)
         end if
c ---    Implies barrier
         call xcsync(flush_lp)
      enddo
c
c --  lgridp not switched off for now, but we might deactivate it later for "empty" tiles
      if (.not.lgridp) then
         write(lp,'(a,i5)') 'GRIDP DE-activated for tile ',ijqr
      end if
      end subroutine gridpinp


c --- -----------------------------------------------------------------
C --- Routine Allocates global (in module) vars. Only needed once
C --- in new version.
C --- Allocates:  gpstore, itime
C --- KAL: gptdim is now based on 31 days of month
c --- -----------------------------------------------------------------
      subroutine gridpalloc()
      use mod_diagnostics
c      use mod_year_info
      implicit none
      allocate(gpstore(gptdim,gpdim))
      allocate(itime(gptdim))
      itime=-1
      end subroutine gridpalloc

c --- -----------------------------------------------------------------
C --- Clears itime array - it couldnt be simpler
c --- -----------------------------------------------------------------
      subroutine gridpclear()
      use mod_diagnostics
c      use mod_year_info
      implicit none
      itime=-1
      end subroutine gridpclear


c --- -----------------------------------------------------------------
C --- Get gp data from hycom fields, once an hour
C --- Put them into the gp variables at correct time record
c --- -----------------------------------------------------------------
      subroutine gridpstore(n,m,rt)
      use mod_year_info, only: year_info
      use mod_forcing_nersc
#if defined (ICE)
      use mod_common_ice
#endif
      implicit none
      integer             , intent(in) :: n  ! New time indexes
      integer             , intent(in) :: m  ! old time indexes
      type(year_info)     , intent(in) :: rt ! Time variable
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   tuwind,tvwind, tslp
      character(len=17) tag17
      real fli,flj
      real up,vp,ud,vd,theta_up,theta_vp
      real dlon,dlat
      integer ilon,ilat
      character(len=10)ctmp
      integer i,j,k,kn,igp,km,l, gpit
      include 'common_blocks.h'

c --- Time record. day 1 of month, 00 hours = 1
c --- Time record. day 1 of month, 02 hours = 2
c --- etc .. etc ..
      gpit=rt%idm*24 + rt%ihh +1
c
c --- Put integer description of time into itime
C     ctmp(1:4)=rt%cyy
C     ctmp(5:6)=rt%cmm
C     ctmp(7:8)=rt%cdm
C     ctmp(9:10)=rt%chh
C     read(ctmp,'(i10.10)')itime(gpit)
      itime(gpit) = rt%iyy * 10**6 + rt%imm * 10**4 + rt%idm * 10**2 +
     &              rt%ihh 
c
c --- Diagnostics
      if (mod(gpit,12)==0) print '(a,i5,a,i10,a,4i4,a,i5)',
     &   'gridpstore: gpit=',gpit,
     &   ' itime=',itime(gpit),
     &   ' rt= ',rt%iyy, rt%imm, rt%idm,  rt%ihh,
     &   ' tile= ',ijqr
c 
c --- This routine can now be called when rt%iss/=0. The closest
c --- hour is chosen.
      if (rt%iss>baclin) then
         if (mnproc==1) 
     &   print *,'(gridpstore should only be called on the hour or'//
     &           ' immediately after the hour has passed)'
         call xcstop('gripstore')
      end if
c
c --- Get slp and winds
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
      do i=1,ii
         tuwind(i,j)=uwind(i,j,l0)*w0+uwind(i,j,l1)*w1
     &              +uwind(i,j,l2)*w2+uwind(i,j,l3)*w3
         tvwind(i,j)=vwind(i,j,l0)*w0+vwind(i,j,l1)*w1
     &              +vwind(i,j,l2)*w2+vwind(i,j,l3)*w3
         tslp  (i,j)=slp  (i,j,l0)*w0+slp  (i,j,l1)*w1
     &              +slp  (i,j,l2)*w2+slp  (i,j,l3)*w3
      end do
      end do
c
c --- Interpolation time
!$OMP PARALLEL DO PRIVATE(igp,i,j,dlon,dlat,theta_up,theta_vp,
!$OMP& k,kn,ud,vd,km) SHARED(gpit,m,n) 
!$OMP& SCHEDULE(STATIC,1)
      do igp=1,gpdim
         i=gp(igp)%i
         j=gp(igp)%j
         dlon=plon(i+1,j)-plon(i-1,j)
         dlat=plat(i+1,j)-plat(i-1,j)
c ---    Grid directions -- "i" rel lon/lat
         if(dlon.LT.180.)dlon=360.0+dlon
         if(dlon.GT.180.)dlon=dlon-360.0
C ---    KAL -- From HYCOM 2.1.34 - radian is pi/180
         theta_up=atan2(dlat,dlon*cos(radian*.5*
     &                               (plat(i-1,j)+plat(i+1,j))))

         ! Grid directions -- "j" rel lon/lat
         dlon=plon(i,j+1)-plon(i,j-1)
         dlat=plat(i,j+1)-plat(i,j-1)
         if(dlon.LT.180.)dlon=360.0+dlon
         if(dlon.GT.180.)dlon=dlon-360.0
C ---    KAL -- From HYCOM 2.1.34 - radian is pi/180
         theta_vp = atan2(dlat,dlon*cos(radian*.5*
     &                                 (plat(i,j-1)+plat(i,j+1))))

         do k=1,kdm
            gp(igp)%z(k)=p(i,j,k+1)/onem 
            gp(igp)%u(k)=(u(i,j,k,n)+ubavg(i,j,n))*cos(theta_up)+
     &                   (v(i,j,k,n)+vbavg(i,j,n))*cos(theta_vp)
            gp(igp)%v(k)=(u(i,j,k,n)+ubavg(i,j,n))*sin(theta_up)+
     &                   (v(i,j,k,n)+vbavg(i,j,n))*sin(theta_vp)
            gp(igp)%t(k)=temp(i,j,k,n) 
            gp(igp)%s(k)=saln(i,j,k,n)
         enddo
c
c ---    Wind rotation introduced AFTER NWAG and WANE projects by GE
         gp(igp)%uwind   =tuwind(i,j)*cos(theta_up)
     &                   +tvwind(i,j)*cos(theta_vp)
         gp(igp)%vwind   =tuwind(i,j)*sin(theta_up)
     &                   +tvwind(i,j)*sin(theta_vp)
         gp(igp)%sealevel=srfhgt(i,j)/g
#if defined (ICE)
         gp(igp)%hice    = hicem(i,j)
         gp(igp)%fice    = ficem(i,j)
         gp(igp)%uice    = iceu(i,j)*cos(theta_up)+ 
     &                     icev(i,j)*cos(theta_vp)
         gp(igp)%vice    = iceu(i,j)*sin(theta_up)+ 
     &                     icev(i,j)*sin(theta_vp)
#endif

         gpstore(gpit,igp)=gp(igp)  
      enddo
C$OMP END PARALLEL DO
      end subroutine gridpstore





c --- -----------------------------------------------------------------
C --- Put grid point data on file. One directory per month, one 
c --- file per point.
c --- TODO: Netcdf output?
c --- -----------------------------------------------------------------
      subroutine gridpdump(n,rt)
      use mod_xc
      use mod_year_info, only: year_info
      use mod_hycom_nersc
      use netcdf
      implicit none
      integer, intent(in) :: n
      type(year_info), intent(in) :: rt

      type(year_info) rtt
      character(len=17) tag17
      real fli,flj
      real up,vp,ud,vd,theta_up,theta_vp
      real dlon,dlat
      integer ilon,ilat,it
      integer ilen
      integer m,i,j
      character(len=10)ctmp
      character(len=12) gpdir
      character(len=31) gpfile
      character(len=80) gpfilenc
      character(len=10) tag10a,tag10b
      character(len=6)  csi
      logical ex, exnc
      integer :: firstrec, lastrec, si
      real*4 :: r4lon, r4lat, r4deep
      type(grid_points) :: tgp ! temporary gp 
      type(grid_points) :: tgpnc(1:31*24) ! temporary gp 
      logical,save :: lfirst=.true.
      include 'common_blocks.h'
c --- Fetch first an last "real" value calculated
      firstrec=gptdim+1
      lastrec=1
      do i=1,gptdim
         if(itime(i)>0) then
            firstrec=min(i,firstrec)
            lastrec =max(i,lastrec)
         end if
      end do

c --- Put itime into tag10
      write(tag10a,'(i10.10)')itime(lastrec)
      write(tag10b,'(i10.10)')itime(lastrec) ! This should work now
      if (mnproc==1) then
         print *,'itime(lastrec)  = ',itime(lastrec)
      end if
c
c --- Name of directory to store GP data
      gpdir(1:3)=rungen
      gpdir(4:5)='GP'
      gpdir(6:9)=tag10b(1:4)
      gpdir(10:10)='_'
      gpdir(11:12)=tag10b(5:6)
      inquire(file=gpdir,exist=ex)
      if (.not.ex.and.mnproc==1) then
         write(*,'(a,a)')'GP: creating the directory ./', gpdir
         call mkdir_wrap( trim(gpdir)//char(0) )
      endif
      call xcsync(flush_lp)
c
c --- dump legend in text file in GP dir
      if (mnproc==1) then
         open(10,file=gpdir//'/'//rungen//'gprecord.asc', 
     &           status='unknown')
         write(10,'(a12,a3)') 'var name    ','dim'
         write(10,'(a12,i3)') 'itime       ',1
         write(10,'(a12,i3)') 'i           ',1
         write(10,'(a12,i3)') 'j           ',1
         write(10,'(a12,i3)') 'lon         ',1
         write(10,'(a12,i3)') 'lat         ',1
         write(10,'(a12,i3)') 'depth       ',1
         write(10,'(a12,i3)') 'intf        ',kdm
         write(10,'(a12,i3)') 'ut          ',kdm
         write(10,'(a12,i3)') 'vt          ',kdm
         write(10,'(a12,i3)') 'temp        ',kdm
         write(10,'(a12,i3)') 'saln        ',kdm
         write(10,'(a12,i3)') 'ssh         ',1
         write(10,'(a12,i3)') 'uwind       ',1
         write(10,'(a12,i3)') 'vwind       ',1
#if defined (ICE)
         write(10,'(a12,i3)') 'fice        ',1
         write(10,'(a12,i3)') 'hice        ',1
         write(10,'(a12,i3)') 'uice        ',1
         write(10,'(a12,i3)') 'vice        ',1
#endif
         close (10)
      end if

c --- Get IO-length for one record
      inquire(iolength=ilen)itime(1),tgp%i,tgp%j, 
     &   r4lon, r4lat, r4deep,    ! lon/lat/depth
     &   tgp%z,tgp%u,tgp%v,tgp%t,tgp%s,
     &   tgp%sealevel, tgp%uwind, tgp%vwind
#if defined (ICE)
     &   ,tgp%fice,tgp%hice,tgp%uice,tgp%vice
#endif

c --- Go through the points and dump data 
      do m=1,gpdim
         i=gp(m)%i
         j=gp(m)%j
         if (plon(i,j) < 0.0) then
            tag17(1:1)='-'
         else
            tag17(1:1)='+'
         endif
c
c ---    Write lon/lat into file name
         ilon=int(abs(plon(i,j)))
         ilat=int(abs(plat(i,j)))
         write(tag17(2:4),'(I3.3)')ilon
         tag17(5:5)='.'
         write(tag17(6:8),'(I3.3)')
!LB     &      nint((abs(plon(i,j))-float(ilon))*1000.0)
     &      min(nint((abs(plon(i,j))-float(ilon))*1000.0),999)
         tag17(9:9)='x'
         if (plat(i,j) < 0.0) then
            tag17(10:10)='-'
         else
            tag17(10:10)='+'
         endif
         write(tag17(11:13),'(I3.3)')ilat
         tag17(14:14)='.'
         !LB Same for lat, avoids stars in file names
         write(tag17(15:17),'(I3.3)')min(nint((abs(plat(i,j))
     &        -float(ilat))*1000.0),999)

c ---    Create final file name
         gpfile( 1: 3)='gp_'
         gpfile( 4: 7)=tag10b(1:4)
         gpfile( 8: 8)='_'
         gpfile( 9:10)=tag10b(5:6)
         gpfile(11:11)='_'
         gpfile(12:28)=tag17
         gpfile(29:31)='.uf'
c
c ---    Alternative (and much simpler ) scheme
         write(csi,'(i6.6)') si
c ---    gpfile='gp_'//csi//'.nc'
c
c ---    Diagnostics for first point on  each tile
         if (m==1) then
            print '(a,2i4,a,i10)',
     &         'gpdump: tile ',mproc,nproc,
     &         ' First  time : ',itime(firstrec)
            print '(a,2i4,a,i10)',
     &         'gpdump: tile ',mproc,nproc,
     &         ' Last  record and time : ', itime(lastrec)
         end if
         if (lfirst) then
            print '(a,2i3,a)',
     &         'gpdump: tile ',mproc,nproc,
     &         ' Dumping gridp data to '//gpdir//'/'//gpfile
         end if
c
         inquire(exist=ex,file=gpdir//'/'//gpfile)
         open(10,file=gpdir//'/'//gpfile,form='unformatted',
     &           access='direct',recl=ilen,status='unknown')
        
c ---    Fill with itime < 0 if new file
         if (.not.ex) then
            do it=1,firstrec-1
               write(10,rec=it) -1
            end do
         end if
c
         do it=firstrec,lastrec
            if (itime(it) >  0) then
               tgp=gpstore(it,m)
               r4lon =plon(tgp%i,tgp%j)
               r4lat =plat(tgp%i,tgp%j)
               r4deep=depths(tgp%i,tgp%j)
               write(10,rec=it) itime(it),tgp%i,tgp%j, 
     &            r4lon, r4lat, r4deep,    
     &            tgp%z,tgp%u,tgp%v,tgp%t,tgp%s,
     &            tgp%sealevel, tgp%uwind, tgp%vwind
#if defined (ICE)
     &            ,tgp%fice,tgp%hice,tgp%uice,tgp%vice
#endif
            else
               if (mnproc==1) then
                  print *,'gridpdump: record number= ',it,
     &                 ' not written, time= ',itime(it),
     &                 ' for tile ',ijqr
               end if
            endif
         enddo
         close(10)
      enddo
      if (mnproc == 1) then
         print *,'gpdump: done...'
         print *
      end if
      lfirst=.false.
      end subroutine gridpdump

      end module mod_gridp
