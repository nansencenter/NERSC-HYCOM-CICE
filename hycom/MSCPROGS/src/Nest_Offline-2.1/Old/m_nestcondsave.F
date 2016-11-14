      module m_nestcondsave
! The stuff below sets up MPI real types... see mod_xc_mp.inc


      contains
C --- Save nesting conditions from this model to inner models,
C --- Specified in infile read below (nestname)
      subroutine nestcondsave(rt,n,nestdto,iens)
      use mod_xc 
      use mod_year_info
      use mod_nesting
#if defined (ICE) && defined (ICE_NEST)
      use mod_common_ice
#endif

      use m_filenesting
      use m_rotate
      use m_initconfmap
      use m_oldtonew
      use m_newtoold
      use m_pivotp
      use m_nearestpoint
      use m_bilincoeff
      implicit none
      type(year_info), intent(in)    :: rt
      integer,         intent(in)    :: n
      integer,         intent(in)    :: nestdto
      integer,         intent(in)    :: iens

      integer, parameter :: max_nrgrids=100
      integer nrgrids
      character(len=80) nestname(max_nrgrids)
      character(len=100) nestingfile
      logical ex,ass,lerr,checknest
      integer irec, irec_offset
      character(len=11) tag7
      character(len=15) char15
      character(len=5) char5
      character(len=3) char3,char3_2
      character(len=2) cgrid,cproc
      real lon,lat
      real ba1,ba2,ba3,ba4
      real lat_n,lon_n
      integer ipiv,jpiv,gipiv,gjpiv
      integer i,j,k,m,l,klevel,tmpint

      logical :: isvelocity
      integer :: ivar, iedge
    
      integer :: inner_idm, inner_jdm, inner_iidm, inner_jjdm,
     &   inner_inest
      integer :: first_i_inner, last_i_inner
      integer :: first_j_inner, last_j_inner

      real, allocatable, dimension(:,:) :: 
     &   inner_tmp1,inner_tmp2,inner_tmp3
      real, allocatable, dimension(:,:) :: 
     &   inner_depths,inner_lon,inner_lat
      real*8, allocatable :: inner_io(:,:)
      real*4, allocatable :: inner_io2(:,:)
      integer, allocatable, dimension(:,:) :: 
     &   inner_ipiv, inner_jpiv
      real, allocatable, dimension(:,:) :: 
     &   inner_a1, inner_a2,
     &   inner_a3, inner_a4
      logical, allocatable, dimension(:,:) :: 
     &   inner_tile_mask

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   fldsum, fld1, fld2

      real, dimension(nx,ny) :: modlon, modlat, gdepths
      logical, dimension(nx,ny) :: gmsk

      integer :: ntecsave


#if defined (MPI)
      integer :: ierr
      include 'mpif.h'
#endif
      include 'common_blocks.h'
      include 'stmt_fns.h'

      ntecsave=1

      ! read dimensions and locations for internal grids   
      inquire(file='nesting.in',exist=ex)
      if (.not.ex) stop 'nestcond_save: nesting.in does not exist'
      open(10,file='nesting.in')
         do m=1,max_nrgrids
            read(10,'(a80)',end=999)nestname(m)
            if (mnproc==1)
     &      print *,'nestname read from file= ',trim(nestname(m)),
     &              ' for grid no ',m
         enddo
  999 nrgrids=m-1
      close(10)

      call xcaget(modlon,plon,0)
      call xcaget(modlat,plat,0)
      call xcaget(gdepths,depths,0)
      gmsk=gdepths>0.
      where (modlon<0.) modlon=modlon+360.

      ! Initialize conformal mapping
      call initconfmap


! Loop through the internal grids specified in nestname files (read above)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do m=1,nrgrids
      inquire(file=trim(nestname(m))//'nestloc.uf',exist=ex)
      if (ex) then

         ! Read nesting positions for this internal grid
         open(10,file=trim(nestname(m))//'nestloc.uf',
     &        form='unformatted')

         read(10)inner_idm,inner_jdm,inner_inest
         inner_iidm=inner_idm-inner_inest+1
         inner_jjdm=inner_jdm-inner_inest+1

         if (mnproc==1) write(lp,'(a,a,5i5)') trim(nestname(m)),
     &      'nestloc.uf  has dimensions: ', inner_idm,inner_jdm,
     &      inner_iidm,inner_jjdm,inner_inest

         allocate(inner_lon(inner_idm,inner_jdm))
         allocate(inner_lat(inner_idm,inner_jdm))
         allocate(inner_io (inner_idm,inner_jdm))
         allocate(inner_io2(inner_idm,inner_jdm))
         read(10) inner_io ; inner_lon=inner_io
         read(10) inner_io ; inner_lat=inner_io
         close(10)

         if (mnproc==1) write(lp,'(a,a)') 
     &      trim(nestname(m)),'nestloc.uf  is read'

         ! KAL - new scheme - saves space ...
         allocate(inner_tmp1  (inner_idm,inner_jdm))
         allocate(inner_tmp2  (inner_idm,inner_jdm))
         allocate(inner_tmp3  (inner_idm,inner_jdm))
         allocate(inner_depths(inner_idm,inner_jdm))




         ! Read the local grid
C --- tag7 actually has 11 chars to accomodate for huge grids in future
         if (inner_idm>999 .or. inner_jdm > 999) then
            write(tag7,'(i5.5,a,i5.5)')inner_idm,'x',inner_jdm
         else
            write(tag7,'(i3.3,a,i3.3)')inner_idm,'x',inner_jdm
         end if
         inquire(file=trim(nestname(m))//'ndepths'//trim(tag7)//'.uf',
     &      exist=ex)
         if (.not.ex) then
            if (mnproc==1) write(lp,'(a)')
     &      'nesting depths file for local grid does not exist:',
     &              'ndepths'//trim(tag7)//'.uf'
            call xcstop('nestcondsave')
            stop '(nestcondsave)'
         else
            open(10,file=trim(nestname(m))//'ndepths'//trim(tag7)//'.uf'
     &              ,form='unformatted')
               read(10) inner_io
               inner_depths=inner_io
            close(10)
         endif

         
         ! Extend local grid to boundary
         where (inner_depths(2,:)>0.) 
     &      inner_depths(1,:)=inner_depths(2,:)
         where (inner_depths(inner_idm-1,:)>0.) 
     &      inner_depths(inner_idm,:)=inner_depths(inner_idm-1,:)
         where (inner_depths(:,2)>0.) 
     &      inner_depths(:,1)=inner_depths(:,2)
         where (inner_depths(:,inner_jdm-1)>0.)
     &      inner_depths(:,inner_jdm)=inner_depths(:,inner_jdm-1)

Cdiag    if (mnproc==1) then
Cdiag    open(15,file='gridtest.tec')
Cdiag    WRITE(15,'(''TITLE= "Nested fields"'')')
Cdiag    write(15,'(a)')'VARIABLES="i" "j" "lon" "lat"'//
Cdiag&       '"depths[m]" '
Cdiag    WRITE(15,'(a,a,i4,a,i4)') 'ZONE T="test"',
Cdiag&         ' I=',inner_idm,' J=', inner_jdm, ' F=BLOCK'
Cdiag    WRITE(15,101)((i,i=1,inner_idm),j=1,inner_jdm)
Cdiag    WRITE(15,101)((j,i=1,inner_idm),j=1,inner_jdm)
Cdiag    WRITE(15,102)((inner_lon(i,j),
Cdiag&                  i=1,inner_idm),j=1,inner_jdm)
Cdiag    WRITE(15,102)((inner_lat(i,j),
Cdiag&                  i=1,inner_idm),j=1,inner_jdm)
Cdiag    WRITE(15,102)((inner_depths(i,j),
Cdiag&                  i=1,inner_idm),j=1,inner_jdm)
Cdiag    close(15)
Cdiag    end if


         !if (mnproc==1) call oldtonew(-4.71,106.22,lat_n,lon_n)
         !call xcstop('nestcondsave test')
         !stop 'nestcondsave test'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign points on inner grid to points on the global grid of this model (if possible)
! NB: ipiv,jpiv refers to points on the Global grid
         allocate(inner_ipiv(inner_idm,inner_jdm))
         allocate(inner_jpiv(inner_idm,inner_jdm))
         allocate(inner_a1  (inner_idm,inner_jdm))
         allocate(inner_a2  (inner_idm,inner_jdm))
         allocate(inner_a3  (inner_idm,inner_jdm))
         allocate(inner_a4  (inner_idm,inner_jdm))
         inquire(file=trim(nestname(m))//'pivots'//trim(tag7)//'.uf',
     &      exist=ex)
         ! In testing fase, always redo the ipivs...
         ex=.false.
         if (ex) then ! read from file
            open(10,file=trim(nestname(m))//'pivots'//trim(tag7)//'.uf',
     &           form='unformatted')
            read(10)inner_ipiv,inner_jpiv
            read(10)inner_a1  ,inner_a2  
            read(10)inner_a3  ,inner_a4  
            close(10)
         else
            if (mnproc==1) print *,'Assigning points on first pass...'
            inner_ipiv=-1
            inner_jpiv=-1
            inner_a1=0.
            inner_a2=0.
            inner_a3=0.
            inner_a4=0.
            do j=1,inner_jdm
            do i=1,inner_idm
            if (inner_depths(i,j)>0.) then 

               ! Get pivot points of internal grid positions
               lon=inner_lon(i,j)
               lat=inner_lat(i,j)
               call oldtonew(lat,lon,lat_n,lon_n)
               call pivotp(lon_n,lat_n,ipiv,jpiv)



               ! Check if this point is on global grid
               if (ipiv>=1.and.ipiv<nx.and.jpiv>=1.and.jpiv<ny) then

C                 call newtoold(lat_n,lon_n,lat,lon)
C                 print '(a,3i5,9f10.2)','lonlat:',mnproc,ipiv,jpiv,
C    &               inner_lon(i,j),inner_lat(i,j),
C    &               modlon(ipiv,jpiv),modlat(ipiv,jpiv),
C    &               lon,lat,lon_n,lat_n
C                 call xcstop('nestcondsave test')
C                 stop 'nestcondsave test'

                  ! Check if surrounded by sea cells on the global  grid
                  if (gmsk(ipiv  ,jpiv  ).and.gmsk(ipiv+1,jpiv  ).and.
     &                gmsk(ipiv  ,jpiv+1).and.gmsk(ipiv+1,jpiv+1))
     &               then
                     call bilincoeff(modlon,modlat,nx,ny,lon,lat,
     &                  ipiv,jpiv,ba1,ba2,ba3,ba4)

                  ! Otherwise search for nearest point
                  else
                     call nearestpoint(modlon,modlat,nx,ny,lon,lat,
     &                                 ipiv,jpiv,ba1,ba2,ba3,ba4,
     &                                 gmsk,ass)
                  endif
   
                  inner_ipiv(i,j)=ipiv
                  inner_jpiv(i,j)=jpiv
                  inner_a1  (i,j)=ba1
                  inner_a2  (i,j)=ba2
                  inner_a3  (i,j)=ba3
                  inner_a4  (i,j)=ba4
C                 print '(3i5,6f10.2)',mnproc,ipiv,jpiv,
C    &                  inner_lon(i,j),inner_lat(i,j),
C    &                  inner_a1(i,j),inner_a2(i,j),
C    &                  inner_a3(i,j),inner_a4(i,j)
               end if
            end if
            end do
            end do

            ! Store for later
            if (mnproc==1) then
               write(tag7,'(i3.3,a,i3.3)')nx,'x',ny
               print *,'Assigning points - done'
               print *, 'Storing to '//trim(nestname(m))
     &                  //'pivots'//trim(tag7)//'.uf'
               open(10,file=trim(nestname(m))//'pivots'//trim(tag7)
     &              //'.uf',  form='unformatted')
               write(10)inner_ipiv,inner_jpiv
               write(10)inner_a1  ,inner_a2  
               write(10)inner_a3  ,inner_a4  
               close(10)

               ! For tecplot
               open(10,file=trim(nestname(m))//'pivots'//trim(tag7)
     &              //'.tec', form='formatted',status='replace')
               do j=1,inner_jdm
               do i=1,inner_idm
               if (inner_depths(i,j)>0) then
                  write(10,'(2i6,2f14.4,2i6,4f14.4)') i,j,
     &               inner_lon(i,j),inner_lat(i,j),
     &               inner_ipiv(i,j),inner_jpiv(i,j),
     &               inner_a1(i,j), inner_a2(i,j),
     &               inner_a3(i,j), inner_a4(i,j)
               end if
               end do
               end do
               close(10)
            end if
         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
! Local mask for this tile/process 
         allocate(inner_tile_mask(inner_idm,inner_jdm))
         inner_tile_mask=.false.
         do j=1,inner_jdm
         do i=1,inner_idm
         if (inner_depths(i,j)>0.) then 

            ipiv=inner_ipiv(i,j)
            jpiv=inner_jpiv(i,j)

            ! Checks if pivot points are on this tile
            inner_tile_mask(i,j)= ipiv > i0 .and. ipiv <= i0+ii .and.
     &                            jpiv > j0 .and. jpiv <= j0+jj

            ! Cheks if in nesting zone of inner model
            inner_tile_mask(i,j)=inner_tile_mask(i,j).and.
     &         min(i,inner_idm-i,j,inner_jdm-j)<=inner_inest


         end if
         end do
         end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
! Record offsets on nesting file

         if (nestdto < 24) then
            irec_offset=(iens-1)*24/nestdto+rt%ihh/nestdto+1
         else
            irec_offset=iens
         endif
         nestingfile=filenesting(rt,trim(nestname(m)))

         !nvar=4*nz+2         ! For loop (velocity treated as one)
         !nvar_tot=5*nz+3     ! For offset, velocity treated as two

         irec=(irec_offset-1)*nvar_tot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!Open files for IO on node 1
         if (mnproc==1) then
            write(lp,'(a,a,a)') 'saving records ', 
     &         ' to files=',trim(nestingfile)//'_xx'
            inquire(iolength=j)inner_io2(1:inner_inest,1:inner_jdm)
            open(10,file=trim(nestingfile)//'_i1',form='unformatted',
     &           access='direct',recl=j)
            inquire(iolength=j)inner_io2(1:inner_inest,1:inner_jdm)
            open(11,file=trim(nestingfile)//'_ii',form='unformatted',
     &           access='direct',recl=j)
            inquire(iolength=j)inner_io2(1:inner_idm,1:inner_inest)
            open(12,file=trim(nestingfile)//'_j1',form='unformatted',
     &           access='direct',recl=j)
            inquire(iolength=j)inner_io2(1:inner_idm,1:inner_inest)
            open(13,file=trim(nestingfile)//'_jj',form='unformatted',
     &           access='direct',recl=j)
            open(16,file=trim(nestingfile)//'.hdr',form='formatted',
     &           access='direct',recl=100,status='unknown')
            open(15,file=trim(nestname(m))//'testfld.tec',
     &           form='formatted',status='replace')
         end if
            


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!Variable loop ...
         fldsum=0.
         do ivar=1,nvar! Ivar needs to be rubbed out

            if     (ivar>=1 .and. ivar <=nz) then ! Salinity
               klevel=ivar
               fld1=saln(:,:,klevel,n)
               isvelocity=.false.
               write(char15,'(a11,i4)') 'saln level ',klevel
               write(char5,'(a3,i2.2)') 'SAL',klevel
               write(char3,'(a3)') 'SAL'
            elseif (ivar>=nz+1 .and. ivar <=2*nz) then ! Temperature
               klevel=ivar-nz
               fld1=temp(:,:,klevel,n)
               isvelocity=.false.
               write(char15,'(a11,i4)') 'temp level ',klevel
               write(char5,'(a3,i2.2)') 'TEM',klevel
               write(char3,'(a3)') 'TEM'
            elseif (ivar>=2*nz+1 .and. ivar <=3*nz) then ! Interface values
               klevel=ivar-2*nz
               fldsum=fldsum+dp(:,:,klevel,n)
               fld1=fldsum
               isvelocity=.false.
               write(char15,'(a11,i4)') 'intf level ',klevel
               write(char5,'(a3,i2.2)') 'INT',klevel
               write(char3,'(a3)') 'INT'
            elseif (ivar>=3*nz+1 .and. ivar <=4*nz) then ! Velocities
               klevel=ivar-3*nz
               fld1=u(:,:,klevel,n)
               fld2=v(:,:,klevel,n)
               isvelocity=.true.
               call rotate(fld1,fld2,
     &                     plat,plon,idm+2*nbdy,jdm+2*nbdy,'m2l')
               write(char15,'(a11,i4)') 'Vel  level ',klevel
               write(char5,'(a3,i2.2)') 'VEL',klevel
               write(char3,'(a3)') 'UT '
               write(char3_2,'(a3)') 'VT '
            elseif (ivar==4*nz+1) then ! Baro Pressure
               fld1=pbavg(:,:,n)
               isvelocity=.false.
               char15='pbavg'
               write(char5,'(a3,i2.2)') 'PBA',1
               write(char3,'(a3)') 'PB '
               klevel=0
            elseif (ivar==4*nz+2) then ! Baro Velocities
               fld1=ubavg(:,:,n)
               fld2=vbavg(:,:,n)
               call rotate(fld1,fld2,
     &                     plat,plon,idm+2*nbdy,jdm+2*nbdy,'m2l')
               isvelocity=.true.
               char15='Bar vel'
               write(char5,'(a3,i2.2)') 'VBA',1
               write(char3,'(a3)') 'UB '
               write(char3_2,'(a3)') 'VB '
               klevel=0
#if defined (ICE) && defined (ICE_NEST)
            elseif (ivar==4*nz+3) then  ! Ice velocities
               fld1=iceU(:,:)
               fld2=iceV(:,:)
               call rotate(fld1,fld2,
     &                     plat,plon,idm+2*nbdy,jdm+2*nbdy,'m2l')
               isvelocity=.true.
               char15='ice vel'
               write(char5,'(a3,i2.2)') 'VI ',1
               write(char3,'(a3)') 'UI '
               write(char3_2,'(a3)') 'VI '
               klevel=0
            elseif (ivar==4*nz+4) then  ! Ice Thickness
               fld1=hicem(:,:)
               isvelocity=.false.
               char15='hicem'
               write(char5,'(a3,i2.2)') 'HI ',1
               write(char3,'(a3)') 'HI '
               klevel=0
            elseif (ivar==4*nz+5) then  ! Ice Concentration
               fld1=ficem(:,:)
               isvelocity=.false.
               char15='hicem'
               write(char5,'(a3,i2.2)') 'FI ',1
               write(char3,'(a3)') 'FI '
               klevel=0
#endif
            else
               if (mnproc==1) then
                  print *,'You shoudnt be here ...'
                  print *,irec
               end if
               call xcstop('nestcondsave')
               stop '(nestcondsave)'
            end if


            ! Other variables could be output like this, examples are
            ! ICE variables, ICESTATE variables, tracer variables etc
            ! etc ... The main layout is here now, all that remains 
            ! is to code it in a slightly more flexible way....
            ! Note that input to inner model relies on the header file, 
            ! which makes input very flexible...

            inner_tmp1=0.
            inner_tmp2=0.


! Do actual interpolation using what we calculated earlier 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!$OMP PARALLEL DO PRIVATE(i,j,lon,lat,lon_n,lat_n,ipiv,jpiv, 
!$OMP&                    ba1,ba2,ba3,ba4,ass,k)
            do j=1,inner_jdm
            do i=1,inner_idm
            if (inner_tile_mask(i,j)) then 
               ipiv=inner_ipiv(i,j)-i0 ! On local grid
               jpiv=inner_jpiv(i,j)-j0 ! On local grid
               ba1 =inner_a1  (i,j)
               ba2 =inner_a2  (i,j)
               ba3 =inner_a3  (i,j)
               ba4 =inner_a4  (i,j)


               ! We should have a representation of this point now...
               inner_tmp1(i,j)=ba1*fld1(ipiv  ,jpiv  )     
     &                        +ba2*fld1(ipiv+1,jpiv  )
     &                        +ba3*fld1(ipiv+1,jpiv+1) 
     &                        +ba4*fld1(ipiv  ,jpiv+1)

Cdiag          if (char5(1:3)=='SAL' .and. klevel==1) then
Cdiag             print '(4i5,5f14.2,i5)',
Cdiag&               i,j,ipiv,jpiv,fld1(ipiv,jpiv),depths(ipiv,jpiv),
Cdiag&               ba1,saln(i,j,klevel,n),inner_tmp1(i,j),mnproc
Cdiag          end if

               if (isvelocity)  then
                  inner_tmp2(i,j)=ba1*fld2(ipiv  ,jpiv  )     
     &                           +ba2*fld2(ipiv+1,jpiv  )
     &                           +ba3*fld2(ipiv+1,jpiv+1) 
     &                           +ba4*fld2(ipiv  ,jpiv+1)
               end if
            end if
            end do
            end do

            !Diagnose local mask
            if (char5(1:3)=='SAL' .and. klevel==1) then
            do j=1,inner_jdm
            do i=1,inner_idm
            if (inner_tile_mask(i,j)) then
               inner_tmp3(i,j)=1.
            else
               inner_tmp3(i,j)=0.
            end if
            end do
            end do
            write(cproc,'(i2.2)') mnproc
            open(98,file='masktest'//cproc//'.tec')
            WRITE(98,'(''TITLE= "Nested fields"'')')
            write(98,'(a)')'VARIABLES="i" "j" "lon" "lat"'//
     &       '"depths[m]" "mask" "fld"'
            WRITE(98,'(a,a,i4,a,i4)') 'ZONE T="test"',
     &            ' I=',inner_idm,' J=', inner_jdm, ' F=BLOCK'
            WRITE(98,101)((i,i=1,inner_idm),j=1,inner_jdm)
            WRITE(98,101)((j,i=1,inner_idm),j=1,inner_jdm)
            WRITE(98,102)((inner_lon(i,j),
     &                     i=1,inner_idm),j=1,inner_jdm)
            WRITE(98,102)((inner_lat(i,j),
     &                     i=1,inner_idm),j=1,inner_jdm)
            WRITE(98,102)((inner_depths(i,j),
     &                     i=1,inner_idm),j=1,inner_jdm)
            WRITE(98,102)((inner_tmp3(i,j),
     &                     i=1,inner_idm),j=1,inner_jdm)
            WRITE(98,102)((inner_tmp1(i,j),
     &                     i=1,inner_idm),j=1,inner_jdm)
            close(98)
            end if

            !print *,minval(inner_tmp1),maxval(inner_tmp1)


            ! Rotate if velociy
            if (isvelocity) then
               call rotate(inner_tmp1,inner_tmp2,inner_lat,inner_lon,
     &                     inner_idm,inner_jdm,'l2m')
            end if

            ! Reduce inner_tmp between tiles. Each node may have only part
            ! of the full inner_tmp. After MPI, it should be on all nodes
#if defined MPI 
Cdiag       write(lp,*) char5,klevel,mnproc
Cdiag       if (char5(1:3)=='SAL')  then
Cdiag          write(lp,'(a,i5,a,2e14.2,a,i5)') 'Writing record ',irec,
Cdiag&            ' '//char15,minval(inner_tmp1),maxval(inner_tmp1),
Cdiag&            ' node ',mnproc
Cdiag       end if
Cdiag       call flush(lp)

           call mpi_allreduce(
     &        inner_tmp1,inner_tmp3,inner_idm*inner_jdm,
     &        MTYPERNEST,MPI_SUM,MPI_COMM_WORLD,ierr)
            inner_tmp1=inner_tmp3
            if (isvelocity) then
CKAL           write(lp,'(a,i5,a,2e14.2,a,i5)') 'Writing record ',irec,
CKAL &            ' '//char15,minval(inner_tmp1),maxval(inner_tmp1),
CKAL &            ' node ',mnproc
              call mpi_allreduce(
     &           inner_tmp2,inner_tmp3,inner_idm*inner_jdm,
     &           MTYPERNEST,MPI_SUM,MPI_COMM_WORLD,ierr)
              inner_tmp2=inner_tmp3
            end if
#endif

            ! We just put the result on all nodes, make mnproc=1 do the
            ! actual output
            irec=irec+1
            if (mnproc==1) then
               inner_io2=inner_tmp1
               write(lp,'(a,i5,a,2e14.2)') 'Writing record ',irec,
     &            ' '//char15,minval(inner_io2),maxval(inner_io2)
               write(10,rec=irec) inner_io2(1:inner_inest,1:inner_jdm)
               write(11,rec=irec) 
     &            inner_io2(inner_iidm:inner_idm,1:inner_jdm)
               write(12,rec=irec) inner_io2(1:inner_idm,1:inner_inest)
               write(13,rec=irec) 
     &            inner_io2(1:inner_idm,inner_jjdm:inner_jdm)
               write(16,103,rec=irec)char3,klevel,irec,irec_offset,
     &               rt%iyy,rt%idd,rt%ihh,kdm,minval(inner_io2),
     &               maxval(inner_io2)
            end if

            if (isvelocity) then
               irec=irec+1
               if (mnproc==1) then
               inner_io2=inner_tmp2
               write(lp,'(a,i5,a,2e14.2)') 'Writing record ',irec,
     &            char15,minval(inner_io2),maxval(inner_io2)
                write(10,rec=irec) inner_io2(1:inner_inest,1:inner_jdm)
                write(11,rec=irec) 
     &             inner_io2(inner_iidm:inner_idm,1:inner_jdm)
                write(12,rec=irec) inner_io2(1:inner_idm,1:inner_inest)
                write(13,rec=irec) 
     &             inner_io2(1:inner_idm,inner_jjdm:inner_jdm)
                write(16,103,rec=irec)char3_2,klevel,irec,irec_offset,
     &               rt%iyy,rt%idd,rt%ihh,kdm,minval(inner_io2),
     &               maxval(inner_io2)
               end if
            end if

            ! Store for later - this is diagnostics
            if (mnproc==1) then
            if (char5(1:3)=='SAL') then
               if (ntecsave==1) then
                  WRITE(15,'(''TITLE= "Nested fields"'')')
                  write(15,'(a)')'VARIABLES="i" "j" "lon" "lat"'//
     &                '"depths[m]" "pivot-i" "pivot-j" '//
     &                '"variable"'
                  WRITE(15,'(a,a15,a,i4,a,i4)') 'ZONE T=',trim(char5),
     &                  ' I=',inner_idm,' J=', inner_jdm, ' F=BLOCK'
                  WRITE(15,101)((i,i=1,inner_idm),j=1,inner_jdm)
                  WRITE(15,101)((j,i=1,inner_idm),j=1,inner_jdm)
                  WRITE(15,102)((inner_lon(i,j),
     &                           i=1,inner_idm),j=1,inner_jdm)
                  WRITE(15,102)((inner_lat(i,j),
     &                           i=1,inner_idm),j=1,inner_jdm)
                  WRITE(15,102)((inner_depths(i,j),
     &                           i=1,inner_idm),j=1,inner_jdm)
                  WRITE(15,101)((inner_ipiv(i,j),
     &                           i=1,inner_idm),j=1,inner_jdm)
                  WRITE(15,101)((inner_jpiv(i,j),
     &                           i=1,inner_idm),j=1,inner_jdm)
               else
                  WRITE(15,'(a,a15,a,i4,a,i4)') 'ZONE T=',trim(char5),
     &                  ' I=',inner_idm,' J=', inner_jdm, ' F=BLOCK'
                  write(15,'(a)')'D=(1,2,3,4,5,6,7)'
               end if
               WRITE(15,102)((inner_tmp1(i,j),
     &                        i=1,inner_idm),j=1,inner_jdm)
               ntecsave=ntecsave+1
            end if
            end if



         end do ! ivar


         deallocate( inner_depths   )
         deallocate( inner_lon      )
         deallocate( inner_lat      )
         deallocate( inner_io       )
         deallocate( inner_io2      )
         deallocate( inner_tmp1     )
         deallocate( inner_tmp2     )
         deallocate( inner_tmp3     )
         deallocate( inner_a1       )
         deallocate( inner_a2       )
         deallocate( inner_a3       )
         deallocate( inner_a4       )
         deallocate( inner_ipiv     )
         deallocate( inner_jpiv     )
         deallocate( inner_tile_mask)

         if (mnproc==1) then
            close(10)
            close(11)
            close(12)
            close(13)
            close(16)
            close(15)
         end if

      else ! Could not find nestloc ...
         if (mnproc==1) then
            write(lp,*) 'Could not find'//
     &         trim(nestname(m))//'nestloc.uf'
            write(lp,*) 'I quit'
         end if
         call xcstop('nestcondsave')
         stop '(nestcondsave)'
      endif
      enddo ! Cycle grids

  101 format(30i5)      
  102 format(10e14.3)
  103 format(a3," level=",i4," record=",i5," offset=",i5," date=",i4,
     &       " ",i3," ",i2," kdm=",i3," min/max:",2e12.2)


      end subroutine nestcondsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
      end module m_nestcondsave
