      program confmap_gmap
      use mod_za
      use mod_confmap
      implicit none
c
c     create an array index map to a diferent-grid subregion from a 
c     full region hycom file.
c
c     subregion grid is arbitrary, but must be based on confmap 
c     Source grid must be a confmap grid
c     Files:  
c        grid.info          : conformal mapping configuration 
c                             file for input grid
c        regional.grid.[ab] : hycom grid files for target grid
c
c     Based on routine isuba_gmapi by Alan J. Wallcraft.
c
c
      character*80         :: chead
      character*256        :: flnm_map,flnm_reg
c
      !real,    allocatable :: plat_in( :,:),plon_in( :,:)
      real,    allocatable :: plat_out(:,:),plon_out(:,:),
     &   x_out(:,:), y_out(:,:)
      integer, allocatable :: m_out(:,:)
      integer :: idm_in, jdm_in, nir, nor, icnt, i, j, ipiv, jpiv, no, l
      real                 :: hmina,hminb,hmaxa,hmaxb,
     &                        deg2rad,dist,dist_max,plat_max,plat_min,
     &                        dx,dy,qdx,xp,yp, lat_n, lon_n,x ,y
      logical :: ongrid

c
      call blkdat(  chead, flnm_map)
      call xcspmd
      call zaiost
      allocate( plat_out(idm,jdm), plon_out(idm,jdm), m_out(idm,jdm) )
      allocate( x_out(idm,jdm), y_out(idm,jdm))
c
c       
c     read the input and output grid locations (no error checking).
c
      nir = 24
      call zaiopf('regional.grid.a','old', nir)
      call zaiord(plon_out,m_out,.false., hmina,hmaxa, nir)
      call zaiord(plat_out,m_out,.false., hmina,hmaxa, nir)
      call zaiocl(nir)

c --- Initialize mapping of input grid
      call initconfmap(idm_in,jdm_in,.true.)

c --- Go through target grid, and map from source grid. Find nearest point
      icnt=1
      do j=1,jdm
      do i=1,idm

         call ll2gind(plon_out(i,j),plat_out(i,j),x,y)
         ipiv=floor(x)
         jpiv=floor(y)

         if (ipiv>=1.and.jpiv>=1.and.ipiv<=idm_in.and.jpiv<=jdm_in) then
            ongrid=.true.
            x_out(i,j) = x
            y_out(i,j) = y
            m_out(i,j)=1
         else 
            ongrid=.false.
            x_out(i,j) = 0
            y_out(i,j) = 0
            m_out(i,j)=0

            ! Uncomment this if you want diagnostic output in gmap ab files
            print *,"Target point is outside of source grid:",i,j
            stop
         end if


         if (mod(icnt,10000)==0) 
C        if (mod(icnt,1)==0) 
     &      print '(2i6, 4f14.4,l12)',i,j,plat_out(i,j),plon_out(i,j),
     &             x,y,ongrid
C
         icnt = icnt + 1
      end do
      end do

      no = 15
      l  = len_trim(flnm_map)
      open (unit=no,file=flnm_map(1:l-2)//'.b',form='formatted',
     .      status='replace',action='write')
      call zaiopf(flnm_map(1:l-2)//'.a','replace', no)
c
      write(no,'(a)') trim(chead)
      call flush(no)
      write(6, *)
      write(6, *)
      write(6, '(a)') trim(chead)
      call flush(6)
c
      call zaiowr(x_out,m_out,.true., hmina,hmaxa, no, .false.)
      write(no,'(a,2f12.4)') 'xmap:  min,max =',hmina,hmaxa
      call flush(no)
      write(6, '(a,2f12.4)') 'xmap:  min,max =',hmina,hmaxa
      call flush(6)
c
      call zaiowr(y_out,m_out,.true., hmina,hmaxa, no, .false.)
      write(no,'(a,2f12.4)') 'ymap:  min,max =',hmina,hmaxa
      call flush(no)
      write(6, '(a,2f12.4)') 'ymap:  min,max =',hmina,hmaxa
      call flush(6)
c
      close(no)
      call zaiocl(no)
      write(6, '(a,2f12.4)') 'xmap:  min,max =',hmina,hmaxa
      call flush(6)

      end program

      subroutine blkdat(  chead, flnm_map)
      use mod_xc  ! HYCOM communication interface
      implicit none
      character*80  :: chead
      character*256 :: flnm_reg,flnm_map
c
c --- read blkdat.input for interpolated subregion.
c
c --- 'flnm_map'  = output grid map   filename
c --- 'chead   '  = single line header for grid map file
c
      read( *,'(a)')      flnm_map
      write(6,'(a)') trim(flnm_map)
      write(6,*)
      read( *,'(a)')      chead
      write(6,'(a)') trim(chead)
      write(6,*)
      call flush(6)
c
c --- 'idm   ' = output longitudinal array size
c --- 'jdm   ' = output latitudinal  array size
c --- 'maxinc' = maximum input array index jump on target grid
c
      !call blkini(idm_out,   'idm   ')
      !call blkini(jdm_out,   'jdm   ')
      !call blkini(maxinc,    'maxinc')
      write(6,*)
      call flush(6)
c
      return
      end subroutine blkdat
