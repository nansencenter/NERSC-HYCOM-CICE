      program isubaregion
      use mod_za  ! HYCOM I/O interface
      use mod_zb  ! HYCOM I/O interface for subregion
      implicit none
c
c     create a diferent-grid subregion from a full region archive file.
c
c     subregion grid is arbitrary, except that any part of this grid
c     that is outside the input grid must be land.
c
c     same resolution grid is allowed.  this is usually handled by
c     subregion.f, but use this program (or isubregion.f) when the
c     bathymetries are not identical over the subregion.
c
c     subregion same as region is allowed.  this provides a way to
c     restart a run with a different coastline and/or bathymetry,
c     i.e. original archive -> modified archive -> modified restart.
c
c     if the output grid is an integer multiple of the original, with
c     co-located p-grid points, then isubregion.f will be faster and
c     it will also produce more accurate velocity fields.
c
c     Alan J. Wallcraft,  NRL,  July 2002.
c
      character*80         :: cline,cline_u,cline_out
      character*128        :: flnm_in, flnm_tin,
     &                        flnm_out,flnm_top,flnm_reg
      integer              :: idm_out,jdm_out,
     &                        iref_out,jref_out,iref_in,jref_in,
     &                        ijgrd
      integer              :: i,ii,ijsafe,ios,ip,itmp,j,jj,jp,jq,l,
     &                        ni,nir,no,nor
      integer              :: k,l0,l1, ibadl,ibads
      integer              :: if_sm,il_sm,jf_sm,jl_sm
      logical              :: laxis(2),laxis2,lrot(2),
     &                        lperiod(2),ldebug,lcheck
      logical              :: smooth,icegln
      real                 :: hmina,hminb,hmaxa,hmaxb,
     &                        dist,dx,dy,xp,yp,rtmp,up,vp
      integer, allocatable :: m_sm(:,:),  iv_sm(:,:)
      integer, allocatable :: m_in(:,:),  m_out(:,:),   m_osm(:,:)
      real,    allocatable :: a_in(:,:),  a_out(:,:)
      real,    allocatable :: u_in(:,:),  u_out(:,:)
      real,    allocatable :: v_in(:,:),  v_out(:,:)
      real,    allocatable :: t_in(:,:),  t_out(:,:,:)
      real,    allocatable :: p_in(:,:,:),p_out(:,:,:)
c
      real,    allocatable :: plat_in( :,:),plon_in( :,:),pang_in( :,:)
      real,    allocatable :: plat_out(:,:),plon_out(:,:),pang_out(:,:)
      integer, allocatable :: i_out(:,:),j_out(:,:)
      real,    allocatable :: x_out(:,:),y_out(:,:)
      real,    allocatable :: plat_in_min(:),plat_in_max(:)
c
      real,    parameter   :: hspval=0.5*2.0**100  ! half spval
      real,    parameter   :: onem=9806.0          ! g/thref
      real,    parameter   :: tenm=10.0*onem
c
      call xcspmd
      call zaiost
      call blkdat(idm_out,jdm_out,
     &            smooth,icegln,
     &            flnm_reg,flnm_in,flnm_tin,flnm_out,flnm_top,cline_out)
      call zbiost(idm_out,jdm_out)
c
      allocate(  iv_sm(jdm,2) )
c
      allocate(   m_sm(idm,jdm),     m_osm(idm_out,jdm_out) )
      allocate(   m_in(idm,jdm),     m_out(idm_out,jdm_out) )
      allocate(   a_in(idm,jdm),     a_out(idm_out,jdm_out) )
      allocate(   u_in(idm,jdm),     u_out(idm_out,jdm_out) )
      allocate(   v_in(idm,jdm),     v_out(idm_out,jdm_out) )
      allocate(   t_in(idm,jdm),     t_out(idm_out,jdm_out,3) )
      allocate(   p_in(idm,jdm,0:1), p_out(idm_out,jdm_out,0:1) )
      allocate( plat_in(idm,jdm), plat_out(idm_out,jdm_out) )
      allocate( plon_in(idm,jdm), plon_out(idm_out,jdm_out) )
      allocate( pang_in(idm,jdm), pang_out(idm_out,jdm_out) )
c                                                            
      allocate( plat_in_max(jdm), plat_in_min(jdm) )         
c                                                   
      allocate( i_out(idm_out,jdm_out), j_out(idm_out,jdm_out) )
      allocate( x_out(idm_out,jdm_out), y_out(idm_out,jdm_out) )
c
c     get the output p-grid mask from the bathymetry.
c
      l  = len_trim(flnm_top)
      open (unit=13,file=flnm_top(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
      write(lp,'(/a)') 'OUTPUT BATHYMETRY:'
      do i= 1,6
        read( 13,'(a)') cline
        write(lp,'(a)') cline(1:len_trim(cline))
      enddo
      l = index(cline,'=')
      read (cline(l+1:),*)  hminb,hmaxb
      close(unit=13)
c
      l  = len_trim(flnm_top)
      call zbiopf(flnm_top(1:l-2)//'.a','old', 13)
      call zbiord(t_out(1,1,1),m_out,.false., hmina,hmaxa, 13)
      call zbiocl(13)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - output bathymetry .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif
c
      do j= 1,jdm_out
        do i= 1,idm_out
          if     (t_out(i,j,1).gt.hspval .or.
     &            t_out(i,j,1).le.0.0        ) then
            t_out(i,j,1) = 0.0
            m_out(i,j)   = 0
          else
            t_out(i,j,1) = t_out(i,j,1)*onem
            m_out(i,j)   = 1
          endif
        enddo
      enddo
c
c     get the input p-grid mask from the bathymetry.
c
      l  = len_trim(flnm_tin)
      open (unit=13,file=flnm_tin(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
      write(lp,'(/a)') ' INPUT BATHYMETRY:'
      do i= 1,6
        read( 13,'(a)') cline
        write(lp,'(a)') cline(1:len_trim(cline))
      enddo
      l = index(cline,'=')
      read (cline(l+1:),*)  hminb,hmaxb
      close(unit=13)
c
      l  = len_trim(flnm_tin)
      call zaiopf(flnm_tin(1:l-2)//'.a','old', 13)
      call zaiord(t_in,m_in,.false., hmina,hmaxa, 13)
      call zaiocl(13)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - bathymetry .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif
c
      do j= 1,jdm
        do i= 1,idm
          if     (t_in(i,j).gt.hspval .or.
     &            t_in(i,j).le.0.0        ) then
            t_in(i,j) = 0.0
            m_in(i,j) = 0
            m_sm(i,j) = 0
          else
            t_in(i,j) = t_in(i,j)*onem
            m_in(i,j) = 1
            m_sm(i,j) = 1
          endif
        enddo
      enddo
c       
c     read the input and output grid locations (no error checking).
c
      nir = 24
      call zaiopf('regional.grid.a','old', nir)
      call zaiord(plon_in,m_in,.false., hmina,hmaxa, nir)
      call zaiord(plat_in,m_in,.false., hmina,hmaxa, nir)
      call zaiosk(nir) !qlon
      call zaiosk(nir) !qlat
      call zaiosk(nir) !ulon
      call zaiosk(nir) !ulat
      call zaiosk(nir) !vlon
      call zaiosk(nir) !vlat
      call zaiord(pang_in,m_in,.false., hmina,hmaxa, nir)
      call zaiocl(nir)
c     
      nor = 25
      call zbiopf(trim(flnm_reg),'old', nor)
      call zbiord(plon_out,m_out,.false., hmina,hmaxa, nor)
      call zbiord(plat_out,m_out,.false., hmina,hmaxa, nor)
      call zbiosk(nor) !qlon
      call zbiosk(nor) !qlat
      call zbiosk(nor) !ulon
      call zbiosk(nor) !ulat
      call zbiosk(nor) !vlon
      call zbiosk(nor) !vlat
      call zbiord(pang_out,m_out,.false., hmina,hmaxa, nor)
      call zbiocl(nor)
c     
c     calculate the output grid location w.r.t. the input grid.
c
      lperiod(1) = maxval(plon_in(:,:))-
     &             minval(plon_in(:,:))  .gt. 350.0
      if     (lperiod(1)) then
        write(6,'(a)') ' input domain assumed to be periodic'
      else
        write(6,'(a)') ' input domain assumed to be non-periodic'
      endif
c
      lperiod(2) = maxval(plon_out(:,:))-
     &             minval(plon_out(:,:))   .gt. 350.0
      if     (lperiod(2)) then
        write(6,'(a)') 'output domain assumed to be periodic'
      else
        write(6,'(a)') 'output domain assumed to be non-periodic'
      endif
c       
      laxis(1) = .true.
      do i= 2,idm
        laxis(1) = laxis(1) .and.
     &             maxval(abs(plat_in(1,:)-plat_in(i,:))).le.1.e-2
      enddo
      do j= 2,jdm
        laxis(1) = laxis(1) .and.
     &             maxval(abs(plon_in(:,1)-plon_in(:,j))).le.1.e-2
      enddo
      if     (laxis(1)) then
        write(6,'(a)') ' input domain has 1-d lat/lon axes'
      else
        write(6,'(a)') ' input domain is curvilinear'
      endif
c
      lrot(1) = minval(pang_in(:,:)).ne.0.0 .or.
     &          maxval(pang_in(:,:)).ne.0.0
      if     (.not.lrot(1)) then
        write(6,'(a)') ' input domain is     aligned E-W,N-S'
      else
        write(6,'(a)') ' input domain is not aligned E-W,N-S'
      endif
c
      lrot(2) = minval(pang_out(:,:)).ne.0.0 .or.
     &          maxval(pang_out(:,:)).ne.0.0
      if     (.not.lrot(2)) then
        write(6,'(a)') 'output domain is     aligned E-W,N-S'
      else
        write(6,'(a)') 'output domain is not aligned E-W,N-S'
      endif
c       
      laxis(2) = .true.
      do i= 2,idm_out
        laxis(2) = laxis(2) .and.
     &             maxval(abs(plat_out(1,:)-plat_out(i,:))).le.1.e-2
      enddo
      do j= 2,jdm_out
        laxis(2) = laxis(2) .and.
     &             maxval(abs(plon_out(:,1)-plon_out(:,j))).le.1.e-2
      enddo
      if     (laxis(2)) then
        write(6,'(a)') 'output domain has 1-d lat/lon axes'
      else
        write(6,'(a)') 'output domain is curvilinear'
      endif
      write(6,*)
      call flush(6)
c     
      laxis2 = laxis(1) .and. laxis(2)
c     
      do j= 1,jdm
        plat_in_min(j) = minval(plat_in(:,j))
        plat_in_max(j) = maxval(plat_in(:,j))
      enddo
c       
      ldebug = .false.
c     
      do jj= 1,jdm_out
        do ii= 1,idm_out
*         ldebug = mod(ii,20).eq.1 .and. mod(jj,20).eq.1
*         ldebug = mod(ii,50).eq.1 .and. mod(jj,50).eq.1
          if     (laxis2 .and. ii.ne.1 .and. jj.ne.1) then
c         
c           shortcut for 1-d axes.
c
            i_out(ii,jj) = i_out(ii, 1)
            j_out(ii,jj) = j_out( 1,jj)
            x_out(ii,jj) = x_out(ii, 1)
            y_out(ii,jj) = y_out( 1,jj)
            cycle
          endif
c           
c         find the nearest point by exhaustive search, but improve
c         efficiency by using plat_*_min/_max to exclude far away rows.
c
          xp   = plon_out(ii,jj)
          yp   = plat_out(ii,jj)
c         
c         start with a nearby point.
c
          if     (ii.eq.1) then
            if     (jj.eq.1) then
              ip = 1
              jp = 1
            else
              ip = i_out(1,jj-1)
              jp = j_out(1,jj-1)
            endif
          endif
          dy =      abs(plat_in(ip,jp) - yp)
          dx = mod( abs(plon_in(ip,jp) - xp), 360.0 )
          if     (dx.gt.180.0) then
            dx = 360.0 - dx
          endif
          dist = dx+dy
          if     (ldebug) then
            write(6,'(a,4i5,3f9.2)')
     &        'ii,jj,ip,jp,dx,dy,dist = ',
     &         ii,jj,ip,jp,dx,dy,dist
            call flush(6)
          endif
c           
          do jq= 0,jdm
            if     (jq.eq.0) then
              j = jp  ! search estimated row location first
            else
              j = jq
            endif
            dy = min( abs(plat_in_min(j)-yp),
     &                abs(plat_in_max(j)-yp) )
            if     (dy.gt.dist) then
              if     (ldebug) then
                write(6,'(a,3i5,2f9.2)')
     &            'ii,jj,j,dy,dist (cycle)',
     &             ii,jj,j,dy,dist
                call flush(6)
              endif
              cycle  ! far away row
            endif
            do i = 1,idm
              dy =      abs(plat_in(i,j) - yp)
              dx = mod( abs(plon_in(i,j) - xp), 360.0 )
              if     (dx.gt.180.0) then
                dx = 360.0 - dx
              endif
              if     (dx+dy.le.dist) then
                ip   = i
                jp   = j
                dist = dx+dy
                if     (ldebug) then
                  write(6,'(a,4i5,3f9.2)')
     &              'ii,jj,ip,jp,dx,dy,dist = ',
     &               ii,jj,ip,jp,dx,dy,dist
                  call flush(6)
                endif
              endif
            enddo
          enddo
c           
c         convert nearest point into bilinear cell and distances.
c
          if     (ldebug) then
            dy =      abs(plat_in(ip,jp)-yp)
            dx = mod( abs(plon_in(ip,jp)-xp), 360.0 )
            if     (dx.gt.180.0) then
              dx = 360.0 - dx
            endif
            write(6,'(a,4i5,3f9.2)')
     &        'ii,jj,i_,j_,dx,dy = ',
     &         ii,jj,ip,jp,dx,dy
            call flush(6)
          endif
          i = ip
          j = jp
          if     (j.eq.jdm) then
            jp = jp-1
          elseif (j.ne.1)   then
            dy =               (yp-plat_in(ip,jp))/
     &           (plat_in(ip,jp+1)-plat_in(ip,jp))
            if     (dy.lt.0.0) then
              jp = jp-1
            endif
          endif
          if     (i.eq.idm) then
            if     (lperiod(1)) then
              dx = mod( plon_in(ip,jp)-xp, 360.0 )
              if     (dx.lt.-180.0) then
                dx = 360.0 + dx
              elseif (dx.gt. 180.0) then
                dx = 360.0 - dx
              endif
              dy = mod( plon_in(ip,jp)-plon_in(   1,jp), 360.0 )
              if     (dy.lt.-180.0) then
                dy = 360.0 + dy
              elseif (dy.gt. 180.0) then
                dy = 360.0 - dy
              endif
              if     (abs(dy).gt.0.001) then
                if (dx/dy.lt.0.0) then
                  ip = ip-1
                endif
              endif  !abs(dy).gt.0.001
            else
              ip = ip-1
            endif  !lperiod:else
          else
            dx = mod( plon_in(ip,jp)-xp, 360.0 )
            if     (dx.lt.-180.0) then
              dx = 360.0 + dx
            elseif (dx.gt. 180.0) then
              dx = 360.0 - dx
            endif
            dy = mod( plon_in(ip,jp)-plon_in(ip+1,jp), 360.0 )
            if     (dy.lt.-180.0) then
              dy = 360.0 + dy
            elseif (dy.gt. 180.0) then
              dy = 360.0 - dy
            endif
            if     (abs(dy).gt.0.001) then
              if     (dx/dy.lt.0.0) then
                if     (ip.ne.1) then
                  ip = ip-1
                elseif (lperiod(1)) then
                  ip = idm
                endif
              endif  !dx/dy.lt.0.0
            endif  !abs(dy).gt.0.001
          endif  !i.eq.idm:else
          i_out(ii,jj) = ip
          j_out(ii,jj) = jp
          y_out(ii,jj) =               (yp-plat_in(ip,jp))/
     &                   (plat_in(ip,jp+1)-plat_in(ip,jp))
          if     (ldebug) then
            write(6,'(a,2i5,4f9.2)')
     &        'ii,jj,yp,plat,plat+,y_',
     &         ii,jj,yp,plat_in(ip,jp),plat_in(ip,jp+1),y_out(ii,jj)
            call flush(6)
          endif
          dx = mod( plon_in(ip,jp)-xp, 360.0 )
          if     (dx.lt.-180.0) then
            dx = 360.0 + dx
          elseif (dx.gt. 180.0) then
            dx = 360.0 - dx
          endif
          dy = mod( plon_in(ip,jp)-plon_in(ip+1,jp), 360.0 )
          if     (dy.lt.-180.0) then
            dy = 360.0 + dy
          elseif (dy.gt. 180.0) then
            dy = 360.0 - dy
          endif
          if     (abs(dy).gt.0.001) then
            x_out(ii,jj) = dx/dy
          else
            x_out(ii,jj) = 0.0
          endif
          if     (ldebug) then
            write(6,'(a,4i5,2f9.2,2f9.3)')
     &        'ii,jj,i_,j_,x_,y_ = ',
     &         ii,jj,i_out(ii,jj),j_out(ii,jj),
     &               x_out(ii,jj),y_out(ii,jj),dx,dy
            call flush(6)
          endif
        enddo
      enddo
