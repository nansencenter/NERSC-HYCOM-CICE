      program isuba_topog
      use mod_za  ! HYCOM I/O interface
      use mod_zb  ! HYCOM I/O interface for subregion
      implicit none
c
c     create a different-grid subregion bathymetry from
c     a full region HYCOM 2.0 bathymetry.
c
c     subregion grid is arbitrary, but any part of this grid
c     that is outside the input grid will be land.
c
c     if the output grid is an integer multiple finer than
c     the input grid, then isub_topog.f will be slightly faster.
c
c     Alan J. Wallcraft,  NRL,  July 2002.
c
      character*79         :: preambl(5)
      character*80         :: cline,cline_out
      character*128        :: flnm_in,flnm_out,flnm_reg
      logical              :: laxis(2),laxis2,lperiod,ldebug
      integer              :: idm_out,jdm_out
      integer              :: i,ii,ip,j,jj,jp,jq,l,ni,nir,no,nor,nsmooth
      integer              :: iidebug,jjdebug
      real                 :: hmina,hminb,hmaxa,hmaxb,dist,dx,dy,xp,yp
      integer, allocatable :: m_in(:,:),m_sm(:,:),m_out(:,:)
      real,    allocatable :: a_in( :,:),plat_in( :,:),plon_in( :,:)
      real,    allocatable :: a_out(:,:),plat_out(:,:),plon_out(:,:)
      integer, allocatable :: i_out(:,:),j_out(:,:)
      real,    allocatable :: x_out(:,:),y_out(:,:)
      real,    allocatable :: plat_in_min(:),plat_in_max(:)
      real,    allocatable :: plat_in_minall(:),plat_in_maxall(:)
c
      real,    parameter   :: hspval=0.5*2.0**100
c
      call xcspmd
      call zaiost
      call blkdat(idm_out,jdm_out,
     &            flnm_reg,flnm_in,flnm_out,cline_out)
      call zbiost(idm_out,jdm_out)
c
      allocate(    m_sm(idm,jdm) )
      allocate(    m_in(idm,jdm),    m_out(idm_out,jdm_out) )
      allocate(    a_in(idm,jdm),    a_out(idm_out,jdm_out) )
      allocate( plat_in(idm,jdm), plat_out(idm_out,jdm_out) )
      allocate( plon_in(idm,jdm), plon_out(idm_out,jdm_out) )
c
      allocate( plat_in_max(jdm), plat_in_min(jdm) )
      allocate( plat_in_maxall(jdm), plat_in_minall(jdm) )
c
      allocate( i_out(idm_out,jdm_out), j_out(idm_out,jdm_out) )
      allocate( x_out(idm_out,jdm_out), y_out(idm_out,jdm_out) )
c
c     open input and output bathymetry files.
c
      ni = 14
      l  = len_trim(flnm_in)
      open (unit=ni,file= flnm_in(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
      call zaiopf( flnm_in(1:l-2)//'.a','old', ni)
c
      no = 15
      l  = len_trim(flnm_out)
      open (unit=no,file=flnm_out(1:l-2)//'.b',form='formatted',
     .      status='new',action='write')
      call zbiopf(flnm_out(1:l-2)//'.a','new', no)
c
c     process the bathymetry header
c
      read( ni,'(a79)') preambl
      preambl(1) = cline_out(1:79)
      write(no,'(a79)') preambl
      call flush(no)
      write(lp,*)
      write(lp,'(a79)') preambl
      write(lp,*)
      call flush(lp)
c
c     input the full domain bathymetry.
c
      read( ni,'(a)') cline
      l = index(cline,'=')
      read (cline(l+1:),*)  hminb,hmaxb
      call zaiord(a_in,m_in,.false., hmina,hmaxa, ni)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - full domain .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif
c
      do j= 1,jdm
        do i= 1,idm
          if     (a_in(i,j).gt.hspval .or.
     &            a_in(i,j).le.0.0        ) then
            a_in(i,j) = 0.0
            m_in(i,j) = 0
          else
            m_in(i,j) = 1
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
      call zaiocl(nir)
c
      nor = 25
      call zbiopf(trim(flnm_reg),'old', nor)
      call zbiord(plon_out,m_out,.false., hmina,hmaxa, nor)
      call zbiord(plat_out,m_out,.false., hmina,hmaxa, nor)
      call zbiocl(nor)
c
c     calculate the output grid location w.r.t. the input grid.
c
      lperiod = maxval(plon_in(:,:))-minval(plon_in(:,:)) .gt. 350.0
      if     (lperiod) then
        write(6,'(a)') ' input domain assumed to be periodic'
      else
        write(6,'(a)') ' input domain assumed to be non-periodic'
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
        plat_in_min(j) = minval(plat_in(:,j))  !min in row
        plat_in_max(j) = maxval(plat_in(:,j))  !max in row
      enddo
      do j= 1,jdm
        plat_in_minall(j) = minval(plat_in_min(j:jdm))  !min row and above
        plat_in_maxall(j) = maxval(plat_in_max(j:jdm))  !max row and above
      enddo
c
      ldebug  = .false.
*
      iidebug = min(idm_out/5, 500)
      jjdebug = min(jdm_out/5, 500)
c
      do jj= 1,jdm_out
        do ii= 1,idm_out
          ldebug = mod(ii,iidebug).eq.1 .and. mod(jj,jjdebug).eq.1
c
          if     (laxis2 .and. ii.ne.1 .and. jj.ne.1) then
c
c           shortcut for 1-d axes.
c
            i_out(ii,jj) = i_out(ii,  jj-1)
            j_out(ii,jj) = j_out(ii-1,jj)
            x_out(ii,jj) = x_out(ii,  jj-1)
            y_out(ii,jj) = y_out(ii-1,jj)
            if     (i_out(ii,jj).ne.1 .and. j_out(ii,jj).ne.1) then
              if     (ldebug) then
                write(6,'(a,4i5,2f9.2)')
     &            'ii,jj,i_,j_,x_,y_ = ',
     &             ii,jj,i_out(ii,jj),j_out(ii,jj),
     &                   x_out(ii,jj),y_out(ii,jj)
                call flush(6)
              endif
              cycle
            endif
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
            if     (dy.le.plat_in_minall(j)-dist .or.
     &              dy.ge.plat_in_maxall(j)+dist     ) then
*             if     (ldebug) then
*               write(6,'(a,3i5,f9.2)')
*    &            'ii,jj,j,dist (exit)',
*    &             ii,jj,j,dist
*               call flush(6)
*             endif
              exit  !already found nearest point
            endif
            if     (dy.le.plat_in_min(j)-dist .or.
     &              dy.ge.plat_in_max(j)+dist     ) then
*             if     (ldebug) then
*               write(6,'(a,3i5,f9.2)')
*    &            'ii,jj,j,dist (cycle)',
*    &             ii,jj,j,dist
*               call flush(6)
*             endif
              cycle  ! far away row
            endif
            do i = 1,idm
              dy =      abs(plat_in(i,j) - yp)
              dx = mod( abs(plon_in(i,j) - xp), 360.0 )
              if     (dx.gt.180.0) then             
                dx = 360.0 - dx                     
              endif           
              if     (dx+dy.lt.dist) then
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
            enddo !i
          enddo !jq
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
            if     (lperiod) then
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
                elseif (lperiod) then
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
        enddo !ii
      enddo !jj
